%% Nicholas Czarnek
% SSPACISS Laboratory, Duke University
% 24 March 2015
%
% runParsonAnalysis.m
% The purpose of this function is to run the Parson analysis on an input
% dataset.  It is assumed that only one device is sent in, so that use is
% the first column, and the device is the second.

function [errorStruct,assignedPower] = runParsonAnalysis(eDS,varargin)

%% Handle varargin.
% Set up the appliance prior models
options.fridgeStateMeans = [2 180 160];
options.fridgeStateCovs = [5 100 100];
options.fridgeTrans = [0.95 0.03 0.02; 0 0 1; 0.2 0 0.8];

options.mwStateMeans = [4 1700];
options.mwStateCovs = [100 1000];
options.mwTrans = [0.99 0.01; 1 0];

options.wdStateMeans = [0 5000];
options.wdStateCovs = [100 5000];
options.wdTrans = [0.9 0.1; 0.5 0.5];

options.dwStateMeans = [0 1400];
options.dwStateCovs = [100 1300];
options.dwTrans = [0.9999 0.0001; 0.3 0.7];

options.acStateMeans = [4 2300];
options.acStateCovs = [100 2300];
options.acTrans = [0.99 0.01; 0.9 0.1];

options.fridgeOn = true;
options.mwOn = false;
options.wdOn = false;
options.dwOn = false;
options.acOn = false;

options.startingPoint = 1;
% Default window lengths for REDD (should not change for others):
% fridge -> 200
% mw -> 6
% wd ->20
% dw -> 30
% ac -> 20
options.fridgeWindowLength = 200; % default for the refrigerator
options.mwWindowLength = 6;
options.wdWindowLength = 20;
options.dwWindowLength = 30;
options.acWindowLength = 20;

% Default training_length for all:
% fridge - 3000
% mw - 4000
% wd - 5000
% dw - 5000
% ac - 2500
options.fridgeTrainingLength = 3000;
options.mwTrainingLength = 4000;
options.wdTrainingLength = 5000;
options.dwTrainingLength  = 5000;
options.acTrainingLength  = 2500;

% Default likelihood thresholds:
% fridge - 0.01
% mw - 0.00001
% wd - 0.001
% dw - 0.001
% ac - 0.001
options.fridgeLikThres = 0.01;
options.mwLikThres = 0.00001;
options.wdLikThres = 0.001;
options.dwLikThres = 0.001;
options.acLikThres = 0.001;

% Default number of training windows
% fridge - 3
% mw - 4
% wd - 3
% dw - 2
% ac - 4
options.fridgeNumOfWindows = 3;
options.mwNumOfWindows = 4;
options.wdNumOfWindows = 3;
options.dwNumOfWindows = 2;
options.acNumOfWindows = 4;

parsedOuts = prtUtilSimpleInputParser(options,varargin);

%%
fridgeStateMeans = parsedOuts.fridgeStateMeans;
fridgeStateCovs = parsedOuts.fridgeStateCovs;
fridgeTrans = parsedOuts.fridgeTrans;
mwStateMeans = parsedOuts.mwStateMeans;
mwStateCovs = parsedOuts.mwStateCovs;
mwTrans = parsedOuts.mwTrans;
wdStateMeans = parsedOuts.wdStateMeans;
wdStateCovs = parsedOuts.wdStateCovs;
wdTrans = parsedOuts.wdTrans;
dwStateMeans = parsedOuts.dwStateMeans;
dwStateCovs = parsedOuts.dwStateCovs;
dwTrans = parsedOuts.dwTrans;
acStateMeans = parsedOuts.acStateMeans;
acStateCovs = parsedOuts.acStateCovs;
acTrans = parsedOuts.acTrans;
fridgeOn = parsedOuts.fridgeOn;
mwOn = parsedOuts.mwOn;
wdOn = parsedOuts.wdOn;
dwOn = parsedOuts.dwOn;
acOn = parsedOuts.acOn;
startingPoint = parsedOuts.startingPoint;
% windowLength = parsedOuts.windowLength;
% trainingLength = parsedOuts.trainingLength;
% likThres = parsedOuts.likThres;
% numOfWindows = parsedOuts.numOfWindows;

%% Define the state means, covariances, and transition matrix.
if fridgeOn
    stateMeans = fridgeStateMeans;
    stateCovs = fridgeStateCovs;
    stateTrans = fridgeTrans;
    windowLength = parsedOuts.fridgeWindowLength;
    trainingLength = parsedOuts.fridgeTrainingLength;
    likThres = parsedOuts.fridgeLikThres;
    numOfWindows = parsedOuts.fridgeNumOfWindows;
elseif mwOn
    stateMeans = mwStateMeans;
    stateCovs = mwStateCovs;
    stateTrans = mwTrans;
    windowLength = parsedOuts.mwWindowLength;
    trainingLength = parsedOuts.mwTrainingLength;
    likThres = parsedOuts.mwLikThres;
    numOfWindows = parsedOuts.mwNumOfWindows;
elseif wdOn
    stateMeans = wdStateMeans;
    stateCovs = wdStateCovs;
    stateTrans = wdTrans;
    windowLength = parsedOuts.wdWindowLength;
    trainingLength = parsedOuts.wdTrainingLength;
    likThres = parsedOuts.wdLikThres;
    numOfWindows = parsedOuts.wdNumOfWindows;
elseif dwOn
    stateMeans = dwStateMeans;
    stateCovs = dwStateCovs;
    stateTrans = dwTrans;
    windowLength = parsedOuts.dwWindowLength;
    trainingLength = parsedOuts.dwTrainingLength;
    likThres = parsedOuts.dwLikThres;
    numOfWindows = parsedOuts.dwNumOfWindows;
elseif acOn
    stateMeans = acStateMeans;
    stateCovs = acStateCovs;
    stateTrans = acTrans;
    windowLength = parsedOuts.acWindowLength;
    trainingLength = parsedOuts.acTrainingLength;
    likThres = parsedOuts.acLikThres;
    numOfWindows = parsedOuts.acNumOfWindows;
end    

% always_on = {0, 3, 0, 0, 0, 4};

residue = eDS.data(:,1);
appliancePower = eDS.data(:,2);

%% Only deal with the on times for the current dataset
% Take care of this before sending it in.
% kL = [eDS.observationInfo.keepLogicals]';
% eDS = eDS.retainObservations(kL);

%% Calculate the difference model parameters.
nStates = numel(stateMeans);
init = ones(1,nStates)/nStates;
idx = repmat(1:nStates,nStates,1);
emitMean = stateMeans(idx) - stateMeans(idx');
emitCov = stateCovs(idx) + stateCovs(idx');

%% Make the difference hmm
priorDhmmBnet = make_dhmm(init,stateMeans,stateCovs,...
    emitMean,emitCov,...
    stateTrans);

priorHmmBnet = make_hmm(init,stateMeans,stateCovs,stateTrans);

observedPower = residue;
diffs = [0 diff(observedPower)'];
evidence = {};
evidence(2,:) = num2cell(diffs);
T = length(evidence);

%% Create smoothing engine
engine = smoother_engine(jtree_2TBN_inf_engine(priorDhmmBnet));

%% Assume trainng type is 2
trainingData = find_training_ranges_generic(...
    residue(startingPoint:startingPoint + trainingLength - 1),...
    windowLength,priorDhmmBnet,numOfWindows);

%% Train the models
[trainedDhmmBnet,loglik] = learn_params_generic(priorDhmmBnet,diff(trainingData));
try
    [trainedHmmBnet,hmmLoglik] = learn_params_generic(priorHmmBnet,trainingData);
catch err
    trainedHmmBnet = priorHmmBnet;
end
hmmEmissions = struct(trainedHmmBnet.CPD{2});
trainedDhmmBnet.CPD{1} = tabular_CPD(trainedDhmmBnet,1,'CPT',ones(1,length(init))/length(init));
trainedDhmmBnet.CPD{2} = trainedHmmBnet.CPD{2};
trainedDhmmBnet.CPD{3} = priorDhmmBnet.CPD{3};

%% Infer power with the Viterbi algorithm
evidence(3,:) = num2cell(observedPower);
[mpe,ll,max_prob,ignored_obs] = my_viterbi_diff2(trainedDhmmBnet,evidence,1,likThres);

ignored_obs = logical(ignored_obs);

%% What power was inferred from this calculation?
[~,stateAssignments] = max(max_prob);

structHmm = struct(trainedHmmBnet.CPD{2});
stateMeans = structHmm.mean;

assignedPower = zeros(T,1);

for stateIdx = unique(stateAssignments)
    assignedPower(stateAssignments == stateIdx) = stateMeans(stateIdx);
end

%% Total disaggregation error
totalDisaggregationErrorInclusive = abs(sum(assignedPower) - sum(appliancePower))/sum(appliancePower);

%% RMS error
rmsErrorInclusive = sqrt(1/T * sum((assignedPower - appliancePower).^2));

%% Account for the ignored observations
assignedPower(ignored_obs) = 0;
appliancePower(ignored_obs) = 0;

totalDisaggregationError = abs(sum(assignedPower) - sum(appliancePower))/sum(appliancePower);

rmsError = sqrt(1/T * sum((assignedPower - appliancePower).^2));

errorStruct.totalDisaggregationErrorInclusive = totalDisaggregationErrorInclusive;
errorStruct.totalDisaggregationError = totalDisaggregationError;
errorStruct.rmsErrorInclusive = rmsErrorInclusive;
errorStruct.rmsError = rmsError;
errorStruct.ignoredObs = ignored_obs;