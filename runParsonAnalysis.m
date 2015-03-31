%% Nicholas Czarnek
% SSPACISS Laboratory, Duke University
% 24 March 2015
%
% runParsonAnalysis.m
% The purpose of this function is to run the Parson analysis on an input
% dataset.  It is assumed that only one 

function runParsonAnalysis(eDS,varargin)

%% Handle varargin.
% Set up the appliance prior models
options.fridge_state_means = [2 180 160];
options.fridge_state_covs = [5 100 100];
options.fridge_trans = [0.95 0.03 0.02; 0 0 1; 0.2 0 0.8];

options.micro_state_means = [4 1700];
options.micro_state_covs = [100 1000];
options.micro_trans = [0.99 0.01; 1 0];

options.wd_state_means = [0 5000];
options.wd_state_covs = [100 5000];
options.wd_trans = [0.9 0.1; 0.5 0.5];

options.dw_state_means = [0 1400];
options.dw_state_covs = [100 1300];
options.dw_trans = [0.9999 0.0001; 0.3 0.7];

options.ac_state_means = [4 2300];
options.ac_state_covs = [100 2300];
options.ac_trans = [0.99 0.01; 0.9 0.1];

options.fridgeOn = true;
options.mwOn = false;
options.wdOn = false;
options.dwOn = false;
options.acOn = false;

options.starting_point = 1;

% Device signals
options.fridge = [];
options.mw = [];
options.wd = [];
options.dw = [];
options.ac = [];


% Default window lengths for REDD (should not change for others):
% fridge -> 200
% mw -> 6
% wd ->20
% dw -> 30
% ac -> 20
options.window_length = 200;

% Default training_length for all:
% fridge - 3000
% mw - 4000
% wd - 5000
% dw - 5000
% ac - 2500
options.training_length = 3000;

% Default likelihood thresholds:
% fridge - 0.01
% mw - 0.00001
% wd - 0.001
% dw - 0.001
% ac - 0.001
options.lik_thres = 0.01;

% Default number of training windows
% fridge - 3
% mw - 4
% wd - 3
% dw - 2
% ac - 4
options.num_of_windows = 3;

parsedOuts = prtUtilSimpleInputParser(options,varargin);

%%
fridge_state_means = parsedOuts.fridge_state_means;
fridge_state_covs = parsedOuts.fridge_state_covs;
fridge_trans = parsedOuts.fridge_trans;
micro_state_means = parsedOuts.micro_state_means;
micro_state_covs = parsedOuts.micro_state_covs;
micro_trans = parsedOuts.micro_trans;
wd_state_means = parsedOuts.wd_state_means;
wd_state_covs = parsedOuts.wd_state_covs;
wd_trans = parsedOuts.wd_trans;
dw_state_means = parsedOuts.dw_state_means;
dw_state_covs = parsedOuts.dw_state_covs;
dw_trans = parsedOuts.dw_trans;
ac_state_means = parsedOuts.ac_state_means;
ac_state_covs = parsedOuts.ac_state_covs;
ac_trans = parsedOuts.ac_trans;
fridgeOn = parsedOuts.fridgeOn;
mwOn = parsedOuts.mwOn;
wdOn = parsedOuts.wdOn;
dwOn = parsedOuts.dwOn;
acOn = parsedOuts.acOn;
starting_point = parsedOuts.starting_point;
window_length = parsedOuts.window_length;
training_length = parsedOuts.training_length;
lik_thres = parsedOuts.lik_thres;
fridge = parsedOuts.fridge;
mw = parsedOuts.mw;
dw = parsedOuts.dw;
wd = parsedOuts.wd;
ac = parsedOuts.ac;

%%
transition_matrix = {fridge_trans, micro_trans, wd_trans, dw_trans, dw_trans, ac_trans};
state_means = {fridge_state_means, micro_state_means, wd_state_means, dw_state_means, dw_trans, ac_state_means};
state_covs = {fridge_state_covs, micro_state_covs, wd_state_covs, dw_state_covs, dw_trans, ac_state_covs};
always_on = {0, 3, 0, 0, 0, 4};

training_type = 2; %1 = no training, 2 = agg training, 3 = submetered training

residue = eDS.data(:,1);

inferred_power = zeros(numel(residue),1);

deviceOn = [fridgeOn;mwOn;wdOn;dwOn;acOn];

for dInc = 1:numel(deviceOn)
    if deviceOn(dInc)
        appliance = 
    end
end
appliance_idx = appliance + 1;
