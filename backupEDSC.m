%% Nicholas Czarnek
% SSPACISS Laboratory, Duke University
% Search this file for FINDME to see places that need to be updated.
%
% energyDataSet.m
% This class will be used for analysis of different datasets for energy
% disaggregation.
% Methods:
%   findErrorRMS          - find the RMS error between assigned and actual
%                           power of each component
%   findErrorAbsolute     - find the percent error between the assigned and
%                           actual energy assigned to each component
%   findEnergy            - 
%   plotDensity           - plot ksdensities for each component
%   visualization tools   - 

classdef energyDataSetClass < prtDataSetClass
  properties
    %% Hart category indicates type from Hart et al, 92.
    % 0 - not categorized yet
    % 1 - on/off
    % 2 - FSM
    % 3 - continuously variable
    % 4 - continuously on
    hartCategory;
    
    %% Events is a cell array indicating the events that were detected for each type of appliance.
    dataEventIdx;
    
  end
  %% Setup methods.
  methods
    %% Constructor for class instantiation.
    function obj = energyDataSetClass(varargin)
      obj = prtUtilAssignStringValuePairs(obj,varargin{:});
    end
    
    
    %% findErrorRmsChance
    % This method finds the RMS error that would occur if zero power were 
    % assigned across all time.
    function [chanceRMS,chanceRMS_struct] = findErrorRmsChance(obj,varargin)
      if isempty(varargin)
        % Output an error for every component
        chanceRMS = rms(obj.data,1)';
        
        featureNames = obj.getFeatureNames';
        
        nC = obj.nFeatures;
      else
        chanceRMS = rms(obj.data(:,varargin{1}(:)))';
        
        featureNames = obj.getFeatureNames(varargin{1}(:))';
        
        nC = max(size(varargin{1}));
      end
      
      rmsCells = cell(nC,1);
      for rInc = 1:nC
        featureNames{rInc}(ismember(featureNames{rInc},' ,.:;!')) = [];
        rmsCells{rInc} = chanceRMS(rInc);
      end
      
      chanceRMS_struct = cell2struct(rmsCells',featureNames',2);
    end
    
    
    
    %% findComponentEnergy
    % This method calculates the energy used by each device
    function [componentEnergy,energyStruct] = findComponentEnergy(obj,varargin)
      if isempty(varargin)
        % Use a trapezoidal approximation of the energy.
        componentEnergy = trapz(obj.data,1)';
        
        featureNames = obj.getFeatureNames';
        
        nC = obj.nFeatures;
      else
        componentEnergy = trapz(obj.data(:,varargin{1}(:)))';
        
        featureNames = obj.getFeatureNames(varargin{1}(:))';
        
        nC = max(size(varargin{1}));
      end
      
      energyCells = cell(nC,1);
      
      for rInc = 1:nC
        featureNames{rInc}(ismember(featureNames{rInc},' ,.:;!')) = [];
        energyCells{rInc} = componentEnergy(rInc);
      end
      
      energyStruct = cell2struct(energyCells',featureNames',2);
    end
    
    
    %% rankComponentsByEnergy
    % Same as findComponentEnergy, execpt with ranking.
    function [componentEnergy,energyStruct] = rankComponentsByEnergy(obj,varargin)
      if isempty(varargin)
        % Use a trapezoidal approximation of the energy.
        componentEnergy = trapz(obj.data,1)';
        
        featureNames = obj.getFeatureNames';
        
        nC = obj.nFeatures;
      else
        componentEnergy = trapz(obj.data(:,varargin{1}(:)))';
        
        featureNames = obj.getFeatureNames(varargin{1}(:))';
        
        nC = max(size(varargin{1}));
      end
      
      [~,eIdx] = sort(componentEnergy,'descend');
      
      featureNames = featureNames(eIdx);
      componentEnergy = componentEnergy(eIdx);
      
      energyCells = cell(nC,1);
      
      for rInc = 1:nC
        featureNames{rInc}(ismember(featureNames{rInc},' ,.:;!')) = [];
        energyCells{rInc} = componentEnergy(rInc);
      end
      
      energyStruct = cell2struct(energyCells',featureNames',2);
    end
    
    
    %% downsampleData - downsamples the data based on the input downsample factor
    function obj = downsampleData(obj,dsFactor)
      dIdx = 1:obj.nObservations;
      
      dIdx = downsample(dIdx,dsFactor);
      
      obj = obj.retainObservations(dIdx);
    end
    
    
    %% detectBinaryEvents - keeps events within the object itself
    function obj = detectEvents(obj,detectorType,detectorDetails,varargin)
      %% FINDME
      % This should be replaced with an improved detector.
      % If nothing is sent in, run event detection on everything.
      % Otherwise, run detection on the input columns.
      if isempty(varargin)
        focusColumns = 1:obj.nFeatures;
      else
        focusColumns = varargin{1};
      end
      
      obj.dataEventIdx = cell(obj.nFeatures,1);
      
      %% Threshold based detector.
      if strcmp(detectorType,'threshold')
        threshold = detectorDetails.threshold;
        
        %% If the threshold is singular, use the same everywhere.
        % Otherwise, check to make sure that the number of thresholds is
        % the same as the number of focusColumns
        if max(size(threshold)) ~= 1
          singleValueCheck = 0;
          if max(size(threshold)) ~= max(size(focusColumns))
            error('Either input one threshold or a threshold for each device')
          end
        else
          singleValueCheck = 1;
        end
        
        threshInc = 1;
        
        if singleValueCheck
          for featureInc = 1:obj.nFeatures
            if any(focusColumns == featureInc)
              aboveThreshold = obj.data(:,featureInc)>threshold;
              
              threshDiff = [0;diff(aboveThreshold)];
              
              obj.dataEventIdx{featureInc} = find(abs(threshDiff));
            end
            
          end
          
        else
          for featureInc = 1:obj.nFeatures
            if any(focusColumns == featureInc)
              aboveThreshold = obj.data(:,featureInc)>threshold(threshInc);
              
              % Increment which threshold will be used next.
              threshInc = threshInc + 1;
              
              threshDiff = [0;diff(aboveThreshold)];
              
              obj.dataEventIdx{featureInc} = find(abs(threshDiff));
            end
            
          end
        end
        
      %% Kyle's detector
      elseif strcmp(detectorType,'bradbury')
        
        %% Set up the event detection module
        ds.windowLength    = 15;
        ds.bufferLength    = 2;
        ds.threshold       = 1;
        ds.smoothFactor    = 0.75;
        ds.timeStamp      = [obj.observationInfo.times];
        
        for featureInc = 1:max(size(focusColumns))
          ds.data = obj.data(:,focusColumns(featureInc));
          
          fEvents = detectEvents(ds);
          
          obj.dataEventIdx{featureInc} = sort(cat(1,fEvents.onEventsIndex,fEvents.offEventsIndex));
          
        end
      end
      
      
      
      
    end
    
    %% 
    function eventFeatures = extractEventFeatures(obj,featureType,featureDetails)
      % This function outputs a prtDataSet based on the input data and
      % details stored in the structure featureDetails.
      
      if isempty(obj.dataEventIdx)
        error('Please first run event detection or input times at which to extract events as a cell array with the same number of cells as there are features');
      end
      
      eventFeatures = prtDataSetClass;
      
      if strcmp(featureType,'raw')
        %% This means that the desired features are simply the data 
        % themselves surrounding the events. Go through each event of each
        % device.
        if mod(featureDetails.windowSize,2) == 0
          error('Please enter an odd windowSize for an even number of features on both sides of the event')
        end
        
        windowSize = featureDetails.windowSize;
        
        halfWin = (windowSize - 1)/2;
        
        for deviceInc = 1:obj.nFeatures
          
          if ~isempty(obj.dataEventIdx{deviceInc})
            
            currentIdx = obj.dataEventIdx{deviceInc};
            
            currentIdx = currentIdx(currentIdx>halfWin);
            
            currentIdx = currentIdx(currentIdx<obj.nObservations - halfWin);
            
            numEvents = max(size(currentIdx));
            
            featureArray = zeros(numEvents,windowSize);
            
            %% Go through each time and extract features
            for eventInc = 1:max(size(currentIdx))
              startIdx = max(currentIdx(eventInc) - halfWin,1);
              endIdx = min(currentIdx(eventInc) + halfWin,obj.nObservations);
              
              % zeroMin is set to 1 if the features extracted are assumed
              % to be from windows of data floored to a minimum of 0.
              if isfield(featureDetails,'zeroMin')
                if featureDetails.zeroMin
                  featureArray(eventInc,:) = obj.data(startIdx:endIdx,deviceInc) - ...
                    min(obj.data(startIdx:endIdx,deviceInc));
                else
                  featureArray(eventInc,:) = obj.data(startIdx:endIdx,deviceInc);
                end
              else
                featureArray(eventInc,:) = obj.data(startIdx:endIdx,deviceInc);
              end
            end
            
            
            %% FINDME
            % Establish targets for the given device.
            % This needs to be updated for consistency.
            featureLabel = obj.userData.goodColumns(deviceInc);
            featureTargets = featureLabel*ones(numEvents,1);
            
            featureDS = prtDataSetClass(featureArray,featureTargets);
            
            eventTimes = struct('times',num2cell([obj.observationInfo(currentIdx).times]));
            
            featureDS.observationInfo = eventTimes;
            
            %% Concatenate onto eventFeatures.
            eventFeatures = catObservations(eventFeatures,featureDS);
            
          end
        end
        
      end
      
    end
    
    %% trainClassifier
    % ds - prtDataSet that contains the features and labels for the set to
    % be 
    function trainedClassifier = trainClassifier(obj,ds,classifier)
      trainedClassifier = classifier.train(ds);
    end
    
    function classOuts = runClassifier(obj,ds,trainedClass)
      if ~trainedClass.isTrained
        error('Please train your classifier')
      end
      
      cOuts = trainedClass.run(ds);
      
      %% Assign the proper labels to the classOuts
      [~,maxLoc] = max(cOuts.data,[],2);
      
      classOuts = cOuts;
      
      classOuts.data = zeros(classOuts.nObservations,1);
      
      uniqueClasses = classOuts.uniqueClasses;
      
      for mInc = 1:cOuts.nClasses
        classOuts.data(maxLoc == mInc) = uniqueClasses(mInc);
      end
      
      
    end
    
    
    
    
    
    
    
    
    function kOuts = crossVal(obj,ds,classifier,numFolds)
      cOuts = classifier.kfolds(ds,numFolds);
      
      %% Assign the proper labels to the classOuts
      [~,maxLoc] = max(cOuts.data,[],2);
      
      kOuts = cOuts;
      
      kOuts.data = zeros(kOuts.nObservations,1);
      
      uniqueClasses = kOuts.uniqueClasses;
      
      for mInc = 1:cOuts.nClasses
        kOuts.data(maxLoc == mInc) = uniqueClasses(mInc);
      end
      
    end
    
    
    
    
    
    
    
    %% Check the performance of supervised classification.
    function scoreConfusionMatrix(obj,kOuts)
      figure;
      prtScoreConfusionMatrix(kOuts.data,kOuts.targets)
    end
    
    
    
    
    
    
    
    
    
    
    %% Find the means of clusters based on a GMM.
    function [dataGmms,dataLls] = gmmClusterData(obj,nClusters)
      
      dataLls = zeros(1,obj.nFeatures);
      
      dataGmms = [];
      
      if numel(nClusters) == 1
        %% Assume the same number of clusters for all devices.
        dataLls = zeros(1,obj.nFeatures);
        
        %% Use the same number of clusters for each of the devices.
        for cInc = 1:obj.nFeatures
          %% Add a bias to ensure that none of the components have 0 covariance
          covBias = 1e-5;
          
          fComponents = repmat(prtRvMvn('covarianceBias',covBias),nClusters,1);
          
          clusterAlgo = prtRvGmm('components',fComponents);
          
          deviceSet = prtDataSetClass(obj.data(:,cInc));
          
          clusterAlgo = clusterAlgo.train(deviceSet);
          
          dataGmms = cat(1,dataGmms,clusterAlgo);
          
          dataLls(cInc) = sum(clusterAlgo.logPdf(deviceSet));
        end
        
        
        
        
      else
        %% If the same number is not going to be used, ensure that the vector 
        % has the same length as the number of features.
        if numel(nClusters)~=obj.nFeatures
          error('Either send in one value for the number of clusters or the same number of clusters as there are features.')
        end
        
        %% Assume the same number of clusters for all devices.
        dataLls = zeros(1,obj.nFeatures);
        
        %% Use the same number of clusters for each of the devices.
        for cInc = 1:obj.nFeatures
          currentNClusters = nClusters(cInc);
          
          if currentNClusters == 0
            clusterAlgo = prtRvGmm;
            
            dataGmms = cat(1,dataGmms,clusterAlgo);
            
            dataLls(cInc) = -1e6;
          else
            
            %% Add a bias to ensure that none of the components have 0 covariance
            covBias = 1e-5;
            
            fComponents = repmat(prtRvMvn('covarianceBias',covBias),currentNClusters,1);
            
            clusterAlgo = prtRvGmm('components',fComponents);
            
            deviceSet = prtDataSetClass(obj.data(:,cInc));
            
            clusterAlgo = clusterAlgo.train(deviceSet);
            
            dataGmms = cat(1,dataGmms,clusterAlgo);
            
            dataLls(cInc) = sum(clusterAlgo.logPdf(deviceSet));
          end
        end
        

      end
    end
    
    
    
    
    
    
    %% How much energy was assigned based on the classified features and associated timestamps?
    function assignedEnergy = assignEnergy(obj,assigned,dataGmms)
      keyboard
      
      assignedEnergy = obj;
      
      assignedEnergy.data = zeros(size(obj.data));
      
      uniqueAssigned = unique(assigned.data);
      
      %% Go through each of the appliances in the assigned object.
      for oInc = 1:max(size(uniqueAssigned))
        % Sort the timestamps for the given class.
        
      end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% Score ROCs for the performance of the event detection.
    
  end
  
  
  
  
  
end
