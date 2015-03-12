%% Nicholas Czarnek
% SSPACISS Laboratory, Duke University
% Search this file for FINDME to see places that need to be updated.
%
% energyDataSet.m
% This class will be used for analysis of different datasets for energy
% disaggregation.  There are several methods that make life easier.
% Methods:
%   findErrorRMSChance      - find the RMS error between assigned and actual
%                           power of each component
%   findComponentEnergy     - find the percent error between the assigned and
%                           actual energy assigned to each component
%   rankComponentsByEnergy  - plot ksdensities for each component
%   downsampleData          - downsample the data by the given factor
%   gmmClusterData          - cluster the data with a gmm
%   markEvents              - given events for each device within a data
%                           set, mark them on top of the data

classdef energyDataSetClass < prtDataSetClass
    
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
            %% Parse everything.
            options.devices = 1:obj.nFeatures;
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            devices = parsedOut.devices;
            
            xT = obj.getTimesFromUTC('timeScale','s');
            
            componentEnergy = trapz(xT,obj.data(:,devices));
            
            nD = numel(devices);
            
            energyCells = cell(nD,1);
            
            featureNames = obj.getFeatureNames(devices);
            
            for rInc = 1:nD
                featureNames{rInc}(ismember(featureNames{rInc},' ,.:;!')) = [];
                energyCells{rInc} = componentEnergy(rInc);
            end
            
            energyStruct = cell2struct(energyCells',featureNames',2);
            
%             
%             if isempty(varargin)
%                 % Use a trapezoidal approximation of the energy.
%                 componentEnergy = trapz(obj.data,1)';
%                 
%                 featureNames = obj.getFeatureNames';
%                 
%                 nC = obj.nFeatures;
%             else
%                 componentEnergy = trapz(obj.data(:,varargin{1}(:)))';
%                 
%                 featureNames = obj.getFeatureNames(varargin{1}(:))';
%                 
%                 nC = max(size(varargin{1}));
%             end
%             
%             energyCells = cell(nC,1);
%             
%             for rInc = 1:nC
%                 featureNames{rInc}(ismember(featureNames{rInc},' ,.:;!')) = [];
%                 energyCells{rInc} = componentEnergy(rInc);
%             end
%             
%             energyStruct = cell2struct(energyCells',featureNames',2);
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
        
        
        %% Find the means of clusters based on a GMM.
        function [dataGmms,dataLls] = gmmClusterData(obj,nClusters)
            
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
        
        %% Visualize the events on top of the aggregate and submetered signals
        % eventTimes is an energyEventClass object
        % Send in the eventTimes object from the device that you want to
        % see, and specify the device with 'focusDevice'.
        function markEventTimes(obj,eventTimes,varargin)
            % Assume that the user wants the events to be marked in
            % aggregate.
            options.focusDevice = 1;%2:obj.nFeatures;
            options.timeScale = 'hrs';
            options.XL = [];
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            focusDevice = parsedOut.focusDevice;
            timeScale = parsedOut.timeScale;
            XL = parsedOut.XL;
            
            xT = obj.getTimesFromUTC('zeroTimes',false,'timeScale','days');
            yT = obj.getTimesFromUTC('timeScale',timeScale);

            
            for dInc = 1:numel(focusDevice)
                focusDevice = focusDevice(dInc);
                
                %% Aggregate plot
                figure;
                plot(yT,obj.data(:,focusDevice))
                hold on
                
                legendStr = {[obj.getFeatureNames{focusDevice},' power']};
                
                pOnColors = ['g' 'b' 'r' 'c' 'm' 'y' 'k'];
                pOffColors = ['r' 'g' 'b' 'c' 'm' 'y' 'k'];
                
                for eInc = 1:numel(eventTimes)
                    [~,onIdx] = intersect(xT,eventTimes(eInc).onEventsTimes);
                    [~,offIdx] = intersect(xT,eventTimes(eInc).offEventsTimes);
                    
                    plot(yT(onIdx),obj.data(onIdx,focusDevice),[pOnColors(eInc),'o'])
                    plot(yT(offIdx),obj.data(offIdx,focusDevice),[pOffColors(eInc),'x'])
                    
                    if ~strcmp(eventTimes(eInc).className,'')
                        legendStr = cat(2,legendStr,[eventTimes(eInc).className{1},' on events'],....
                            [eventTimes(eInc).className{1},' off events']);
                    else
                        legendStr = cat(2,legendStr,['Device ',num2str(eInc),' on events'],....
                            ['Device ',num2str(eInc),' off events']);
                    end
                end
                title([obj.getFeatureNames{focusDevice},' with marked events'],'Interpreter','None')
                l = legend(legendStr);%[obj.getFeatureNames{focusDevice},' power'],'On events','Off events');
                set(l,'Interpreter','None')
                hold off
                
                xlabel(['Time (',timeScale,')'])
                ylabel('Power (W)')
                
                if ~isempty(XL)
                    xlim([XL(1) XL(2)])
                end
            end
            
            
            
        end
        
        
        function xT = getTimesFromUTC(obj,varargin)
            %% Note that this function assumes that times are formatted as UTC.
            % By default, the time scale is set to minutes, and it is assumed
            % that the first measurement occurs at time 0.
            options.timeScale = 'min';
            options.zeroTimes = true;
            
            pOut = prtUtilSimpleInputParser(options,varargin);
            
            timeScale = pOut.timeScale;
            zeroTimes = pOut.zeroTimes;
            
            xT = [obj.observationInfo.times]';
            
            if zeroTimes
                xT = xT - min(xT);
            end
            
            switch timeScale
                case 'hrs'
                    xT = xT*24;
                case 'min'
                    xT = xT*1440;
                case 's'
                    xT = xT*86400;
                case 'days'
                    % do nothing.  already stored as fractions of days
                    
            end
            
        end
        
        
        function plotDevices(obj,varargin)
            options.device = 1;
            options.combinedPlot = false;
            options.timeScale = 'hrs';
            options.zeroTimes = true;
            options.XL = [];
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            device = parsedOut.device;
            combinedPlot = parsedOut.combinedPlot;
            timeScale = parsedOut.timeScale;
            zeroTimes = parsedOut.zeroTimes;
            XL = parsedOut.XL;
            
            xT = obj.getTimesFromUTC('timeScale',timeScale,'zeroTimes',zeroTimes);
            
            if combinedPlot
                figure;
                plot(xT,obj.data(:,device))
                legend(obj.getFeatureNames(device))
                xlabel(['Time (',timeScale,')'])
                ylabel('Power (W)')
                title('Combined plot of input devices')
                if ~isempty(XL)
                    xlim([XL(1) XL(2)])
                end
            else
                for dInc = 1:numel(device)
                    figure;
                    plot(xT,obj.data(:,device(dInc)))
                    legend(obj.getFeatureNames(device(dInc)))
                    xlabel(['Time (',timeScale,')'])
                    ylabel('Power (W)')
                    title(obj.getFeatureNames(device(dInc)))
                    if ~isempty(XL)
                        xlim([XL(1) XL(2)])
                    end
                    
                end
            end
            
            
            
        end
        
        %% Get an eventTimes object based on the events in the userData of the dataset
        function bluedEvents = getBluedEvents(obj,varargin)
            options.phase = 'both';
            options.device = [];
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            bluedPhase = lower(parsedOut.phase);
            device = parsedOut.device;
            
            yT = obj.getTimesFromUTC('timeScale','days','zeroTimes',false);
            
            obj.userData.eventIdx = getEventIdx(yT,obj.userData.eventTimes);
            
            if strcmp(bluedPhase,'a')
                bluedPhase = 1;
            elseif strcmp(bluedPhase,'b')
                bluedPhase = 2;
            end
            
            bluedEvents = energyEventClass;
            bluedEvents.className = 'combinedEvents';
            bluedEvents.house = obj.name;
            bluedEvents.houseNumber = 1;
            
            if isscalar(device)
                bluedEvents.classNumber = device;
                bluedEvents.className = obj.userData.deviceNames(device);
            else
                bluedEvents.classNumber = 1;
            end
            
            uD = obj.userData;
            
            nEvents = numel(uD.eventIdx);
            
            
            if isempty(device)
                includeIndex = true(nEvents,1);
            else
                includeIndex = false(nEvents,1);
                for eInc = 1:numel(uD.eventIdx)
                    includeIndex(eInc) = device == uD.eventTypes(eInc);
                end
            end
            
            
            
            switch bluedPhase
                case 'both'
                    deviceOn = uD.onEvents(includeIndex);
                    deviceIdx = uD.eventIdx(includeIndex);
                    
                    onEventsIdx = deviceIdx(deviceOn);
                    
                    offEventsIdx = deviceIdx(~deviceOn);
                    
                case 1 % Just phase A
                    keepIdx = includeIndex' & strcmp(uD.phase(includeIndex),'A');
                    
                    deviceOn = uD.onEvents(keepIdx);
                    deviceIdx = uD.eventIdx(keepIdx);
                    
                    onEventsIdx = deviceIdx(deviceOn);
                    offEventsIdx = deviceIdx(~deviceOn);
                    
                    if any(strcmp(uD.phase(keepIdx),'B'))
                        warning('You included devices from phase B that were ignored here');
                    end
                    
                case 2
                    keepIdx = includeIndex' & strcmp(uD.phase(includeIndex),'B');
                    
                    deviceOn = uD.onEvents(keepIdx);
                    deviceIdx = uD.eventIdx(keepIdx);
                    
                    onEventsIdx = deviceIdx(deviceOn);
                    offEventsIdx = deviceIdx(~deviceOn);
                    
                    
                    if any(strcmp(uD.phase(keepIdx),'A'))
                        warning('You included devices from phase A that were ignored here');
                    end
                    
                otherwise
                    %% Assume both are desired
                    % Get events from both phases.
                    deviceOn = uD.onEvents(includeIndex);
                    deviceIdx = uD.eventIdx(includeIndex);
                    
                    onEventsIdx = deviceIdx(deviceOn);
                    
                    offEventsIdx = deviceIdx(~deviceOn);
                    
            end
            
            bluedEvents.onEventsIndex = onEventsIdx;
            bluedEvents.offEventsIndex = offEventsIdx;
            
            bluedEvents.onEventsTimes = yT(onEventsIdx);
            bluedEvents.offEventsTimes = yT(offEventsIdx);
                    
        end
        
        
        %% Extract features based on the input times
        % eventTimes must be an eventTimes object to make life easier
        function energyFeats = extractEventData(obj,eventTimes,varargin)
            %% Deal with varargin
            options.windowInS = 61;% take data from 1 minute on each side
%             options.secondIncrement = 1;% adjust the resolution of the
%             features, e.g. take data at 1:1:60 or 1:8:60

            % default to extract features from the aggregate load
            options.devices = 1;
            % Allow on, off, or both to be extracted
            options.featureType = 'both';
            % allow the class to be specified
            options.classNumber = 1;% as
            % allow the class name to be provided
            options.className = {''};
            % subtract out the minimum from each feature such that
            % everything starts at the baseline
            options.zeroFeatures = true;
            
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            
            windowInS = parsedOut.windowInS;
            devices = parsedOut.devices;
            classNumber = parsedOut.classNumber;
            featureType = parsedOut.featureType;
            className = parsedOut.className;
            zeroFeatures = parsedOut.zeroFeatures;
            
            %% Only keep the with proper measurements.
            if isfield(obj.observationInfo,'keepLogicals')
                kL = [obj.observationInfo.keepLogicals]';
            else
                kL = ones(obj.nObservations,1);
            end
            
            xT = obj.getTimesFromUTC('timeScale','days');
            yT = obj.getTimesFromUTC('timeScale','days','zeroTimes',false);
            
            tRes = xT(2);
            
            obsPerDay = round(1/tRes);
            
            %% How many features should there be per observation?
            obsPerHalf = round(windowInS/(tRes*86400));
            
            %% Round everything to the nearest unit of time.
            yT = round(yT*obsPerDay)/obsPerDay;
            
            eventTimes.onEventsTimes = round(eventTimes.onEventsTimes*obsPerDay)/obsPerDay;
            
            eventTimes.offEventsTimes = round(eventTimes.offEventsTimes*obsPerDay)/obsPerDay;
            
            energyFeats = repmat(energyFeatureClass,numel(devices),1);
            
            %% Extract features from each event
            for dInc = 1:numel(devices)
%                 currentDevice = devices(dInc);
%                 
                %% Go through the different devices and create a new energyFeatureClass for each
                currentFeats = energyFeatureClass;
                
                %% Set up the observation info.
                if strcmp(featureType,'both')
                    timestamp = cat(1,eventTimes.onEventsTimes,...
                        eventTimes.offEventsTimes);
                    
                    eventType = cat(1,ones(numel(eventTimes.onEventsTimes),1),...
                        zeros(numel(eventTimes.offEventsTimes),1));
                    
                    removeIdx = [];
                    
                    featureMat = [];
                    
                    for eventInc = 1:numel(eventTimes.onEventsTimes)
                        %% Find the event index
                        [~,eventIdx] = intersect(yT,eventTimes.onEventsTimes(eventInc));
                        
                        %% Make the data vector.
                        featureIdx = eventIdx - obsPerHalf:eventIdx + obsPerHalf;
                        
                        if ~any(featureIdx<1)&&~any(featureIdx>obj.nObservations) && all(kL(featureIdx) == 1)
                            featureVector = obj.data(featureIdx,dInc)';
                            
                            if zeroFeatures
                                featureVector = featureVector - min(featureVector);
                            end
                            
                            featureMat = cat(1,featureMat,featureVector);
                        else
                            removeIdx = cat(1,removeIdx,eventInc);
                        end
                    end
                    
                    nOn = numel(eventTimes.onEventsTimes);
                    
                    for eventInc = 1:numel(eventTimes.offEventsTimes)
                        %% Find the event index
                        [~,eventIdx] = intersect(yT,eventTimes.offEventsTimes(eventInc));
                        
                        %% Make the data vector.
                        featureIdx = eventIdx - obsPerHalf:eventIdx + obsPerHalf;
                        
                        if ~any(featureIdx<1)||~any(featureIdx>obj.nObservations)
                            featureVector = obj.data(featureIdx,dInc)';
                            
                            if zeroFeatures
                                featureVector = featureVector - min(featureVector);
                            end
                            
                            featureMat = cat(1,featureMat,featureVector);
                        else
                            removeIdx = cat(1,removeIdx,eventInc+nOn);
                        end
                    end
                    
                elseif strcmp(featureType,'on')
                    timestamp = eventTimes.onEventsTimes;
                    
                    eventType = ones(numel(eventTimes.onEventsTimes),1);
                    
                    removeIdx = [];
                    
                    featureMat = [];
                    
                    for eventInc = 1:numel(eventTimes.onEventsTimes)
                        %% Find the event index
                        [~,eventIdx] = intersect(yT,eventTimes.onEventsTimes(eventInc));
                        
                        %% Make the data vector.
                        featureIdx = eventIdx - obsPerHalf:eventIdx + obsPerHalf;
                        
                        if ~any(featureIdx<1)&&~any(featureIdx>obj.nObservations) && all(kL(featureIdx) == 1)
                            featureVector = obj.data(featureIdx,dInc)';
                            
                            if zeroFeatures
                                featureVector = featureVector - min(featureVector);
                            end
                            
                            featureMat = cat(1,featureMat,featureVector);
                        else
                            removeIdx = cat(1,removeIdx,eventInc);
                        end
                    end
                    
                    
                    
                elseif strcmp(featureType,'off')
                    timestamp = eventTimes.offEventsTimes;
                    
                    eventType = ones(numel(eventTimes.offEventsTimes),1);
                    
                    removeIdx = [];
                    
                    featureMat = [];
                    
                    for eventInc = 1:numel(eventTimes.offEventsTimes)
                        %% Find the event index
                        [~,eventIdx] = intersect(yT,eventTimes.offEventsTimes(eventInc));
                        
                        %% Make the data vector.
                        featureIdx = eventIdx - obsPerHalf:eventIdx + obsPerHalf;
                        
                        if ~any(featureIdx<1)&&~any(featureIdx>obj.nObservations) && all(kL(featureIdx) == 1)
                            featureVector = obj.data(featureIdx,dInc)';
                            
                            if zeroFeatures
                                featureVector = featureVector - min(featureVector);
                            end
                            
                            featureMat = cat(1,featureMat,featureVector);
                        else
                            removeIdx = cat(1,removeIdx,eventInc);
                        end
                    end
                    
                else
                    error('Please choose on, off, or both for the eventType\n')
                end
                
                currentFeats.data = featureMat;
                timestamp(removeIdx) = [];
                eventType(removeIdx) = [];
                
                timestamp = num2cell(timestamp);
                eventType = num2cell(eventType);
                
                obsInfo = struct('timestamp',timestamp,'eventType',eventType);
                
                currentFeats.observationInfo = obsInfo;
                
                currentFeats.targets = classNumber*ones(currentFeats.nObservations,1);
                
                currentFeats.classNames = className;
                
                energyFeats(dInc) = currentFeats;
                
            end
            
        end
        
        
        %% Determine the average power consumed based on events
        function powerOut = findAveragePower(obj,eventTimes,varargin)
            %% Handle varargin
            options.features = [obj.featureInfo(2:end).pecanClass]';
            options.fromAggregate = true;
            options.zeroPower = true;
            options.removeTails = true;
            options.tailPercent = 5;
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            features = parsedOut.features;
            fromAggregate = parsedOut.fromAggregate;
            zeroPower = parsedOut.zeroPower;
            removeTails = parsedOut.removeTails;
            tailPercent = parsedOut.tailPercent;
            
            pecanClasses = [obj.featureInfo.pecanClass]';
            
            %% Get the time in seconds.
            xT = obj.getTimesFromUTC('timeScale','s');
            yT = obj.getTimesFromUTC('timeScale','days','zeroTimes',false);
            
            %% Assume that we're extracting from the aggregate load
            extractionFeature = pecanClasses(1);
            
            powerOut = repmat(struct('averagePower',[],'averageTimeInS',[],...
                'classNumber',[],'fromAggregate',[],'numberInstances',[]),numel(features),1);
            
            %% Go through the events.
            for dInc = 1:numel(features)
                currentFeature = features(dInc);
                
                fIdx = pecanClasses == currentFeature;
                
                eIdx = false(numel(eventTimes),1);
                
                %% Find the eventTimes structure that corresponds to the current device.
                for eInc = 1:numel(eventTimes)
                    if (eventTimes(eInc).classNumber == currentFeature)
                        eIdx(eInc) = true;
                    end
                end
                
                eventIndex = find(eIdx);
                
                if ~isempty(eventIndex)
                    
                    if numel(eventIndex)>1
                        error('You have multiple devices with the same class')
                    end
                    
                    if ~fromAggregate
                        %% Find the matching datastream in the energyDataSet
                        extractionFeature = find(pecanClasses == currentFeature);
                    end
                    
                    energyUsed = zeros(numel(eventTimes(eventIndex).onEventsTimes),1);
                    onTime = zeros(numel(eventTimes(eventIndex).onEventsTimes),1);
                    
                    %% Go through all of the on events, find the energy used
                    % and the time taken.
                    for eInc = 1:numel(eventTimes(eventIndex).onEventsTimes)
                        %% Find the next off time that is larger than the current time.
                        currentOn = eventTimes(eventIndex).onEventsTimes(eInc);
                        
                        offIdx = find(eventTimes(eventIndex).offEventsTimes>currentOn,1,'first');
                        
                        if ~isempty(offIdx)
                            currentOff = eventTimes(eventIndex).offEventsTimes(offIdx);
                            
                            focusIndices = yT>=currentOn & yT<=currentOff;
                            
                            startInS = xT(find(focusIndices,1,'first'));
                            endInS = xT(find(focusIndices,1,'last'));
                            
                            powerSignature = obj.data(focusIndices,extractionFeature);
                            
                            if zeroPower
                                powerSignature = powerSignature - min(powerSignature);
                            end
                            
                            energyUsed(eInc) = trapz(xT(focusIndices),powerSignature);
                            
                            onTime(eInc)= endInS - startInS + 1;
                        end
                        
                    end
                    
                    %% Remove the times that were not covered.
                    zeroTimes = find(energyUsed == 0);
                    
                    energyUsed(zeroTimes) = [];
                    onTime(zeroTimes) = [];
                    
                    numInstances = numel(energyUsed);
                    
                    if numInstances ~= 0
                        newIdx = 1:numel(energyUsed);
                        
                        %% Handle outliers
                        if removeTails
                            [energyUsed,newIdx] = sort(energyUsed);
                            
                            
                            startKeep = round(tailPercent/100*numInstances);
                            endKeep = round((100-tailPercent)/100*numInstances);
                            
                            focusIdx = max(startKeep,1):min(endKeep,numInstances);
                            
                            energyUsed = energyUsed(focusIdx);
                            newIdx = newIdx(focusIdx);
                            onTime = onTime(newIdx);
                        end
                        
                    end
                    
                    %% Find the average power for the current device.
                    powerOut(dInc).averagePower = sum(energyUsed)/sum(onTime);
                    powerOut(dInc).averageTimeInS = mean(onTime);
                    powerOut(dInc).classNumber = currentFeature;
                    powerOut(dInc).fromAggregate = fromAggregate;
                    powerOut(dInc).className = obj.getFeatureNames(fIdx);
                    powerOut(dInc).numberInstances = numel(energyUsed);
                    
                end
            end
        end
        
        %% Assign power to specific devices based on the current 
        function assignedPower = assignPower(obj,eventTimes,averagePower,varargin)
            % eventTimes -  energyEventClass object or array of objects
            % averagePower - output structure or structure array from findAveragePower
            
            %% Work with varargin.
            options.useOff = false;
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            useOff = parsedOut.useOff;
            
            %% Establish the new dataset, assignedPower
            assignedPower = obj;
            
            assignedPower.data = zeros(obj.nObservations,obj.nFeatures);
            
            pecanClasses = [obj.featureInfo.pecanClass]';
            
            %% Get the times
            xT = obj.getTimesFromUTC('timeScale','s');
            yT = obj.getTimesFromUTC('timeScale','days','zeroTimes',false);
            
            %% Go through the eventTimes array to figure out which devices 
            % to assign power to
            for dInc = 1:numel(eventTimes)
                %% Find the corresponding feature vector.
                fIdx = pecanClasses == eventTimes(dInc).classNumber;
                
                aPIdx = false(numel(averagePower),1);
                
                %% Find the averagePower structure corresponding to the current class
                for aInc = 1:numel(averagePower)
                    if averagePower(aInc).classNumber == eventTimes(dInc).classNumber
                        aPIdx(aInc) = true;
                    end
                end
                
                if ~useOff
                    %% Make assignments based on the eventTimes
                    for eInc = 1:numel(eventTimes(dInc).onEventsTimes)
                        %% Work with the proper indices
                        startIdx = find(yT>=eventTimes(dInc).onEventsTimes(eInc),1,'first');
                        
                        startXT = xT(startIdx);
                        endXT = startXT + averagePower(aPIdx).averageTimeInS;
                        
                        focusIdx = xT>=startXT & xT<=endXT;
                        
                        %% Assign power to the current indices.
                        assignedPower.data(focusIdx,fIdx) = averagePower(aPIdx).averagePower;
                        
                        
                    end
                else
                    %% Make assignments based both on the on and off times.
                    for eInc = 1:numel(eventTimes(dInc).onEventsTimes)
                        %% Work with the proper indices
                        % Find the next off.
                        nextOffIdx = find(eventTimes(dInc).offEventsTimes>=eventTimes(dInc).onEventsTimes(eInc),1,'first');
                        
                        if ~isempty(nextOffIdx)
                            focusIdx = yT>=eventTimes(dInc).onEventsTimes(eInc) &...
                                yT<=eventTimes(dInc).offEventsTimes(nextOffIdx);
                            
                            assignedPower.data(focusIdx,fIdx) = averagePower(aPIdx).averagePower;
                        end
                        
                        
                    end
                end
            end
            
        end
        
        
        
        %% Determine the error when comparing assigned to true power.
        function errorMetrics = calculateAssignmentErrors(obj,assignedPower,varargin)
            options.focusDevice = 1:obj.nFeatures;
            options.subset = false;
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            focusDevice = parsedOut.focusDevice;
            subset = parsedOut.subset;
            
            errorMetrics = prtDataSetClass('data',zeros(4,numel(focusDevice)));
            
            errorNames = {'RMS';'percentEnergyExplained';'meanNormError';'chanceRMS'};
            
            obsInfo = struct('errorType',errorNames);
            
            aggregateEnergyUsed = obj.findComponentEnergy('device',1);
                
            if ~subset
                errorMetrics = errorMetrics.setFeatureNames(obj.getFeatureNames);
                errorMetrics.observationInfo = obsInfo;
                
                for dInc = 1:numel(focusDevice)
                    %% RMS
                    errorMetrics.data(1,dInc) = sqrt(1/obj.nObservations * ...
                        sum((obj.data(:,dInc) - assignedPower.data(:,dInc)).^2));
                    
                    %% Percent energy explained.  Note that since this is a percent,
                    % we don't need to worry about the timescale.
                    assignedEnergy = assignedPower.findComponentEnergy('device',dInc);
                    
                    errorMetrics.data(2,dInc) = assignedEnergy/aggregateEnergyUsed;
                    
                    %% Mean normalized error
                    actualEnergyUsed = obj.findComponentEnergy('device',dInc);
                    
                    errorMetrics.data(3,dInc) = abs(actualEnergyUsed - ...
                        assignedEnergy)/actualEnergyUsed;
                    
                    %% Chance RMS. Assign average energy across all time.
                    actualPower = obj.data(:,dInc);
                    
                    meanPower = ones(obj.nObservations,1) * mean(actualPower);
                    
                    errorMetrics.data(4,dInc) = sqrt(1/obj.nObservations * ...
                        sum((obj.data(:,dInc) - meanPower).^2));
                    
                end
            else
                errorMetrics = errorMetrics.setFeatureNames(obj.getFeatureNames(focusDevice));
                errorMetrics.observationInfo = obsInfo;
                for dInc = 1:numel(focusDevice)
                    %% RMS
                    errorMetrics.data(1,dInc) = sqrt(1/obj.nObservations * ...
                        sum((obj.data(:,focusDevice(dInc)) - assignedPower.data(:,dInc)).^2));
                    
                    %% Percent energy explained.  Note that since this is a percent,
                    % we don't need to worry about the timescale.
                    assignedEnergy = assignedPower.findComponentEnergy('device',dInc);
                    
                    errorMetrics.data(2,dInc) = assignedEnergy/aggregateEnergyUsed;
                    
                    %% Mean normalized error
                    actualEnergyUsed = obj.findComponentEnergy('device',focusDevice(dInc));
                    
                    errorMetrics.data(3,dInc) = abs(actualEnergyUsed - ...
                        assignedEnergy)/actualEnergyUsed;
                    
                    %% Chance RMS. Assign average energy across all time.
                    actualPower = obj.data(:,focusDevice(dInc));
                    
                    meanPower = ones(obj.nObservations,1) * mean(actualPower);
                    
                    errorMetrics.data(4,dInc) = sqrt(1/obj.nObservations * ...
                        sum((obj.data(:,focusDevice(dInc)) - meanPower).^2));
                    
                end
            end
        end
        
        
        %% Run the entire supervised disaggregation system 
        function [assignedPowerPerformance] = runSupervisedSystem(obj,eventTimes,varargin)
            % How many cross val folds do you want for the full
            % disaggregation system?
            options.nXFolds = 5;
            % What detection window should be used for the detector?
            options.detectionHalfWinInS = 61;
            % When is an event considered to be detected?
            options.eventHaloInS = 61;
            % What is the threshold selection criterion based on the FAR?
            options.farThreshold = 1;
            % How big should the features be?
            options.extractionWindow = 61;
            % Should the principal components be used?
            options.usePca = true;
            % If so, how many?
            options.nPcaComponents = 10;
            % Should the true event times be used for the classification
            % step?
            options.useTrueTimes = false;
            % Do you only want to focus on detection?
            options.onlyRunDetection = false;
            % Do you want to zmuv the components?
            options.zmuvFeatures = true;
            % What classifier do you want?
            options.classifier = 'knn';
            % If KNN, what is k?
            options.k = 5;
            % If random forrest, how many trees?
            options.nTrees = 20;
            % If SVM, rbf or linear?
            options.useRbf = true;
            
            parsedOuts = prtUtilSimpleInputParser(options,varargin);
            
            
            nXFolds = parsedOuts.nXFolds;
            detectionHalfWinInS = parsedOuts.detectionHalfWinInS;
            eventHaloInS = parsedOuts.eventHaloInS;
            farThreshold = parsedOuts.farThreshold;
            extractionWindow = parsedOuts.extractionWindow;
            usePca = parsedOuts.usePca;
            nPcaComponents = parsedOuts.nPcaComponents;
            useTrueTimes = parsedOuts.useTrueTimes;
            onlyRunDetection = parsedOuts.onlyRunDetection;
            zmuvFeatures = parsedOuts.zmuvFeatures;
            classifier = parsedOuts.classifier;
            k = parsedOuts.k;
            nTrees = parsedOuts.nTrees;
            useRbf = parsedOuts.svmRbf;
            
            %% Split up the data evenly into the number of folds
            nObsPerFold = round(obj.nObservations/nXFolds);
            
            xT = obj.getTimesFromUTC('timeScale','days','zeroTimes',false);
            
            %% Set up the aggregate energyEventClass
            detectedTestingEvents = energyEventClass;
            assignedEventsTotal = repmat(energyEventClass,numel(eventTimes),1);
            

            %% Run through each fold separately.
            for fInc = 1:nXFolds
                testIdx = false(obj.nObservations,1);
                trainIdx = true(obj.nObservations,1);
                
                if fInc == nXFolds
                    startIdx = nObsPerFold*(fInc-1)+1;
                    endIdx = obj.nObservations;
                else
                    startIdx = nObsPerFold*(fInc-1)+1;
                    endIdx = nObsPerFold*fInc - 1;
                end
                
                testIdx(startIdx:endIdx) = true;
                
                trainIdx(startIdx:endIdx) = false;
                
                
                if fInc == 1
                    trainTimes = [xT(endIdx+1) xT(end)];
                elseif fInc == nXFolds
                    trainTimes = [xT(1) xT(startIdx - 1)];
                else
                    trainTimes = [xT(1) xT(startIdx - 1);xT(endIdx+1) xT(end)];
                end
                
                trainData = obj.retainObservations(trainIdx);
                testData = obj.retainObservations(testIdx);
                
                %% Combine the energy events from all classes into one for the aggregate
                eventTimes(1) = combineEnergyEvents(eventTimes,'devices',2:numel(eventTimes));
                
                
                %% Score the event detector based on the eventTimes
                trainingTimes = eventTimes;
                testingTimes = eventTimes;
                
                %% Figure out which on and off times to keep for training
                for tInc = 1:numel(eventTimes)
                    %% Training times.
                    % Find the off times that are within the training
                    % times.
                    keepOff = false(numel(eventTimes(tInc).offEventsTimes),1);
                    keepOn = false(numel(eventTimes(tInc).onEventsTimes),1);
                    for eInc = 1:size(trainTimes,1)
                        keepOff(eventTimes(tInc).offEventsTimes>=trainTimes(eInc,1)...
                            & eventTimes(tInc).offEventsTimes<=trainTimes(eInc,2)) = true;
                        keepOn(eventTimes(tInc).onEventsTimes>=trainTimes(eInc,1)...
                            & eventTimes(tInc).onEventsTimes<=trainTimes(eInc,2)) = true;
                        
                        
                    end
                    
                    
                    % onEvents and offEvents are not always established since they
                    % are never used
                    if ~isempty(trainingTimes(tInc).offEvents)
                        trainingTimes(tInc).offEvents = trainingTimes(tInc).offEvents(keepOff);
                    end
                    if ~isempty(trainingTimes(tInc).onEvents)
                        trainingTimes(tInc).onEvents = trainingTimes(tInc).onEvents(keepOn);
                    end
                    
                    trainingTimes(tInc).offEventsIndex = trainingTimes(tInc).offEventsIndex(keepOff);
                    trainingTimes(tInc).offEventsTimes = trainingTimes(tInc).offEventsTimes(keepOff);
                    
                    trainingTimes(tInc).onEventsIndex = trainingTimes(tInc).onEventsIndex(keepOn);
                    trainingTimes(tInc).onEventsTimes = trainingTimes(tInc).onEventsTimes(keepOn);
                    
                    %% Testing times
                    % onEvents and offEvents are not always established since they
                    % are never used
                    if ~isempty(testingTimes(tInc).offEvents)
                        testingTimes(tInc).offEvents = testingTimes(tInc).offEvents(~keepOff);
                    end
                    if ~isempty(testingTimes(tInc).onEvents)
                        testingTimes(tInc).onEvents = testingTimes(tInc).onEvents(~keepOn);
                    end
                    
                    testingTimes(tInc).onEventsIndex = testingTimes(tInc).onEventsIndex(~keepOn);
                    testingTimes(tInc).onEventsTimes = testingTimes(tInc).onEventsTimes(~keepOn);
                    
                    testingTimes(tInc).offEventsIndex = testingTimes(tInc).offEventsIndex(~keepOff);
                    testingTimes(tInc).offEventsTimes = testingTimes(tInc).offEventsTimes(~keepOff);
                end
                
                
                if ~useTrueTimes
                    %% Run the event detector on the training data
                    detectedTrainingEvents = detectEnergyEvents(trainData,...
                        'halfWindowInS',detectionHalfWinInS,'device',1);
                    
                    %% Score the performance.
                    trainingEventDetectionPerformance = scoreEventDetectionQuickly(...
                        detectedTrainingEvents,trainingTimes(1),eventHaloInS);
                    
                    onThreshIdx = find(trainingEventDetectionPerformance.onFa>=farThreshold,1,'first');
                    
                    onEventThreshold = trainingEventDetectionPerformance.onThresholds(onThreshIdx);
                end
                
                %% Extract features for the classifier.
                energyFeatures = energyFeatureClass;
                for tInc = 2:numel(trainingTimes)
                    if ~isempty(trainingTimes(tInc).onEventsTimes)
                        currentFeatures = obj.extractEventData(trainingTimes(tInc),'className',...
                            trainingTimes(tInc).className,'classNumber',...
                            trainingTimes(tInc).classNumber,'featureType','on',...
                            'windowInS',extractionWindow);
                        
                        energyFeatures = catObservations(energyFeatures,currentFeatures);
                    end
                end
                
                %% Train the classifier.
                switch lower(classifier)
                    case 'knn'
                        classifier = prtClassKnn('k',k) + prtDecisionMap;
                    case 'svm'
                        if useRbf
                            classifier = prtClassLibSvm + prtDecisionMap;
                        else
                            classifier = prtClassLibSvm('kernel',0) + prtDecisionMap;
                        end
                    case 'rf'
                        classifier = prtClassTreeBaggingCap('nTrees',nTrees) + prtDecisionMap;
                end
                
                %% Modify the features based on input options to take the 
                % principal components and the zscore
                if usePca
                    [inputFeats,inputPca] = energyFeatures.getPca('nComponents',nPcaComponents,...
                        'plotPca',false,'zmuv',zmuvFeatures);
                else
                    inputFeats = energyFeatures;
                    if zmuvFeatures
                        zM = prtPreProcZmuv;
                        zM = zM.train(energyFeatures);
                        inputFeats = zM.run(energyFeatures);
                    end
                end
                
                classifier = classifier.train(inputFeats);
                
                %% Detect events in the testing data.
                if ~useTrueTimes
                    testingDetectedEvents = detectEnergyEvents(testData,...
                        'halfWindowInS',detectionHalfWinInS,'device',1,'threshold',...
                        onEventThreshold);
                else
                    %% Get the confidences, but use the true times.
                    eventConfidences = detectEnergyEvents(testData,...
                        'halfWindowInS',detectionHalfWinInS,'device',1);
                    testingDetectedEvents = testingTimes(1);
                    testingDetectedEvents.confidences = eventConfidences.confidences;
                    testingDetectedEvents.timeStamps = eventConfidences.timeStamps;
                end
                
                %% Add on to the existing energyEventClass for later scoring
                detectedTestingEvents.onEventsIndex = cat(1,...
                    detectedTestingEvents.onEventsIndex,...
                    testingDetectedEvents.onEventsIndex);
                detectedTestingEvents.onEventsTimes = cat(1,...
                    detectedTestingEvents.onEventsTimes,...
                    testingDetectedEvents.onEventsTimes);
                detectedTestingEvents.offEventsIndex = cat(1,...
                    detectedTestingEvents.offEventsIndex,...
                    testingDetectedEvents.offEventsIndex);
                detectedTestingEvents.offEventsTimes = cat(1,...
                    detectedTestingEvents.offEventsTimes,...
                    testingDetectedEvents.offEventsTimes);
                
                %% Extract features from the current detected events.
                detectedEnergyFeatures = obj.extractEventData(testingDetectedEvents,...
                    'featureType','on','windowInS',extractionWindow);
                
                if usePca
                    detectedEnergyPca = inputPca.run(detectedEnergyFeatures);
                else
                    detectedEnergyPca = detectedEnergyFeatures;
                end
                
                %% Run the classifier.
                testingClassOuts = classifier.run(detectedEnergyPca);
                
                %% Make the assignments.
                % THIS NEEDS TO BE CHANGED IN THE FUTURE TO ACCOUNT FOR
                % DIFFERENT CLASSES!!!!  THIS CURRENTLY REQUIRES THE SAME
                % NUMBER OF EVENTTIMES AS FEATURES, AND IT WOULD BE NICE
                % FOR THIS TO BE AGNOSTIC!!!
                assignedEvents = repmat(energyEventClass,numel(eventTimes),1);
                for eventInc = 1:testingClassOuts.nObservations
                    currentClass = testingClassOuts.data(eventInc);
                    
                    assignedEvents(currentClass).onEventsTimes = ...
                        cat(1,assignedEvents(currentClass).onEventsTimes,...
                        testingClassOuts.observationInfo(eventInc).timestamp);
                end
                
                %% Add on to the aggregate events
                for typeInc = 1:numel(eventTimes)
                    assignedEvents(typeInc).classNumber = typeInc;
                    
                    if ~isempty(assignedEvents(typeInc).onEventsTimes)
                        assignedEventsTotal(typeInc).onEventsTimes = ...
                            cat(1,assignedEventsTotal(typeInc).onEventsTimes,....
                            assignedEvents(typeInc).onEventsTimes);
                    end
                end
                
                
            end
            
            %% THIS SHOULD BE CHANGED SO THAT THE AVERAGE ACCOUNTS FOR ALL 
            % FOLDS, NOT JUST THE LAST FOLD.
            averagePowers = trainData.findAveragePower(trainingTimes(2:end));
            
            
            %% Now, we have the aggregated event times and the average power
            % details, so we can find both the device ROCs and the energy 
            % assignment metrics.
            
            for typeInc = 1:numel(eventTimes)
                assignedEventsTotal(typeInc).classNumber = typeInc;
            end
            
            %% Assign power to each device
            newData = obj.assignPower(assignedEventsTotal,averagePowers);
            
            
            %% Evaluate the assignment performance.
            assignedPowerPerformance = obj.calculateAssignmentErrors(newData);
            
            %% Evaluate the detection performance of the aggregated events.
%             detectionPerformance = scoreEventDetectionQuickly(detectedTestingEvents,...
%                 eventTimes(1),eventHaloInS);
        end
    end
    
    
    
    
    
end














%
%
% %% detectBinaryEvents - keeps events within the object itself
% function obj = detectEvents(obj,detectorType,detectorDetails,varargin)
% %% FINDME
% % This should be replaced with an improved detector.
% % If nothing is sent in, run event detection on everything.
% % Otherwise, run detection on the input columns.
% if isempty(varargin)
%   focusColumns = 1:obj.nFeatures;
% else
%   focusColumns = varargin{1};
% end
%
% obj.dataEventIdx = cell(obj.nFeatures,1);
%
% %% Threshold based detector.
% if strcmp(detectorType,'threshold')
%   threshold = detectorDetails.threshold;
%
%   %% If the threshold is singular, use the same everywhere.
%   % Otherwise, check to make sure that the number of thresholds is
%   % the same as the number of focusColumns
%   if max(size(threshold)) ~= 1
%     singleValueCheck = 0;
%     if max(size(threshold)) ~= max(size(focusColumns))
%       error('Either input one threshold or a threshold for each device')
%     end
%   else
%     singleValueCheck = 1;
%   end
%
%   threshInc = 1;
%
%   if singleValueCheck
%     for featureInc = 1:obj.nFeatures
%       if any(focusColumns == featureInc)
%         aboveThreshold = obj.data(:,featureInc)>threshold;
%
%         threshDiff = [0;diff(aboveThreshold)];
%
%         obj.dataEventIdx{featureInc} = find(abs(threshDiff));
%       end
%
%     end
%
%   else
%     for featureInc = 1:obj.nFeatures
%       if any(focusColumns == featureInc)
%         aboveThreshold = obj.data(:,featureInc)>threshold(threshInc);
%
%         % Increment which threshold will be used next.
%         threshInc = threshInc + 1;
%
%         threshDiff = [0;diff(aboveThreshold)];
%
%         obj.dataEventIdx{featureInc} = find(abs(threshDiff));
%       end
%
%     end
%   end
%
%   %% Kyle's detector
% elseif strcmp(detectorType,'bradbury')
%
%   %% Set up the event detection module
%   ds.windowLength    = 15;
%   ds.bufferLength    = 2;
%   ds.threshold       = 1;
%   ds.smoothFactor    = 0.75;
%   ds.timeStamp      = [obj.observationInfo.times];
%
%   for featureInc = 1:max(size(focusColumns))
%     ds.data = obj.data(:,focusColumns(featureInc));
%
%     fEvents = detectEvents(ds);
%
%     obj.dataEventIdx{featureInc} = sort(cat(1,fEvents.onEventsIndex,fEvents.offEventsIndex));
%
%   end
% end
%
%
%
%
% end
%
% %%
% function eventFeatures = extractEventFeatures(obj,featureType,featureDetails)
% % This function outputs a prtDataSet based on the input data and
% % details stored in the structure featureDetails.
%
% if isempty(obj.dataEventIdx)
%   error('Please first run event detection or input times at which to extract events as a cell array with the same number of cells as there are features');
% end
%
% eventFeatures = prtDataSetClass;
%
% if strcmp(featureType,'raw')
%   %% This means that the desired features are simply the data
%   % themselves surrounding the events. Go through each event of each
%   % device.
%   if mod(featureDetails.windowSize,2) == 0
%     error('Please enter an odd windowSize for an even number of features on both sides of the event')
%   end
%
%   windowSize = featureDetails.windowSize;
%
%   halfWin = (windowSize - 1)/2;
%
%   for deviceInc = 1:obj.nFeatures
%
%     if ~isempty(obj.dataEventIdx{deviceInc})
%
%       currentIdx = obj.dataEventIdx{deviceInc};
%
%       currentIdx = currentIdx(currentIdx>halfWin);
%
%       currentIdx = currentIdx(currentIdx<obj.nObservations - halfWin);
%
%       numEvents = max(size(currentIdx));
%
%       featureArray = zeros(numEvents,windowSize);
%
%       %% Go through each time and extract features
%       for eventInc = 1:max(size(currentIdx))
%         startIdx = max(currentIdx(eventInc) - halfWin,1);
%         endIdx = min(currentIdx(eventInc) + halfWin,obj.nObservations);
%
%         % zeroMin is set to 1 if the features extracted are assumed
%         % to be from windows of data floored to a minimum of 0.
%         if isfield(featureDetails,'zeroMin')
%           if featureDetails.zeroMin
%             featureArray(eventInc,:) = obj.data(startIdx:endIdx,deviceInc) - ...
%               min(obj.data(startIdx:endIdx,deviceInc));
%           else
%             featureArray(eventInc,:) = obj.data(startIdx:endIdx,deviceInc);
%           end
%         else
%           featureArray(eventInc,:) = obj.data(startIdx:endIdx,deviceInc);
%         end
%       end
%
%
%       %% FINDME
%       % Establish targets for the given device.
%       % This needs to be updated for consistency.
%       featureLabel = obj.userData.goodColumns(deviceInc);
%       featureTargets = featureLabel*ones(numEvents,1);
%
%       featureDS = prtDataSetClass(featureArray,featureTargets);
%
%       eventTimes = struct('times',num2cell([obj.observationInfo(currentIdx).times]));
%
%       featureDS.observationInfo = eventTimes;
%
%       %% Concatenate onto eventFeatures.
%       eventFeatures = catObservations(eventFeatures,featureDS);
%
%     end
%   end
%
% end
%
% end
%
% %% trainClassifier
% % ds - prtDataSet that contains the features and labels for the set to
% % be
% function trainedClassifier = trainClassifier(obj,ds,classifier)
% trainedClassifier = classifier.train(ds);
% end
%
% function classOuts = runClassifier(obj,ds,trainedClass)
% if ~trainedClass.isTrained
%   error('Please train your classifier')
% end
%
% cOuts = trainedClass.run(ds);
%
% %% Assign the proper labels to the classOuts
% [~,maxLoc] = max(cOuts.data,[],2);
%
% classOuts = cOuts;
%
% classOuts.data = zeros(classOuts.nObservations,1);
%
% uniqueClasses = classOuts.uniqueClasses;
%
% for mInc = 1:cOuts.nClasses
%   classOuts.data(maxLoc == mInc) = uniqueClasses(mInc);
% end
%
%
% end
%
%
%
%
%
%
%
%
% function kOuts = crossVal(obj,ds,classifier,numFolds)
% cOuts = classifier.kfolds(ds,numFolds);
%
% %% Assign the proper labels to the classOuts
% [~,maxLoc] = max(cOuts.data,[],2);
%
% kOuts = cOuts;
%
% kOuts.data = zeros(kOuts.nObservations,1);
%
% uniqueClasses = kOuts.uniqueClasses;
%
% for mInc = 1:cOuts.nClasses
%   kOuts.data(maxLoc == mInc) = uniqueClasses(mInc);
% end
%
% end
%
%
%
%
%
%
%
% %% Check the performance of supervised classification.
% function scoreConfusionMatrix(obj,kOuts)
% figure;
% prtScoreConfusionMatrix(kOuts.data,kOuts.targets)
% end
%
%
%
%
%
%
%
%
% %% How much energy was assigned based on the classified features and associated timestamps?
% function assignedEnergy = assignEnergy(obj,assigned,dataGmms)
% keyboard
%
% assignedEnergy = obj;
%
% assignedEnergy.data = zeros(size(obj.data));
%
% uniqueAssigned = unique(assigned.data);
%
% %% Go through each of the appliances in the assigned object.
% for oInc = 1:max(size(uniqueAssigned))
%   % Sort the timestamps for the given class.
%
% end
% end
%
%
%
