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
        function [componentEnergy,energyStruct,energyIdx] = rankComponentsByEnergy(obj,varargin)
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
            
            componentEnergy = componentEnergy/trapz(obj.data(:,1));
            
            [~,eIdx] = sort(componentEnergy,'descend');
            
            featureNames = featureNames(eIdx);
            componentEnergy = componentEnergy(eIdx);
            
            energyCells = cell(nC,1);
            
            for rInc = 1:nC
                featureNames{rInc}(ismember(featureNames{rInc},' ,.:;!')) = [];
                energyCells{rInc} = componentEnergy(rInc);
            end
            
            energyStruct = cell2struct(energyCells',featureNames',2);
            
            energyIdx = eIdx;
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
                figure
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
                        try
                            legendStr = cat(2,legendStr,[eventTimes(eInc).className{1},' on events'],....
                                [eventTimes(eInc).className{1},' off events']);
                        catch
                            legendStr = cat(2,legendStr,[eventTimes(eInc).className,' on events'],....
                                [eventTimes(eInc).className,' off events']);
                        end
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
            options.useCurrentFigure = false;
            options.lineStyle = [];
            options.focusIndices = [];
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            device = parsedOut.device;
            combinedPlot = parsedOut.combinedPlot;
            timeScale = parsedOut.timeScale;
            zeroTimes = parsedOut.zeroTimes;
            XL = parsedOut.XL;
            useCurrentFigure = parsedOut.useCurrentFigure;
            lineStyle = parsedOut.lineStyle;
            focusIndices = parsedOut.focusIndices;
            
            xT = obj.getTimesFromUTC('timeScale',timeScale,'zeroTimes',zeroTimes);
            
            if combinedPlot
                if useCurrentFigure
                    fH = gcf;
                    figure(fH)
                    hold on
                else
                    figure;
                end
                if isempty(lineStyle)
                    if isempty(focusIndices)
                        plot(xT,obj.data(:,device))
                    else
                        plot(xT(focusIndices),obj.data(focusIndices,device))
                    end
                else
                    if isempty(focusIndices)
                        plot(xT,obj.data(:,device),lineStyle)
                    else
                        plot(xT(focusIndices),obj.data(focusIndices,device),lineStyle)
                    end
                end
                legend(obj.getFeatureNames(device))
                xlabel(['Time (',timeScale,')'])
                ylabel('Power (W)')
                title('Combined plot of input devices')
                if ~isempty(XL)
                    xlim([XL(1) XL(2)])
                end
                
                if useCurrentFigure
                    %% Modify the legend.
                    lH = findobj(gcf,'tag','legend');
                    legendString = lH.String;
                    legendString = cat(2,legendString,obj.getFeatureNames(device));
                    legend(legendString);
                end
                
            else
                for dInc = 1:numel(device)
                    if useCurrentFigure
                        fH = gcf;
                        figure(fH)
                        hold on
                    else
                        figure;
                    end
                    %                     plot(xT,obj.data(:,device(dInc)))
                    if isempty(lineStyle)
                        if isempty(focusIndices)
                            plot(xT,obj.data(:,device(dInc)))
                        else
                            plot(xT(focusIndices),obj.data(focusIndices,device(dInc)))
                        end
                    else
                        if isempty(focusIndices)
                            plot(xT,obj.data(:,device(dInc)),lineStyle)
                        else
                            plot(xT(focusIndices),obj.data(focusIndices,device(dInc)),lineStyle)
                        end
                    end
                    
                    lH = findobj(gcf,'tag','legend');
                    noCat = false;
                    if isempty(lH)
                        legend(obj.getFeatureNames(device(dInc)),'interpreter','none')
                        lH = findobj(gcf,'tag','legend');
                        noCat = true;
                    end
                    xlabel(['Time (',timeScale,')'])
                    ylabel('Power (W)')
                    title(obj.getFeatureNames(device(dInc)))
                    if ~isempty(XL)
                        xlim([XL(1) XL(2)])
                    end
                    
                    if useCurrentFigure
                        %% Modify the legend.
                        legendString = lH.String;
                        if ~noCat
                            legendString = cat(2,legendString,obj.getFeatureNames(device(dInc)));
                            legend(legendString,'interpreter','none');
                        end
                    end
                    
                end
            end
            
            
            
        end
        
        %% This function allows us to see a breakdown of the component 
        % energy as a percent of the total energy consumed.
        function fH = barEnergy(obj,varargin)
            options.newFigure = true;
            options.retainFeatures = [obj.getFeatureInfo.pecanClass]';
            options.newSum = true;
            options.combineFeatures = [];
            options.percent = false;
            options.ascend = true;
            options.vertical = false;
            options.upsideDown = false;
            parsedOuts = prtUtilSimpleInputParser(options,varargin);
            newFigure = parsedOuts.newFigure;
            retainFeatures = parsedOuts.retainFeatures(:);
            newSum = parsedOuts.newSum;
            combineFeatures = parsedOuts.combineFeatures;
            percent = parsedOuts.percent;
            ascend = parsedOuts.ascend;
            vertical = parsedOuts.vertical;
            upsideDown = parsedOuts.upsideDown;
            
            %% Only retain the features that were sent in.
            pCs = [obj.getFeatureInfo.pecanClass]';
            plotSet = obj.retainFeatures(ismember(pCs,cat(1,pCs(1),retainFeatures)));
            
            if newSum
                plotSet.data(:,1) = sum(plotSet.data(:,2:end),2);
            end
            
            [a,~,c] = plotSet.rankComponentsByEnergy;
            
            if percent
                a = a*100;
            end
            
            componentNames = plotSet.getFeatureNames(c(end:-1:2));
            
            if newFigure
                figure;
            end
            
            if ascend
                fH = bar(a(plotSet.nFeatures:-1:2));
                set(gca,'xtick',1:plotSet.nFeatures - 1)
                set(gca,'xticklabel',componentNames,'FontSize',13)
            else
                fH = bar(a(2:plotSet.nFeatures));
                set(gca,'xtick',1:plotSet.nFeatures - 1)
                set(gca,'xticklabel',componentNames(end:-1:1),'FontSize',13)
            end
            
            xlim([0 plotSet.nFeatures])
            if vertical
                if upsideDown
                    view(180,90)
                end
            else
                view(90,-90)
            end
            
            title('Fraction of energy consumed per appliance')
            
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
            % concatenated features
            options.catFeatures = false;
            
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            
            windowInS = parsedOut.windowInS;
            devices = parsedOut.devices;
            classNumber = parsedOut.classNumber;
            featureType = parsedOut.featureType;
            className = parsedOut.className;
            zeroFeatures = parsedOut.zeroFeatures;
            catFeatures = parsedOut.catFeatures;
            
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
            
            if catFeatures
                energyFeats = energyFeatureClass;
            else
                energyFeats = repmat(energyFeatureClass,numel(devices),1);
            end
            
            %% Extract features from each event
            for dInc = 1:numel(devices)
                currentDevice = devices(dInc);
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
                            featureVector = obj.data(featureIdx,currentDevice)';
                            
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
                            featureVector = obj.data(featureIdx,currentDevice)';
                            
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
                            featureVector = obj.data(featureIdx,currentDevice)';
                            
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
                            featureVector = obj.data(featureIdx,currentDevice)';
                            
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
                
                if catFeatures
                    energyFeats = energyFeats.catFeatures(currentFeats);
                    energyFeats.targets = currentFeats.targets;
                    energyFeats.observationInfo = currentFeats.observationInfo;
                    energyFeats.classNames = className;
                else
                    energyFeats(dInc) = currentFeats;
                end
                
            end
            
        end
        
        %% Determine the average power consumed for BLUED
        function powerOut = findAveragePowerBlued(obj,eventTimes,allEnergyFeatures,varargin)
            options.zeroWindow = 61;
            options.removeTails = true;
            options.tailPercent = 5;
            parsedOuts = prtUtilSimpleInputParser(options,varargin);
            zeroWindow = parsedOuts.zeroWindow;
            removeTails = parsedOuts.removeTails;
            tailPercent = parsedOuts.tailPercent;
            %% Get the time in seconds.
            xT = obj.getTimesFromUTC('timeScale','s');
            yT = obj.getTimesFromUTC('timeScale','days','zeroTimes',false);
            
            %% Assume we're extracting from the aggregate load
            powerOut = repmat(struct('averagePower',[],'averageTimeInS',[],...
                'classNumber',[],'fromAggregate',[],'numberInstances',[]),...
                allEnergyFeatures.nClasses,1);
            
            
            %% Ensure that the event times are all sorted properly.
            for eInc = 1:numel(eventTimes)
                [eventTimes(eInc).onEventsTimes,sortIdx] = unique(eventTimes(eInc).onEventsTimes);
                eventTimes(eInc).onEventsIndex = eventTimes(eInc).onEventsIndex(sortIdx);
                eventTimes(eInc).onClass = eventTimes(eInc).onClass(sortIdx);
                
                [eventTimes(eInc).offEventsTimes,sortIdx] = unique(eventTimes(eInc).offEventsTimes);
                eventTimes(eInc).offEventsIndex = eventTimes(eInc).offEventsIndex(sortIdx);
                eventTimes(eInc).offClass = eventTimes(eInc).offClass(sortIdx);
            end
            
            
            for classInc = 1:allEnergyFeatures.nClasses
                %% Find the events corresponding to the current class label.
                currentClass = false(numel(eventTimes),1);
                for eInc = 1:numel(eventTimes)
                    if eventTimes(eInc).classNumber == allEnergyFeatures.uniqueClasses(classInc);
                        currentClass(eInc) = true;
                    end
                end
                
                eventIdx = find(currentClass);
                currentFeature = eventTimes(currentClass).classNumber;
                %% Go through each of the events.
                energyUsed = zeros(numel(eventTimes(eventIdx).onEventsTimes),1);
                onTime = energyUsed;
                
                for eInc = 1:numel(eventTimes(eventIdx).onEventsTimes)
                    %% Find the next off time that is larger than the current time
                    currentOn = eventTimes(eventIdx).onEventsTimes(eInc);
                    offIdx = find(eventTimes(eventIdx).offEventsTimes>currentOn,1,'first');
                    
                    if ~isempty(offIdx)
                        currentOff = eventTimes(eventIdx).offEventsTimes(offIdx);
                        focusIndices = yT>=currentOn & yT<=currentOff;
                        
                        if all(focusIndices == false)
                            onIdx = find(yT>=currentOn,1,'first');
                            focusIndices(onIdx) = true;
                        end
                        
                        startInS = xT(find(focusIndices,1,'first'));
                        endInS = xT(find(focusIndices,1,'last'));
                        
                        powerSignature = obj.data(focusIndices,1);
                        
                        %% Zero out the power
                        zeroIndices = yT>=currentOn - zeroWindow & yT<=currentOff + zeroWindow;
                        
                        if all(zeroIndices == false)
                            onIdx = find(yT>=currentOn,1,'first');
                            zeroIndices(onIdx) = true;
                        end
                        
                        windowSignature = obj.data(zeroIndices,1);
                        powerSignature = powerSignature - min(windowSignature);
                        
                        try
                            energyUsed(eInc) = trapz(xT(focusIndices),powerSignature);
                        catch
                            %% Approximate the energy used by on for one minute.
                            % This is appropriate for the downsampled data.
                            energyUsed(eInc) = powerSignature(1)*60;
                        end
                        
                        onTime(eInc)= endInS - startInS + 1;
                    end
                end
                %% Remove the times that were not covered.
                zeroTimes = find(energyUsed == 0);
                
                energyUsed(zeroTimes) = [];
                onTime(zeroTimes) = [];
                
                numInstances = numel(energyUsed);
                
                if numInstances ~= 0
                    %% Handle outliers
                    if removeTails
                        [energyUsed,newIdx] = sort(energyUsed);
                        onTime = onTime(newIdx);
                        
%                         
%                         startKeep = round(tailPercent/100*numInstances);
%                         endKeep = round((100 - tailPercent)/100*numInstances);
%                         

                        %% Use the median and one point on both sides for the power calculation
                        if ~mod(numInstances,2)
                            %% Even
                            energyUsed = energyUsed(numInstances/2:numInstances/2+1);
                            newIdx = newIdx(numInstances/2:numInstances/2+1);
                            onTime = onTime(numInstances/2:numInstances/2+1);
                        else
                            %% odd
                            medianEnergy = median(energyUsed);
                        
                            medianIdx = find(energyUsed == medianEnergy,1,'first');
                            
                            try
                                energyUsed = energyUsed(medianIdx - 1:medianIdx + 1);
                                newIdx = newIdx(medianIdx - 1:medianIdx + 1);
                                onTime = onTime(medianIdx - 1:medianIdx + 1);
                            catch
                                energyUsed = energyUsed(medianIdx);
                                newIdx = newIdx(medianIdx);
                                onTime = onTime(medianIdx);
                            end
                        end
                        
%                         
%                         focusIdx = max(startKeep,1):min(endKeep,numInstances);
%                         
%                         energyUsed = energyUsed(focusIdx);
%                         newIdx = newIdx(focusIdx);
%                         onTime = onTime(newIdx);
                    end
                    
                end
                
                %% Find the average power for the current device
                powerOut(classInc).averagePower = sum(energyUsed)/sum(onTime);
                powerOut(classInc).averageTimeInS = mean(onTime);
                powerOut(classInc).classNumber = currentFeature;
                powerOut(classInc).fromAggregate = true;
                powerOut(classInc).className = eventTimes(currentClass).className;
                powerOut(classInc).numberInstances = numel(energyUsed);
                
                
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
            options.zeroWindow = 61;
            options.bluedDevices = [];
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            features = parsedOut.features;
            fromAggregate = parsedOut.fromAggregate;
            zeroPower = parsedOut.zeroPower;
            removeTails = parsedOut.removeTails;
            tailPercent = parsedOut.tailPercent;
            zeroWindow = parsedOut.zeroWindow;
            zeroWindow = zeroWindow/86400;
            pecanClasses = [obj.featureInfo.pecanClass]';
            
            %% Ensure that the event times are all sorted properly.
            for eInc = 1:numel(eventTimes)
                [eventTimes(eInc).onEventsTimes,sortIdx] = unique(eventTimes(eInc).onEventsTimes);
                eventTimes(eInc).onEventsIndex = eventTimes(eInc).onEventsIndex(sortIdx);
                eventTimes(eInc).onClass = eventTimes(eInc).onClass(sortIdx);
                
                [eventTimes(eInc).offEventsTimes,sortIdx] = unique(eventTimes(eInc).offEventsTimes);
                eventTimes(eInc).offEventsIndex = eventTimes(eInc).offEventsIndex(sortIdx);
                eventTimes(eInc).offClass = eventTimes(eInc).offClass(sortIdx);
            end
            
            %% Get the time in seconds.
            xT = obj.getTimesFromUTC('timeScale','s');
            yT = obj.getTimesFromUTC('timeScale','days','zeroTimes',false);
            
            %% Assume that we're extracting from the aggregate load
            extractionFeature = 1;
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
                    aveT = tic;
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

                            if all(focusIndices == false)
                                onIdx = find(yT>=currentOn,1,'first');
                                focusIndices(onIdx) = true;
                            end
                            
                            startInS = xT(find(focusIndices,1,'first'));
                            endInS = xT(find(focusIndices,1,'last'));
                            
                            powerSignature = obj.data(focusIndices,extractionFeature);
                            
                            if zeroPower
                                zeroIndices = yT>=currentOn - zeroWindow & yT<=currentOff + zeroWindow;
                                
                                if all(zeroIndices == false)
                                    onIdx = find(yT>=currentOn,1,'first');
                                    zeroIndices(onIdx) = true;
                                end
                                
                                windowSignature = obj.data(zeroIndices,extractionFeature);
                                
                                powerSignature = powerSignature - min(windowSignature);
                            end
                            
                            try
                                energyUsed(eInc) = trapz(xT(focusIndices),powerSignature);
                            catch
                                if numel(find(focusIndices)) == 1
                                    %% This for downsampled datasets, this 
                                    % means that the original event lasted
                                    % less than 1 minute.  Therefore, just
                                    % consider the energy used to be the
                                    % power over 60 s.
                                    energyUsed(eInc) = powerSignature * 60;
                                end
                            end
                            
                            onTime(eInc)= endInS - startInS + 1;
                        end
                        
                    end
                    
                    %% Remove the times that were not covered.
                    zeroTimes = find(energyUsed == 0);
                    
                    energyUsed(zeroTimes) = [];
                    onTime(zeroTimes) = [];
                    
                    numInstances = numel(energyUsed);
                    
                    if numInstances ~= 0
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
                    
                    aveStop = toc(aveT);
                    
                    try
                        fprintf(1,[powerOut(dInc).className,' average found, ',num2str(aveStop),'s\n'])
                    catch
                        fprintf(1,[powerOut(dInc).className{1},' average found, ',num2str(aveStop),'s\n'])
                    end
                    
                end
            end
            
% goes right after if isempty(offIdx)            
%                             if zeroPower
%                                 focusIndices = yT>=currentOn - zeroWindow & yT<=currentOff + zeroWindow;
%                                 
%                             else
%                                 focusIndices = yT>=currentOn & yT<=currentOff;
%                             end
        end
        
        
        %% Assign power to the BLUED data set
        function assignedPower = assignPowerBlued(obj,eventTimes,averagePowers,varargin)
            options.useOff = false;
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            useOff = parsedOut.useOff;
            
            %% Establish the new dataset and zero it out
            assignedPower = energyDataSetClass;
            assignedPower.data = zeros(obj.nObservations,numel(averagePowers)+1);
            
            %% Get the times.
            xT = obj.getTimesFromUTC('timeScale','s');
            yT = obj.getTimesFromUTC('timeScale','days','zeroTimes',false);
            
            featureNames = {'use'};
            featureInfo = repmat(struct('className',[],'classNumber',[]),numel(averagePowers)+1,1);
            featureInfo(1).className = 'use';
            featureInfo(1).classNumber = 0;
            
            %% Go through the event times and assign power based on the event labels
            for dInc = 1:numel(averagePowers)
                %% Find the event times that corresponds to the current class
                currentClass = false(numel(eventTimes),1);
                for eInc = 1:numel(eventTimes)
                    if eventTimes(eInc).classNumber == averagePowers(dInc).classNumber
                        currentClass(eInc) = true;
                    end
                end
                
                eIdx = find(currentClass);
                
                
                if ~useOff
                    focusIdx = false(assignedPower.nObservations,1);
                    
                    for eInc = 1:numel(eventTimes(eIdx).onEventsTimes)
                        startIdx = find(yT>=eventTimes(eIdx).onEventsTimes(eInc),1,'first');
                        
                        startXT = xT(startIdx);
                        endXT = startXT + averagePowers(dInc).averageTimeInS;
                        
                        focusIdx(xT>=startXT&xT<=endXT) = true;
                    end
                    assignedPower.data(focusIdx,dInc + 1) = averagePowers(dInc).averagePower;
                    
                    featureNames = cat(1,featureNames,averagePowers(dInc).className);
                    featureInfo(dInc+1).className = featureNames{dInc+1};
                    featureInfo(dInc+1).classNumber = averagePowers(dInc).classNumber;
                end
                
                
            end
            
            assignedPower = assignedPower.setFeatureNames(featureNames);
            assignedPower = assignedPower.setFeatureInfo(featureInfo);
            
            assignedPower.observationInfo = obj.observationInfo;
            
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
                if ~isempty(eventTimes(dInc).classNumber)
                    assignStart = tic;
                    
                    fIdx = pecanClasses == eventTimes(dInc).classNumber;
                    
                    aPIdx = false(numel(averagePower),1);
                    
                    %% Find the averagePower structure corresponding to the current class
                    for aInc = 1:numel(averagePower)
                        if averagePower(aInc).classNumber == eventTimes(dInc).classNumber
                            aPIdx(aInc) = true;
                        end
                    end
                    
                    if any(aPIdx)
                        if ~useOff
                            %% Make assignments based on the eventTimes
                            % First get the focusIdx, then make the assignments.
                            %                         tS2 = tic;%
                            focusIdx = false(assignedPower.nObservations,1);
                            for eInc = 1:numel(eventTimes(dInc).onEventsTimes)
                                startIdx = find(yT>=eventTimes(dInc).onEventsTimes(eInc),1,'first');
                                
                                startXT = xT(startIdx);
                                endXT = startXT + averagePower(aPIdx).averageTimeInS;
                                
                                focusIdx(xT>=startXT&xT<=endXT) = true;
                                
                            end
                            assignedPower.data(focusIdx,fIdx) = averagePower(aPIdx).averagePower;
                            %                         tStop2 = toc(tS2);% 2 s
                            
                            
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
                    
                    assignStop = toc(assignStart);
                    try
                        fprintf(1,[eventTimes(dInc).className,' assignment complete, ',num2str(assignStop),'s\n'])
                    catch
                        fprintf(1,[eventTimes(dInc).className{1},' assignment complete, ',num2str(assignStop),'s\n'])
                    end
                end
            end
            
        end
        
        function devicePercentExplained = findDevicePercentExplained(obj,assignedPower)
            
            devicePercentExplained = zeros(obj.nFeatures,1);
            
            for dInc = 1:obj.nFeatures
                actualEnergyUsed = obj.findComponentEnergy('device',dInc);
                assignedEnergy = assignedPower.findComponentEnergy('device',dInc);
                
                devicePercentExplained(dInc,1) = assignedEnergy/actualEnergyUsed;
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
        function [systemResults,systemParameters] = runSupervisedSystem(obj,eventTimes,varargin)
            % Note that this assumes that the aggregate load is sent in as
            % the first eventTimes array component.
            
            %%
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
            % Should the true classes be used for the assignment step?
            options.useTrueClasses = false;
            % Do you only want to focus on detection?
            %             options.onlyRunDetection = false;
            % Do you want to zmuv the components?
            options.zmuvFeatures = true;
            % What classifier do you want? Default is knn with k = 5
            options.classifier = 'knn';
            % If KNN, what is k?
            options.k = 5;
            % If random forrest, how many trees?
            options.nTrees = 20;
            % If SVM, rbf or linear?
            options.useRbf = true;
            % Should the power assignments be made from the submetered data or
            % the aggregate load?
            options.fromAggregate = true;
            % Do you want to find the detection performance?
            options.findDetectionPerformance = false;
            % Do you want to find the classification performance?
            options.findClassificationPerformance = true;
            % Do you want to assign energy to devices?
            options.findAssignmentPerformance = true;
            % Should any of the classes be combined into one?
            options.otherClasses = [];
            % Should we score both on and off times or just on?
            options.onlyOn = true;
            % Allow the features to be saved
            options.saveFeatureDir = [];
            % Allow the features to be loaded.
            options.loadFeatureDir = [];
            % Which detector do you want to use?
            options.detectorType = 'glr';
            % Do you want to combine any events?
            options.combinedClasses = {};
            % Do you want to name the combined classes?
            options.combinedClassNames = {};
            % Which devices do you want to extract data from?
            options.devices = 1;
            % Do you want to concatenate faetures?
            options.catFeatures = true;
            % What do you want your min detection threshold to be for the
            % detector?
            options.minDetectorThreshold = 0.01;
            % Do you want the cross val keys to be set?
            options.xValKeys = [];
            % Allow allEnergyFeatures to be sent in
            options.allEnergyFeatures = energyFeatureClass;
            options.findAssignmentErrors = true;
            options.extraSmooth = false;
            options.extraSmoothWindowInS = 30;
            options.comparisonData = obj;
            
            parsedOuts = prtUtilSimpleInputParser(options,varargin);
            systemParameters = parsedOuts;
            
            nXFolds = parsedOuts.nXFolds;
            detectionHalfWinInS = parsedOuts.detectionHalfWinInS;
            eventHaloInS = parsedOuts.eventHaloInS;
            farThreshold = parsedOuts.farThreshold;
            extractionWindow = parsedOuts.extractionWindow;
            usePca = parsedOuts.usePca;
            nPcaComponents = parsedOuts.nPcaComponents;
            useTrueTimes = parsedOuts.useTrueTimes;
            %             onlyRunDetection = parsedOuts.onlyRunDetection;
            zmuvFeatures = parsedOuts.zmuvFeatures;
            classifier = parsedOuts.classifier;
            k = parsedOuts.k;
            nTrees = parsedOuts.nTrees;
            useRbf = parsedOuts.useRbf;
            useTrueClasses = parsedOuts.useTrueClasses;
            fromAggregate = parsedOuts.fromAggregate;
            findDetectionPerformance = parsedOuts.findDetectionPerformance;
            findClassificationPerformance = parsedOuts.findClassificationPerformance;
            otherClasses = parsedOuts.otherClasses;
            findAssignmentPerformance = parsedOuts.findAssignmentPerformance;
            onlyOn = parsedOuts.onlyOn;
            saveFeatureDir = parsedOuts.saveFeatureDir;
            loadFeatureDir = parsedOuts.loadFeatureDir;
            detectorType = parsedOuts.detectorType;
            combinedClasses = parsedOuts.combinedClasses;
            combinedClassNames = parsedOuts.combinedClassNames;
            devices = parsedOuts.devices;
            catFeatures = parsedOuts.catFeatures;
            minDetectorThreshold = parsedOuts.minDetectorThreshold;
            xValKeys = parsedOuts.xValKeys;
            allEnergyFeatures = parsedOuts.allEnergyFeatures;
            findAssignmentErrors = parsedOuts.findAssignmentErrors;
            extraSmooth = parsedOuts.extraSmooth;
            extraSmoothWindowInS = parsedOuts.extraSmoothWindowInS;
            comparisonData = parsedOuts.comparisonData;
            
            if isempty(combinedClassNames)
                combinedClassNames = cell(numel(combinedClasses),1);
            end
            
            %% Combine classes into one
            if ~isempty(otherClasses)
                otherIdx = false(numel(eventTimes),1);
                for dInc = 1:numel(eventTimes)
                    if ~isempty(eventTimes(dInc).classNumber)
                        if any(otherClasses == eventTimes(dInc).classNumber)
                            otherIdx(dInc) = true;
                        end
                    end
                end
                
                firstDevice = find(otherIdx,1,'first');
                
                eventTimes(firstDevice) = combineEnergyEvents(eventTimes(otherIdx),...
                    'classNumber',1000);
                
                otherIdx(firstDevice) = false;
            end
            
            if ~isempty(combinedClasses)
                combinedIdx = false(numel(eventTimes),1);
                for cellInc = 1:numel(combinedClasses)
                    cellIdx = false(numel(eventTimes),1);
                    
                    for dInc = 1:numel(eventTimes)
                        if ~isempty(eventTimes(dInc).classNumber)
                            if any(combinedClasses{cellInc} == eventTimes(dInc).classNumber)
                                cellIdx(dInc) = true;
                            end
                        end
                    end
                    
                    firstDevice = find(cellIdx,1,'first');
                    %                     firstDevice = find(combinedClasses{cellInc},1,'first');
                    
                    eventTimes(firstDevice) = combineEnergyEvents(eventTimes(cellIdx),...
                        'classNumber',1000 + cellInc);
                    
                    try
                        if ~isempty(combinedClassNames{cellInc})
                            eventTimes(firstDevice).className = combinedClassNames{cellInc};
                        end
                    catch
                        error('Make sure that names are sent in as cells')
                    end
                    
                    cellIdx(firstDevice) = false;
                    
                    combinedIdx = combinedIdx | cellIdx;
                end
            end
            
            if ~isempty(combinedClasses)&&~isempty(otherClasses)
                otherAndCombinedIdx = otherIdx | combinedIdx;
            elseif ~isempty(combinedClasses)
                otherAndCombinedIdx = combinedIdx;
            elseif ~isempty(otherClasses)
                otherAndCombinedIdx = otherIdx;
            else
                otherAndCombinedIdx = false(numel(eventTimes),1);
            end
            
            eventTimes = eventTimes(~otherAndCombinedIdx);
            
            
            %% Split up the data evenly into the number of folds
            nObsPerFold = round(obj.nObservations/nXFolds);
            
            xT = obj.getTimesFromUTC('timeScale','days','zeroTimes',false);
            
            %% Set up the aggregate energyEventClass
            detectedTestingEvents = energyEventClass;
            %             assignedEventsTotal = repmat(energyEventClass,numel(eventTimes),1);
            
            %% Extract features for the classifier.
%             allEnergyFeatures = energyFeatureClass;
            
            
            
            featureStart = tic;
            
            %% Check if the feature dataset was already sent in.
            if allEnergyFeatures.nFeatures == 0
                if ~exist(fullfile(loadFeatureDir,'allEnergyFeatures.mat'))
                    for tInc = 2:numel(eventTimes)
                        if ~isempty(eventTimes(tInc).onEventsTimes)
                            tStart = tic;
                            
                            currentFeatures = obj.extractEventData(eventTimes(tInc),'className',...
                                eventTimes(tInc).className,'classNumber',...
                                eventTimes(tInc).classNumber,'featureType','on',...
                                'windowInS',extractionWindow,...
                                'devices',devices,...
                                'catFeatures',catFeatures);
                            
                            allEnergyFeatures = catObservations(allEnergyFeatures,currentFeatures);
                            tStop = toc(tStart);
                            
                            try
                                className = currentFeatures.classNames{1};
                            catch
                                className = currentFeatures.classNames;
                            end
                            
                            fprintf(1,[num2str(tStop),'s: ',num2str(currentFeatures.nObservations),' ',className,' features extracted\n'])
                        end
                    end
                    
                    if ~isempty(saveFeatureDir)
                        save(fullfile(saveFeatureDir,'allEnergyFeatures.mat'),'allEnergyFeatures')
                    end
                else
                    load(fullfile(loadFeatureDir,'allEnergyFeatures.mat'))
                end
            end
            featureStop = toc(featureStart);
            fprintf(1,['Features extracted, ',num2str(featureStop),'s\n'])
            featureTimeStamps = [allEnergyFeatures.observationInfo.timestamp]';
            
            
            if findClassificationPerformance
                
                %% Set up the classifier
                switch lower(classifier)
                    case 'knn'
                        prtClassifier = prtClassKnn('k',k) + prtDecisionMap;
                    case 'svm'
                        if useRbf
                            prtClassifier = prtClassBinaryToMaryOneVsAll + prtDecisionMap;
                            prtClassifier.actionCell{1}.baseClassifier = prtClassLibSvm;
                        else
                            prtClassifier = prtClassBinaryToMaryOneVsAll + prtDecisionMap;
                            prtClassifier.actionCell{1}.baseClassifier = prtClassLibSvm('kernelType',0);
                        end
                    case 'rf'
                        prtClassifier = prtClassTreeBaggingCap('nTrees',nTrees) + prtDecisionMap;
                end
                
                classStart = tic;

                %% Take the top PCs if desired after zmuving, if desired
                if usePca
                    if zmuvFeatures
                        pcedFeatures = allEnergyFeatures.getPca('onlyOn',1,...
                            'zmuv',true,'plotPca',false);
                    else
                        pcedFeatures = allEnergyFeatures.getPca('onlyOn',1,...
                            'plotPca',false);
                    end
                    
                else
                    pcedFeatures = allEnergyFeatures;
                end
            end
            
            %% Perform 5 fold cross val if desired.
            if findClassificationPerformance % ~40-50 s for 18k 3 feature observations with knn and k = 5
                tStart = tic;
                
                if isempty(xValKeys)
                    classificationPerformance = prtClassifier.kfolds(pcedFeatures,nXFolds);
                else
                    try
                        %% Allow multiple cross val folds to be sent in
                        for xvInc = 1:size(xValKeys,2)
                            classificationPerformance(xvInc) = ...
                                prtClassifier.crossValidate(pcedFeatures,xValKeys(:,xvInc));
                        end
                    catch
                        warning('Check your cross val keys against the number of observations. \n Defaulting back to random keys.');
                        classificationPerformance = prtClassifier.kfolds(pcedFeatures,nXFolds);
                    end
                end
                tStop = toc(tStart);
                fprintf(1,[num2str(tStop),'s for ',classifier,' ',num2str(nXFolds),'-fold cross val\n'])
                classStop = toc(classStart);
                fprintf(1,[num2str(nXFolds),' cross val complete, ',num2str(classStop),'s\n'])
            else
                classificationPerformance = [];
            end
            
            %% Set up the testing feature object to store all of the extracted features.
            allTestingFeatures = energyFeatureClass;
            
            %% Run through each fold separately.
            for fInc = 1:nXFolds
                foldStart = tic;
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
                    trainingTimes(tInc).offClass = trainingTimes(tInc).offClass(keepOff);
                    
                    trainingTimes(tInc).onEventsIndex = trainingTimes(tInc).onEventsIndex(keepOn);
                    trainingTimes(tInc).onEventsTimes = trainingTimes(tInc).onEventsTimes(keepOn);
                    trainingTimes(tInc).onClass = trainingTimes(tInc).onClass(keepOn);
                    
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
                    testingTimes(tInc).onClass = testingTimes(tInc).onClass(~keepOn);
                    
                    testingTimes(tInc).offEventsIndex = testingTimes(tInc).offEventsIndex(~keepOff);
                    testingTimes(tInc).offEventsTimes = testingTimes(tInc).offEventsTimes(~keepOff);
                    testingTimes(tInc).offClass = testingTimes(tInc).offClass(~keepOff);
                end
                
                
                if ~useTrueTimes
                    if findClassificationPerformance
                        %% Run the event detector on the training data
                        detectedTrainingEvents = detectEnergyEvents(trainData,...
                            'halfWindowInS',detectionHalfWinInS,'device',1,...
                            'detectorType',detectorType,'extraSmooth',extraSmooth,...
                            'extraSmoothWindowInS',extraSmoothWindowInS);
                        
                        
                        tStart = tic;
                        %% Score the performance.
                        
                        trainingEventDetectionPerformance = scoreEventDetectionQuickly(...
                            detectedTrainingEvents,trainingTimes(1),eventHaloInS,'onlyOn',onlyOn,...
                            'minThreshold',minDetectorThreshold);
                        
                        
                        tStop = toc(tStart);
                        fprintf(1,[num2str(tStop),'s for event detection\n'])
                        
                        
                        onThreshIdx = find(trainingEventDetectionPerformance.onFa>=farThreshold,1,'first');
                        
                        onEventThreshold = trainingEventDetectionPerformance.onThresholds(onThreshIdx);
                    end
                end
                
                %% Get the training features.
                startTime = xT(startIdx);
                endTime = xT(endIdx);
                
                if findClassificationPerformance
                    if fInc == 1
                        energyFeatures = allEnergyFeatures.retainObservations(featureTimeStamps>=endTime);
                    elseif fInc == nXFolds
                        energyFeatures = allEnergyFeatures.retainObservations(featureTimeStamps<=startTime);
                    else
                        firstFeatures = allEnergyFeatures.retainObservations(featureTimeStamps<=startTime);
                        lastFeatures = allEnergyFeatures.retainObservations(featureTimeStamps>=endTime);
                        energyFeatures = catObservations(firstFeatures,lastFeatures);
                    end
                    
                    
                    %% Modify the features based on input options to take the
                    % principal components and the zscore
                    if usePca
                        [inputFeats,inputPca,zmuv] = energyFeatures.getPca('nComponents',nPcaComponents,...
                            'plotPca',false,'zmuv',zmuvFeatures);
                    else
                        inputFeats = energyFeatures;
                        if zmuvFeatures
                            zM = prtPreProcZmuv;
                            zM = zM.train(energyFeatures);
                            inputFeats = zM.run(energyFeatures);
                        end
                    end
                    
                    prtClassifier = prtClassifier.train(inputFeats);
                end
                
                if findClassificationPerformance
                    %% Are events detected or just based on truth?
                    if ~useTrueTimes
                        %% Detect events in the testing data.
                        testingDetectedEvents = detectEnergyEvents(testData,...
                            'halfWindowInS',detectionHalfWinInS,'device',1,'threshold',...
                            onEventThreshold,'detectorType',detectorType,'extraSmooth',extraSmooth,...
                            'extraSmoothWindowInS',extraSmoothWindowInS);
                    else
                        %% Get the confidences, but use the true times.
                        eventConfidences = detectEnergyEvents(testData,...
                            'halfWindowInS',detectionHalfWinInS,'device',1,...
                            'detectorType',detectorType,'extraSmooth',extraSmooth,...
                            'extraSmoothWindowInS',extraSmoothWindowInS);
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
                    detectedTestingEvents.timeStamps = cat(1,....
                        detectedTestingEvents.timeStamps,...
                        testingDetectedEvents.timeStamps);
                    detectedTestingEvents.confidences = cat(1,...
                        detectedTestingEvents.confidences,...
                        testingDetectedEvents.confidences);
                    detectedTestingEvents.onClass = cat(1,...
                        detectedTestingEvents.onClass,...
                        testingDetectedEvents.onClass);
                    detectedTestingEvents.offClass = cat(1,...
                        detectedTestingEvents.offClass,...
                        testingDetectedEvents.offClass);
                end
                
                
                
                if findClassificationPerformance
                    %% Extract features from the current detected events.
                    testingEvents = allEnergyFeatures.retainObservations(featureTimeStamps>=startTime & ...
                        featureTimeStamps<=endTime);
                    if ~useTrueTimes
                        detectedEnergyFeatures = obj.extractEventData(testingDetectedEvents,...
                            'featureType','on','windowInS',extractionWindow,...
                            'devices',devices,'catFeatures',catFeatures);
                    else
                        detectedEnergyFeatures = testingEvents;
                    end
                    
                    if zmuvFeatures
                        detectedEnergyFeatures = zmuv.run(detectedEnergyFeatures);
                    end
                    
                    
                    if usePca
                        detectedEnergyPca = inputPca.run(detectedEnergyFeatures);
                    else
                        detectedEnergyPca = detectedEnergyFeatures;
                    end
                    
                    %% Add on to the existing features
                    allTestingFeatures = catObservations(allTestingFeatures,...
                        detectedEnergyPca);
                    
                    
                    %% Run the classifier.
                    testingClassOuts = prtClassifier.run(detectedEnergyPca);
                    
                    %% Make the assignments.
                    possibleClasses = inputFeats.uniqueClasses;
                    assignedEvents = repmat(energyEventClass,numel(possibleClasses),1);
                    
                    for eInc = 1:numel(possibleClasses)
                        cNum = possibleClasses(eInc);
                        assignedEvents(eInc).classNumber = cNum;
                        assignedEvents(eInc).className = inputFeats.getClassNames(cNum);
                        assignedEvents(eInc).house = eventTimes(1).house;
                        assignedEvents(eInc).houseNumber = eventTimes(1).houseNumber;
                    end
                    
                    for eventInc = 1:testingClassOuts.nObservations
                        currentClass = testingClassOuts.data(eventInc);
                        
                        %% Find the corresponding class
                        currentLogicals = possibleClasses == currentClass;
                        assignedEvents(currentLogicals).onEventsTimes = ...
                            cat(1,assignedEvents(currentLogicals).onEventsTimes,...
                            testingClassOuts.observationInfo(eventInc).timestamp);
                    end
                    
                    
                    %% Add on to the aggregate events
                    if ~exist('assignedEventsTotal')% initialize the total events
                        assignedEventsTotal = assignedEvents;
                    else
                        for typeInc = 1:numel(assignedEvents)
                            currentClass = assignedEvents(typeInc).classNumber;
                            
                            focusClass = false(numel(assignedEvents),1);
                            
                            for focusInc = 1:numel(focusClass)
                                if assignedEventsTotal(focusInc).classNumber == currentClass
                                    focusClass(focusInc) = true;
                                end
                            end
                            
                            assignedEventsTotal(focusClass).onEventsTimes = ...
                                cat(1,assignedEventsTotal(focusClass).onEventsTimes,...
                                assignedEvents(typeInc).onEventsTimes);
                            
                            
                        end
                    end
                    
                    
                    foldStop = toc(foldStart);
                    fprintf(1,['Fold ',num2str(fInc),' complete, ',num2str(foldStop),'s\n'])
                end
            end
            
            if findAssignmentPerformance
                if strcmp(obj.name,'BLUED')
                    averagePowers = obj.findAveragePowerBlued(eventTimes,allEnergyFeatures);
                    
                    %% Assign power to each device
                    if ~useTrueClasses
                        newData = obj.assignPowerBlued(assignedEventsTotal,averagePowers);
                    else
                        newData = obj.assignPowerBlued(eventTimes,averagePowers);
                    end
                    
                    if findAssignmentErrors
                        %% Remember that we're now only focusing on the aggregate data
                        energySubset = comparisonData.retainFeatures(1);
                        
                        sumData = energySubset;
                        sumData.data = sum(newData.data(:,2:end),2);
                        
                        %% correct for times where the sum is greater than aggregate.
                        greaterThan = sumData.data>energySubset.data;
                        sumData.data(greaterThan) = energySubset.data(greaterThan);
                        
                        assignedPowerPerformance = energySubset.calculateAssignmentErrors(sumData);
                    else
                        assignedPowerPerformance = [];
                    end
                    
                else
                    %% Determine how to assign power.
                    averagePowers = obj.findAveragePower(eventTimes(2:end),...
                        'fromAggregate',fromAggregate);
                    
                    %% Now, we have the aggregated event times and the average power
                    % details, so we can find both the device ROCs and the energy
                    % assignment metrics.
                    %                 for typeInc = 1:numel(eventTimes)
                    %                     assignedEventsTotal(typeInc).classNumber = eventTimes(typeInc).classNumber;
                    %                 end
                    
                    %% Assign power to each device
                    if ~useTrueClasses
                        newData = obj.assignPower(assignedEventsTotal,averagePowers);
                    else
                        newData = obj.assignPower(eventTimes,averagePowers);
                    end
                    
                    if findAssignmentErrors
                        %% Evaluate the assignment performance.
                        assignedPowerPerformance = comparisonData.calculateAssignmentErrors(newData);
                    else
                        assignedPowerPerformance = [];
                    end
                end
            else
                assignedPowerPerformance = [];
                newData = [];
                averagePowers = [];
            end
            
            %% Evaluate the detection performance of the aggregated events.
            if findDetectionPerformance
                detectionPerformance = scoreEventDetectionQuickly(detectedTestingEvents,...
                    eventTimes(1),eventHaloInS,'onlyOn',onlyOn,...
                    'minThreshold',minDetectorThreshold);
            else
                detectionPerformance = [];%scoreEventDetectionQuickly(eventTimes(1),...
                %eventTimes(1),eventHaloInS,'onlyOn',onlyOn);
            end
            
            
            systemResults.assignedPowerPerformance = assignedPowerPerformance;
            systemResults.detectionPerformance = detectionPerformance;
            systemResults.classificationPerformance = classificationPerformance;
            systemResults.newData = newData;
            systemResults.allEnergyFeatures = allEnergyFeatures;
            systemResults.detectedTestingEvents = detectedTestingEvents;
            systemResults.aggregateEvents = eventTimes(1);
            systemResults.averagePowers = averagePowers;
            
        end
    end
    
    
    
    
    
end




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