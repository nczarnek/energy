%% Nicholas Czarnek
% SSPACISS Laboratory, Duke university
% 27 February 2015
%
% energyFeatureClass.m
% This prtDataSetClass subclass is used to store energy features extracted
% from energy events.

classdef energyFeatureClass < prtDataSetClass
    methods
        %% Constructor for class instantiation.
        function obj = energyDataSetClass(varargin)
            obj = prtUtilAssignStringValuePairs(obj,varargin{:}); %#ok<NODEF>
        end
        
        %% Include a function to only keep on events based on the observationInfo
        function obj = keepOn(obj,varargin)
            options.classes = [];
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            classes = parsedOut.classes;
            
            if ~isempty(classes)
                obj = obj.retainClasses(classes);
            end
            
            obj = obj.retainObservations([obj.observationInfo.eventType] == 1);
                
        end
        
        %% Include a function to only keep off events based on the observationInfo
        function obj = keepOff(obj,varargin)
            options.classes = [];
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            classes = parsedOut.classes;
            
            if ~isempty(classes)
                obj = obj.retainClasses(classes);
            end
            
            obj = obj.retainObservations([obj.observationInfo.eventType] == 0);

        end
        
        %% Visualize the pca components of the features
        function pcaOuts = getPca(obj,varargin)
            options.nComponents = 3;
            options.classes = [];
            options.onlyOn = 0;
            options.onlyOff = 0;
            options.individual = false;
            options.zmuv = false;
            options.plotPca = true;
            
            parsedOut = prtUtilSimpleInputParser(options,varargin);
            
            nComponents = parsedOut.nComponents;
            classes = parsedOut.classes;
            onlyOn = parsedOut.onlyOn;
            onlyOff = parsedOut.onlyOff;
            individual = parsedOut.individual;
            zmuv = parsedOut.zmuv;
            plotPca = parsedOut.plotPca;
            
            %% Check for stupid
            if onlyOn == 1 && onlyOff == 1
                warning('You can''t have both onlyOn and onlyOff be 1.  onlyOff has been set to 0')
                onlyOff = ~onlyOn;
            end
            
            %% Check if they only want a certain class
            if ~isempty(classes)
                obj = obj.retainClasses(classes);
            end
            
            if onlyOn
                obj = obj.keepOn;
            end
            
            if onlyOff
                obj = obj.keepOff;
            end
            
            if zmuv
                zM = prtPreProcZmuv;
                zM = zM.train(obj);
                obj = zM.run(obj);
            end
            
            if ~individual
                pca = prtPreProcPca('nComponents',nComponents);
                pca = pca.train(obj);
                
                pcaOuts = pca.run(obj);
            else
                for cInc = 1:obj.nClasses
                    currentClass = obj.retainClasses(obj.uniqueClasses(cInc));
                    
                    pca = prtPreProcPca('nComponents',nComponents);
                    pca = pca.train(currentClass);
                    
                    pcaOuts(cInc) = pca.run(currentClass);
                end
            end
            
            
            if plotPca
                for pInc = 1:numel(pcaOuts)
                    figure;
                    imagesc(pcaOuts(pInc))
                    title(['PCA components for ',pcaOuts(pInc).classNames{1}])
                    
                    figure;
                    pcaOuts(pInc).plot
                end
            end
        end
    end
end