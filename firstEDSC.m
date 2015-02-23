%% Nicholas Czarnek
% SSPACISS Laboratory, Duke University
%
% energyDataSet.m
% This class will be used for analysis of different datasets for energy
% disaggregation.
% z
% THIS NEEDS TO BE MADE MORE ROBUST:
% 1) WHEN SETTING COMPONENTS OR TIMESTAMPS, THE CODE SHOULD CHECK IF THE
% LENGTH OF BOTH IS EQUAL IN A BETTER WAY THAN I DO NOW!
%
%     cEvents               - sparse matrix of detected component events
%     mEvents               - sparse matrix of mains events
%
% Properties:
%   componentStructure    - object containing all component information
%     components            - nObsC x nComponents matrix with each component's data
%     nComponents           - number of components
%     nObsC                 - number of component observations
%     cTS                   - component timestamps
%     cLabels               - nComponents x 1 cell array of labels for each component
%     cType                 - "real" or "apparent"
%     cOverlap              - nObsC x 1 logical array indicating where components
%                             overlap with mains
%     cClass                - nComponents x 1 array indicating type of
%                             appliance
%                           - based on Hart et al 92 categorization
%                           - 1: on/off
%                           - 2: FSM
%                           - 3: continuously variable
%                           - 4: continuously on
%
%   mainStructure         - object containing all mains information
%     mains                 - tLx1 array of mains data
%     nObsM                 - same as mainsLength
%     mTS                   - mains timestamps
%     mType                 - "real" or "apparent"
%     mOverlap              - nObsM x 1 logical array indicating where mains
%                             overlaps with components
%     mType                 - "real" or "apparent"
%
%   assignedPower         - nObsC x nComponents matrix of assigned power for each 
%                           component
%
% Methods:
%   findErrorRMS          - find the RMS error between assigned and actual
%                           power of each component
%   findErrorAbsolute     - find the percent error between the assigned and
%                           actual energy assigned to each component
%   findEnergy            - 
%   plotDensity           - plot ksdensities for each component
%   visualization tools   - 

classdef energyDataSet < prtDataSetClass
  
  %% Components setup.
  properties
    components;
  end
  
  properties (Dependent,SetAccess = protected)
    % These properties cannot be set by the user.  They are dependent on
    % the values from components.
    nComponents;
    nObsC;
  end
  
  properties
    cOverlap;
    cTS;
    cLabels;
    cType;
    cClass;
  end
  
  %% Mains setup
  properties
    mains;
  end
  
  properties (Dependent,SetAccess = protected)
    nObsM;
  end
  
  properties
    mOverlap;
    mTS;
    mType;
  end
  
  %% Setup methods.
  methods
    %% Constructor for class instantiation.
    function obj = energyDataSet(varargin)
      obj = prtUtilAssignStringValuePairs(obj,varargin{:});
    end
    
    %% Set up get functions for the dependent variables.
    function out = get.nComponents(self)
      out = size(self.components,2);
    end
    
    function out = get.nObsC(self)
      out = size(self.components,1);
    end

    function out = get.nObsM(self)
      out = size(self.mains,1);
    end
    
    %% Timestamp length must be consistent with components or mains.
    function obj = set.mTS(obj,mTimestamps)
      if ~isempty(obj.mains)
        if size(obj.mains,1)~=size(mTimestamps,1)
          error('Your timestamps must have the same length as your mains')
        end
      end
      
      
      if max(size(mTimestamps)) == numel(mTimestamps)
        obj.mTS = mTimestamps(:);
      end
    end
    
    function obj = set.cTS(obj,cTimestamps)
      if ~isempty(obj.components)
        if size(obj.components,1)~=size(cTimestamps,1)
          error('Your timestamps must have the same length as your components')
        end
      end
      
      if max(size(cTimestamps)) == numel(cTimestamps)
        obj.cTS = cTimestamps(:);
      end
    end
    
    function obj = set.mains(obj,Mains)
      if ~isempty(obj.mTS)
        if size(obj.mTS,1)~=size(Mains,1)
          error('Your timestamps must have the same length as your mains')
        end
      end
      
      obj.mains = Mains;
    end
    
    
    
    function obj = set.components(obj,Components)
      if ~isempty(obj.cTS)
        if size(obj.cTS,1)~=size(Components,1)
          error('Your timestamps must have the same length as your mains')
        end
      end
      
      obj.components = Components;
    end
    
    
    %% Check the overlap of components and mains.
    function out = get.cOverlap(self)
      if isempty(self.mTS)
        out = zeros(self.nObsC,1);
      else
        out = ismember(self.cTS,self.mTS);
      end
    end
    
    function out = get.mOverlap(self)
      if isempty(self.cTS)
        out = zeros(self.nObsM,1);
      else
        out = ismember(self.mTS,self.cTS);
      end
    end
    
  end
  
  
  
  
  
end