%% Nicholas Czarnek
% 8 January 2014
% SSPACISS Laboratory, Duke University
%
% energyEventClass.m
% The purpose of this class is to establish a consistent event storing
% module such that we can easily perform event detection and comparisons.

classdef energyEventClass
   properties
       className;
       classNumber;
       confidences;
       timeStamps;
       keepLogicals;
       offEvents;
       offEventsIndex;
       offEventsTimes;
       offClass;
       onEvents;
       onEventsIndex;
       onEventsTimes;
       onClass;
       house;
       houseNumber;
   end
   
   methods
       %% Constructor for class instantiation.
        function obj = energyEventClass(varargin)
            obj = prtUtilAssignStringValuePairs(obj,varargin{:});
        end
        
   end
end