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
       onEvents;
       onEventsIndex;
       onEventsTimes;
       house;
       houseNumber;
   end
   
   methods
       
   end
end