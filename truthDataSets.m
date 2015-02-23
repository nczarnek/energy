%% Nicholas Czarnek
% SSPACISS Laboratory, Duke University
% 12 July 2014
% 
% truthDataSet.m
% The purpose of this function is to truth events from submetered data in
% order to be able to evaluate the performance of event detection modules.
%
% This function takes an input folder containing
% 
% dsFolder: folder which contains datasets.
% dsFile:   specific file to truth

function truthDataSets(dsFolder,dsFile)

folderFiles = dir(fullfile(dsFolder,'*.mat'));

nFiles = size(folderFiles,1);

for fInc = 1:nFiles
  if nargin == 1
    load(fullfile(dsFolder,folderFiles(fInc).name))
    fileName = strrep(folderFiles(fInc).name,'.mat','');
  else
    load(fullfile(dsFolder,dsFile))
    fileName = strrep(dsFile,'.mat','');
  end
  
  %% Make a directory for the current house.
  if ~exist(fullfile(dsFolder,fileName),'dir')
    mkdir(fullfile(dsFolder,fileName))
  end
  
  timeStamps = [energyDataSet.observationInfo.times];
  
  %% Truth each component separately
  % Don't score use since we'll be doing detection on that later.
  for cInc = 2:energyDataSet.nFeatures
    currentData = energyDataSet.data(:,cInc)';
    fName = energyDataSet.getFeatureNames(cInc);

    xTimes = (timeStamps - min(timeStamps))*1440;
    
    figure;
    plot(xTimes,currentData)
    title(fName,'Interpreter','none')
    xlabel('Time (min)')
    xlim([0 max(xTimes)])
    
    fprintf(1,'Check the current figure to see if you want to skip it or not. \n');
    
    %% Check if the current file already exists
    
    fprintf(1,['Current component: ',fName{1},'\n'])
    
    if exist(fullfile(dsFolder,fileName,[fName{1},'.mat']))
      fprintf(1,'A time file already exists for the current component\n')
      keepGoing = input('Enter 0 to skip, 1 to replace, or 2 to append to the current file: ');
    else
      keepGoing = input('Enter 0 to skip current device or 1 to truth the current device: ');
    end
    
    if keepGoing == 0
      % Continue to the next component
      continue
    else
      %% Determine the number of points to include per 2 hour window
      pointsPerDay = round(1/(timeStamps(2) - timeStamps(1)));
      nObsPerWindow = 120; % 2 hours per window
      
      % Put everything in terms of minutes
      sRate = pointsPerDay/1440;
      
      appendOn = 0;
      
      truthedEvents = {};
      
      if keepGoing == 2
        appendOn = 1;
        load(fullfile(dsFolder,fileName,fName{1}))
        
      end
      
      repeatMarking = 1;
      
      while repeatMarking
        try
          
          if appendOn == 0
            truthedEvents = markEvents(currentData,{'event'},nObsPerWindow,'srate',sRate);
          else
            truthedEvents = cat(1,truthedEvents,markEvents(currentData,{'event'},nObsPerWindow,'srate',sRate));
          end
          removeRow = input('\n\nInput 1 to remove the last marking or 0 to continue to the next device.\n If you input 1, you can append to the previous markings: ');
          
          if removeRow
            truthedEvents = truthedEvents(1:end-1,:);
            appendOn = 1;
            repeatMarking = 1;
          else
            repeatMarking = 0;
            appendOn = 0;
          end
          
        catch
          repeatMarking = input('\n\nInput 0 to go to the next device or 1 to repeat this device: ');
        end
      end
      close all
      keyboard
      
      onEvents = cell2mat(truthedEvents(:,2));
      offEvents = cell2mat(truthedEvents(:,3));
      
      if keepGoing == 2
        trueTimes.onTimes = cat(1,trueTimes.onTimes,timeStamps(round(onEvents))');
        trueTimes.offTimes = cat(1,trueTimes.offTimes,timeStamps(round(offEvents))');
        trueTimes.onIdx = cat(1,trueTimes.onIdx,round(onEvents));
        trueTimes.offIdx = cat(1,trueTimes.offIdx,round(offEvents));
      else
        trueTimes.onTimes = timeStamps(round(onEvents))';
        trueTimes.offTimes = timeStamps(round(offEvents))';
        trueTimes.onIdx = round(onEvents);
        trueTimes.offIdx = round(offEvents);
      end
      
      %% Save it.
%       save(fullfile(dsFolder,fileName,fName{1}),'trueTimes')
      
    end
    
    close all
    
  end
  
  
  if nargin == 1
    break
  end
end