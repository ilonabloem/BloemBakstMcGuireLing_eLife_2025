% Behavior/eye tracking analysis

clearvars
close all

%setup paths - projectRootPath script uses the code directory as the main project folder
projectDir 	= projectRootPath;
figureDir   = fullfile(projectDir, 'Figures');
dataDir     = fullfile(projectDir, 'Data', 'eyeTracking');

%load in data
load('eyeData.mat')

disp('load')
%% Convert cells to strings/numbers

[allData] 	= behConvertCell(allData,nFiles,nSubs);
disp('convert')
%% Convert to degrees

[allData] 	= behDegree(allData,nFiles,ppd);
disp('degree')
%% Performance and cue width
saveFig 	= false;
allData 	= behAccuracy(allData,nFiles,nWidths,widths,figureDir,saveFig);
disp('accuracy')

%% Remove blink artifact from eye tracking data
showFigs 	= false;
allData 	= behBlink(allData,nSubs,nFiles,showFigs);
%% Eye position check

maxScreen   = sp/ppd;
nLetters    = 20;
eccThresh   = 1.2;

condColors 	= ...
[0.3622    0.1804    0.5000; ...
0.8346    0.4157    0.5000; ...
1.0000    0.5333    0.5000; ...
1.0000    0.7490    0.5000]; 
saveFig   	= false;

behGazePos(allData,nFiles,nSubs,nLetters,widths,eccThresh,condColors,figureDir,saveFig);
disp('Gaze check')