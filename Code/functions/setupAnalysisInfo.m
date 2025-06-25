% All parameters underlying the experiments

function p = setupAnalysisInfo(p)


switch p.Task
    case 'CM' % Continuous monitoring
        p.SubjNames       = {'013','004','019','022','012','021','003','025'};
        p.numRuns         = [11 12 8 11 12 10 10 9];
        p.numAttWindows   = 4;
        p.WidthAttWindow  = [1 3 5 9];
        p.numTRreps       = 2;
        p.numTRblock      = 10;
        p.numLocations    = 20; % 20 letters in display
        p.attLocations    = linspace(0, 2*pi-(2*pi/p.numLocations), p.numLocations);
        p.savestr         = 'attention';
    case 'PCM' % physical contrast manipulation
        p.SubjNames       = {'022','012','021','003','025'};
        p.numRuns         = [2 2 2 2 2];
        p.numAttWindows   = 5;
        p.WidthAttWindow  = [1 3 5 7 9];   
        p.numTRreps       = 10;
        p.numLocations    = 20; % 4 letters in display
        p.attLocations    = linspace(0, 2*pi-(2*pi/p.numLocations), p.numLocations);
        p.savestr         = 'perception';
end

p.ROInames        = {'V1' 'V2' 'V3' };%'hV4'
p.numROIs         = numel(p.ROInames);
p.Voxel_cutoff    = 0.1; % percentage of voxels to use for analysis
p.pFiles          = 'BehavioralData/run_';
p.nSubs           = numel(p.SubjNames);

%% Experimental params
p.ScreenRes       = [1024 768];
p.newPixPerDegree = 38; %% Updated scan viewing conditions: 43
p.PixPerDegree    = 43; %% Updated scan viewing conditions: 43

% visual space in radians (0-360deg)
p.Eccen_ring      = 6;
p.ecc_inner       = 0.7;      % Inner eccen bounds
p.ecc_outer       = 9.1;      % Outer eccen bounds, 8.9 for newer subjects..
p.letter_size     = 2.1;
p.nBins           = 60;
p.binCenters      = linspace(0,2*pi -(2*pi/p.nBins), p.nBins);
p.binEdges        = (p.binCenters(1:end)+[p.binCenters(2:end) 2*pi])/2 ;
% binEdges        = linspace(0,2*pi,nBins+1);
% binCenters      = (binEdges(1:end-1)+binEdges(2:end))/2;

p.numRepeatLoc    = 5;
p.locBlock        = 10; %10 TRs/contrast stimulus in physical contrast localizer

p.CondColors      = parula(numel(p.WidthAttWindow)+1);

% pRF param thresholds
p.R2_thres        = 10;
p.RF_thres        = 0.01;

%% Timing parameters
p.TR              = 1.550;                    % TR duration (in secs)
p.Fixation        = 0.5;
p.StimDisplay     = 0.5;
p.TargetDisplay   = 0.35;                     % event duration (in secs)
p.refresh         = 0.05;
p.Response        = 1.75;                     % Response + feedback time
p.initialPeriod   = 10;                       % num TRs
p.finalPeriod     = 10;                       % num TRs
p.TrialDur        = p.Fixation + p.StimDisplay + p.TargetDisplay + p.Response;
p.numAttLocations = p.numLocations * p.numRepeatLoc; % 4 repeats for all locations
p.attBlockDur     = p.TrialDur*p.numAttLocations;
p.numTotalEvents  = p.numAttLocations;          %(numAttRep * numAttLocations * numAttWindows) * numRepBlock;
switch p.Task
    case {'CM', 'PCM'}
        p.initialPeriod   = 10;               % num TRs
        p.finalPeriod     = 10;
    case 'TD'
        p.initialPeriod = 0;                  % num TRs
        p.finalPeriod = 10;                   % num TRs
end
p.totalTR         = (p.numTotalEvents * p.TrialDur)/p.TR + p.initialPeriod + p.finalPeriod;
p.timeLag         = 3; % in TRs (1.55 s)
p.doSmooth        = true;   

