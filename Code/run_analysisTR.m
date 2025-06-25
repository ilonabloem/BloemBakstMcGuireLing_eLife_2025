function out = run_analysisTR(Task, recompute, numCores)

%% Run attentional window analysis
%
% function out = s_runAnalysis(task, recompute)
%
% Task can be: 'CM' (Continuous Monitoring) or 'PCM' (Physical Contrast Manipulation)
% If recompute is set to 1 - model analysis will be rerun even if a resuls file already exists
% IB - 2025

if ~exist('Task','var') || isempty(Task)
    p.Task      = 'CM'; % Default is 'CM'
else
    p.Task      = Task;
end

if ~exist('recompute','var') || isempty(recompute)
    p.recompute = 0; % Default is to not rerun the analyses
else
    p.recompute = recompute;
end

% Update based on requested amount of cores:
if ~exist('numCores','var') || isempty(numCores)

    num_cores   = 4; %maxNumCompThreads;
else
    num_cores   = numCores;
end
maxNumCompThreads(num_cores); % 

%% Setup directories
p.projectDir    = projectRootPath;
p.dataDir       = p.projectDir; % Folder 'Data' should reside inside this directory 
if ~exist(fullfile(p.dataDir, 'Data'), 'dir') > 0
    error('Data folder not found within the project directory')
end
p.codeDir     	= fileparts(which('projectRootPath'));
p.figureDir   	= fullfile(p.projectDir, 'Figures'); % Figures will be saved inside main project directory

p.saveDir     	= fullfile(p.dataDir, 'Data', 'modelOutput');
p.saveName      = sprintf('attWindow_smooth_TRmodelResults_%s.mat', p.Task);

if ~exist(p.figureDir, 'dir'), mkdir(p.figureDir); end

%% Get experimental parameters
p               = setupAnalysisInfo(p);

%% Preprocess data
showFigures     = false;
saveFig         = false;
Data            = preprocessData(p, showFigures, saveFig);

%% Extract task data
Design          = attWindow_taskData(p);

%% Fit models
showFigures     = false;
saveFig         = false;

out             = struct('avgResults', cell(1,5), ...
                        'avgSummary', cell(1,5));
count = 0;
for n = [1 2 3 5 10]
    p.nTR = n;
    
    count = count + 1;    
    [out(count).avgResults, out(count).avgSummary]  = analysisAttWindowTR(p, Data, Design, showFigures, saveFig);

    out(count).numTRs = n;
end

out(1).Data        = Data;
out(1).params      = p;
out(1).Design      = Design;

if ~exist(fullfile(p.saveDir, p.saveName), 'file') > 0 || p.recompute ~= 0
    save(fullfile(p.saveDir, p.saveName), 'out');
end
