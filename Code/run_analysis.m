function out = run_analysis(Task, recompute, numCores)

%% Run attentional window analysis
%
% function out = s_runAnalysis(task, recompute, numCores)
%
% Task can be: 'CM' (Continuous Monitoring) or 'PCM' (Physical Contrast Manipulation)
% If recompute is set to 1 - model analysis will be overwritten even if a results file already exists  
% numCores constraints the available cores in matlab 
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

    num_cores   = 2; %maxNumCompThreads;
else
    num_cores   = numCores;
end
maxNumCompThreads(num_cores); % 

%% Setup directories
p.projectDir  	= projectRootPath;
% p.dataDir       = p.projectDir; % Folder 'Data' should reside inside data directory 
p.dataDir       = fullfile('~', 'Documents', 'BloemBakstMcGuireLing_eLife_2025');
p.codeDir     	= fileparts(which('projectRootPath'));
p.figureDir   	= fullfile(p.projectDir, 'Figures'); % Figures will be saved inside main project directory

p.saveDir     	= fullfile(p.dataDir, 'Data', 'modelOutput');
p.saveName      = sprintf('attWindow_smooth_modelResults_%s.mat', p.Task);


%% Get experimental parameters
p               = setupAnalysisInfo(p);

%% Preprocess data
showFigures     = false;
saveFig         = false;
Data            = preprocessData(p, showFigures, saveFig);

%% Extract task data
Design          = attWindow_taskData(p);

%% Create 2D representations
showFigures     = false;
saveFig         = false;

p.analysisType  = 'blockAvg';
attWindow_2D_visualization(p, Data, Design, showFigures, saveFig)

%% Fit model
p.analysisType  = 'blockAvg';
[out.avgResults, out.avgSummary]  = analysisAttWindow(p, Data, Design);

out.Data        = Data;
out.params      = p;
out.Design      = Design;

if ~exist(fullfile(p.saveDir, p.saveName), 'file') > 0 || p.recompute > 0
    save(fullfile(p.saveDir, p.saveName), 'out');
end

