%% Attentional window figures
function run_createFiguresTR(Task, figureDir, dataDir, saveFig)

if ~exist('Task','var') || isempty(Task)
    Task            = 'CM'; 
end
if ~exist('figureDir','var') || isempty(figureDir)
    figureDir       = fullfile(projectRootPath, 'Figures', 'variableTR'); % find project root dir based on where this code is located 
end
if ~exist('dataDir', 'var') || isempty(dataDir)
    dataDir         = projectRootPath;
    if ~exist(fullfile(dataDir, 'Data'), 'dir') > 0
        error('Data folder not found within the project directory')
    end
end
if ~exist('saveFig', 'var') || isempty(saveFig)
    saveFig         = true; 
end
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end
% Add paths
projectRootPath;

% make sure necessary toolbox exist
if ~exist('Violin', 'file') > 0
    error('Toolbox ''Violinplot-Matlab'' is necessary to reproduce the figures')
end

%% Load data and model results
resultsName         = fullfile(dataDir, 'Data', 'modelOutput', sprintf('attWindow_smooth_TRmodelResults_%s.mat', Task));

if exist(resultsName, 'file') > 0 
    load(resultsName, 'out');
    data            = out;
else
    error('Analysis results do not exist ..')
end


%% setup how to bin parameter estimates
p                   = data.params;

p.histEdgesAng      = linspace(-pi,pi,21); % 9 degrees bins 
p.histCentersAng    = (p.histEdgesAng(1:end-1)+p.histEdgesAng(2:end))/2;

p.histEdgesWidth    = linspace(0, 2*pi,21); % Spaces in 4.5 degrees bins - matches distance between stimuli 
p.histCentersWidth  = (p.histEdgesWidth(1:end-1)+p.histEdgesWidth(2:end))/2;

p.histEdgesResp     = linspace(-1,4,21); % 
p.histCentersResp   = (p.histEdgesResp(1:end-1)+p.histEdgesResp(2:end))/2;

p.histEdgesAmpl     = linspace(0,10,21); % 
p.histCentersAmpl   = (p.histEdgesAmpl(1:end-1)+p.histEdgesAmpl(2:end))/2;

p.histEdgesBase     = linspace(-5,5,21);
p.histCentersBase   = (p.histEdgesBase(1:end-1)+p.histEdgesBase(2:end))/2;

p.histEdgesBeta     = linspace(1.8,50,21);
p.histCentersBeta   = (p.histEdgesBeta(1:end-1)+p.histEdgesBeta(2:end))/2;

p.histEdgesStd      = linspace(deg2rad(6),pi,21);
p.histCentersStd    = (p.histEdgesStd(1:end-1)+p.histEdgesStd(2:end))/2;

p.R2_cutoff         = 0.8; % proportion of R2 to use


customColor         = [ones(1,256/2) linspace(1,0,256/2); ...
                       linspace(1,0,256); ...
                       0.5 * ones(1,256)]';

Colors              = customColor([65 90 120 150 210],:);
p.CondColors        = Colors([1:3 5],:);
exmpleSubj          = '019';


p.exampleSubj       = find(strcmp(p.SubjNames, exmpleSubj)); %{'022','012','021','003','025'}

%% Preallocate variable for results:
numTRs              = numel(data);  
p.numTRs            = numTRs; 
results             = struct('pooledR2', cell(p.nSubs, numTRs), ...
                     'R2idx', cell(p.nSubs, numTRs), ...
                     'WhichWidth', cell(p.nSubs, numTRs), ...
                     'tmp', cell(p.nSubs, numTRs));

%% R2 summary 
numR2trials         = NaN(numTRs,p.nSubs,p.numAttWindows);
propR2trials        = NaN(numTRs,p.nSubs,p.numAttWindows);

% figure('Color', [1 1 1])
for tr = 1:numel(data)
    for sub = 1:p.nSubs

        whichWidth  = false(size(data(tr).avgResults(1).widthCondition{sub},1),p.numAttWindows);
        medR2       = NaN(p.numAttWindows, p.numROIs);

        % Sum R2 of model fit for max num TRs over all ROIs to find which timepoints to exclude
        allR2       = [data(tr).avgResults(1).R2{sub}, ...
                           data(tr).avgResults(2).R2{sub}, ...
                           data(tr).avgResults(3).R2{sub}];

        pooledR2CM  = sum(allR2,2);

        % Make sure not to use blank periods
        pooledR2CM(isnan(data(end).avgResults(1).widthCondition{sub})) = 0;

        %{
        subplot(2,4,sub)
        scatter(pooledR2CM, data.indvResults(3).R2{sub})
        xlabel('Pooled R2'), ylabel('V3 R2')
        title(p.SubjNames{sub}), xlim([0 3]), ylim([0 1]), axis square
        %}

        % Sort based on largest value - best fit
        [~,R2sort]  = sort(pooledR2CM, 'descend');
        cutoffIdx   = R2sort(1:sum(~isnan(data(tr).avgResults(1).widthCondition{sub}))*p.R2_cutoff);
        R2idx       = false(size(data(tr).avgResults(1).R2{sub}));
        R2idx(cutoffIdx) = true;

        % Loop through the different width conditions
        for ii = 1:p.numAttWindows
            whichWidth(:,ii) = data(tr).avgResults(1).widthCondition{sub} == p.WidthAttWindow(ii);
            propR2trials(tr,sub,ii) = sum(whichWidth(:,ii) & R2idx) / sum(whichWidth(:,ii));
            numR2trials(tr,sub,ii) = sum(whichWidth(:,ii) & R2idx);

            medR2(ii,:)  = median(allR2(whichWidth(:,ii) & R2idx, :));
        end

        results(sub, tr).pooledR2   = pooledR2CM;
        results(sub, tr).R2idx      = R2idx;
        results(sub, tr).WhichWidth = whichWidth;
        results(sub, tr).medR2      = medR2;

    end
end
% sgtitle('Individual trial R2 summary')

%% Variable TR ---
figureVariableTR(p, data, results, fullfile(figureDir), saveFig);




