%% Attentional window figures
function run_createFigures(Task, figureDir, dataDir, saveFig)

if ~exist('Task','var') || isempty(Task)
    Task            = 'CM'; % 'CM' or 'PCM'. Default is 'CM' 
end
if ~exist('figureDir','var') || isempty(figureDir)
    figureDir       = fullfile(projectRootPath, 'Figures'); % find project root dir based on where this code is located 
end
if ~exist('dataDir', 'var') || isempty(dataDir)
    dataDir         = fullfile(projectRootPath);
end
if ~exist('saveFig', 'var') || isempty(saveFig)
    saveFig         = true; 
end

% Add paths
projectRootPath;

% make sure necessary toolbox exist
if ~exist('Violin', 'file') > 0
    error('Toolbox ''Violinplot-Matlab'' is necessary to reproduce the figures')
end

%% Load data and model results

resultsName         = fullfile(dataDir, 'Data', 'modelOutput', sprintf('attWindow_smooth_modelResults_%s.mat', Task));

if exist(resultsName, 'file') > 0 
    load(resultsName, 'out');
    data            = out;
else
    error('Analysis results do not exist ..')
end


%% setup how to bin parameter estimates
p               = data.params;
p               = setupFigureInfo(p);

%% Preallocate variable for results:
results         = struct('avgpooledR2', cell(1,p.nSubs), ...
                         'avgR2idx', cell(1,p.nSubs), ...
                         'avgWhichWidth', cell(1,p.nSubs), ...
                         'tmp', cell(1,p.nSubs));

%% R2 summary 
propR2blocks    = NaN(p.nSubs,p.numAttWindows);
numR2blocks     = NaN(p.nSubs,p.numAttWindows);

% figure('Color', [1 1 1])
for sub = 1:p.nSubs
    
    whichWidth  = false(size(data.avgResults(1).widthCondition{sub},1),p.numAttWindows);
    medR2       = NaN(p.numAttWindows, p.numROIs);
    
    % Sum R2 of model fit for max num TRs over all ROIs to find which timepoints to exclude
    allR2       = [data(end).avgResults(1).R2{sub}, ...
                   data(end).avgResults(2).R2{sub}, ...
                   data(end).avgResults(3).R2{sub}];
                   
    pooledR2CM  = sum(allR2,2);  
    
    % Make sure to exclude blank periods
    pooledR2CM(isnan(data.avgResults(1).widthCondition{sub})) = 0;
    
    %{
    subplot(2,4,sub)
    scatter(pooledR2CM, data.avgResults(3).R2{sub})
    xlabel('Pooled R2'), ylabel('V3 R2')
    title(p.SubjNames{sub}), xlim([0 3]), ylim([0 1]), axis square
    %}
    
    % Sort based on largest value - best fit
    [~,R2sort]  = sort(pooledR2CM, 'descend');
    cutoffIdx   = R2sort(1:sum(~isnan(data.avgResults(1).widthCondition{sub}))*p.R2_cutoff);
    R2idx       = false(size(data.avgResults(1).R2{sub}));
    R2idx(cutoffIdx) = true;
    
    % Loop through the different width conditions
    for ii = 1:p.numAttWindows
        whichWidth(:,ii) = data.avgResults(1).widthCondition{sub} == p.WidthAttWindow(ii);
        propR2blocks(sub,ii) = sum(whichWidth(:,ii) & R2idx) / sum(whichWidth(:,ii));
        numR2blocks(sub,ii) = sum(whichWidth(:,ii) & R2idx);
        
        medR2(ii,:)  = median(allR2(whichWidth(:,ii) & R2idx, :));
    end
    
    results(sub).avgpooledR2   = pooledR2CM;
    results(sub).avgR2idx      = R2idx;
    results(sub).avgWhichWidth = whichWidth;
    results(sub).medR2         = medR2;
end

%% Compute avg R2 across ROI and width condition
if ~exist(fullfile(figureDir, Task),'dir'), mkdir(fullfile(figureDir, Task)); end

if ~exist(fullfile(figureDir, Task, sprintf('%s_SEMR2.csv', p.savestr)), 'file') > 0
    tmp         = [results(:).medR2];
    tmp         = reshape(tmp, [p.numAttWindows, p.numROIs, p.nSubs]);

    %-- Mean
    avgR2       = mean(tmp,3);

    attWidths   = p.WidthAttWindow(:);
    avgR2_V1    = avgR2(:,1);
    avgR2_V2    = avgR2(:,2);
    avgR2_V3    = avgR2(:,3);
    T           = table(attWidths, avgR2_V1, avgR2_V2, avgR2_V3);
    writetable(T, fullfile(figureDir, Task, sprintf('%s_avgR2.csv', p.savestr)))

    %-- SEM
    SEMR2       = std(tmp, [], 3) / sqrt(p.nSubs);

    SEMR2_V1    = SEMR2(:,1);
    SEMR2_V2    = SEMR2(:,2);
    SEMR2_V3    = SEMR2(:,3);
    T           = table(attWidths, SEMR2_V1, SEMR2_V2, SEMR2_V3);
    writetable(T, fullfile(figureDir, Task, sprintf('%s_SEMR2.csv', p.savestr)))
end

%% Activity profiles
figureActivityProfilesFWHM(p, data.avgSummary, results, 'avg', fullfile(figureDir, Task), saveFig);

%% Angular error ---
figureAngularError(p, data.avgResults, results, 'avg', fullfile(figureDir, Task, 'paramEst'), saveFig);

%% FWHM ---
figureFWHMEstimate(p, data.avgResults, results, 'avg', fullfile(figureDir, Task, 'paramEst'), saveFig);

%% Gain --- 
figureGainEstimate(p, data.avgResults, results, 'avg', fullfile(figureDir, Task, 'paramEst'), saveFig);

%% Amplitude --- 
figureAmplitudeEstimate(p, data.avgResults, results, 'avg', fullfile(figureDir, Task, 'paramEst'), saveFig);

%% Baseline ---
figureBaselineEstimate(p, data.avgResults, results, 'avg', fullfile(figureDir, Task, 'paramEst'), saveFig);


%% Miscellaneous figures
% Model schematic
if strcmp(Task, 'CM')
    
    figureModelDemo(p, data.avgSummary, results, 'avg', fullfile(figureDir, Task), saveFig);

end

