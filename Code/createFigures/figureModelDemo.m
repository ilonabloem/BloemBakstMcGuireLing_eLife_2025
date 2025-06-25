

function figureModelDemo(p, Data, results, type, figureDir, saveFig)


names               = fieldnames(results);

%% Extract data Attention
avgDataTM_R2        = NaN(p.numAttWindows, p.numROIs, p.nSubs, p.nBins); %trueMean

for sub = 1:p.nSubs
    
    R2idx = results(sub).(names{contains(names, sprintf('%sR2idx', type))});
    whichWidth = results(sub).(names{contains(names, sprintf('%sWhichWidth', type))});
    
    
    for roi = 1:p.numROIs
%         FWHMoutliers = results(sub).(names{contains(names, 'FWHMoutliers')})(roi,:);
        
        for ii = 1:p.numAttWindows
            
            % Recentered data profiles
            avgDataTM_R2(ii,roi,sub,:)  = nanmean(Data(roi).TrueMeanRecenteredData{sub}(:, whichWidth(:,ii) & R2idx),2);

        end
    end
end



%% Example spatial profiles V3 narrow and wide cue conditions

avgResponseWide     = squeeze(avgDataTM_R2(end,3,4,:));
avgResponseNarrow   = squeeze(avgDataTM_R2(2,3,4,:));
% shift curve so that peak is not at pi
shiftIndx           = cat(2, 16:60, 1:15);
avgResponseNarrow   = avgResponseNarrow(shiftIndx);

%et example fits
genGVals            = fitGenGauss([], {'initialize'});

%-- Wide
fixedParams         = cell(3,1);
fixedParams{1}      = 'optimize';
fixedParams{2}      = linspace(-pi, 3*pi - (4*pi/(p.nBins*2)), p.nBins*2)';
fixedParams{3}      = avgResponseWide(:);

paramsWide          = ...
                    fmincon(@(x) fitGenGauss(x, fixedParams), genGVals.init, [], [], [], [], genGVals.lb, genGVals.ub, [],  genGVals.opt);
                
fixedParams{1}      = 'prediction';
wideResults         = fitGenGauss(paramsWide, fixedParams);

%-- Narrow
fixedParams         = cell(3,1);
fixedParams{1}      = 'optimize';
fixedParams{2}      = linspace(-pi, 3*pi - (4*pi/(p.nBins*2)), p.nBins*2)';
fixedParams{3}      = avgResponseNarrow(:);

paramsNarrow        = ...
                    fmincon(@(x) fitGenGauss(x, fixedParams), genGVals.init, [], [], [], [], genGVals.lb, genGVals.ub, [],  genGVals.opt);
                
fixedParams{1}      = 'prediction';
narrowResults       = fitGenGauss(paramsNarrow, fixedParams);

%% Create various gen gauss fits for demo purposes

% param order: 'mu', 'SD', 'Beta', 'ampl', 'offset'
exampleParams       = [pi 2 log10(2) 1 0; ... 
                       pi/2 1 log10(2) 1 0; ...  
                       1.3*pi 1 log10(8) 1 0; ...  
                       ];
                  
fixedParams         = cell(3,1);
fixedParams{1}      = 'prediction';
fixedParams{2}      = linspace(-pi, 3*pi - (4*pi/(p.nBins*2)), p.nBins*2)';

demoGenGauss        = NaN(size(exampleParams,1), p.nBins);
for ii = 1:size(exampleParams,1)
    fits                = fitGenGauss(exampleParams(ii,:), fixedParams);
    demoGenGauss(ii,:)  = fits.y_est;
end

%% Visualize

figure('Color', [1 1 1], 'Position', [100 100 1000 1000])
subplot(2,2,1)
plot(ones(2,1) * exampleParams(:,1)', repmat([0; 1], [1 size(exampleParams,1)]), 'k')
hold on, 
plot(p.binCenters, demoGenGauss)
box off, xlim([0 2*pi])
set(gca, 'XTick', [0 pi 2*pi], 'XTickLabel', [-180 0 180], 'YTick', [], 'TickDir', 'out')
xlabel('Polar angle (deg)')
ylabel('Response (a.u.)')
title([{sprintf('GG1: mu:%.2f, std:%.2f, beta:%.2f', exampleParams(1,1), exampleParams(1,2), 10.^exampleParams(1,3))}; ...
    {sprintf('GG2: mu:%.2f, std:%.2f, beta:%.2f', exampleParams(2,1), exampleParams(2,2), 10.^exampleParams(2,3))}; ...
    {sprintf('GG3: mu:%.2f, std:%.2f, beta:%.2f', exampleParams(3,1), exampleParams(3,2), 10.^exampleParams(3,3))}])

subplot(2,2,3)
plot(p.binCenters, avgResponseWide, 'k.', 'MarkerSize', 30)
hold on, 
plot(p.binCenters, wideResults.y_est, 'r', 'LineWidth', 3)
box off, xlim([0 2*pi])
set(gca, 'XTick', [0 pi 2*pi], 'XTickLabel', [-180 0 180], 'YTick', [], 'TickDir', 'out')
xlabel('Polar angle (deg)')
ylabel('Response (a.u.)'), ylim([-0.5 0.8])
title(sprintf('GG: mu:%.2f, std:%.2f, beta:%.2f', paramsWide(1,1), paramsWide(1,2), 10.^paramsWide(1,3)))


subplot(2,2,4)
plot(p.binCenters, avgResponseNarrow, 'k.', 'MarkerSize', 30)
hold on, 
plot(p.binCenters, narrowResults.y_est, 'r', 'LineWidth', 3)
box off, xlim([0 2*pi])
set(gca, 'XTick', [0 pi 2*pi], 'XTickLabel', [-180 0 180], 'YTick', [], 'TickDir', 'out')
xlabel('Polar angle (deg)')
ylabel('Response (a.u.)'), ylim([-0.5 0.8])
title(sprintf('GG: mu:%.2f, std:%.2f, beta:%.2f', (paramsNarrow(1,1)), paramsNarrow(1,2), 10.^paramsNarrow(1,3)))


if saveFig
    if ~exist(figureDir, 'dir'), mkdir(figureDir); end
    print(fullfile(figureDir, sprintf('%s_%sTrials_%s', mfilename, type, p.savestr)), '-dpdf', '-painters', '-bestfit')
end

