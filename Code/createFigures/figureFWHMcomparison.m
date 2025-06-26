
function figureFWHMcomparison(csvDir, figureDir, saveFig)

if ~exist('csvDir','var') || isempty(csvDir)
    error('Missing path to csv folder')
end
if ~exist('figureDir','var') || isempty(figureDir)
    figureDir       = fullfile(projectRootPath, 'Figures'); % find project root dir based on where this code is located 
end
if ~exist('saveFig', 'var') || isempty(saveFig)
    saveFig         = true; 
end


%-- setup paths
projectRootPath;
conditions  = {'attention', 'perception'};
ROIs        = {'V1', 'V2', 'V3'};
attWidths   = {[1 3 5 9]; [1 3 5 7 9]};
fwhmAtt     = table(attWidths{1}(:), 'VariableNames', {'attWidths'});
fwhmPerc    = table(attWidths{2}(:), 'VariableNames', {'attWidths'});

%-- load the FWHM csv files for both attention and perception
csvList     = dir(fullfile(csvDir, 'CM', 'attention_FWHM_*.csv'));
csvList     = cat(1, csvList, dir(fullfile(csvDir, 'PCM', 'perception_FWHM_V*.csv')));



if numel(csvList) ~= 6 % 3 ROIs x 2 conditions

    error('attention_FWHM_V*.csv files not found')
else

    for ii = 1:numel(conditions) % attention and perception

        for jj = 1:numel(ROIs) % 3 ROIs
            
            if strcmp(conditions{ii}, 'attention')

                fwhmData = readtable(fullfile(csvDir, 'CM', sprintf('attention_FWHM_%s.csv', ROIs{jj})));
                
                fwhmAtt.(sprintf('FWHM_V%d', jj)) = table2array(mean(fwhmData))';

            elseif strcmp(conditions{ii}, 'perception')
                
                fwhmData = readtable(fullfile(csvDir, 'PCM', sprintf('perception_FWHM_%s.csv', ROIs{jj})));
                fwhmPerc.(sprintf('FWHM_V%d', jj)) = table2array(mean(fwhmData))';

            end
        end


    end

end

%% Save csvs
writetable(fwhmAtt, fullfile(csvDir, 'CM', 'attention_avgFWHM.csv'))
writetable(fwhmPerc, fullfile(csvDir, 'PCM', 'perception_avgFWHM.csv'))

%% Visualize
figure('Color', [1 1 1], 'Position', [50 200 800 500])
dotColors = [0 0.4 0.6];

plot([50 200], [50 200], 'k', 'HandleVisibility','off')
for ii = 1:numel(ROIs)

    hold on, 
    scatter(fwhmAtt{:,sprintf('FWHM_V%d', ii)}, ...
        fwhmPerc{ismember(attWidths{2}, attWidths{1}),sprintf('FWHM_V%d', ii)}, ...
        100, dotColors(ii)*ones(4,3), 'filled')

end
box off; axis square
xlabel('FWHM Attention field', 'FontSize', 14)
ylabel('FWHM Perception field', 'FontSize', 14)
legend(ROIs, 'Location', 'northwest')
set(gca, 'TickDir', 'out', 'XColor', 'k', 'ycolor', 'k')
title('Figure 7c')

if saveFig
    if ~exist(figureDir, 'dir'), mkdir(figureDir); end
    print(fullfile(figureDir, sprintf('fig7c_FWHMcomparison')), '-dpdf', '-vector','-bestfit')
end

end
