function figureFWHMEstimate(p, Data, results, type, figureDir, saveFig)

plotType        = 'violin'; % 'violin' or 'histogram'

names           = fieldnames(results);


FWHMR2          = NaN(p.numROIs,p.nSubs,p.numAttWindows,length(p.histCentersWidth));
WidthMedians    = NaN(p.numROIs,p.nSubs,p.numAttWindows);
WidthMedianExcl = NaN(p.numROIs,p.nSubs,p.numAttWindows);

varFWHM        = cell(p.numROIs, p.numAttWindows);

for sub=1:p.nSubs
    
    R2idx = results(sub).(names{contains(names, sprintf('%sR2idx', type))});
    whichWidth = results(sub).(names{contains(names, sprintf('%sWhichWidth', type))});
    
    for roi=1:p.numROIs
        
        outliers = Data(roi).FWHMEstimates{sub} > rad2deg(p.histCentersWidth(end));
        
        for ii=1:p.numAttWindows
            FWHMR2(roi,sub,ii,:)        = histcounts(Data(roi).FWHMEstimates{sub}(whichWidth(:,ii) & R2idx), rad2deg(p.histEdgesWidth));
            WidthMedians(roi,sub,ii)    = nanmedian(Data(roi).FWHMEstimates{sub}(whichWidth(:,ii) & R2idx));
            WidthMedianExcl(roi,sub,ii) = nanmedian(Data(roi).FWHMEstimates{sub}(whichWidth(:,ii) & R2idx & ~outliers));
            
            if sub == p.exampleSubj
                varFWHM{roi,ii} = cat(1, varFWHM{roi,ii}, Data(roi).FWHMEstimates{sub}(whichWidth(:,ii) & R2idx));
            end
        end
        
    end
end

%% visualize 
figure('Color', [1 1 1], 'Position', [200 20 1000 800]);
clf
for roi = 1:p.numROIs
    switch plotType
        case 'histogram' 
            subplot(2,3,roi)
            for ii = p.numAttWindows:-1:1
                hold on
%                 plot(rad2deg(p.histCentersWidth), (squeeze(FWHMR2(roi,:,ii,:)), 'Color', p.CondColors(ii,:), 'LineWidth', 1.5)
                histogram('BinCounts', squeeze(FWHMR2(roi,p.exampleSubj,ii,:)), 'BinEdges', rad2deg(p.histEdgesWidth), 'FaceColor', p.CondColors(ii,:), 'EdgeColor', p.CondColors(ii,:), 'FaceAlpha', 0.5)

            end
            box off; title(p.ROInames{roi}); ylabel('avg count'); xlabel('FWHM (rad)')
            set(gca, 'XTick', 1:180:360, 'XTickLabel',{'0','180','360'})
            %ylim([0 7])
        case 'violin'
            subplot(2,3,roi)
            for ii = 1:p.numAttWindows
                hold on,
                h               = Violin(varFWHM(roi,ii), p.WidthAttWindow(ii), 'ViolinColor', {p.CondColors(ii,:)}, 'HalfViolin','full', 'QuartileStyle','boxplot', 'DataStyle', 'none', 'Width', 0.7);
                h.BoxColor      = p.CondColors(ii,:);
            end
            box off
            set(gca, 'XTick', p.WidthAttWindow, 'XTickLabel', p.WidthAttWindow, 'YTick', [0 90 180 270 360], 'YTickLabel', [0 90 180 270 360])
            xlabel('Width condition'); ylabel('FWHM (deg)')
            ylim([0 360])
            title(p.ROInames{roi})
    end
    subplot(2,3,roi+3)
    hold on,
    for s = 1:p.nSubs
        if s == p.exampleSubj
            plot(p.WidthAttWindow,squeeze(WidthMedians(roi,s,:))','^--','color',[.5 .5 .5],'MarkerSize',5, 'MarkerFaceColor', [.5 .5 .5])
        else
            plot(p.WidthAttWindow,squeeze(WidthMedians(roi,s,:))','.--','color',[.5 .5 .5],'MarkerSize',12)
        end
    end
    
    errorbar(p.WidthAttWindow, mean(squeeze(WidthMedians(roi,:,:)),1), std(squeeze(WidthMedians(roi,:,:)))/sqrt(p.nSubs), '.-k', 'MarkerSize', 20, 'CapSize', 0, 'LineWidth', 2);
    box off;
    set(gca, 'XTick', p.WidthAttWindow, 'XTickLabel', p.WidthAttWindow, 'YTick', [50 210]), 
    xlim([0.5 9.5]), ylim([50 210])
    title(p.ROInames{roi})
    if roi == 1
        ylabel('Width FWHM (deg)');
        %legend(h,SubjNamesCM)
    end
    xlabel('Width condition')

end
sgtitle(sprintf('%s: %s trial Median FWHM estimates', p.savestr, type))


if saveFig
    if ~exist(figureDir, 'dir'), mkdir(figureDir); end
    print(fullfile(figureDir, sprintf('%s_%sTrials_%s_%s', mfilename, type, p.savestr, plotType)), '-dpdf', '-painters', '-bestfit')

end

%-- save table with mean widths
%% Save csv files for stats
for roi = 1:p.numROIs
    
    FWHM_width = squeeze(WidthMedians(roi,:,:));
    
    T = array2table(FWHM_width);
    T.Properties.RowNames = p.SubjNames;
    T.Properties.Description = 'Model results: FWHM of spatial profile fits';
    T.Properties.VariableUnits = repmat({'deg'}, [1, p.numAttWindows]);
    writetable(T, fullfile(p.csvDir, sprintf('%s_FWHM_%s.csv', p.savestr, p.ROInames{roi})))
    
end


