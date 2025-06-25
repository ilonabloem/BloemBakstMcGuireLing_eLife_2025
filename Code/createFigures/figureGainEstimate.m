function figureGainEstimate(p, modelEstimates, results, type, figureDir, saveFig)

plotType        = 'violin'; % 'violin' or 'histogram'

names           = fieldnames(results);

GainEstimateR2  = NaN(p.numROIs,p.nSubs,p.numAttWindows,length(p.histCentersAmpl));
medGainR2       = NaN(p.numROIs,p.nSubs,p.numAttWindows);

varGain         = cell(p.numROIs, p.numAttWindows);


for sub=1:p.nSubs
    
    R2idx = results(sub).(names{contains(names, sprintf('%sR2idx', type))});
    whichWidth = results(sub).(names{contains(names, sprintf('%sWhichWidth', type))});
    
    for roi=1:p.numROIs
        

        for ii=1:p.numAttWindows
            
            % extract amplitude
            gain    = (modelEstimates(roi).amplEstimates{sub}(whichWidth(:,ii) & R2idx));

            GainEstimateR2(roi,sub,ii,:) = histcounts(gain, p.histEdgesAmpl);
            medGainR2(roi,sub,ii)        = nanmedian(gain);
            
            if sub == p.exampleSubj
                varGain{roi,ii} = cat(1, varGain{roi,ii}, gain);
            end
            
        end

    end
end


%% visualize results
figure('Color', [1 1 1], 'Position', [200 20 1000 800]);
clf
for roi = 1:p.numROIs
    switch plotType
        case 'histogram' 
            subplot(2,3,roi)
            for ii = p.numAttWindows:-1:1
                hold on
                histogram('BinCounts', squeeze(medGainR2(roi,p.exampleSubj,ii,:)), 'BinEdges', rad2deg(p.histEdgesAmpl), 'FaceColor', p.CondColors(ii,:), 'EdgeColor', p.CondColors(ii,:), 'FaceAlpha', 0.5)

            end
            box off; title(p.ROInames{roi}); ylabel('avg count'); xlabel('Gain (%PSC)')
            set(gca, 'XTick', p.histCentersAmpl([1 9 end]), 'XTickLabel', (p.histCentersAmpl([1 9 end])))
            %ylim([0 7])
        case 'violin'
            subplot(2,3,roi)
            for ii = 1:p.numAttWindows
                hold on,
                h               = Violin(varGain(roi,ii), p.WidthAttWindow(ii), 'ViolinColor', {p.CondColors(ii,:)}, 'HalfViolin','full', 'QuartileStyle','boxplot', 'DataStyle', 'none', 'Width', 0.7);
                h.BoxColor      = p.CondColors(ii,:);
            end
            box off
            set(gca, 'XTick', p.WidthAttWindow, 'XTickLabel', p.WidthAttWindow, 'YTick', [-1 0 4], 'YTickLabel', [-1 0 4])
            xlabel('Width condition'); ylabel('Gain (%SC)')
            ylim([-1 4])
            title(p.ROInames{roi})
    end
    subplot(2,3,roi+p.numROIs)
    hold on,
    for s = 1:p.nSubs
        if s == p.exampleSubj
            plot(p.WidthAttWindow,squeeze(medGainR2(roi,s,:))','^--','color',[.5 .5 .5],'MarkerSize',5, 'MarkerFaceColor', [.5 .5 .5])
        else
            plot(p.WidthAttWindow,squeeze(medGainR2(roi,s,:))','.--','color',[.5 .5 .5],'MarkerSize',12)
        end
    end
    
    errorbar(p.WidthAttWindow, mean(squeeze(medGainR2(roi,:,:)),1), std(squeeze(medGainR2(roi,:,:)))/sqrt(p.nSubs), ...
        'k.-', 'MarkerSize', 20, 'CapSize', 0, 'LineWidth', 2);
    box off;
    yl = ylim; yl = [floor(yl(1)/0.5)*0.5  ceil((max(squeeze(medGainR2(roi,:,:)), [], 'all')+0.1)/0.5)*0.5];
    xlim([0.5 9.5]), ylim(yl);
    set(gca, 'XTick', p.WidthAttWindow, 'XTickLabel', p.WidthAttWindow), xlim([0.5 9.5])
    if roi == 3
        ylabel('Gain (%SC)');
        legend(p.ROInames(:))
        xlabel('Width condition')
    end
    

end
sgtitle(sprintf('%s: %s trial Median Gain estimates', p.savestr, type))

if saveFig
    if ~exist(figureDir, 'dir'), mkdir(figureDir); end
    print(fullfile(figureDir, sprintf('%s_%sTrials_%s_%s', mfilename, type, p.savestr, plotType)), '-dpdf', '-painters', '-bestfit')

end

%% Save csv files for stats
for roi = 1:p.numROIs
    
    Gain_width = squeeze(medGainR2(roi,:,:));
    
    T = array2table(Gain_width);
    T.Properties.RowNames = p.SubjNames;
    T.Properties.Description = 'Model results: Gain of spatial profile fits';
    T.Properties.VariableUnits = repmat({'%SC'}, [1, p.numAttWindows]); 
    writetable(T, fullfile(figureDir, sprintf('%s_Gain_%s.csv', p.savestr, p.ROInames{roi})))
    
end


