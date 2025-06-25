function counts = figureStdEstimate(p, Data, results, type, figureDir, saveFig)

names           = fieldnames(results);


StdR2           = NaN(p.numROIs,p.nSubs,p.numAttWindows,length(p.histCentersStd));
StdMedians      = NaN(p.numROIs,p.nSubs,p.numAttWindows);
StdMedianExcl   = NaN(p.numROIs,p.nSubs,p.numAttWindows);

for sub=1:p.nSubs
    
    R2idx = results(sub).(names{contains(names, sprintf('%sR2idx', type))});
    whichWidth = results(sub).(names{contains(names, sprintf('%sWhichWidth', type))});
    
    for roi=1:p.numROIs
        
        outliers = Data(roi).FWHMEstimates{sub} > rad2deg(p.histCentersStd(end));
        
        for ii=1:p.numAttWindows
            StdR2(roi,sub,ii,:)         = histcounts(Data(roi).sdEstimates{sub}(whichWidth(:,ii) & R2idx), rad2deg(p.histEdgesStd));
            StdMedians(roi,sub,ii)      = nanmedian(Data(roi).sdEstimates{sub}(whichWidth(:,ii) & R2idx));
            StdMedianExcl(roi,sub,ii)   = nanmedian(Data(roi).sdEstimates{sub}(whichWidth(:,ii) & R2idx & ~outliers));
        end
        
    end
end

%% visualize 
figure('Color', [1 1 1], 'Position', [200 20 1000 800]);
clf
for roi = 1:p.numROIs
    subplot(2,3,roi)
    for ii = 1:p.numAttWindows %p.numAttWindows:-1:1
        hold on
        plot(rad2deg(p.histCentersStd), mean(squeeze(StdR2(roi,:,ii,:)),1), 'Color', p.CondColors(ii,:), 'LineWidth', 1.5)
%         histogram('BinCounts', mean(squeeze(FWHMR2(roi,:,ii,:)),1), 'BinEdges', rad2deg(p.histEdgesWidth), 'FaceColor', p.CondColorsCM(ii,:), 'EdgeColor', p.CondColorsCM(ii,:))
        if roi == 3
            counts(ii,:) = mean(squeeze(StdR2(roi,:,ii,:)),1);
        end
    end
    box off; title(p.ROInames{roi}); ylabel('avg count'); xlabel('Std (deg)')
    set(gca, 'XTick', 1:90:180, 'XTickLabel',{'0','90','180'})
    %ylim([0 7])
    
    subplot(2,3,roi+3)
    hold on, 
    for s = 1:p.nSubs
        plot(p.WidthAttWindow, squeeze(StdMedians(roi,s,:))', '.--', 'color',[.5 .5 .5],'MarkerSize', 12);
    end
    errorbar(p.WidthAttWindow, mean(squeeze(StdMedians(roi,:,:)),1), std(squeeze(StdMedians(roi,:,:)))/sqrt(p.nSubs), '.-k', 'MarkerSize', 20, 'CapSize', 0, 'LineWidth', 2);
    box off;
    set(gca, 'XTick', p.WidthAttWindow, 'XTickLabel', p.WidthAttWindow), xlim([0.5 9.5]), ylim([0 180])
    title(p.ROInames{roi})
    if roi == 1
        ylabel('Width Std (deg)');
        %legend(h,SubjNamesCM)
    end
    xlabel('Width condition')

end
sgtitle(sprintf('%s: %s trial Median Std estimates', p.savestr, type))


if saveFig
    if ~exist(figureDir, 'dir'), mkdir(figureDir); end
    print(fullfile(figureDir, sprintf('%s_%sTrials_%s', mfilename, type, p.savestr)), '-dpdf', '-painters', '-bestfit')

end


