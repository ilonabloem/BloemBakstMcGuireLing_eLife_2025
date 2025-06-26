function figureBaselineEstimate(p, modelEstimates, results, type, figureDir, saveFig)

plotType        = 'violin'; % 'violin' or 'histogram'

names           = fieldnames(results);

BaseEstimateR2  = NaN(p.numROIs,p.nSubs,p.numAttWindows,length(p.histCentersBase));
medBaseR2       = NaN(p.nSubs,p.numAttWindows,p.numROIs);

varBase         = cell(p.numROIs, p.numAttWindows);

for sub=1:p.nSubs
    
    R2idx = results(sub).(names{contains(names, sprintf('%sR2idx', type))});
    whichWidth = results(sub).(names{contains(names, sprintf('%sWhichWidth', type))});
    
    for roi=1:p.numROIs
        
        for ii=1:p.numAttWindows
            
            % extract amplitude
            base    = (modelEstimates(roi).baseEstimates{sub}(whichWidth(:,ii) & R2idx));
            
            BaseEstimateR2(roi,sub,ii,:) = histcounts(base, p.histEdgesBase);
            medBaseR2(sub,ii,roi) = nanmedian((base));
            
            if sub == p.exampleSubj
                varBase{roi,ii} = cat(1, varBase{roi,ii}, base);
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
                histogram('BinCounts', (squeeze(BaseEstimateR2(roi,p.exampleSubj,ii,:))), ...
                    'BinEdges', p.histEdgesBase, 'FaceColor', p.CondColors(ii,:), ...
                    'EdgeColor', p.CondColors(ii,:), 'FaceAlpha', 0.5)
            end
            box off; title(p.ROInames{roi}); ylabel('avg count'); xlabel('Est baseline')
            set(gca, 'XTick', p.histEdgesBase(1):1:p.histEdgesBase(end))
            if strcmpi(type, 'indv')
                if strcmpi(p.Task, 'CM')
                    ylim([0 200])
                else
                    ylim([0 40])
                end
            else
                if strcmpi(p.Task, 'CM')
                    ylim([0 30])
                else
                    ylim([0 5])
                end
            end
        case 'violin'
            subplot(2,3,roi)
            plot([0 10], [0 0], 'k')
            for ii = 1:p.numAttWindows
                hold on,
                h               = Violin(varBase(roi,ii), p.WidthAttWindow(ii), 'ViolinColor',{p.CondColors(ii,:)}, 'HalfViolin','full', 'QuartileStyle','boxplot', 'DataStyle', 'none', 'Width', 0.7);
                h.BoxColor      = p.CondColors(ii,:);
            end
            box off
            set(gca, 'XTick', p.WidthAttWindow, 'XTickLabel', p.WidthAttWindow, 'YTick', -5:1, 'YTickLabel', -5:1)
            xlabel('Width condition'); ylabel('Baseline (%SC)')
            title(p.ROInames{roi})
            ylim([-5 1])
    end
    
    subplot(2,3,roi+p.numROIs)
    hold on
    plot([0 9.5], [0 0], 'k')
    for s = 1:p.nSubs
        if s == p.exampleSubj
            plot(p.WidthAttWindow,medBaseR2(s,:,roi),'^--','color',[.5 .5 .5],'MarkerSize',5, 'MarkerFaceColor', [.5 .5 .5])
        else
            plot(p.WidthAttWindow,medBaseR2(s,:,roi),'.--','color',[.5 .5 .5],'MarkerSize',12)
        end

    end
    errorbar(p.WidthAttWindow,mean(medBaseR2(:,:,roi)),std(medBaseR2(:,:,roi))/sqrt(p.nSubs),'.-k','MarkerSize',20,'CapSize',0,'LineWidth',2)
    set(gca, 'XTick', p.WidthAttWindow, 'XTickLabel', p.WidthAttWindow, 'YTick', -0.75:0.25:0.25, 'YTickLabel', -0.75:0.25:0.25), 
    yl = ylim; yl = [-0.75  0.25];
    xlim([0.5 9.5]), ylim(yl);    title(p.ROInames{roi})
    if roi==1
        ylabel('baseline (%SC)')
    end
    xlabel('Width condition')

end

sgtitle(sprintf('%s: %s trial Median Baseline estimates', p.savestr, type))

if saveFig
    if ~exist(figureDir, 'dir'), mkdir(figureDir); end
    print(fullfile(figureDir, sprintf('%s_%sTrials_%s_%s', mfilename, type, p.savestr, plotType)), '-dpdf', '-painters', '-bestfit')
end

%% Save csv files for stats
for roi = 1:p.numROIs
    
    Base_width = squeeze(medBaseR2(:,:,roi));
    
    T = array2table(Base_width);
    T.Properties.RowNames = p.SubjNames;
    T.Properties.Description = 'Model results: Baseline of spatial profile fits';
    T.Properties.VariableUnits = repmat({'%SC'}, [1, p.numAttWindows]); 
    writetable(T, fullfile(p.csvDir, sprintf('%s_Base_%s.csv', p.savestr, p.ROInames{roi})))
    
end

