function figureAngularError(p, modelEstimates, results, type, figureDir, saveFig)

plotType        = 'violin'; % 'violin' or 'histogram'

names           = fieldnames(results);

AngularErrorR2  = NaN(p.numROIs,p.nSubs,p.numAttWindows,length(p.histCentersAng));
medErrR2        = NaN(p.nSubs,p.numAttWindows,p.numROIs);

varError        = cell(p.numROIs, p.numAttWindows);
for sub=1:p.nSubs
    
    R2idx = results(sub).(names{contains(names, sprintf('%sR2idx', type))});
    whichWidth = results(sub).(names{contains(names, sprintf('%sWhichWidth', type))});
    
    for roi=1:p.numROIs

        for ii=1:p.numAttWindows
            
            AngularErrorR2(roi,sub,ii,:) = histcounts(modelEstimates(roi).diffAng{sub}(whichWidth(:,ii) & R2idx), p.histEdgesAng);
            medErrR2(sub,ii,roi) = nanmedian(abs(modelEstimates(roi).diffAng{sub}(whichWidth(:,ii) & R2idx)));
            
            if sub == p.exampleSubj
                varError{roi,ii} = cat(1, varError{roi,ii}, modelEstimates(roi).diffAng{sub}(whichWidth(:,ii) & R2idx));
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
                histogram('BinCounts', (squeeze(AngularErrorR2(roi,p.exampleSubj,ii,:))), ...
                    'BinEdges', p.histEdgesAng, 'FaceColor', p.CondColors(ii,:), ...
                    'EdgeColor', p.CondColors(ii,:), 'FaceAlpha', 0.5)
            end
            box off; title(p.ROInames{roi}); ylabel('avg count'); xlabel('Est angle error')
            set(gca, 'XTick', [-pi 0 pi], 'XTickLabel', {'-180', '0', '180'})
            if strcmpi(type, 'indv')
                if strcmpi(p.Task, 'CM')
                    ylim([0 90])
                else
                    ylim([0 10])
                end
            else
                if strcmpi(p.Task, 'CM')
                    ylim([0 20])
                else
                    ylim([0 3])
                end
            end
        case 'violin'
            subplot(2,3,roi)
            plot([0 10], [0 0], 'k')
            for ii = 1:p.numAttWindows
                hold on,
                h               = Violin(varError(roi,ii), p.WidthAttWindow(ii), 'ViolinColor',{p.CondColors(ii,:)}, 'HalfViolin','full', 'QuartileStyle','boxplot', 'DataStyle', 'none', 'Width', 0.7);
                h.BoxColor      = p.CondColors(ii,:);
            end
            box off
            set(gca, 'XTick', p.WidthAttWindow, 'XTickLabel', p.WidthAttWindow, 'YTick', [-pi -pi/2 0 pi/2 pi], 'YTickLabel', -180:90:180)
            xlabel('Width condition'); ylabel('Angular error (deg)')
            title(p.ROInames{roi})
            ylim([-pi pi])
    end
    
    subplot(2,3,roi+p.numROIs)
    hold on
    for s = 1:p.nSubs
        if s == p.exampleSubj
            plot(p.WidthAttWindow,medErrR2(s,:,roi),'^--','color',[.5 .5 .5],'MarkerSize',5, 'MarkerFaceColor', [.5 .5 .5])
        else
            plot(p.WidthAttWindow,medErrR2(s,:,roi),'.--','color',[.5 .5 .5],'MarkerSize',12)
        end

    end
    errorbar(p.WidthAttWindow,mean(medErrR2(:,:,roi)),std(medErrR2(:,:,roi))/sqrt(p.nSubs),'.-k','MarkerSize',20,'CapSize',0,'LineWidth',2)
    set(gca, 'XTick', p.WidthAttWindow, 'XTickLabel', p.WidthAttWindow, 'YTick', 0:pi/4:(1/2*pi), 'YTickLabel', {'0','45','90'}), 
    xlim([0.5 9.5]), ylim([0 1/2*pi])
    title(p.ROInames{roi})
    if roi==1
        ylabel('absolute error (deg)')
    end
    xlabel('Width condition')

end
sgtitle(sprintf('%s: Angular error in %s trial width location estimate', p.savestr, type))

if saveFig
    if ~exist(figureDir, 'dir'), mkdir(figureDir); end
    print(fullfile(figureDir, sprintf('%s_%sTrials_%s_%s', mfilename, type, p.savestr, plotType)), '-dpdf', '-painters', '-bestfit')
end

%% Save csv files for stats
for roi = 1:p.numROIs
    
    angErr_width = rad2deg(medErrR2(:,:,roi));
    
    T = array2table(angErr_width);
    T.Properties.RowNames = p.SubjNames;
    T.Properties.Description = 'Model results: absolute angular error of spatial profile fits';
    T.Properties.VariableUnits = repmat({'deg'}, [1, p.numAttWindows]);
    writetable(T, fullfile(p.csvDir, sprintf('%s_locError_%s.csv', p.savestr, p.ROInames{roi})))
    
end



