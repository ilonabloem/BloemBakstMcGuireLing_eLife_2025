
function figureActivityProfilesFWHM(p, Data, results, type, figureDir, saveFig)


names           = fieldnames(results);

%% Extract data Attention
avgData_R2          = NaN(p.numAttWindows, p.numROIs, p.nSubs, p.nBins);
avgDataTM_R2        = NaN(p.numAttWindows, p.numROIs, p.nSubs, p.nBins); %trueMean

for sub = 1:p.nSubs
    
    R2idx = results(sub).(names{contains(names, sprintf('%sR2idx', type))});
    whichWidth = results(sub).(names{contains(names, sprintf('%sWhichWidth', type))});
    
    
    for roi = 1:p.numROIs
%         FWHMoutliers = results(sub).(names{contains(names, 'FWHMoutliers')})(roi,:);
        
        for ii = 1:p.numAttWindows
            
            % Recentered data profiles
            avgData_R2(ii,roi,sub,:)    = nanmean(Data(roi).RecenteredData{sub}(:, whichWidth(:,ii) & R2idx),2);
            avgDataTM_R2(ii,roi,sub,:)  = nanmean(Data(roi).TrueMeanRecenteredData{sub}(:, whichWidth(:,ii) & R2idx),2);

        end
                 
    end

end

%% True mean recentered data %%%%
bounds  = @(x) [floor((min(x(:)))*10)/10 ceil((max(x(:)))*10)/10];
yBounds = bounds(squeeze(nanmean(avgDataTM_R2,3)));

figure('Color', [1 1 1], 'Position', [50 200 1500 500])
for roi = 1:p.numROIs
    
    subplot(1, p.numROIs, roi)
    for ii = p.numAttWindows:-1:1
        hold on,
        
        hP = patch([p.binCenters fliplr(p.binCenters)], [nanmedian(squeeze(avgDataTM_R2(ii,roi,:,:)),1)+nanstd(squeeze(avgDataTM_R2(ii,roi,:,:)),1)/sqrt(p.nSubs) ...
            fliplr(nanmedian(squeeze(avgDataTM_R2(ii,roi,:,:)),1)-nanstd(squeeze(avgDataTM_R2(ii,roi,:,:)),1)/sqrt(p.nSubs))], 'g');
        set(hP, 'FaceColor', p.CondColors(ii,:), 'FaceAlpha', 0.7, 'EdgeColor', p.CondColors(ii,:))
        h(ii) = plot(p.binCenters, nanmedian(squeeze(avgDataTM_R2(ii,roi,:,:)),1), 'Color', p.CondColors(ii,:),'LineWidth',2);
        
    end
    box off
    xlabel('Angle (Rad)'), ylabel('BOLD response'); title('Recentered true mean'),
    set(gca, 'XTick', [0 pi 2*pi], 'XTickLabel', {'0' 'pi' '2pi'}); xlim([0 2*pi])
    ylim(yBounds), axis square
end
legend(h, num2str(p.WidthAttWindow'));
sgtitle(sprintf('fig2_%s_%s_trueMean_activityProfile', p.savestr, type), 'interpreter', 'none');

if saveFig
    if ~exist(figureDir, 'dir'), mkdir(figureDir); end
    print(fullfile(figureDir, sprintf('fig2_%s_%s_trueMean_activityProfile', p.savestr, type)), '-dpdf', '-vector','-bestfit')
end

