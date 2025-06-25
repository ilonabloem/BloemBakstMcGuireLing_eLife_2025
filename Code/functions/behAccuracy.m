function [allData]=behAccuracy(allData,nFiles,nWidths,widths,figureDir,saveFig)

nSubs       = length(allData);
nTrials     = length(allData{1}.data{1}.width);
allPctCorr  = NaN(nSubs,nWidths+1);
nanwidth    = zeros(nSubs,nWidths+1);

for sub = 1:nSubs
    perf    = zeros(1,nWidths+1);
    tot     = zeros(1,nWidths+1);
    for run = 1:nFiles(sub)
        idxResponse = ~isnan(allData{sub}.data{run}.keyEO);
        for ww = 1:nWidths
            idxWidth    = allData{sub}.data{run}.width==widths(ww);
            perf(ww)    = perf(ww) + sum(allData{sub}.data{run}.Correct(idxWidth&idxResponse));
            tot(ww)     = tot(ww) + sum(idxWidth&idxResponse);
        end
    end
    allData{sub}.pctCorr    = perf./tot;
    allPctCorr(sub,:)       = perf./tot;
end

testvChance = NaN(1,nWidths);

figure('Color', [1 1 1], 'Position', [200 20 600 500]);
clf
hold on
bar(nanmean(allPctCorr(:,1:nWidths)),'EdgeColor','k','FaceColor',[.8 .8 .8])
errorbar(1:nWidths,nanmean(allPctCorr(:,1:nWidths)),nanstd(allPctCorr(:,1:nWidths))/sqrt(nSubs),'.k','LineWidth',2)
for wid=1:nWidths
    scatter(ones(1,nSubs)*wid,allPctCorr(:,wid),'markerfacecolor',[.3 .3 .3],'markeredgecolor','k','jitter','on','markerfacealpha',0.5)
    [~,testvChance(wid)] = ttest(allPctCorr(:,wid)-0.5);
end
set(gca,'XTick',1:nWidths,'XTickLabel',widths)
set(gca,'YTick',.5:.1:1)
set(gca,'FontSize',16)
xlabel('Attention cue')
ylabel('Proportion correct')
ylim([.5 1])
title('Task performance')
if saveFig
    print(fullfile(figureDir, sprintf('behavioralPerformance')), '-dpdf', '-painters', '-bestfit')
end

% regression to test for linear decreasing relationship
allX = repmat([1 3 5 9]*18,[8,1]);
allY = allPctCorr(:,1:nWidths);
disp('Test for linear decreasing relationship between width and behavioral accuracy...')
fitlm(allX(:),allY(:))

% ttest for better performance in widest condition
[~,p] = ttest(allPctCorr(:,3),allPctCorr(:,4));
disp('Test for better performance in width 9 compared to width 5...')
disp(['p = ' num2str(p)])
