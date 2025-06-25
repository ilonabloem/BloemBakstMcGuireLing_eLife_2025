function [] = behGazePos(allData,nFiles,nSubs,nLetters,widths,eccThresh,condColors,figureDir,saveFig)
%% look at gaze position by width

blockLen    = 5;
nWidths     = length(widths);
rhoSummary  = NaN(nWidths,nLetters+1,nSubs);
overallRho  = NaN(2,nSubs);

for sub=1:nSubs

    rhoWidth    = cell(nWidths,nLetters);
    rhoAll      = [];
    
    for run=1:nFiles(sub)
        nTrials = length(allData{sub}.data{run}.Correct);
        [allData{sub}.eye{run}.theta,allData{sub}.eye{run}.rho] = cart2pol(allData{sub}.eye{run}.xposTB,allData{sub}.eye{run}.yposTB);
        rhoAll = [rhoAll,allData{sub}.eye{run}.rho];
        
        for trial=1:blockLen:nTrials
            widthIdx                    = find(widths==allData{sub}.data{run}.width(trial));
            locIdx                      = allData{sub}.data{run}.winCtrLetterPos(trial)+1;
            startIdx                    = find(allData{sub}.eye{run}.sample==allData{sub}.msg{run}.cueOn(trial));
            endIdx                      = find(allData{sub}.eye{run}.sample==allData{sub}.msg{run}.key(trial+(blockLen-1)));
            rhoWidth{widthIdx,locIdx}   = horzcat(rhoWidth{widthIdx,locIdx},allData{sub}.eye{run}.rho(startIdx:endIdx));
        end
    end
    for ww=1:nWidths
        for loc=1:nLetters
            rhoSummary(ww,loc,sub) = nanmean(rhoWidth{ww,loc});
        end
        rhoSummary(ww,end,sub) = rhoSummary(ww,1,sub);
   
    end
    overallRho(1,sub) = nanmean(rhoAll);
end

figure('Color', [1 1 1], 'Position', [200 20 600 500]);
avgRho = squeeze(nanmean(rhoSummary,3));
p1 = polarplot([linspace(0,(2*pi)-((2*pi)/nLetters),nLetters),0],avgRho(1,:),'LineWidth',2,'color',condColors(1,:));
hold on
p2 = polarplot([linspace(0,(2*pi)-((2*pi)/nLetters),nLetters),0],avgRho(2,:),'LineWidth',2,'color',condColors(2,:));
p3 = polarplot([linspace(0,(2*pi)-((2*pi)/nLetters),nLetters),0],avgRho(3,:),'LineWidth',2,'color',condColors(3,:));
p4 = polarplot([linspace(0,(2*pi)-((2*pi)/nLetters),nLetters),0],avgRho(4,:),'LineWidth',2,'color',condColors(4,:));
set(gca,'FontSize',16)
%title(num2str(widths(width)))
legend([p1 p2 p3 p4],num2str(widths'))
thetaticks(0:45:359)
ax = ancestor(p1, 'polaraxes');
ax.RAxisLocation = 60;
title('Gaze eccentricity')

if saveFig
    print(fullfile(figureDir, sprintf('gazeEccentricity')), '-dpdf', '-painters', '-bestfit')
end

disp('Check whether the average gaze eccentricity for each participant is less than threshold...')
disp(['Number of participants exceeding ' num2str(eccThresh) ' = ' num2str(sum(overallRho(1,:)>eccThresh))])