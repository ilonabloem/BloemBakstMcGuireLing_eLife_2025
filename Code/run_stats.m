
function run_stats(dataDir)
% Replicate stats for Bloem, Bakst, McGuire & Ling (2025) eLife
% Scripts reads CSV files that are created in the scripts run_createFigures and run_createFiguresTR
% The location of the CSV files is fullfile(dataDir, 'Data', 'csvs').
%
% LRB May 2024


%% Setup paths
prm.projectDir  = projectRootPath;
if ~exist('dataDir', 'var') || isempty(dataDir)
    prm.dataDir     = fullfile(projectRootPath);
end
if ~exist(fullfile(prm.dataDir, 'Data'), 'dir') > 0
    error('Data folder not found within the project directory')
end
csvDir          = fullfile(prm.dataDir, 'Data', 'csvs', 'CM');
variableTRdir   = fullfile(csvDir,'variableTR');

%% Setup
prm.Task          = 'CM';
prm               = setupAnalysisInfo(prm);

%% variable TR analysis
%%%% Angular error vs TR/cue width; FWHM vs TR %%%%
disp('****VARIABLE TR****')

%-- find number of TR analysis:
TRs             = [1 2 3 5 10];
nTRs            = length(TRs);

xTR             = repmat(TRs,[1 prm.numAttWindows])';
xWidth          = repmat(prm.WidthAttWindow, [nTRs 1]);
xWidth          = xWidth(:);

angBetasTR      = NaN(prm.numROIs, prm.nSubs, 2);
ampBetasTR      = NaN(prm.numROIs, prm.nSubs, 2);
widBetasTR      = NaN(prm.numROIs, prm.nSubs, 2);
basBetasTR      = NaN(prm.numROIs, prm.nSubs, 2);
R2BetasTR       = NaN(prm.numROIs, prm.nSubs, 2);

for roi = 1:prm.numROIs
    
    angBOLD     = NaN(prm.nSubs,prm.numAttWindows*nTRs,1);
    ampBOLD     = NaN(prm.nSubs,prm.numAttWindows*nTRs,1);
    widBOLD     = NaN(prm.nSubs,prm.numAttWindows*nTRs,1);
    basBOLD     = NaN(prm.nSubs,prm.numAttWindows*nTRs,1);
    R2BOLD      = NaN(prm.nSubs,prm.numAttWindows*nTRs,1);
    
    counter = 1;
    for wd = 1:prm.numAttWindows 
        angData     = readtable(fullfile(variableTRdir, sprintf('ang_variableTR_Width%d_%s.csv', wd, prm.ROInames{roi})));
        widData     = readtable(fullfile(variableTRdir, sprintf('width_variableTR_Width%d_%s.csv', wd, prm.ROInames{roi})));
        R2Data      = readtable(fullfile(variableTRdir, sprintf('R2_variableTR_Width%d_%s.csv', wd, prm.ROInames{roi})));
        ampData     = readtable(fullfile(variableTRdir, sprintf('gain_variableTR_Width%d_%s.csv', wd, prm.ROInames{roi})));
        basData     = readtable(fullfile(variableTRdir, sprintf('base_variableTR_Width%d_%s.csv', wd, prm.ROInames{roi})));

        for sub = 1:prm.nSubs
            angBOLD(sub, counter:counter+(nTRs-1))   = angData{sub,1:nTRs};
            ampBOLD(sub, counter:counter+(nTRs-1))   = ampData{sub,1:nTRs};
            widBOLD(sub, counter:counter+(nTRs-1))   = widData{sub,1:nTRs};
            basBOLD(sub, counter:counter+(nTRs-1))   = basData{sub,1:nTRs};
            R2BOLD(sub, counter:counter+(nTRs-1))    = R2Data{sub,1:nTRs};
        end

        counter = counter + nTRs;
    end
    
    for sub = 1:prm.nSubs
        angLM                   = fitlm([xTR, xWidth],angBOLD(sub,:)');
        angBetasTR(roi,sub,:)   = angLM.Coefficients.Estimate(2:3);
        
        ampLM                   = fitlm([xTR, xWidth],ampBOLD(sub,:)');
        ampBetasTR(roi,sub,:)   = ampLM.Coefficients.Estimate(2:3);
        
        widLM                   = fitlm([xTR, xWidth],widBOLD(sub,:)');
        widBetasTR(roi,sub,:)   = widLM.Coefficients.Estimate(2:3);
        
        basLM                   = fitlm([xTR, xWidth],basBOLD(sub,:)');
        basBetasTR(roi,sub,:)   = basLM.Coefficients.Estimate(2:3);
        
        R2LM                    = fitlm([xTR, xWidth],R2BOLD(sub,:)');
        R2BetasTR(roi,sub,:)    = R2LM.Coefficients.Estimate(2:3);
    end
    
    disp(['~~~~~~' prm.ROInames{roi} '~~~~~~'])
    
    disp('---angular error---')
    [~,angRes,~,stats] = ttest(angBetasTR(roi,:,1));
    disp(['TR: t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(angRes,3)) ', Bonferroni corrected p = ' num2str(round(angRes*prm.numROIs,3))])
    [~,angRes,~,stats] = ttest(angBetasTR(roi,:,2));
    disp(['Cue width: t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(angRes,3)) ', Bonferroni corrected p = ' num2str(round(angRes*prm.numROIs,3))])
    
    disp('---FWHM---')
    [~,widRes,~,stats] = ttest(widBetasTR(roi,:,1));
    disp(['TR: t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(widRes,3)) ', Bonferroni corrected p = ' num2str(round(widRes*prm.numROIs,3))])
    [~,widRes,~,stats] = ttest(widBetasTR(roi,:,2));
    disp(['Cue width: t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(widRes,3)) ', Bonferroni corrected p = ' num2str(round(widRes*prm.numROIs,3))])
    
    disp('---Gain---')
    [~,ampRes,~,stats] = ttest(ampBetasTR(roi,:,1));
    disp(['TR: t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(ampRes,3)) ', Bonferroni corrected p = ' num2str(round(ampRes*prm.numROIs,3))])
    [~,ampRes,~,stats] = ttest(ampBetasTR(roi,:,2));
    disp(['Cue width: t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(ampRes,3)) ', Bonferroni corrected p = ' num2str(round(ampRes*prm.numROIs,3))])
    
    disp('---Baseline---')
    [~,basRes,~,stats] = ttest(basBetasTR(roi,:,1));
    disp(['TR: t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(basRes,3)) ', Bonferroni corrected p = ' num2str(round(basRes*prm.numROIs,3))])
    [~,basRes,~,stats] = ttest(basBetasTR(roi,:,2));
    disp(['Cue width: t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(basRes,3)) ', Bonferroni corrected p = ' num2str(round(basRes*prm.numROIs,3))])
    
    disp('---R2---')
    [~,R2Res,~,stats] = ttest(R2BetasTR(roi,:,1));
    disp(['TR: t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(R2Res,3)) ', Bonferroni corrected p = ' num2str(round(R2Res*prm.numROIs,3))])
    [~,R2Res,~,stats] = ttest(R2BetasTR(roi,:,2));
    disp(['Cue width: t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(R2Res,3)) ', Bonferroni corrected p = ' num2str(round(R2Res*prm.numROIs,3))])
    
    pause
    
end

%% 
%%%% FWHM vs cue width %%%%
disp('---FWHM vs cue width, by TR---')
widBetasbyTR = NaN(prm.numROIs, prm.nSubs);
for roi = 1:prm.numROIs
    disp(['~~~~~' prm.ROInames{roi} '~~~~~'])
    disp('TR      t(DF)                p-value')
    for tr = 1:nTRs
        widTR = NaN(prm.nSubs, prm.numAttWindows);
        for wd = 1:prm.numAttWindows
            widData     = readtable(fullfile(variableTRdir, sprintf('width_variableTR_Width%d_%s.csv', wd, prm.ROInames{roi})));
            widTR(:,wd) = widData{:,tr};
        end
        for sub = 1:prm.nSubs
            widLM                = fitlm(prm.WidthAttWindow, widTR(sub,:));
            widBetasbyTR(roi,sub)= widLM.Coefficients.Estimate(2);
        end
        [~,p,~,stats]       = ttest(widBetasbyTR(roi,:));
        disp([num2str(TRs(tr)) '      t(' num2str(stats.df) ') = ' num2str(stats.tstat) ...
            ', p = ' num2str(p)])
    end
    pause
end

%% Summary stats: Attention
disp('****ATTENTION SUMMARY****')

%%%% cue width
angBetas = NaN(prm.numROIs, prm.nSubs);
ampBetas = NaN(prm.numROIs, prm.nSubs);
widBetas = NaN(prm.numROIs, prm.nSubs);
basBetas = NaN(prm.numROIs, prm.nSubs);

for roi = 1:prm.numROIs
    angData = readtable(fullfile(csvDir, sprintf('attention_locError_%s.csv', prm.ROInames{roi})));
    ampData = readtable(fullfile(csvDir, sprintf('attention_Ampl_%s.csv', prm.ROInames{roi})));
    widData = readtable(fullfile(csvDir, sprintf('attention_FWHM_%s.csv', prm.ROInames{roi})));
    basData = readtable(fullfile(csvDir, sprintf('attention_Base_%s.csv', prm.ROInames{roi})));
    
    if roi==1
        V1 = ampData;
    elseif roi==2
        V2 = ampData;
    elseif roi==3
        V3 = ampData;
    end
    
    for sub = 1:prm.nSubs
        angLM = fitlm(prm.WidthAttWindow,angData{sub,:});
        angBetas(roi,sub) = angLM.Coefficients.Estimate(2);
        
        ampLM = fitlm(prm.WidthAttWindow,ampData{sub,:});
        ampBetas(roi,sub) = ampLM.Coefficients.Estimate(2);
        
        widLM = fitlm(prm.WidthAttWindow,widData{sub,:});
        widBetas(roi,sub) = widLM.Coefficients.Estimate(2);
        
        basLM = fitlm(prm.WidthAttWindow,basData{sub,:});
        basBetas(roi,sub) = basLM.Coefficients.Estimate(2);
    end
        
    disp(prm.ROInames{roi})
    disp('---angular error---')
    [~,p,~,stats] = ttest(angBetas(roi,:));
    disp(['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(p,3)), ', Bonferroni corrected p = ' num2str(round(p*prm.numROIs,3))])
      
    disp('---gain---')
    [~,p,~,stats] = ttest(ampBetas(roi,:));
    disp(['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(p,3)), ', Bonferroni corrected p = ' num2str(round(p*prm.numROIs,3))])
    
    disp('---FWHM---')
    [~,p,~,stats] = ttest(widBetas(roi,:));
    disp(['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(p,3)), ', Bonferroni corrected p = ' num2str(round(p*prm.numROIs,3))])
    
    disp('---Baseline---')
    [~,p,~,stats] = ttest(basBetas(roi,:));
    disp(['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(p,3)), ', Bonferroni corrected p = ' num2str(round(p*prm.numROIs,3))])
    
    pause
end
V1 = V1{:,:};
V2 = V2{:,:};
V3 = V3{:,:};
disp('~~~~ gain across ROIs ~~~~')
[~,p] = ttest(V1(:),V2(:));
disp(['V1 vs V2, p = ' num2str(round(p,3))])
[~,p] = ttest(V1(:),V3(:));
disp(['V1 vs V3, p = ' num2str(round(p,3))])
[~,p] = ttest(V2(:),V3(:));
disp(['V2 vs V3, p = ' num2str(round(p,3))])


%% V3 AMP CHECK %%
figure;
hold on
for sub = 1:prm.nSubs
    ampLM = fitlm(prm.WidthAttWindow,ampData{sub,:});
    plot(prm.WidthAttWindow,ampData{sub,:},'-')
    plot(prm.WidthAttWindow,polyval(flipud(ampLM.Coefficients.Estimate)',prm.WidthAttWindow),':','color',[.5 .5 .5])
    xlabel('width')
    xticks(prm.WidthAttWindow)
    ylabel('amplitude')
    pause
end

%% Summary stats: Perception
disp('****PERCEPTION SUMMARY****')
prm.Task        = 'PCM';
csvDir          = fullfile(prm.dataDir, 'Data', 'csvs', prm.Task);
prm             = setupAnalysisInfo(prm);

angBetasPCM     = NaN(prm.numROIs, prm.nSubs);
ampBetasPCM     = NaN(prm.numROIs, prm.nSubs);
widBetasPCM     = NaN(prm.numROIs, prm.nSubs);
basBetasPCM     = NaN(prm.numROIs, prm.nSubs);

%%%% cue width
for roi = 1:prm.numROIs
    angData = readtable(fullfile(csvDir, sprintf('perception_locError_%s.csv', prm.ROInames{roi})));
    ampData = readtable(fullfile(csvDir, sprintf('perception_Ampl_%s.csv', prm.ROInames{roi})));
    widData = readtable(fullfile(csvDir, sprintf('perception_FWHM_%s.csv', prm.ROInames{roi})));
    basData = readtable(fullfile(csvDir, sprintf('perception_Base_%s.csv', prm.ROInames{roi})));
    
    for sub = 1:prm.nSubs 
        angLM = fitlm(prm.WidthAttWindow, angData{sub,:});
        angBetasPCM(roi,sub) = angLM.Coefficients.Estimate(2);
        
        ampLM = fitlm(prm.WidthAttWindow,ampData{sub,:});
        ampBetasPCM(roi,sub) = ampLM.Coefficients.Estimate(2);
        
        widLM = fitlm(prm.WidthAttWindow,widData{sub,:});
        widBetasPCM(roi,sub) = widLM.Coefficients.Estimate(2);
        
        basLM = fitlm(prm.WidthAttWindow,basData{sub,:});
        basBetasPCM(roi,sub) = basLM.Coefficients.Estimate(2);
    end
    
    disp(prm.ROInames{roi})
    disp('---angular error---')
    [~,p,~,stats] = ttest(angBetasPCM(roi,:));
    disp(['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(p,3))])
    
    disp('---gain---')
    [~,p,~,stats] = ttest(ampBetasPCM(roi,:));
    disp(['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(p,3))])
    
    disp('---FWHM---')
    [~,p,~,stats] = ttest(widBetasPCM(roi,:));
    disp(['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(p,3))])
    
    disp('---Baseline---')
    [~,p,~,stats] = ttest(basBetasPCM(roi,:));
    disp(['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(p,3))])
    
    pause
end

%% Perception & Attention comparison
disp('****PERCEPTION vs ATTENTION ****')
attData     = readtable(fullfile(csvDir, '..', 'CM', 'attention_avgFWHM.csv'));
attData     = attData{:,2:4};
pcptData    = readtable(fullfile(csvDir, '..', 'PCM', 'perception_avgFWHM.csv'));
pcptData    = pcptData{[1:3 5],2:4};
disp(['Pearson correlation = ' num2str(corr(pcptData(:),attData(:)))])

[~, p, ~, stats]       = ttest2(pcptData(:),attData(:));
disp(['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat,2)) ', p = ' num2str(round(p,3))]) 
    