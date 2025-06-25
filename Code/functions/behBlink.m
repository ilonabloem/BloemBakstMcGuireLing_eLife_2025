function [allData]= behBlink(allData,nSubs,nRuns,viz)
%% Remove blink artifacts using pupil trace

for sub=1:nSubs
    for run=1:nRuns(sub)
        sWin        = 601;
        diffThresh  = 7;
        start       = find(allData{sub}.eye{run}.sample==allData{sub}.msg{run}.cueOn(1));
        stop        = find(allData{sub}.eye{run}.sample==allData{sub}.msg{run}.key(end));
        
        if isempty(start)
            start   = 1;      
        end
        
        if isempty(stop)
            stop    = length(allData{sub}.eye{run}.sample);
        end
        
        grandMean                                   = nanmean(allData{sub}.eye{run}.psize(start:stop));
        allData{sub}.eye{run}.pupilRaw              = (allData{sub}.eye{run}.psize-grandMean)./grandMean*100;
        allData{sub}.eye{run}.pupilSmooth           = movmedian(allData{sub}.eye{run}.pupilRaw,sWin,'omitnan');
        absDiff                                     = abs(allData{sub}.eye{run}.pupilRaw-allData{sub}.eye{run}.pupilSmooth);
        diffIdx                                     = absDiff<diffThresh;
        allData{sub}.eye{run}.pupilDiff             = allData{sub}.eye{run}.pupilRaw;
        allData{sub}.eye{run}.pupilDiff(~diffIdx)   = NaN;
        
        
        allData{sub}.eye{run}.xposTB                = allData{sub}.eye{run}.xposT;
        allData{sub}.eye{run}.yposTB                = allData{sub}.eye{run}.yposT;
        allData{sub}.eye{run}.xposTB(~diffIdx)      = NaN;
        allData{sub}.eye{run}.yposTB(~diffIdx)      = NaN;

        if viz
            figure(1)
            clf
            hold on
            plot(allData{sub}.eye{run}.pupilRaw)
            plot(allData{sub}.eye{run}.pupilDiff)

            figure(2)
            clf
            hold on
            plot(allData{sub}.eye{run}.xposT,allData{sub}.eye{run}.yposT)
            plot(allData{sub}.eye{run}.xposTB,allData{sub}.eye{run}.yposTB)

            figure(3)
            clf
            hold on
            plot(sort(absDiff))
            plot([0 length(absDiff)], [diffThresh diffThresh],'--k')


            disp(num2str(sub))
            disp(num2str(run))
            pause
        end
    end
end
disp('blink')