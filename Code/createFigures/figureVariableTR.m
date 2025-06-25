function figureVariableTR(p, Data, results, figureDir, saveFig)

% Visualize FWHM and angular error as a function of #TRs used to average the
% data

names           = fieldnames(results);
alphas          = linspace(0.3,1,numel(Data));  
  
%-- setup figure
figure('Color', [1 1 1], 'Position', [200 20 700 700]);

for roi=1:p.numROIs
      
    %-- preallocate variables
    AngMedians      = NaN(p.numTRs,p.nSubs,p.numAttWindows);    
    WidthMedians    = NaN(p.numTRs,p.nSubs,p.numAttWindows);
    AmplMedians     = NaN(p.numTRs,p.nSubs,p.numAttWindows);    
    R2Medians       = NaN(p.numTRs,p.nSubs,p.numAttWindows); 
    BaseMedians     = NaN(p.numTRs,p.nSubs,p.numAttWindows);   
    
    for tr = 1:numel(Data)
        
        for sub=1:p.nSubs
            
            R2idx = results(sub, tr).(names{contains(names, 'R2idx')});
            whichWidth = results(sub, tr).(names{contains(names, 'WhichWidth')});
           
            for ii=1:p.numAttWindows
                WidthMedians(tr,sub,ii)    = nanmedian(Data(tr).avgResults(roi).FWHMEstimates{sub}(whichWidth(:,ii) & R2idx));
                
                AngMedians(tr,sub,ii)       = nanmedian(abs(Data(tr).avgResults(roi).diffAng{sub}(whichWidth(:,ii) & R2idx)));
                
                ampl                        = Data(tr).avgResults(roi).amplEstimates{sub}(whichWidth(:,ii) & R2idx);
                AmplMedians(tr,sub,ii)      = nanmedian(ampl);
                
                BaseMedians(tr,sub,ii)      = nanmedian(Data(tr).avgResults(roi).baseEstimates{sub}(whichWidth(:,ii) & R2idx));

                R2Medians(tr,sub,ii)        = nanmedian(Data(tr).avgResults(roi).R2{sub}(whichWidth(:,ii) & R2idx));
            end
            
            
        end
        
        axIndx = 0:7:7*numel(Data)-1;
        for wd = 1:p.numAttWindows
            
            subplot(5,p.numROIs,roi)
            hold on,
            h = errorbar(tr+axIndx(wd), mean(squeeze(WidthMedians(tr,:,wd))), std(squeeze(WidthMedians(tr,:,wd)))/sqrt(p.nSubs), ...
                '-', 'Color', p.CondColors(wd,:), 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2);

            % Set transparency
            set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alphas(tr)])
            sc1 = scatter(tr+axIndx(wd), mean(squeeze(WidthMedians(tr,:,wd))), [], p.CondColors(wd,:), 'filled');
            alpha(sc1 , alphas(tr));
            
            if wd == p.numAttWindows
                ylim([50 200])
                
                set(gca, 'XTick', axIndx+numel(Data)/2, 'XTickLabel', p.WidthAttWindow')
                if roi == 1
                    ylabel('FWHM (deg)')
                end
            end

            subplot(5,p.numROIs,roi+p.numROIs)
            hold on,
            h = errorbar(tr+axIndx(wd), mean(squeeze(AngMedians(tr,:,wd))), std(squeeze(AngMedians(tr,:,wd)))/sqrt(p.nSubs), ...
                '-', 'Color', p.CondColors(wd,:), 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2);
 
            % Set transparency
            set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alphas(tr)])
            sc1 = scatter(tr+axIndx(wd), mean(squeeze(AngMedians(tr,:,wd))), [], p.CondColors(wd,:), 'filled');
            alpha(sc1 , alphas(tr));
            
            if wd == p.numAttWindows
                title(p.ROInames{roi})
                ylim([0 pi/2]), set(gca, 'YTick', linspace(0,pi/2, 3), 'YTickLabel', {'0', '45' '90'})
                set(gca, 'XTick', axIndx+numel(Data)/2, 'XTickLabel', p.WidthAttWindow')
                if roi == 1
                    ylabel('Absolute angular error (deg)')
                end
            end
            
            subplot(5,p.numROIs,roi+p.numROIs*2)
            hold on,
            h = errorbar(tr+axIndx(wd), mean(squeeze(AmplMedians(tr,:,wd))), std(squeeze(AmplMedians(tr,:,wd)))/sqrt(p.nSubs), ...
                '-', 'Color', p.CondColors(wd,:), 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2);

            % Set transparency
            set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alphas(tr)])
            sc1 = scatter(tr+axIndx(wd), mean(squeeze(AmplMedians(tr,:,wd))), [], p.CondColors(wd,:), 'filled');
            alpha(sc1 , alphas(tr));
            
            if tr == numel(Data) && wd == p.numAttWindows
                yl = ylim; yl = [0 ceil(yl(2)/0.5)*0.5];
                ylim(yl)
                
                set(gca, 'XTick', axIndx+numel(Data)/2, 'XTickLabel', p.WidthAttWindow')
                if roi == 1
                    ylabel('Gain (%SC)')
                end
            end
            
            subplot(5,p.numROIs,roi+p.numROIs*3)
            hold on,
            h = errorbar(tr+axIndx(wd), mean(squeeze(BaseMedians(tr,:,wd))), std(squeeze(BaseMedians(tr,:,wd)))/sqrt(p.nSubs), ...
                '-', 'Color', p.CondColors(wd,:), 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2);

            % Set transparency
            set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alphas(tr)])
            sc1 = scatter(tr+axIndx(wd), mean(squeeze(BaseMedians(tr,:,wd))), [], p.CondColors(wd,:), 'filled');
            alpha(sc1 , alphas(tr));
            
            if tr == numel(Data) && wd == p.numAttWindows
                ylim([-0.6 0])
                
                set(gca, 'XTick', axIndx+numel(Data)/2, 'XTickLabel', p.WidthAttWindow')
                if roi == 1
                    ylabel('Baseline (%SC)')
                end
            end
            
            subplot(5,p.numROIs,roi+p.numROIs*4)
            hold on,
            h = errorbar(tr+axIndx(wd), mean(squeeze(R2Medians(tr,:,wd))), std(squeeze(R2Medians(tr,:,wd)))/sqrt(p.nSubs), ...
                '-', 'Color', p.CondColors(wd,:), 'MarkerSize', 10, 'CapSize', 0, 'LineWidth', 2);
            % Set transparency
            set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alphas(tr)])
            
            sc1 = scatter(tr+axIndx(wd), mean(squeeze(R2Medians(tr,:,wd))), [], p.CondColors(wd,:), 'filled');
            alpha(sc1 , alphas(tr));
            
            if tr == numel(Data) && wd == p.numAttWindows
                
                ylim([0.3 0.6])
                
                set(gca, 'XTick', axIndx+numel(Data)/2, 'XTickLabel', p.WidthAttWindow')
                if roi == 1
                    ylabel('R2')
                end
            end
            
        end
        
    end
    
    %-- Save CSV files 
    for ii = 1:p.numAttWindows
        ang_TR      = AngMedians(: ,:, ii)';
        saveName    = fullfile(figureDir, sprintf('ang_variableTR_Width%i_%s.csv', ii, p.ROInames{roi}));
        T           = table(ang_TR(:,1), ang_TR(:,2), ang_TR(:,3), ang_TR(:,4), ang_TR(:,5));
        T.subjID    = p.SubjNames(:);
        writetable(T, saveName)
        
        width_TR     = WidthMedians(: ,:, ii)';
        saveName    = fullfile(figureDir, sprintf('width_variableTR_Width%i_%s.csv', ii, p.ROInames{roi}));
        T           = table(width_TR(:,1), width_TR(:,2), width_TR(:,3), width_TR(:,4), width_TR(:,5));
        T.subjID    = p.SubjNames(:);
        writetable(T, saveName)
        
        ampl_TR     = AmplMedians(: ,:, ii)';
        saveName    = fullfile(figureDir, sprintf('gain_variableTR_Width%i_%s.csv', ii, p.ROInames{roi}));
        T           = table(ampl_TR(:,1), ampl_TR(:,2), ampl_TR(:,3), ampl_TR(:,4), ampl_TR(:,5));
        T.subjID    = p.SubjNames(:);
        writetable(T, saveName)
        
        base_TR     = BaseMedians(: ,:, ii)';
        saveName    = fullfile(figureDir, sprintf('base_variableTR_Width%i_%s.csv', ii, p.ROInames{roi}));
        T           = table(base_TR(:,1), base_TR(:,2), base_TR(:,3), base_TR(:,4), base_TR(:,5));
        T.subjID    = p.SubjNames(:);
        writetable(T, saveName)
        
        r2_TR       = R2Medians(: ,:, ii)';
        saveName    = fullfile(figureDir, sprintf('R2_variableTR_Width%i_%s.csv', ii, p.ROInames{roi}));
        T           = table(r2_TR(:,1), r2_TR(:,2), r2_TR(:,3), r2_TR(:,4), r2_TR(:,5));
        T.subjID    = p.SubjNames(:);
        writetable(T, saveName)
        
    end
        
end

sgtitle(sprintf('Variable TR analysis'))
   
if saveFig
    if ~exist(figureDir, 'dir'), mkdir(figureDir); end
    print(fullfile(figureDir, sprintf('%s_%s', mfilename, p.savestr)), '-dpdf', '-painters', '-bestfit')

end



