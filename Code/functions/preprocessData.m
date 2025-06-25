
function out = preprocessData(p, plotLocalizer, saveFig)

%% Prepare data for attentional window analysis
%
% function out = preprocessData(p, plotLocalizer, saveFig)
%
% p = set of params
% plotLocalizer = true or false
% saveFig = true or false
%
% IB - 2025

if ~exist(fullfile(p.figureDir, 'voxelSelection'), 'dir'), mkdir(p.figureDir, 'voxelSelection'), end

fprintf('Loading fMRI data... \n')

%% setup figure
switch plotLocalizer
    case true
        for ii= 1:p.numROIs
            figHandle{ii}           = length(findobj('type', 'figure'))+1;
            figure(figHandle{ii})
            set(figHandle{ii},'color', [1 1 1], 'pos', [10 300 1800 600] ...
                ,'Units', 'Pixels', 'PaperPositionMode','Auto', ...
                'PaperUnits','points','PaperSize',[1800 600])


            colormap(parula);
        end
end
    
%% Preprocess fMRI data and select voxels in ROIs
out     = struct('final_tSeries', cell(1,numel(p.SubjNames)), ...
                  'selectedVoxels', cell(1,numel(p.SubjNames)), ...
                  'pRF_maps', cell(1,numel(p.SubjNames)), ...
                  'ang_pRFestimates', cell(1,numel(p.SubjNames)), ...
                  'ecc_pRFestimates', cell(1,numel(p.SubjNames)), ...
                  'rf_pRFestimates', cell(1,numel(p.SubjNames)), ...
                  'R2_pRFestimates', cell(1,numel(p.SubjNames)));

             
for s   = 1:numel(p.SubjNames)
    
    fileName    = sprintf('%s_finaldata_nolocSig.mat', p.SubjNames{s});
    
    if exist(fullfile(p.dataDir, 'Data', 'fmriData', fileName), 'file') == 0
        if exist(fullfile(p.dataDir, 'Data', 'fmriData'), 'dir') == 0; mkdir(fullfile(p.dataDir, 'Data', 'fmriData')); end
        
        subjpath    = fullfile(p.maindir, p.SubjNames{s}, sprintf( '%s.sess', p.SubjNames{s}), 'bold');
        extractTimeSeries(subjpath, fileName, p.totalTR, false);
        
    end
    
    load(fullfile(p.dataDir, 'Data', 'fmriData', fileName), 'pRF_maps', 'final_tSeries')
    
    out(s).final_tSeries   = final_tSeries;
    out(s).pRF_maps        = pRF_maps;
    
    switch p.SubjNames{s}
        case {'013', '004', '019'}% For subjects 013 004 019
            
            % Correct pixels per degree used to generate stimuli in the experiment
            newEccen_ring   = (p.newPixPerDegree / p.PixPerDegree) * p.Eccen_ring;
            newLetter_size  = (p.newPixPerDegree / p.PixPerDegree) * p.letter_size;
            
            innerRadius     = newEccen_ring - newLetter_size/2;
            outerRadius     = newEccen_ring + newLetter_size/2;
            PixPerDegree    = p.PixPerDegree;
            
        otherwise
            newEccen_ring   = p.Eccen_ring;
            newLetter_size  = p.letter_size;
            
            innerRadius     = newEccen_ring - newLetter_size/2;
            outerRadius     = newEccen_ring + newLetter_size/2;
            PixPerDegree    = p.newPixPerDegree;
    end
    
    for roi    = 1:numel(p.ROInames)
 
        % Extract pRF information for voxel selection
        ang_pRFestimates    = pRF_maps(roi).ang;
        ecc_pRFestimates    = pRF_maps(roi).ecc;
        rf_pRFestimates     = pRF_maps(roi).rfsize;
        R2_pRFestimates     = pRF_maps(roi).R2;

        out(s).ang_pRFestimates{roi} = ang_pRFestimates;
        out(s).ecc_pRFestimates{roi} = ecc_pRFestimates;
        out(s).rf_pRFestimates{roi}  = rf_pRFestimates;
        out(s).R2_pRFestimates{roi}  = R2_pRFestimates;
        
        % Include voxels whose RF falls within the attended region
        selected_RF         = (abs(ecc_pRFestimates - newEccen_ring) - newLetter_size/1.5) <= rf_pRFestimates;
        selectedpRFVoxels   = selected_RF & ...
            ecc_pRFestimates > p.ecc_inner & ...
            ecc_pRFestimates < p.ecc_outer & ...
            rf_pRFestimates > p.RF_thres & ...
            R2_pRFestimates > p.R2_thres;

        selectedVoxels      = selectedpRFVoxels;
        
        out(s).selectedVoxels{roi} = selectedVoxels;
        
        % Plot visual field coverage of voxel selection
        switch plotLocalizer
            case true
               
                % Set up some params
                gridEdge            = (p.ScreenRes(2)/2);
                gridDownSample      = gridEdge/20;
                xCoords             = -gridEdge:gridDownSample:gridEdge;
                yCoords             = -gridEdge:gridDownSample:gridEdge;
                
                % Compute cartesian coordinates based on pRF estimates
                [cartX, cartY]  = pol2cart(ang_pRFestimates(selectedVoxels), ecc_pRFestimates(selectedVoxels));
                
                % transform to pixels   
                cartX = cartX * PixPerDegree;
                cartY = cartY * PixPerDegree;

                nVoxSel = histcounts2(cartX, cartY, xCoords, yCoords, ...
                                                'Normalization', 'probability'); 
                
                x_innerRadius = (innerRadius * PixPerDegree) * cos(0:pi/50:2*pi);% + paddedRes/2;
                y_innerRadius = (innerRadius * PixPerDegree)* sin(0:pi/50:2*pi);% + paddedRes/2;
                x_outerRadius = (outerRadius * PixPerDegree) * cos(0:pi/50:2*pi);% + paddedRes/2;
                y_outerRadius = (outerRadius * PixPerDegree)* sin(0:pi/50:2*pi);% + paddedRes/2;
                 
                figure(figHandle{roi});               
                nexttile
                imagesc(xCoords, yCoords, nVoxSel)
                hold on,
                plot(x_innerRadius,y_innerRadius, 'w', 'LineWidth', 1)
                plot(x_outerRadius,y_outerRadius, 'w', 'LineWidth', 1)
                set(gca, 'Xtick', linspace(-gridEdge,gridEdge,5), 'XtickLabel', [-9 -4.5 0 4.5 9],...
                    'Ytick', linspace(-gridEdge,gridEdge,5), 'YtickLabel',  [-9 -4.5 0 4.5 9], 'YDir', 'normal')
                box off; axis square; colorbar
                
                title(sprintf('%s ',p.SubjNames{s}))
                
                if saveFig > 0 && s == numel(p.SubjNames)
                    sgtitle(sprintf('Voxel selection ROI %s ', p.ROInames{roi}))

                    print(gcf, fullfile(p.figureDir, 'voxelSelection', sprintf('VisualFieldCoverage_%s', p.ROInames{roi})), '-dpdf'); 
                end
        end

    end
    
end


