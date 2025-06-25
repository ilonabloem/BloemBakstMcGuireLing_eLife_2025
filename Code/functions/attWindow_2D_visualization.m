
function attWindow_2D_visualization(p, data, design, showFigures, saveFig)


saveName            = sprintf('%s_AttWindow_2Drepresentation_%s_visualization.mat', p.Task, p.analysisType);
saveDir             = fullfile(p.dataDir, 'Data', '2Drepresentations');

if ~exist(saveDir, 'dir'), mkdir(saveDir); end
if ~exist(fullfile(p.figureDir, '2Drepresentations'),'dir'), mkdir(fullfile(p.figureDir, '2Drepresentations')); end

% Set up some params
gridEdge            = (p.ScreenRes(2)/2);
gridDownSample      = 2;
xCoords             = -gridEdge:gridDownSample:gridEdge;
yCoords             = -gridEdge:gridDownSample:gridEdge;
    
[gridX,gridY]       = meshgrid(xCoords,yCoords);

%% Compute 2D maps
if exist(fullfile(saveDir, saveName), 'file') > 0

    fprintf('Loading 2D representations ... \n')
    load(fullfile(saveDir, saveName), 'recenteredAverages');

else
    
    fprintf('Computing 2D representations ... \n')

    if (~exist('data', 'var') > 0 || isempty(data)) || ...
        (~exist('design', 'var') > 0 || isempty(design))
        error('Missing data and design structures')
    end

    output              = struct('all_RecenteredData', cell(p.nSubs, p.numROIs), ...
                                 'all_unRecenteredData', cell(p.nSubs, p.numROIs), ...
                                 'TrialEvents', cell(p.nSubs, p.numROIs));
    
    recenteredAverages  = cell(p.numROIs);

    % Loop through all observers
    for s = 1:p.nSubs
    
        fprintf('Participant %s ... \n', p.SubjNames{s})
        
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
        
        % Task design:
        trueMean            = design(s).trueMean;
        trueWidth           = design(s).trueWidth;

        % Loop through regions of interest
        for roi = 1:p.numROIs
            
            switch p.Task
                case 'CM'
                    % First 2 runs are physical contrast manipulation - rest is task data
                    indx                = 1:p.totalTR:(p.numRuns(s)+2)*p.totalTR+1;
                    TSeries             = data(s).final_tSeries{roi}(:,indx(3):indx(end)-1);
                    
                case 'PCM'
                    % First 2 runs are physical contrast manipulation - rest is task data
                    indx                = 1:p.totalTR:(p.numRuns(s))*p.totalTR+1;
                    TSeries             = data(s).final_tSeries{roi}(:,indx(1):indx(3)-1);
            end
            
            switch p.analysisType
                case 'indvTrials'
                    whichTimePoints     = 1:size(TSeries,2);
                case 'blockAvg'
                    whichTimePoints     = (1:10:size(TSeries,2)-10)+p.timeLag;
            end
            
            if s == 1
                recenteredAverages{roi} = NaN(p.numAttWindows,p.nSubs,size(gridX,1),size(gridX,2));
            end
            
            % Selected voxels (whole visual field)
            ang_pRFestimates    = data(s).ang_pRFestimates{roi};
            ecc_pRFestimates    = data(s).ecc_pRFestimates{roi};
            rf_pRFestimates     = data(s).rf_pRFestimates{roi};
            R2_pRFestimates     = data(s).R2_pRFestimates{roi};
            selectedVoxels      = ecc_pRFestimates > p.ecc_inner & ...
                                  ecc_pRFestimates < p.ecc_outer & ...
                                  rf_pRFestimates > p.RF_thres & ...
                                  R2_pRFestimates > p.R2_thres; 
            
            % Compute cartesian coordinates based on pRF estimates
            [cartX, cartY] = pol2cart(ang_pRFestimates(selectedVoxels), ecc_pRFestimates(selectedVoxels));

            % Add uniform noise to x and y coordinates so that there are only unique combinations        
            cartX = (cartX * PixPerDegree) + (-0.1 + (.1--.1) .* rand(size(cartX)));
            cartY = (cartY * PixPerDegree) + (-0.1 + (.1--.1) .* rand(size(cartY)));

            % Preallocate some variables
            all_RecenteredImages = NaN(size(gridX,1), size(gridX,2), length(whichTimePoints));
            all_UnrecenteredImages = NaN(size(gridX,1), size(gridX,2), length(whichTimePoints));

            tmpMean = [trueMean(whichTimePoints) trueWidth(whichTimePoints)];
            counter = 0;
            
            for tp = whichTimePoints
            
                counter = counter + 1;

                switch p.analysisType
                    case 'indvTrials'
                        select_TSeries = TSeries(selectedVoxels,tp);
                    case 'blockAvg'
                        select_TSeries = mean(TSeries(selectedVoxels,tp:tp+9),2);
                end
 
                % Interpolate activity at grid points based on BOLD activity
                gridFit = griddata(cartX,cartY,select_TSeries,gridX,gridY, 'linear');
                
                % Z-score representation
                normGridFit     =  (gridFit - nanmean(gridFit(:))) ./ nanstd(gridFit(:));

                all_UnrecenteredImages(:,:, counter) = normGridFit;

                if ~isnan(tmpMean(counter,1))
                    if tmpMean(counter,1) == 0
                        all_RecenteredImages(:,:,counter)  = normGridFit;
                    else
                        all_RecenteredImages(:,:,counter)  = imrotate(normGridFit, rad2deg(tmpMean(counter,1)), 'nearest', 'crop');
                    end
                end

            end
            
            % Average for each width condition:
            for ii=1:p.numAttWindows
               
                whichTRs                        = tmpMean(:,2) == p.WidthAttWindow(ii);
                recenteredAverages{roi}(ii,s,:,:)    = squeeze(nanmean(all_RecenteredImages(:,:,whichTRs),3));

            end
            
            
            output(s,roi).all_RecenteredData = all_RecenteredImages;
            output(s,roi).all_unRecenteredData = all_UnrecenteredImages;
            output(s,roi).TrialEvents        = tmpMean;

        end
        
    end
        
    save(fullfile(saveDir, saveName), 'output', 'recenteredAverages', '-v7.3');

end


%% Visualize 2D maps
if showFigures
        
    innerRadius     = p.Eccen_ring - p.letter_size/2;
    outerRadius     = p.Eccen_ring + p.letter_size/2;
    PixPerDegree    = p.newPixPerDegree;

    x_innerRadius   = (innerRadius * PixPerDegree) * cos(0:pi/50:2*pi);
    y_innerRadius   = (innerRadius * PixPerDegree)* sin(0:pi/50:2*pi);
    x_outerRadius   = (outerRadius * PixPerDegree) * cos(0:pi/50:2*pi);
    y_outerRadius   = (outerRadius * PixPerDegree)* sin(0:pi/50:2*pi);
       
    eccenMap        = sqrt(gridX.^2 + gridY.^2);
    
    mask            =  eccenMap <= (p.ecc_outer * PixPerDegree);
    
    colorBounds     = [-.5 0.9; ...
                        -.5 0.9; ...
                        -.5 0.9];

    idx             = 1:p.numAttWindows:p.numAttWindows*p.numROIs;
    figGroupavg = figure('Color', [1 1 1]);
    for roi = 1:p.numROIs

        % Group figures
        for ii = 1:p.numAttWindows
            figure(figGroupavg)
            subplot(p.numROIs,p.numAttWindows,ii+(idx(roi)-1))
            img = squeeze(nanmean(recenteredAverages{roi}(ii,:,:,:),2)) .* mask;
            img(isnan(img)) = 0;
            
            imagesc(xCoords/PixPerDegree,xCoords/PixPerDegree, ...
                    img)
            hold on, 
            plot(x_innerRadius/PixPerDegree,y_innerRadius/PixPerDegree, 'k', 'LineWidth', 1)
            plot(x_outerRadius/PixPerDegree,y_outerRadius/PixPerDegree, 'k', 'LineWidth', 1)
            box off
            axis square
            caxis(colorBounds(roi,:))
            if roi == 1
               title(sprintf('Attentional window %i', ii))               
            end
            
        end
        subplot(p.numROIs,p.numAttWindows,1)
        ylabel('Eccentricity (°)'); xlabel('Eccentricity (°)')
       
    end
       
    cbh             = colorbar;
    cbh.Position(3) = cbh.Position(3)*2;
    cbh.Position(4) = cbh.Position(4)*2;
    cbh.Position(1) = .95-cbh.Position(3);
    cbh.Position(2) = 0.5-cbh.Position(4)/2;
    
    sgtitle(sprintf('%s recentered group average, %s estimates', p.Task, p.analysisType))
    if saveFig
       
        figure(figGroupavg)
        print(fullfile(p.figureDir, '2Drepresentations', sprintf('%s_groupAvg_2Drepresentation_%s', p.Task, p.analysisType)), '-dpdf', '-vector', '-bestfit')
        
    end
    
    % Individual participant figures
    %{
    for s = 1:p.nSubs
         
        figure('Color', [1 1 1])
        
        for roi = 1:p.numROIs
            for ii   = 1:p.numAttWindows

                img = squeeze(recenteredAverages{roi}(ii,s,:,:)) .* mask;
                img(isnan(img)) = 0;
                
                subplot(p.numROIs,p.numAttWindows,ii+(idx(roi)-1))
                imagesc(xCoords/PixPerDegree,xCoords/PixPerDegree, ...
                    img)
                hold on,
                plot(x_innerRadius/PixPerDegree, y_innerRadius/PixPerDegree, 'k', 'LineWidth', 1)
                plot(x_outerRadius/PixPerDegree, y_outerRadius/PixPerDegree, 'k', 'LineWidth', 1)

                %set(gca, 'Xtick', linspace(1,ScreenRes(1),5), 'XtickLabel', [-16 -8 0 8 16],...
                %    'Ytick', linspace(1,ScreenRes(1),5), 'YtickLabel', [-16 -8 0 8 16])
                box off; axis square
                caxis(colorBounds(roi,:))
                if roi == 1
                   title(sprintf('Attentional window %i', ii))

                end
            end

        end
        subplot(p.numROIs,p.numAttWindows,1)
        ylabel('Eccentricity (°)'); xlabel('Eccentricity (°)')
            
        sgtitle(sprintf('%s %s recentered average, %s estimates', p.Task, p.SubjNames{s}, p.analysisType))

        if saveFig

            print(fullfile(p.figureDir, '2Drepresentations', sprintf('%s_%s_2Drepresentation_%s', p.Task,  p.SubjNames{s},p.analysisType)), '-dpdf', '-painters', '-bestfit')

        end
        
    end
    %}
    
end

