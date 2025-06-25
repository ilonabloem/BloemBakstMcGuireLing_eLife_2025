
function [ModelEstimates, Summary] = analysisAttWindowTR(p, out, design, showFigures, saveFig)

%% Attentional window analysis
%
% function out = analysisAttWindowTR(p, out, design, showFigures, saveFig)
%
% IB - 2025


%% Set up model fittings
grid_Locparams      = linspace(0,2*pi-(2*pi/6), 6); %attLocations;
grid_Widthparams    = linspace(deg2rad(360/p.numLocations/2), deg2rad((360/p.numLocations)*9), 6);
grid_Widthparams    = grid_Widthparams(1:end-1);

genGVals            = fitGenGauss([], {'initialize'});
zscoreData          = false;

%%
switch p.doSmooth
    case true
        saveName    = sprintf('%s_%s_genGauss_%iTR.mat', p.Task, mfilename, p.nTR);
    case false
        saveName    = sprintf('%s_%s_nosmooth_genGauss_%iTR.mat', p.Task, mfilename, p.nTR);
end

saveDir         = fullfile(p.saveDir, p.Task);
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

%% Only compute if file doesn't exist
if exist(fullfile(saveDir, saveName), 'file') > 0 && p.recompute == 0
    
    results         = load(fullfile(saveDir, saveName), 'ModelEstimates', 'Summary');
    ModelEstimates  = results.ModelEstimates;
    Summary         = results.Summary;
    
    fprintf('Load model results: %s %s \n', p.Task, num2str(p.nTR))
    
else
    
    % Preallocate variables
    Summary             = struct('all_RecenteredData', cell(1,numel(p.ROInames)), ...
        'all_trueMean_RecenteredData', cell(1,numel(p.ROInames)), ... 
        'angleEst', cell(1,numel(p.ROInames)), ...
        'RecenteredData', cell(1,numel(p.ROInames)), ...
        'TrueMeanRecenteredData', cell(1,numel(p.ROInames)), ...
        'avg_errorByLocation', cell(1,numel(p.ROInames)), ...
        'ErrorVerticalMeridian', cell(1,numel(p.ROInames)), ...
        'ErrorHorizontalMeridian', cell(1,numel(p.ROInames)));
    ModelEstimates      = struct('cuedLocation', cell(1,numel(p.ROInames)), ...
        'widthCondition', cell(1,numel(p.ROInames)), ...
        'trialAccuracy', cell(1,numel(p.ROInames)), ...
        'angleEstimates', cell(1,numel(p.ROInames)),...
        'FWHMEstimates', cell(1,numel(p.ROInames)), ...
        'sdEstimates', cell(1,numel(p.ROInames)), ...
        'betaEstimates', cell(1,numel(p.ROInames)), ...
        'amplEstimates', cell(1,numel(p.ROInames)), ...
        'baseEstimates', cell(1,numel(p.ROInames)), ...
        'diffAng', cell(1,numel(p.ROInames)), ...
        'SSE', cell(1,numel(p.ROInames)), ...
        'R2', cell(1,numel(p.ROInames)), ...
        'startVals', cell(1,numel(p.ROInames)), ...
        'blockAccuracy', cell(1,numel(p.ROInames)), ...
        'notOutliers', cell(1,numel(p.ROInames)));
    
    % Loop through all participants
    for subj = 1:numel(p.SubjNames)
        
        trueMean            = design(subj).trueMean;
        trueWidth           = design(subj).trueWidth;
        trialAccuracy       = design(subj).trialAccuracy;
        
        ModelEstimates(1).trialAccuracy{subj} = trialAccuracy;
        
        % loop through ROIs
        for roi = 1:numel(p.ROInames)
            
            fprintf('Run model for participant: %s, ROI: %s \n', p.SubjNames{subj}, p.ROInames{roi})
            
            if subj == 1
                Summary(roi).all_RecenteredData             = NaN(numel(p.SubjNames), p.numAttWindows, p.nBins);
                Summary(roi).all_trueMean_RecenteredData    = NaN(numel(p.SubjNames), p.numAttWindows, p.nBins);
                Summary(roi).angleEst                       = NaN(numel(p.SubjNames), p.numAttWindows, 20);
            end
            
            switch p.Task
                case 'CM'
                    % First 2 runs are physical contrast manipulation - rest is task data
                    indx                = 1:p.totalTR:(p.numRuns(subj)+2)*p.totalTR+1;
                    allTSeries          = out(subj).final_tSeries{roi}(:,indx(3):indx(end)-1);
                    
                case 'PCM'
                    % First 2 runs are physical contrast manipulation - rest is task data
                    indx                = 1:p.totalTR:(p.numRuns(subj))*p.totalTR+1;
                    allTSeries          = out(subj).final_tSeries{roi}(:,indx(1):indx(3)-1);
            end
            
            
            % order the voxels based on polar angle
            selectedVoxels              = out(subj).selectedVoxels{roi};
            coordinates                 = NaN(size(selectedVoxels));
            coordinates(selectedVoxels) = out(subj).ang_pRFestimates{roi}(selectedVoxels);
            
            numVoxels                   = sum(selectedVoxels);
            [location_ang, order_ang]   = sort(coordinates);
            
            %% Fit the data
            ang_data                = location_ang(1:numVoxels);
            x_data                  = p.binCenters(:);
            whichVoxels             = order_ang(1:numVoxels);
            sp_data                 = linspace(-pi, 3*pi - (4*pi/(p.nBins*2)), p.nBins*2)';
            counter                 = 0;
            cond_TSeries            = cell(p.numLocations, p.numAttWindows);
            R2                      = [];
            widthEst_deg            = [];
            widthEst_FWHM_deg       = [];
            SSE                     = [];
            exitflag                = [];
            est_params              = [];
            startVals               = [];
            all_RecenteredData      = [];
            all_UnrecenteredData    = [];
            all_trueMean_RecenteredData     = [];
            trueMean_plot           = [];
            trueWidth_plot          = [];
            blockAccuracy           = [];
            
            for loc = 1:p.numLocations
                for wid = 1:p.numAttWindows
                    
                    whichTRs                = (trueMean == p.attLocations(loc)) & (trueWidth == p.WidthAttWindow(wid));
                    cond_TSeries{loc,wid}   = allTSeries(whichVoxels,whichTRs);
                    cond_Accuracy           = trialAccuracy(whichTRs);
                    maxTR                   = size(cond_TSeries{loc,wid},2);
                    nTR                     = min([maxTR,p.nTR]);
                    nBlockTRs               = maxTR/p.numTRblock;
                    whichTimePoints         = 1:nBlockTRs;
 
                    for tp = whichTimePoints

                        counter     = counter + 1;
                        
                        %choose nTR number from the condition block
                        if nTR==p.numTRblock
                            selectTR = (1:p.numTRblock)+((tp-1)*p.numTRblock);
                        else
%                             selectTR = randperm(maxTR,nTR);
                            selectTR = (1:p.numTRblock)+((tp-1)*p.numTRblock);
                            selectTR = selectTR(randperm(p.numTRblock, nTR));
                        end
                        
                        switch zscoreData
                            case true
                                
                                % Z-score responses for each voxel across trials  - suggested by Ruben
                                % Then average over window
                                  
                                tmpData                 = NaN(size(whichVoxels,1), nTR);
                                for n = 1:nTR
                                    tmpData(:,n)        = normalize(cond_TSeries{loc,wid}(:,selectTR(n)), 'zscore', 'std');
                                end
                                select_TSeries          = nanmean(tmpData,2);     
                                
                            case false

                                if numel(selectTR) > 1
                                    select_TSeries              = mean(cond_TSeries{loc,wid}(:,selectTR),2);
                                else
                                    select_TSeries              = cond_TSeries{loc,wid}(:,selectTR);
                                end
                                
                        end

                        blockAccuracy      = cat(1, blockAccuracy, nanmean(cond_Accuracy(selectTR)));
                        trueMean_plot      = cat(1, trueMean_plot, p.attLocations(loc));
                        trueWidth_plot     = cat(1, trueWidth_plot, p.WidthAttWindow(wid));
                        
                        % Bin data
                        tmp_edges = cat(2, p.binEdges(end), p.binEdges);
                        binned_data = NaN(1, p.nBins);
                        for ii = 1:p.nBins
                            if ii == 1
                                tmpIndx = ang_data >= tmp_edges(ii) | ang_data < tmp_edges(ii+1);

                                binned_data(ii) = nanmedian(select_TSeries(tmpIndx,1));
                            else
                                if sum(ang_data >= tmp_edges(ii) & ang_data < tmp_edges(ii+1)) >= 1
                                    binned_data(ii) = nanmedian(select_TSeries(ang_data >= tmp_edges(ii) & ang_data < tmp_edges(ii+1)));
                                end
                            end
                        end

                        switch p.doSmooth
                            case true
                                smoothWindow    = 3; 
                                tmp_data        = cat(2, binned_data(end-smoothWindow+1:end), binned_data, binned_data(1:smoothWindow));
                                tmp_data        = smoothdata(tmp_data, "movmean",smoothWindow);
                                smooth_TSeries  = tmp_data(smoothWindow+1:length(binned_data)+smoothWindow);
                                TSeries         = smooth_TSeries(:);
                        case false
                                TSeries         = binned_data(:);
                        end
                        
                        %grid search for mu
                        sseGrid         = zeros(numel(grid_Widthparams),numel(grid_Locparams));
                        ampGrid         = zeros(numel(grid_Widthparams),numel(grid_Locparams));
                        baseGrid        = zeros(numel(grid_Widthparams),numel(grid_Locparams));
                        fixedParams     = cell(3,1);
                        fixedParams{1}  = 'prediction';
                        fixedParams{2}  = sp_data;
                        fixedParams{3}  = TSeries(:);

                        for grdW = 1:size(sseGrid,1)
                            for grdA = 1:size(sseGrid,2)

                                modelInit   = [grid_Locparams(grdA), ...
                                    grid_Widthparams(grdW), ...
                                    genGVals.init(3), ...
                                    genGVals.init(4), ...
                                    genGVals.init(5)];

                                fit         = fitGenGauss(modelInit, fixedParams);

                                % Scale model (baseline and amplitude)
                                model               = cat(2, fit.y_est, ones(size(fit.y_est)));
                                modelInv            = pinv(model);
                                notNaN              = ~isnan(fixedParams{3});
                                b_grd               = modelInv(:,notNaN) * fixedParams{3}(notNaN);

                                yEst                = (model * b_grd);
                                sseGrid(grdW, grdA) = nansum((yEst-fixedParams{3}).^2);
                                ampGrid(grdW, grdA) = b_grd(1);
                                baseGrid(grdW, grdA)= b_grd(2);

                            end
                        end

                        % find which start param to seed for optimization
                        sseGrid(ampGrid < 0)    = Inf; % Amplitude should be postive
                        [~, J]                  = sort(sseGrid(:), 'ascend');
                        [x, y]                  = ind2sub(size(sseGrid), J(1));
                        modelInit(1)            = grid_Locparams(y);
                        modelInit(2)            = grid_Widthparams(x);
                        modelInit(3)            = genGVals.init(3);
                        modelInit(4)            = ampGrid(x,y);
                        modelInit(5)            = baseGrid(x,y);

                        %%% Now use the best mu as start param for real fitting procedure
                        startVals       = cat(1, startVals, modelInit);

                        fixedParams{1}  = 'optimize';
                        [est_params(counter,:), SSE(counter), exitflag(counter)] = ...
                            fmincon(@(x) fitGenGauss(x, fixedParams), startVals(counter,:), [], [], [], [], genGVals.lb, genGVals.ub, [],  genGVals.opt);

                        fixedParams{1}  = 'prediction';
                        genGausResults  = fitGenGauss(est_params(counter,:), fixedParams);

                        R2              = cat(1, R2, genGausResults.R2);

                        %% recentering data

                        % List all the edges that are greater than the estimated angle
                        binsGreater             = find(wrapTo2Pi(est_params(counter,1)) <= [p.binEdges 2*pi]);
                        % find the bin just before the smallest edge in our list
                        meanBin                 = binsGreater(1) - 1;
                        % Find the difference between the bin in which the estimated
                        % angle falls and pi (edges(30:31))
                        recenter_mean           = meanBin - round(p.nBins/2);
                        % Make sure indx is wrapped to nBins
                        indx                    = (1:p.nBins) + recenter_mean(1);
                        positiveInput           = (indx >= 0);
                        indx                    = mod(indx, p.nBins);
                        indx((indx == 0) & positiveInput) = p.nBins;

                        all_RecenteredData      = cat(1, all_RecenteredData, fixedParams{3}(indx)');
                        all_UnrecenteredData    = cat(1, all_UnrecenteredData, fixedParams{3}');


                        recenter_trueMean       = (find(wrapTo2Pi(p.attLocations(loc)) < p.binEdges)-1) - round(p.nBins/2);
                        indx_trueMean           = (1:p.nBins) + recenter_trueMean(1);
                        positiveInput           = (indx_trueMean >= 0);
                        indx_trueMean           = mod(indx_trueMean, p.nBins);
                        indx_trueMean((indx_trueMean == 0) & positiveInput) = p.nBins;
                        all_trueMean_RecenteredData = cat(1, all_trueMean_RecenteredData, fixedParams{3}(indx_trueMean)');

                        % Convert width estimate from radians to deg
                        widthEst_deg            = cat(1, widthEst_deg, rad2deg(est_params(counter,2)));

                        % Generate von mises centered on pi to compute FWHM
                        %             tmp_x           = linspace(0,2*pi-(2*pi/359), 359)';
                        %             tmp_sp          = round(linspace(-pi + (4*pi/(359*4)), 3*pi - (4*pi/(359*4)), 359*2),4);
                        tmpParams               = est_params(counter,:);
                        tmpParams(1)            = pi;

                        fixedParams{1}          = 'prediction';
                        fit                     = fitGenGauss(tmpParams, fixedParams);

                        find_halfMax            = fit.y_est >= (max(fit.y_est)-min(fit.y_est))/2 + min(fit.y_est);
                        halfMax_Index           = find(diff(find_halfMax)~=0);
                        if length(halfMax_Index) < 2
                            widthEst_FWHM_deg   = cat(1, widthEst_FWHM_deg, NaN);
                        else
                            widthEst_FWHM_deg   = cat(1, widthEst_FWHM_deg, rad2deg(x_data(halfMax_Index(2)) - x_data(halfMax_Index(1))));
                        end
                       
                        
                    end
            
                end
            end
            %% Plot the estimated parameters
            notOutliers             = exitflag ~= 0;      
            diffAng                 = circ_dist(trueMean_plot, est_params(:,1));
            
            
            ModelEstimates(roi).angleEstimates{subj}    = est_params(:,1);
            ModelEstimates(roi).FWHMEstimates{subj}     = widthEst_FWHM_deg(:);
            ModelEstimates(roi).sdEstimates{subj}       = widthEst_deg(:);
            ModelEstimates(roi).betaEstimates{subj}     = 10.^est_params(:,3);
            ModelEstimates(roi).amplEstimates{subj}     = est_params(:,4);
            ModelEstimates(roi).baseEstimates{subj}     = est_params(:,5);
            ModelEstimates(roi).SSE{subj}               = SSE(:);
            ModelEstimates(roi).R2{subj}                = R2(:);
            
            ModelEstimates(roi).cuedLocation{subj}      = trueMean_plot(:);
            ModelEstimates(roi).widthCondition{subj}    = trueWidth_plot(:);
            ModelEstimates(roi).diffAng{subj}           = diffAng(:);
            ModelEstimates(roi).notOutliers{subj}       = notOutliers;
            ModelEstimates(roi).blockAccuracy{subj}     = blockAccuracy(:);
            ModelEstimates(roi).startVals{subj}         = startVals;
            
            
            bin_edges               = linspace(-pi,pi,21);
            for ii = 1:p.numAttWindows
                whichWidth                                              = trueWidth_plot == p.WidthAttWindow(ii);
                Summary(roi).all_RecenteredData(subj,ii,:)              = nanmean(all_RecenteredData(whichWidth,:));
                Summary(roi).all_trueMean_RecenteredData(subj,ii,:)     = nanmean(all_trueMean_RecenteredData(whichWidth,:));
                Summary(roi).angleEst(subj,ii,:)                        = histcounts(diffAng(whichWidth), bin_edges);
                Summary(roi).blockAccuracy(subj,ii,:)                   = nanmean(blockAccuracy(whichWidth));

                
            end
            
            Summary(roi).RecenteredData{subj}           = all_RecenteredData';
            Summary(roi).TrueMeanRecenteredData{subj}   = all_trueMean_RecenteredData';
            
            switch showFigures
                case true
                    
                    figure('color', [1 1 1], 'pos', [400 10 600 1000]),
                    errorRange          = cell(4,1);
                    
                    for ii = 1:numel(p.WidthAttWindow)
                        subplot(2,1,1)
                        whichWidth = find(trueWidth_plot(:) == p.WidthAttWindow(ii));
                        which_TRs = whichWidth(:);
                        hold on,
                        plot(p.binCenters, squeeze(Summary(roi).smooth_RecenteredData(subj,ii,:))),
                        
                        errorRange{ii} = [num2str(round(min(diffAng(which_TRs)),3)), ' - ', num2str(round(max(diffAng(which_TRs)),3))];
                        
                    end
                    subplot(2,1,1)
                    xlabel('Angle (Rad)'), ylabel('BOLD response'); title([{'Continuous monitoring'};{'Recentered to est mean'}]), legend({'1' '3' '5' '9'})
                    set(gca, 'XTick', [0 pi 2*pi], 'XTickLabel', {'0' 'pi' '2pi'}); xlim([0 2*pi])
                    
                    subplot(2,1,2)
                    for ii = 1:numel(p.WidthAttWindow)
                        hold on,
                        plot(p.binCenters, squeeze(Summary(roi).smooth_trueMean_RecenteredData(subj,ii,:))),
                    end
                    %         hold on, plot([ones(1,2).*attLocations']', [ones(20,1).*[-0.4 1.4]]', 'k', 'HandleVisibility', 'off')
                    xlabel('Angle (Rad)'), ylabel('BOLD response'); title('Recentered to true mean')
                    set(gca, 'XTick', [0 pi 2*pi], 'XTickLabel', {'0' 'pi' '2pi'}); xlim([0 2*pi])
                    
                    if saveFig
                        saveas(gcf, fullfile(figureDir, sprintf('%s_ActivityProfile_S%s_%s.png', p.Task, p.SubjNames{subj}, p.ROInames{roi})), 'png'); 
                    end
                    close;
            end
            
            avg_errorByLocation = NaN(1, p.numLocations);
            
            for kk = 1:p.numLocations
                whichLocation = trueMean_plot == p.attLocations(kk);
                if sum(whichLocation) > 0
                    avg_errorByLocation(kk) = circ_mean(diffAng(whichLocation));
                end
            end
            
            
            Summary(roi).avg_errorByLocation(subj,:) = avg_errorByLocation;
            Summary(roi).ErrorVerticalMeridian(subj,:) = nanmean([avg_errorByLocation(1,4:8); avg_errorByLocation(1,14:18)]);
            Summary(roi).ErrorHorizontalMeridian(subj,:) = nanmean([avg_errorByLocation(1,9:13); [avg_errorByLocation(1,19:end) avg_errorByLocation(1,1:3)]]);
            
        end
        
    end
    p.date = date;
    save(fullfile(saveDir, saveName), 'ModelEstimates', 'Summary', 'p')
end


