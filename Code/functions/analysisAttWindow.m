
function [ModelEstimates, Summary] = analysisAttWindow(p, out, design)

%% Attentional window analysis
%
% function out = analysisAttWindow(p, out, design)
%

%% Set up model fitting
grid_Locparams      = linspace(0,2*pi-(2*pi/6), 6); %attLocations;
grid_Widthparams    = linspace(deg2rad(360/p.numLocations/2), deg2rad((360/p.numLocations)*9), 6);
grid_Widthparams    = grid_Widthparams(1:end-1);

genGVals            = fitGenGauss([], {'initialize'});
numParams           = numel(genGVals.labels);
zscoreData          = false;

%%
switch p.doSmooth
    case true
        saveName    = sprintf('%s_smoothAttWindow_genGauss_%s.mat', p.Task, p.analysisType);
    case false
        saveName    = sprintf('%s_nosmoothAttWindow_genGauss_%s.mat', p.Task, p.analysisType);
end

saveDir         = fullfile(p.saveDir, p.Task);
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

%% Only compute if file doesn't exist
if exist(fullfile(saveDir, saveName), 'file') > 0 && p.recompute == 0
    
    results         = load(fullfile(saveDir, saveName), 'ModelEstimates', 'Summary');
    ModelEstimates  = results.ModelEstimates;
    Summary         = results.Summary;
    
    fprintf('Load model results: %s %s \n', p.Task, p.analysisType)
    
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
            ang_data        = location_ang(1:numVoxels);
            x_data          = p.binCenters(:);
            whichVoxels     = order_ang(1:numVoxels);
            sp_data         = linspace(-pi, 3*pi - (4*pi/(p.nBins*2)), p.nBins*2)';
            
            counter     = 0;
            switch p.analysisType
                case 'indvTrials'
                    whichTimePoints     = 1:size(allTSeries,2);
                case 'blockAvg'
                    whichTimePoints     = (1:10:size(allTSeries,2)-10)+p.timeLag;
                    blockAccuracy       = NaN(size(whichTimePoints));
            end
            
            R2                      = NaN(length(whichTimePoints), 1);
            widthEst_deg            = NaN(length(whichTimePoints), 1);
            widthEst_FWHM_deg       = NaN(length(whichTimePoints), 1);
            SSE                     = NaN(length(whichTimePoints), 1);
            exitflag                = NaN(length(whichTimePoints), 1);
            est_params              = NaN(length(whichTimePoints), numParams);
            startVals               = NaN(length(whichTimePoints), numParams);
            all_RecenteredData      = NaN(length(whichTimePoints), p.nBins);
            all_UnrecenteredData    = NaN(length(whichTimePoints), p.nBins);
            all_trueMean_RecenteredData     = NaN(length(whichTimePoints), p.nBins);
            
            for tp = whichTimePoints
                
                counter     = counter + 1;
                
                switch zscoreData
                    case true
                        
                        % Z-score responses for each voxel across trials
                        % Then average over window
                        switch p.analysisType
                            case 'indvTrials'
                                
                                select_TSeries          = normalize(allTSeries(whichVoxels,tp), 'zscore', 'std');
                                blockAccuracy(counter)  = trialAccuracy(tp,1);
                            case 'blockAvg'
                                
                                tmpData                 = NaN(size(whichVoxels,1), 10);
                                % Select TRs to average over for all events
                                for n = 1:10
                                    tmpData(:,n)        = normalize(allTSeries(whichVoxels,(tp+n)-1), 'zscore', 'std');
                                end
                                
                                select_TSeries          = nanmean(tmpData,2);
                                blockAccuracy(counter)  = nanmean(trialAccuracy(tp:2:tp+9,1));
                                
                        end
                    case false
                        
                        switch p.analysisType
                            case 'indvTrials'
                                select_TSeries          = allTSeries(whichVoxels,tp);
                                blockAccuracy(counter)  = trialAccuracy(tp,1);
                            case 'blockAvg'
                                select_TSeries          = mean(allTSeries(whichVoxels,tp:tp+9),2);
                                blockAccuracy(counter)  = nanmean(trialAccuracy(tp:2:tp+9,1));
                        end
                        
                end
                
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
                        smoothTSeries   = tmp_data(smoothWindow+1:length(binned_data)+smoothWindow);
                        TSeries         = smoothTSeries(:);
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
                startVals(counter,:)    = modelInit;
                
                fixedParams{1} = 'optimize';
                [est_params(counter,:), SSE(counter), exitflag(counter)] = ...
                    fmincon(@(x) fitGenGauss(x, fixedParams), startVals(counter,:), [], [], [], [], genGVals.lb, genGVals.ub, [],  genGVals.opt);
                
                fixedParams{1} = 'prediction';
                genGausResults = fitGenGauss(est_params(counter,:), fixedParams);
                
                R2(counter,1)  = genGausResults.R2;
                
                %% recentering data
                
                % List all the edges that are greater than the estimated angle
                binsGreater = find(wrapTo2Pi(est_params(counter,1)) <= [p.binEdges 2*pi]);
                % find the bin just before the smallest edge in our list
                meanBin     = binsGreater(1) - 1;
                % Find the difference between the bin in which the estimated
                % angle falls and pi (edges(30:31))
                recenter_mean = meanBin - round(p.nBins/2);
                % Make sure indx is wrapped to nBins
                indx = (1:p.nBins) + recenter_mean(1);
                positiveInput = (indx >= 0);
                indx = mod(indx, p.nBins);
                indx((indx == 0) & positiveInput) = p.nBins;
                
                all_RecenteredData(counter,:) = fixedParams{3}(indx);
                all_UnrecenteredData(counter,:) = fixedParams{3};
                
                if ~isnan(trueMean(tp))
                    recenter_trueMean = (find(wrapTo2Pi(trueMean(tp)) < p.binEdges)-1) - round(p.nBins/2);
                    indx_trueMean = (1:p.nBins) + recenter_trueMean(1);
                    positiveInput = (indx_trueMean >= 0);
                    indx_trueMean = mod(indx_trueMean, p.nBins);
                    indx_trueMean((indx_trueMean == 0) & positiveInput) = p.nBins;
                    all_trueMean_RecenteredData(counter,:) = fixedParams{3}(indx_trueMean);
                end
                
                % Convert width estimate from radians to deg
                widthEst_deg(counter,:) = rad2deg(est_params(counter,2));
                
                % Generate von mises centered on pi to compute FWHM
                %             tmp_x           = linspace(0,2*pi-(2*pi/359), 359)';
                %             tmp_sp          = round(linspace(-pi + (4*pi/(359*4)), 3*pi - (4*pi/(359*4)), 359*2),4);
                tmpParams       = est_params(counter,:);
                tmpParams(1)    = pi;
                
                fixedParams{1}  = 'prediction';
                fit             = fitGenGauss(tmpParams, fixedParams);
                
                find_halfMax = fit.y_est >= (max(fit.y_est)-min(fit.y_est))/2 + min(fit.y_est);
                halfMax_Index = find(diff(find_halfMax)~=0);
                if length(halfMax_Index) < 2
                    widthEst_FWHM_deg(counter,:) = NaN;
                else
                    widthEst_FWHM_deg(counter,:) = rad2deg(x_data(halfMax_Index(2)) - x_data(halfMax_Index(1)));
                end
                
            end
            
            
            %% Plot the estimated parameters
            notOutliers             = exitflag ~= 0;
            trueMean_plot           = trueMean(whichTimePoints);
            trueWidth_plot          = trueWidth(whichTimePoints);
            
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
                Summary(roi).AngleEst(subj,ii,:)                        = histcounts(diffAng(whichWidth), bin_edges);
                
                if strcmp(p.analysisType, 'blockAvg')
                    Summary(roi).blockAccuracy(subj,ii,:)                   = nanmean(blockAccuracy(whichWidth));
                end
                
            end
            
            Summary(roi).RecenteredData{subj}           = all_RecenteredData';
            Summary(roi).TrueMeanRecenteredData{subj}   = all_trueMean_RecenteredData';
            
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


