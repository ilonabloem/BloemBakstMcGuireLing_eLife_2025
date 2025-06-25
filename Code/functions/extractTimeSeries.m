
function extractTimeSeries(subjpath, fileName, numTRS, runLoc)

PlotFigures = 0;

scanNames = dir(fullfile(subjpath, '0*'));
% load ROI volume masks
% Retinotopy labels
masks_dir = dir(fullfile(subjpath, 'pRF_VOLs', 'pRF.v*'));
v1_mask = MRIread(fullfile(subjpath, 'pRF_VOLs', masks_dir(1).name));
v2_mask = MRIread(fullfile(subjpath, 'pRF_VOLs', masks_dir(2).name));
v3_mask = MRIread(fullfile(subjpath, 'pRF_VOLs', masks_dir(3).name));
roi_maskNames = {'v1_mask' 'v2_mask' 'v3_mask'};
%v4_dir = dir(fullfile(subjpath, 'pRF_VOLs', 'pRF.h*'));
%if ~isempty(v4_dir)
%    v4_mask = MRIread(fullfile(subjpath, 'pRF_VOLs' v4_dir(1).name));
%    roi_maskNames{4} = 'v4_mask';
%end


final_tSeries = cell(1,numel(roi_maskNames));
linearCoordinates = cell(1,numel(roi_maskNames));
localizer_sig = cell(1,numel(roi_maskNames));

% Localizer sig volume
if PlotFigures > 0
    figure('color', [1 1 1])
    lVF = colormap('autumn');
    rVF = colormap('winter');
end
if runLoc == true
    % load Significance maps computed using localizer runs
    sig_dir = dir([subjpath 'pRF_VOLs/sig*']);
    if ~isempty(sig_dir)
        loc_mask = MRIread([subjpath 'pRF_VOLs/sig.VOL.nii.gz']);
    end
end

% load pRF maps
rh_ecc = MRIread([subjpath 'pRF_VOLs/rh.pRF_ecc.VOL.nii.gz']);
lh_ecc = MRIread([subjpath 'pRF_VOLs/lh.pRF_ecc.VOL.nii.gz']);

rh_ang = MRIread([subjpath 'pRF_VOLs/rh.pRF_ang.VOL.nii.gz']);
lh_ang = MRIread([subjpath 'pRF_VOLs/lh.pRF_ang.VOL.nii.gz']);

rh_R2 = MRIread([subjpath 'pRF_VOLs/rh.pRF_R2.VOL.nii.gz']);
lh_R2 = MRIread([subjpath 'pRF_VOLs/lh.pRF_R2.VOL.nii.gz']);

rh_rfsize = MRIread([subjpath 'pRF_VOLs/rh.pRF_rfsize.VOL.nii.gz']);
lh_rfsize = MRIread([subjpath 'pRF_VOLs/lh.pRF_rfsize.VOL.nii.gz']);

%if ~isempty(v4_dir)
%    structSize = cell(1,3);
%else
%    structSize = cell(1,4);
%end
structSize = cell(1,3);
pRF_maps = struct('ecc', structSize, 'ang', structSize, 'R2', structSize,'rfsize', structSize);

% Create high-pass filter
dispfig = 0;
filterCutoff = 0.005;
totalTRs = numTRS;
TR = 1.550;
freqdelta = 1/(totalTRs*TR);
freqs = 0:freqdelta:(freqdelta*(totalTRs-1));
times = 0:TR:(TR*(totalTRs-1));
hipassfilter = ones(1,totalTRs);
hipassfilter(freqs<filterCutoff) = 0;
% smooth with filter edge with gaussian
smoothedge = 1-(1 * exp(-(((freqs-filterCutoff).^2)/(2*(filterCutoff/2).^2))));
% add smooth edge to square filter
hipassfilter(freqs>filterCutoff) = smoothedge(freqs>filterCutoff);
hipassfilter(totalTRs:-1:round((totalTRs/2)+1)) = hipassfilter(2:round(totalTRs/2)+1);
% give back the DC
hipassfilter(1) = 1;

%{
figure('Color', [1 1 1]),
subplot(1,2,1),
plot(freqs,hipassfilter); hold on,
vline(filterCutoff, 'r-'), vline(freqs(first(find(hipassfilter > 0.999))), 'g-');
xlabel('freq (Hz)'); ylabel('Magnitude');
subplot(1,2,2),
plot(times, abs(fft(hipassfilter)));
xlabel('times (sec)'); ylabel('filter magnitude'); xaxis(0,60);
%}
for n = 1:numel(scanNames)
    
    whichScan = scanNames(n).name;
    % grab the timeseries
    tSeries = MRIread([subjpath whichScan '/fmc.siemens.nii.gz']);

    for roi = 1:numel(roi_maskNames)
        
        mask = eval(['find(' roi_maskNames{roi} '.vol == 1);']);
        
        %% pull pRF data one time
        if n == 1
            % find pRF maps
            % First combine hemispheres
            lh_ecc_map = lh_ecc.vol(mask);
            rh_ecc_map = rh_ecc.vol(mask);
            ecc_map = lh_ecc_map + rh_ecc_map;
            ecc_map(lh_ecc.vol(mask) ~= 0 & rh_ecc.vol(mask) ~= 0) = NaN;
            
            lh_r2_map = lh_R2.vol(mask);
            rh_r2_map = rh_R2.vol(mask);
            r2_map = lh_r2_map + rh_r2_map;
            r2_map(lh_R2.vol(mask) ~= 0 & rh_R2.vol(mask) ~= 0) = NaN;
            
            lh_ang_map = lh_ang.vol(mask) * (pi/180);
            rh_ang_map = rh_ang.vol(mask) * (pi/180);
            ang_map = lh_ang_map + rh_ang_map;
            ang_map(lh_ang.vol(mask) ~= 0 & rh_ang.vol(mask) ~= 0) = NaN;
            
            lh_rfsize_map = lh_rfsize.vol(mask);
            rh_rfsize_map = rh_rfsize.vol(mask);
            rfsize_map = lh_rfsize_map + rh_rfsize_map;
            rfsize_map(lh_rfsize.vol(mask) ~= 0 & rh_rfsize.vol(mask) ~= 0) = NaN;
            
            pRF_maps(roi).ecc = ecc_map;
            pRF_maps(roi).R2 = r2_map;
            pRF_maps(roi).ang = ang_map;
            pRF_maps(roi).rfsize = rfsize_map;
            if PlotFigures
                thresholded_mask = r2_map > 0;
                rVF_colormap = abs([interp(rVF(:,1), ceil(sum(thresholded_mask)/size(rVF,1))) ...
                    interp(rVF(:,2), ceil(sum(thresholded_mask)/size(rVF,1))) ...
                    interp(rVF(:,3), ceil(sum(thresholded_mask)/size(rVF,1)))]);
                
                lVF_colormap = abs([interp(lVF(:,1), ceil(sum(thresholded_mask)/size(lVF,1))) ...
                    interp(lVF(:,2), ceil(sum(thresholded_mask)/size(lVF,1))) ...
                    interp(lVF(:,3), ceil(sum(thresholded_mask)/size(lVF,1)))]);
                
                figure(1)
                subplot(2,numel(roi_maskNames),roi)
                [~, rfsize_order] = sort(lh_rfsize_map(thresholded_mask));
                polarscatter(lh_ang_map(thresholded_mask), lh_ecc_map(thresholded_mask), [], rVF_colormap(rfsize_order,:), 'filled', 'MarkerFaceAlpha', .5)
                [~, rfsize_order] = sort(rh_rfsize_map(thresholded_mask));
                hold on,
                polarscatter(rh_ang_map(thresholded_mask), rh_ecc_map(thresholded_mask), [], lVF_colormap(rfsize_order,:), 'filled', 'MarkerFaceAlpha', .5)
                title(roi_maskNames{roi})
                
                subplot(2,numel(roi_maskNames),numel(roi_maskNames)+roi)
                scatter(lh_ecc_map, lh_rfsize_map, [], [0 0 1], 'filled', 'MarkerFaceAlpha', .5)
                hold on,
                scatter(rh_ecc_map, rh_rfsize_map, [], [1 0 0],'filled', 'MarkerFaceAlpha', .5)
                box off; ylabel('RF size (deg)'), xlabel('Eccentricity bins (deg)')
                xlim([0 max(ecc_map)]); axis square
            end
            linearCoordinates{roi} = mask;

            if runLoc == true
                % Find voxels based on localizer sig maps
                if exist('loc_mask', 'var')
                    Loc_sig = loc_mask.vol(mask);
                    localizer_sig{roi} = Loc_sig;
                end
            end

        end
        
        %% TASK
        % reshape 3d volume to voxels x time
        tmp_tSeries = reshape(tSeries.vol,[(size(tSeries.vol,1)*size(tSeries.vol,2) *size(tSeries.vol,3)),size(tSeries.vol,4)]);
        % remove fixation and adaptation periods
        cut_tSeries = tmp_tSeries(linearCoordinates{roi}(:,1),:);
        %         % transform BOLD signal to Percent Signal Change
        %         mean_tSeries = mean(cut_tSeries,2);
        %         cut_tSeries = (cut_tSeries./repmat(mean_tSeries, [1, size(cut_tSeries,2)])-1)*100;
        %remove slow trend (drift) and high-pass filter data
        detrend_tSeries = NaN(size(cut_tSeries,1), size(cut_tSeries,2));
        X = [linspace(0,1,size(cut_tSeries,2))' ones(size(cut_tSeries,2),1)];
        for ii = 1:size(cut_tSeries,1)

            % detrend
            b = pinv(X)*cut_tSeries(ii,:)';
            detrend_timecourse = cut_tSeries(ii,:)' - X(:,1)*b(1,:);
            % high-pass filter
            detrend_timecourse = ifft(fft(detrend_timecourse) .* hipassfilter');
            detrend_timecourse = real(detrend_timecourse);
            % transform BOLD signal to Percent Signal Change
            dc = mean(detrend_timecourse(1:end));
            %             dc = mean(detrend_timecourse);
            dc(dc == 0) = NaN;
            detrend_tSeries(ii,:) = ((detrend_timecourse(1:end)./dc)-1)*100;

        end
        
        final_tSeries{roi} = [final_tSeries{roi} detrend_tSeries];
        clear detrend_tSeries mask tmp_tSeries cut_tSeries
    end
    
end
if runLoc == true
    save([subjpath 'fMRI_data/' fileName], 'final_tSeries', 'linearCoordinates', 'localizer_sig', 'pRF_maps', '-v7.3');
else
    save([subjpath 'fMRI_data/' fileName], 'final_tSeries', 'linearCoordinates', 'pRF_maps', '-v7.3');
end