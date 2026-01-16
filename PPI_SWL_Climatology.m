% PPI_SWL_Climatology.m
% Description: Main driver for identifying and analyzing SWLs in KASPR PPI
% data. This script takes ~1h to complete.
% Author: Erin Leghart; erin.leghart@stonybrook.edu
% Date: January 14, 2026


%% User-Defined Input Directories

kasprPPIDataDir = '/path/to/kasprdata/';
data_save_dir = '/path/to/save/directory/';

%% Load Data & Set Constants

dates = readtable('/path/to/kaspr_dates.csv');
datalist = dir([kasprPPIDataDir, 'KASPR_PPI_SWL_MOMENTS_', '*.nc']);
disp([num2str(length(datalist)), ' KASPR PPI files']);

% Melting Layer Dataset (removed at end)
% melting_layers = '/path/to/kaspr_melting_layers.csv'
melting_layers = readtable('/path/to/identified_melting_layers_KASPR_PPI.csv');

%% Preallocate lists to store SWL characteristics

layerHeight_km = [];
layerThickness_m = [];
layerAzimuth_deg = [];
layerSquare_km = [];
layerMagnitude = [];
layerVelGrad = [];
layerFitVelGrad = [];
layerResidVel = [];
scanDuration_s = [];
cloudTopHeight = [];
cloudBaseHeight = [];
cloudDepth = [];
profileDateTime = [];
stormDate = [];
stormNum = [];
stormCategory = [];

%% SWL Detection Algorithm

% This main section of the SWL detection algorithm loops through all swl
% files contained within the .zip files swl_files_2017_2018.zip,
% swl_files_2018_2019.zip, swl_files_2019_2020.zip,
% swl_files_2020_2021.zip. There are XX main sections within this code

% 1. Load KASPR_PPI_SWL_moments_YYYYMMDD-HHmmss.nc file & perform
% preliminary QC metrics
% 2. VAD Fitting
% This consistists of the function "compute_fourier_coeffs_by_range" which
% performs a harmonic decomposition of the observed doppler velocity (vel)
% to be used for the VAD fitting.
% 3. Cloud Base, Top, and Depth estimation
% 4. SWL Detection & QC
% 5. Remove SWLs in the Melting Layer
% 6. Save SWL dataset as PPI_SWL_Climatology.csv

% Within Section 4. SWL Detection & QC, there are multiple steps:
% 4.1 Data Prep
% 4.2 Calculate Moving Percentile Windows
% 4.3 Compare Moving Percentile Windows to SW Field
% 4.4 Binary Label Operator & Identification of SWLs
% 4.5 SWL Characteristics & Noise Removal
% 4.6 Append SWLs to Growing Lists


tic;
for i = 1:length(datalist)

    % --- Print statement to check on status while running
    if mod(i, 500) == 0 % print ever 500 files
        disp(['File ', num2str(i), '/', num2str(length(datalist))]);
    end

    % =====================================================================
    % 1. Load KASPR_PPI_SWL_moments_YYYYMMDD-HHmmss.nc file & perform
    % preliminary QC metrics
    % =====================================================================

    kasprdata = [datalist(i).folder,'/', datalist(i).name];
    [timeh, times, ref, zdr, vel, spw, rhohv, kdp, snr, rangekm, xkm, ykm, zkm,...
        elev_deg, az_deg, az_rad, file_duration_s] = PPI_kaspr_variables(kasprdata);
    profileDate_file = kasprdata(end-17:end-10); % scan date
    profileDateTime_file = kasprdata(end-17:end-3); % scan date and time
    profileTime_file = kasprdata(end-8:end-3); % scan date and time
    profileTime_file = duration(str2double(extractBetween(profileTime_file,1,2)), ... % seconds
        str2double(extractBetween(profileTime_file,3,4)), ... % minutes
        str2double(extractBetween(profileTime_file,5,6)));   % seconds

    % --- Calculate relavent gradients
    [~,dz] = gradient(zkm'); % vertical gradient of altitude in km
    verticalRes = mode(dz', 'all') * 1000; % vertical gradient of altitude m. This is not range resolution.
    azimuthal_deg_gradient = mean(gradient(az_deg(:,1)));
    isnan_spw = isnan(spw);

    % --- Compute square kilometer per pixel, which is used for QC (see
    % Leghart et al. 2026)
    [nr, nc] = size(xkm);
    cellArea = nan(nr-1, nc-1); % one area per cell (between 4 corners)
    for m = 1:nr-1
        for n = 1:nc-1
            % Four corners of the cell
            xCorners = [xkm(m,n)   xkm(m,n+1)   xkm(m+1,n+1)   xkm(m+1,n)];
            yCorners = [ykm(m,n)   ykm(m,n+1)   ykm(m+1,n+1)   ykm(m+1,n)];
    
            % Shoelace formula for polygon area
            cellArea(m,n) = 0.5 * abs( sum(xCorners .* circshift(yCorners, -1)) ...
                                     - sum(yCorners .* circshift(xCorners, -1)) );
        end
    end

    % Assign each pixel the mean area of its surrounding cells
    pixelArea_full = nan(nr, nc);
    pixelArea_full(1:end-1, 1:end-1) = cellArea; % upper-left assignment
    pixelArea = pixelArea_full;
    pixelArea(isnan_spw) = nan;

    % =====================================================================
    % 2. VAD Fitting
    % =====================================================================
    
    % --- Compute Fourier coefficients needed for VAD Fitting
    coeffs = compute_fourier_coeffs_by_range(vel, az_rad, 100);
    a0 = coeffs.a0;
    a1 = coeffs.a1;
    b1 = coeffs.b1;
    a2 = coeffs.a2;
    b2 = coeffs.b2;

    % --- Calculate Vr_fit 
    vr_fit = NaN(size(vel)); % Preallocate a Vr_fit array
    for r = 1:size(vel, 2) % Loop through each range gate
        beta = az_rad(:, r);  % 1D azimuth angles [rad] at this range gate
        % Fit Vr using full Fourier series
        vr_fit_range = (a0(r) + ...
                        a1(r) * cos(beta) + ...
                        b1(r) * sin(beta) + ...
                        a2(r) * cos(2 * beta) + ...
                        b2(r) * sin(2 * beta));
        % Store Vr_fit
        vr_fit(:, r) = vr_fit_range;    
    end

    % --- Compute Residual Velocity (vr_resid)
    vr_resid = vel - vr_fit; % vr_resid = vel (observed velocity) - vr_fit (fitted velocity) (Eqn. 1 Leghart et al. 2026)
    abs_vr_resid = abs(vr_resid); % Apply absolute value to quantify vr_resid magnitude, not sign
    
    % --- Mask vr_resid where either Vr or Vr_fit is NaN (QC)
    nan_mask = isnan(vel) | isnan(vr_fit);
    abs_vr_resid(nan_mask) = NaN;

    % --- Calculate Fitted Velocity Radial Gradient (abs_fit_vel_grad) -
    % proxy of the shear of the resolved flow
    fitted_vel_grad = gradient(vr_fit,1, 2);
    abs_fit_vel_grad = abs(fitted_vel_grad);

    % =====================================================================
    % 3. Cloud Base, Top, and Depth estimation
    % =====================================================================

    nan_threshold = 0.75;
    [cloud_top_height_file, cloud_base_height_file, cloud_depth_file] = cloud_depth_PPI(ref', zkm', nan_threshold);

    % =====================================================================
    % 4. SWL Detection & QC
    % ===================================================================== 

    % === 4.1 Data Prep

    % --- Define the percentile needed to be overtaken to identify as a SWL
    percentiles = 75;
    percentiles = percentiles(:)'; % row

    % --- Define Minimum Valid Fraction
    % This is the fraction of data within each radar resolution volume that
    % is required to perform the SWL detection algorithm. This was done to
    % ensure only high-quality data are retained. Please see Leghart et al.
    % 2026 for more details.
    MinValidFrac  = 0.75;

    % --- Define vertical windows for SWL
    layerThicknesses = [100, 250, 500, 750, 1000]; % thicknesses in m of the desired PPI SWLs

    % --- Compute Window Lengths
    gate_sizes = layerThicknesses/verticalRes; % convert vertical thickness to number of range gates
    g100 = fix(gate_sizes(1));
    g250 = fix(gate_sizes(2));
    g500 = fix(gate_sizes(3));
    g750 = fix(gate_sizes(4));
    g1000 = fix(gate_sizes(5));
    gateCounts  = [g100 g250 g500 g750 g1000];   % window lengths (integers)
    gateCounts = gateCounts(:)'; % row
    gateCounts = 2*floor(gateCounts/2) + 1; % force odd lengths so windows are centered
    
    % --- Preallocate Output Percentile maps
    [nAz, nR] = size(spw); %nAz = number of azimuth observations, nR = number of range observations
    nP = numel(percentiles); % number of percentiles to compute (in this case, only one percentile)
    nT = numel(gateCounts); % number of thicknesses to parse through (in this case, 5 thicknesses)
    values = nan(nAz, nR, nP, nT, 'like', spw);   % output percentile maps

    % === 4.2 Calculate Moving Percentile Windows
    
    for k = 1:nT % Main loop over window sizes (thickness scales)
        T = gateCounts(k); % T = thickness (in gate counts, not true thickness)
        for ip = 1:nP
            pctl = percentiles(ip);
            half = floor(T/2);
            q = pctl/100; % Fractional form of percentile - 75% vs 0.75
    
            % --- REPLICATE EDGE HANDLING (no colfilt, no helper needed) ---
            tmp = nan(nAz, nR, 'like', spw);

            for r = 1:nR
                % Build replicate-padded column indices for centered window
                jj = (r-half):(r+half);
                jj = max(1, min(nR, jj));              % clamp -> replicate edges

                seg = spw(:, jj);                      % [nAz x winLen]
                segS = sort(seg, 2, 'ascend');         % NaNs sort to end
                validCnt = sum(~isnan(segS), 2);       % [nAz x 1], per row

                % Fractional rank (1-based), linear interp between order stats
                rr = 1 + (validCnt - 1).*q;
                lo = floor(rr); 
                hi = ceil(rr);
                w  = rr - lo;

                % Safe index guards
                lo(lo < 1) = 1; 
                hi(hi < 1) = 1;

                % Gather lower/upper values per row
                idxL = sub2ind(size(segS), (1:nAz)', lo);
                idxU = sub2ind(size(segS), (1:nAz)', hi);
                yL   = segS(idxL);
                yU   = segS(idxU); % if your editor dislikes [], switch back to ()

                % Interpolate percentile; empty windows -> NaN
                y = (1 - w).*yL + w.*yU;
                y(validCnt == 0) = NaN;

                % Optional: enforce minimum valid fraction in the window
                % (we do enforce this)
                if MinValidFrac > 0
                    needed = ceil(MinValidFrac * numel(jj));
                    y(validCnt < needed) = NaN;
                end

                tmp(:, r) = y;
            end
            % Store result: [n_az x n_range] map for this percentile & thickness
            values(:,:,ip,k) = tmp;
        end
    end

    % === 4.3 Compare Moving Percentile Windows to SW Field
    
    % "values" contains the 75th percentile of SW for each of the SWL
    % thicknesses prescribed in "layerThicknesses". Thus, the true SW value
    % within "spw" must exceed that in "values" to be considered large
    % enough of a SW enhancement to potentially be within a SWL.

    spw_comparison = zeros(nAz, nR);
    for k = 1:nT
        comparison = spw > values(:,:,1,k);
        spw_comparison = spw_comparison + comparison;
    end
    spw_comparison(spw_comparison > 0) = 1;
    spw_comparison = imfill(spw_comparison, 8, 'holes');
    
    % === 4.4 Binary Label Operator & Identification of SWLs

    % The minPix threshold is determined by taking into account the minimum
    % number of range gates required to meet azimuth thresholds. We will
    % automatically remove any SWLs which have a azimuth of 10 degrees or
    % less. Thus, if a SWL were to be one range gate radially, we need to
    % determine how many range gates will account for 10 degrees in
    % azimuth. This value will be our minPix value. The same threshold will
    % be passed for the vertical thickness of the layer, which will be the
    % smallest gate windoe (in range gates) we pass.

    az_thresh_deg = 10;
    minPix = floor(az_thresh_deg/azimuthal_deg_gradient);
    
    % --- bwlabel joins connected elements into one object
    layerNumber = bwlabel(spw_comparison);
    L = layerNumber; 
    L(~isfinite(L)) = 0; % background (non-SWL elements) must be 0 (not NaN/Inf)
    
    % --- Build a mask of labeled pixels, then connected components
    mask = L > 0;
    CC = bwconncomp(mask, 8);

    % If a connected component (CC) meets the azimuthal criteria set above,
    % it will be deemed a SWL. These SWLs will be further QC'd below to
    % remove additional noise.
   
    stats = regionprops(CC, 'Area'); % Measure area and MinorAxisLength per CC
    keep = false(CC.NumObjects,1);
    for idx = 1:CC.NumObjects
        A = stats(idx).Area;
        if A >= minPix % Keep CC's that meet the azimuthal criteria
            keep(idx) = true;
        end
    end
    
    maskKeep = false(size(mask)); % Rebuild a kept mask
    for idx = find(keep)'
        maskKeep(CC.PixelIdxList{idx}) = true;
    end
    
    layerNumber = bwlabel(maskKeep, 8); % Relabel CC's starting over at "1" using an 8 point connection

    % Because we are identifying SWLs within a PPI scan, we have to account
    % for SWLs which may cross the 0-360 degree line. This snippet below
    % "jumps" that 0-360 degree line, and will relabel to CC's as one if
    % they are connected across the 0-360 line. 

    az_comparison = vertcat(layerNumber(1,:), layerNumber(end, :))';
    rowsBothPositive = all(az_comparison > 0, 2);   % logical mask
    validRows = az_comparison(rowsBothPositive, :);
    col1 = az_comparison(rowsBothPositive, 1);
    col2 = az_comparison(rowsBothPositive, 2);
    for layer_num = 1:length(col1)
        locs = find(layerNumber == col1(layer_num));
        layerNumber(locs) = col2(layer_num);
    end

    % === 4.5 SWL Characteristics & Noise Removal
    
    % --- Preallocate SWL characteristics to save from this file 
    layerHeight_list = [];
    layerThickness_list = [];
    layerAzimuth_list = [];
    layerSquare_km_list = [];
    layerMagnitude_list = [];
    layerVelGrad_list = [];
    layerFitVelGrad_list = [];
    layerResidVel_list = [];
    unique_layer_nums = unique(layerNumber);
    unique_layer_nums = unique_layer_nums(2:end); % removes 0 as a number

    % --- Loop through each unique SWL (CC from previous step) and
    % calculate spatial and temporal characteristics
    for layer = 1:length(unique_layer_nums)
        [row, col] = find(layerNumber == unique_layer_nums(layer));

        % --- Height
        layerHeight_gate = round(median(col));
        layerHeight_swl = zkm(1, layerHeight_gate);
        layerHeight_list(end+1) = layerHeight_swl;

        % --- Thickness (avg across all azimuths)
        uniqueRows = unique(row); % unique rows (azimuths) in this object
        thicknessPerRow = nan(size(uniqueRows));
        for r = 1:numel(uniqueRows)
            thisRow = uniqueRows(r);
            colsInRow = col(row == thisRow);
            thicknessPerRow(r) = (max(colsInRow) - min(colsInRow)) * verticalRes;
        end
        layerThickness_swl = mean(thicknessPerRow, 'omitnan'); % Average thickness across rows
        layerThickness_list(end+1) = layerThickness_swl;

        % --- Azimuth [deg] (wrap-safe across 0/360 degree line)
        layerAzimuth_list(end+1) = numel(uniqueRows) * azimuthal_deg_gradient;

        % --- Area [square km]
        layerSquare_km_swl = sum(pixelArea(layerNumber == unique_layer_nums(layer)), 'omitnan');
        layerSquare_km_list(end+1) = layerSquare_km_swl;

        % --- Layer Magnitude [m/s]
        layerMagnitude_swl = mean(spw(layerNumber == unique_layer_nums(layer)), 'omitnan');
        layerMagnitude_list(end+1) = layerMagnitude_swl;

        % --- Layer Fitted Velocity Gradient [m/s]
        layerFitVelGrad_swl = mean(abs_fit_vel_grad(layerNumber == unique_layer_nums(layer)), 'omitnan');
        layerFitVelGrad_list(end+1) = layerFitVelGrad_swl;

        % --- Layer Residual Velocity [m/s]
        layerResidVel_swl = mean(abs_vr_resid(layerNumber == unique_layer_nums(layer)), 'omitnan');
        layerResidVel_list(end+1) = layerResidVel_swl;

    end
    
    % --- Removal of noisy SWLs by spatial thresholds

    % Height
    height_thresh = (layerHeight_list > 1);
    layerHeight_list = layerHeight_list(height_thresh)';
    layerThickness_list = layerThickness_list(height_thresh)';
    layerAzimuth_list = layerAzimuth_list(height_thresh)';
    layerSquare_km_list = layerSquare_km_list(height_thresh)';
    layerMagnitude_list = layerMagnitude_list(height_thresh)';
    layerFitVelGrad_list = layerFitVelGrad_list(height_thresh)';
    layerResidVel_list = layerResidVel_list(height_thresh)';
    unique_layer_nums = unique_layer_nums(height_thresh);

    % Azimuth
    az_thresh = (layerAzimuth_list > az_thresh_deg); % threshold set earlier
    layerHeight_list = layerHeight_list(az_thresh);
    layerThickness_list = layerThickness_list(az_thresh);
    layerAzimuth_list = layerAzimuth_list(az_thresh);
    layerSquare_km_list = layerSquare_km_list(az_thresh);
    layerMagnitude_list = layerMagnitude_list(az_thresh);
    layerFitVelGrad_list = layerFitVelGrad_list(az_thresh);
    layerResidVel_list = layerResidVel_list(az_thresh);
    unique_layer_nums = unique_layer_nums(az_thresh);

    % Area - Square Kilometer
    area_thresh_km = 2;
    area_thresh = (layerSquare_km_list > area_thresh_km);
    layerHeight_list = layerHeight_list(area_thresh);
    layerThickness_list = layerThickness_list(area_thresh);
    layerAzimuth_list = layerAzimuth_list(area_thresh);
    layerSquare_km_list = layerSquare_km_list(area_thresh);
    layerMagnitude_list = layerMagnitude_list(area_thresh);
    layerFitVelGrad_list = layerFitVelGrad_list(area_thresh);
    layerResidVel_list = layerResidVel_list(area_thresh);
    unique_layer_nums = unique_layer_nums(area_thresh);

    overlap = ismember(layerNumber, unique_layer_nums);
    layerNumber(~overlap) = 0;

    % ---  4.6 Append SWLs to Growing Lists

    % SWL characteristics are calculated within one file at a time. All SWL
    % properties are assigned to the preallocated lists set up at the
    % beginning of this code. Once all files have been processed, the lists
    % will be joined together into a single table, and named
    % "PPI_SWL_Climatology.csv".

    layerHeight_km = cat(1, layerHeight_km, layerHeight_list);  
    layerThickness_m = cat(1, layerThickness_m, layerThickness_list); 
    layerAzimuth_deg = cat(1, layerAzimuth_deg, layerAzimuth_list);
    layerSquare_km = cat(1, layerSquare_km, layerSquare_km_list);
    layerMagnitude = cat(1, layerMagnitude, layerMagnitude_list); 
    layerFitVelGrad = cat(1, layerFitVelGrad, layerFitVelGrad_list);
    layerResidVel = cat(1, layerResidVel, layerResidVel_list);
    cloudTopHeight_list = repmat(cloud_top_height_file, [length(layerHeight_list),1]);
    cloudTopHeight = cat(1, cloudTopHeight, cloudTopHeight_list);
    cloudBaseHeight_list = repmat(cloud_base_height_file, [length(layerHeight_list),1]);
    cloudBaseHeight = cat(1, cloudBaseHeight, cloudBaseHeight_list);
    cloudDepth_list = repmat(cloud_depth_file, [length(layerHeight_list),1]);
    cloudDepth = cat(1, cloudDepth, cloudDepth_list);
    file_duration_s_list = repmat(file_duration_s, [length(layerHeight_list),1]);
    scanDuration_s = cat(1, scanDuration_s, file_duration_s_list);
    dateTime_list = repmat(profileDateTime_file, [length(layerHeight_list),1]);
    profileDateTime = cat(1, profileDateTime, dateTime_list);
    stormDate_list = repmat(profileDate_file, [length(layerHeight_list),1]);
    stormDate = cat(1, stormDate, stormDate_list);

    % --- Some dates have more than one storm number, so careful slicing is
    % needed to ensure each SWL is assigned the correct storm number.
    loc = find(string(dates.Dates) == profileDate_file);
    split_times = dates.SplitTimeStart(loc, :);
    if numel(split_times) >= 2
        split_times_dur = timeofday(split_times); % duration array
        later_mask = profileTime_file > split_times_dur; % compare file time to split time
        later_idx = find(later_mask);
        if isempty(later_idx) % means the PPI scan is before any split time
            stormNum_profile = dates.StormNum(loc(1), :);
            stormNum_list = repmat(stormNum_profile, [length(layerHeight_list),1]);
            stormNum = cat(1, stormNum, stormNum_list);
            stormCat_profile = dates.StormTypeNum(loc(1), :);
            stormCat_list = repmat(stormCat_profile, [length(layerHeight_list),1]);
            stormCategory = cat(1, stormCategory, stormCat_list);
        elseif max(later_idx) == 2 % means the PPI file is after the first split time
            stormNum_profile = dates.StormNum(loc(2), :);
            stormNum_list = repmat(stormNum_profile, [length(layerHeight_list),1]);
            stormNum = cat(1, stormNum, stormNum_list);
            stormCat_profile = dates.StormTypeNum(loc(2), :);
            stormCat_list = repmat(stormCat_profile, [length(layerHeight_list),1]);
            stormCategory= cat(1, stormCategory, stormCat_list);
        elseif max(later_idx) > 2 % means the PPI file is after the second split time (only one case has >2 splits)
            stormNum_profile = dates.StormNum(loc(3), :);
            stormNum_list = repmat(stormNum_profile, [length(layerHeight_list),1]);
            stormNum = cat(1, stormNum, stormNum_list);
            stormCat_profile = dates.StormTypeNum(loc(3), :);
            stormCat_list = repmat(stormCat_profile, [length(layerHeight_list),1]);
            stormCategory= cat(1, stormCategory, stormCat_list);
        end
    else
        stormNum_profile = dates.StormNum(loc, :);
        stormNum_list = repmat(stormNum_profile, [length(layerHeight_list),1]);
        stormNum = cat(1, stormNum, stormNum_list);
        stormCat_profile = dates.StormTypeNum(loc, :);
        stormCat_list = repmat(stormCat_profile, [length(layerHeight_list),1]);
        stormCategory = cat(1, stormCategory, stormCat_list);
    end
    if length(stormNum) ~= length(layerHeight_km)
        disp(['File # ', num2str(i), ' incorrect data length']);
    end

end

% --- Create table out of all SWL characteristics list
T_SWLs =  table(profileDateTime, stormDate, stormNum, stormCategory, layerHeight_km, ...
    layerThickness_m, layerMagnitude, layerFitVelGrad, layerResidVel, ...
    layerAzimuth_deg, layerSquare_km, scanDuration_s, cloudTopHeight, cloudBaseHeight, ...
    cloudDepth);

% =========================================================================
% 5. Remove SWLs in the Melting Layer
% =========================================================================

% --- Quality control PPI MLs
ppi_mls = melting_layers;
ppi_mls.profileDatetime_datetime = datetime(ppi_mls.profileDatetime, 'InputFormat', 'yyyyMMdd-HHmmss');

% ZDR Mean Maximum (anomolous maximum)
locs = find(ppi_mls.ZDR_Mean_Anom > 0);
ppi_mls = ppi_mls(locs,:);

% Top/Bottom ZDR Gradient
ppi_mls = ppi_mls(ppi_mls.ZDR_Grad_Top > 0 &...
    ppi_mls.ZDR_Grad_Bot < 0, :);

% Bottom ML altitude, ML Thickness - Texture
ppi_mls = ppi_mls(ppi_mls.Bottom_km < 5, :);
ppi_mls = ppi_mls(ppi_mls.Top_km > 1, :);
ppi_mls = ppi_mls(ppi_mls.MLThickness_m < 2000, :);

% Continuous subcritical CC in deg
ppi_mls = ppi_mls(ppi_mls.CC_RadialDeg > 9,:);

% --- Removal of SWLs that extend into or are fully within a ML
T_SWLs = removeMLfromPPI(T_SWLs, ppi_mls);

% =========================================================================
% 6. Save SWL dataset as PPI_SWL_Climatology.csv
% =========================================================================

writetable(T_SWLs, [data_save_dir, 'PPI_SWL_climatology.csv']);

disp('Done - PPI SWL');
toc;

%% Functions

function T_SWLs = removeMLfromPPI(T_SWLs, ppi_mls)

    % Initialize
    T_SWLs.within_ml = false(size(T_SWLs.layerHeight_km));
    T_SWLs.below_ml = false(size(T_SWLs.layerHeight_km));

    % Loop over MLs
    for j = 1:numel(ppi_mls.Top_km)
        ml_top = ppi_mls.Top_km(j);
        ml_bottom = ppi_mls.Bottom_km(j);
        ml_time = ppi_mls.profileDatetime(j);

        % Match by time
        idx = strcmp(T_SWLs.profileDateTime, ml_time);

        % Within ML
        T_SWLs.within_ml(idx) = ...
            T_SWLs.layerHeight_km(idx) >= ml_bottom & ...
            T_SWLs.layerHeight_km(idx) <= ml_top;

        % Below ML top
        T_SWLs.below_ml(idx) = ...
            T_SWLs.layerHeight_km(idx) <= ml_top;
    end

    % Keep only rows NOT within ML
    T_SWLs = T_SWLs(~T_SWLs.within_ml, :);
end

