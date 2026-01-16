function [cloud_top_height_km, cloud_base_height_km, cloud_depth_km] = cloud_depth_PPI(ref, zkm, nan_threshold)
% cloud_depth_PPI calculates cloud top, base, and depth from a KASPR PPI scan.
% A vertical level is considered cloudy if fewer than `nan_threshold` fraction
% of values are NaN along the azimuth dimension.
%
% Inputs:
%   zhh_2D         - 2D reflectivity array [height x azimuth]
%   heights        - 1D array of height values (same length as size(zhh_2D, 1))
%   nan_threshold  - Scalar (e.g., 0.75); max fraction of NaNs allowed to be considered "cloud"
%
% Outputs:
%   cloud_top_height   - Highest altitude with valid cloud signal
%   cloud_base_height  - Lowest altitude with valid cloud signal
%   cloud_depth        - Depth between cloud base and top

    % Extract 1D height
    heights = zkm(:,1);

    % Mask reflectivity values below threshold
    ref(ref < -30) = NaN;

    % Compute fraction of NaNs at each height level
    frac_nans = sum(isnan(ref), 2) / size(ref, 2);

    % Identify height indices that meet the nan_threshold criterion
    valid_rows = find(frac_nans <= nan_threshold);

    if isempty(valid_rows)
        cloud_top_height_km = NaN;
        cloud_base_height_km = NaN;
        cloud_depth_km = NaN;
        return;
    end

    % Extract top and base heights
    cloud_base_height_km = heights(min(valid_rows));
    cloud_top_height_km = heights(max(valid_rows));
    cloud_depth_km = cloud_top_height_km - cloud_base_height_km;
end
