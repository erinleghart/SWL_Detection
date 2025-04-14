% PPI_SWL_properties.m
% Description: Compute the spatial and magnitude characteristics of PPI SWLs
% identified within layerNumber.
% Author: Erin Leghart; erin.leghart@stonybrook.edu
% Last Updated: April 14, 2025

function [layerHeight_list, layerThickness_list, layerAzimuth_list,  layerMagnitude_list] = ...
    PPI_SWL_properties(layerNumber, spw, zkm, az_deg, verticalRes)
    
    stats = regionprops(layerNumber,'BoundingBox');
    BoundingBox = cat(1, stats.BoundingBox);
    rangeResolution = median(verticalRes, 'all', 'omitnan');
    if ~isempty(BoundingBox) % Will skip any files that have no layers - usually files with bad data
    
        % Use properties computed by 'regionprops' 'BoundingBox' to determine
        % SWL height, thickness, and azimuth.
        layerThicknessRange_list = BoundingBox(:,4); % Layer Thickness - in range gates
        lowestRow = BoundingBox(:,2);
        upperRow = lowestRow + layerThicknessRange_list;
        layerRangeGate = (upperRow+lowestRow)/2; % central range gate of the SWL will be the SWL height
        layerHeight_list = zkm(fix(layerRangeGate),1); % convert range gate to altitude using zkm
        layerThickness_list = layerThicknessRange_list * rangeResolution; % in meters
        
        % Determine the SWL magnitude and azimuth by using the bounding box
        % properties. SWL magnitude is the average SW within a SWL. SWL azimuth
        % is computed by locating the lowest and highest azimuth contained
        % within the SWL.
        layerMagnitude_list = nan(size(layerThickness_list));
        layerAzimuth_list = nan(size(layerThickness_list));
        for i=1:length(layerThickness_list)
            layerMagnitude_list(i) = mean(spw(layerNumber == i),'omitnan');
            [~, col] = find(layerNumber == i);
            cols = unique(col);
            azimuth = az_deg(1, numel(cols));
            layerAzimuth_list(i) = azimuth;
    
        end
    else % Any files without SWLs (usually bad data) will record no SWL properties
        layerHeight_list = [];
        layerThickness_list = [];
        layerAzimuth_list = [];
        layerMagnitude_list = [];
    end 
end
