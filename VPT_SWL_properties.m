% VPT_SWL_properties.m
% Description: Compute the spatial and magnitude characteristics of VPT SWLs
% identified within layerNumber.
% Author: Erin Leghart; erin.leghart@stonybrook.edu
% Last Updated: April 14, 2025

function [layerHeight_list, layerThickness_list, layerMagnitude_list, layerDuration_list] = ...
    VPT_SWL_properties(layerNumber, spw, rangekm, verticalRes)

    stats = regionprops(layerNumber,'BoundingBox');
    BoundingBox = cat(1, stats.BoundingBox);
    if ~isempty(BoundingBox) % Will skip any files that have no layers - usually files with bad data
    
        % Use properties computed by 'regionprops' 'BoundingBox' to determine
        % SWL height, thickness, and duration. Note: VPT SWL duration, in
        % seconds, is equal to the number of SW observations identified
        % within the VPT SWL, as KASPR sampling frequency is 1 Hz.
        layerDuration_list = BoundingBox(:,3); % layer duration in seconds (1 sampling volume = 1 second)
        layerThicknessRange_list = BoundingBox(:,4); % layer thickness in range gates
        lowestRow = BoundingBox(:,2);
        upperRow = lowestRow + layerThicknessRange_list;
        layerRangeGate = (upperRow+lowestRow)/2; % central range gate of the SWL will be the SWL height
        layerHeight_list = rangekm(fix(layerRangeGate),1); % Use midpoint gate of bounding box to estimate height (convert from range gate index to km)
        layerThickness_list = layerThicknessRange_list * verticalRes; % in meters
        
        % Determine the SWL magnitude. SWL magnitude is the average SW within a
        % SWL.
        layerMagnitude_list = nan(size(layerThickness_list));
        for i=1:length(layerThickness_list)
            layerMagnitude_list(i) = mean(spw(layerNumber == i), 'omitnan');
        end
    else % Any files without SWLs (usually bad data) will record no SWL properties
        layerHeight_list = [];
        layerThickness_list = [];
        layerMagnitude_list = [];
        layerDuration_list = [];
    end 
    
end
