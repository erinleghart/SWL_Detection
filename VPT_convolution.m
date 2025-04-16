% VPT_convolution.m
% Description: This function produces the convolution output. The convolution 
% is performed X-times, with X being the length of 'layerThicknesses'. The
% convolution is performed on 'spw', with the threshold of
% 'SWL_threshold'. The output is saved to the cell array 'convCell' of
% dimension {1,X}.
% Author: Erin Leghart; erin.leghart@stonybrook.edu
% Last Updated: April 14, 2025

function [layerNumber] = VPT_convolution(layerThicknesses, spw, verticalRes, SWL_threshold)

    % Create a cell 'convCell' that will contain the convolution outputs
    convCell = cell(size(layerThicknesses));
    
    % We aim to map SWLs of vertical thickess as prescribed within
    % 'layerThicknesses'. To build the convolution kernel, we need to know how
    % many range gates map to the prescribed vertical thicknesses. This section
    % below will build our convolution kernel.
    rangeResolution = verticalRes; 
    k_pos_rows = zeros(size(layerThicknesses));
    % The number of gates MUST be odd to ensure a central row within the
    % kernel. The for loop below tests whether the number of needed gates is
    % odd or even. If numGates is even, +1. If numGates is odd, +2.
    for i=1:length(layerThicknesses)
        result = fix(layerThicknesses(i)/rangeResolution);
        if rem(result, 2) == 0
            numGates = result+1;
        else
            numGates = result+2;
        end
        k_pos_rows(i)=numGates; % number of positive kernel rows
    end
    k_neg_rows=k_pos_rows*2; % number of negative kernel rows
    k_cols = k_pos_rows + k_neg_rows; % our kernel will have the same number of rows and columns
    k_rows = k_pos_rows + k_neg_rows; % our kernel will have the same number of rows and columns
    k_thresh = SWL_threshold; % our convolution threshold

        % Perform the convolution(s) for the number of prescribed layerThicknesses
    for i=1:length(layerThicknesses)
        dist_cols = ((0:k_cols(i)-1)-fix(k_cols(i)/2))*1.0; % distance from kernel center along column dimension
        dist_rows = transpose((0:k_rows(i)-1)-fix(k_rows(i)/2)); % distance from kernel center along row dimension
        r = sqrt((dist_cols .^ 2) .* transpose(ones(k_rows(i))) + (ones(k_cols(i))) .* (dist_rows .^ 2)); % distance of each kernel element from kernel center
    
        % binary mask indicating kernel elements within a radius of
        % (k_active_cols)/2 of center. This will produce a rounded kernel by
        % making the corners of the kernel null.
        active_region = (r <= fix(k_cols(i)/2));
    
        % Build the positive and negative portions of the kernel separately.
        % This ensures we normalize their respective contributions to the
        % kernel individually. Then, merge the positive and negative portions
        % of the kernel.
    
        % Positive portion - place ones in the rows of the positive portion of
        % the kernel, equally split around the center row of the kernel. Set
        % any kernel elements outside of the "active_region" to zero, as kernel
        % elements with a value of zero will not influence the convolution
        % output. Normalize the positive portion of the kernel,
        % 'kernel_positive'.
        kernel_positive = zeros(k_rows(i),k_cols(i));
        kernel_positive(1+fix(k_rows(i)/2)-fix(k_pos_rows(i)/2):1+fix(k_rows(i)/2)+fix(k_pos_rows(i)/2),:) = 1.0;
        kernel_positive = kernel_positive.*active_region;
        kern_pos_sum = sum(kernel_positive, 'all');
        kernel_positive = kernel_positive/kern_pos_sum; %normalize
    
        % Negative portion - place negative ones in the rows of the negative
        % portion of the kernel. The negative portions of the kernel will be on
        % the top and bottom of the kernel. For the rows where the positive
        % portion of the kernel will sit, fill in with ones. Set any kernel 
        % elements outside of the "active_region" to zero, as kernel elements
        % with a value of zero will not influence the convolution output.
        % Normalize the negative portion of the kernel, 'kernel_negative'.
        % Apply the SWL threshold to the negative portion of the kernel. This
        % ensures that any convolution output > 0 will be indicative of a "hit"
        % for a SWL.
        kernel_negative = zeros(k_rows(i),k_cols(i));
        kernel_negative(1+fix(k_rows(i)/2)-fix(k_pos_rows(i)/2)-fix(k_neg_rows(i)/2):1+fix(k_rows(i)/2)+fix(k_pos_rows(i)/2)+fix(k_neg_rows(i)/2),:) = 1.0;
        kernel_negative(1+fix(k_rows(i)/2)-fix(k_pos_rows(i)/2):1+fix(k_rows(i)/2)+fix(k_pos_rows(i)/2),:) = 0.0;
        kernel_negative = kernel_negative.*active_region;
        kern_neg_sum = sum(kernel_negative, 'all');
        kernel_negative = kernel_negative/kern_neg_sum; % normalize
        kernel_negative = kernel_negative * (1.0 + k_thresh);
        kernel_negative = -1.0 * kernel_negative;
    
        kernel = kernel_positive + kernel_negative; % merge into the final kernel
    
        % Peform the convolution
        result = conv2(spw, kernel, 'same');
        % The below line sets any convolution output <= 0 to 0, as these SW
        % elements do not meet the threshold for a "hit" within the SWL
        % algorithm. Any nan values are also set to 0.
        result(result <= 0) = 0; TF = isnan(result); result(TF) = 0;
        convCell{1,i} = result; % assign the convolution output to its assigned cell
    end

    % X-number of convolutions were completed based on the prescribed
    % "layerThicknesses". All convolution outputs are combined into a single
    % array which will identify which SW elements were identified as a "hit"
    % within the SWL convolution. Any element which received at least one "hit"
    % from the convolution will be identified with a SWL. Thus, if an element
    % only received one "hit" from one convolution, or X "hits", it is still
    % retained.
    stacked = cat(4, convCell{:}); % Stack along a new 4th dimension
    isLayer = sum(stacked, 4); % Sum along that new dimension
    isLayer(isLayer > 0) = 1; % binary array containing locations where a SWL is present = 1

    % Now, give each unique SWL a identifying number which is stored in
    % layerNumber.
    layerNumber = bwlabel(isLayer);

end

