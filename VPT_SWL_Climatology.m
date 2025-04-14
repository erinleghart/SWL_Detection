%VPT_SWL_Climatology.m
% Description: Main driver for identifying and analyzing SWLs in KASPR VPT
% data.
% Author: Erin Leghart; erin.leghart@stonybrook.edu
% Last Updated: April 14, 2025

% Load Data
filename = 'kaspr_dates.csv';
dates = readtable(filename);
kasprVPTDataDir = '/path/to/kasprdata';
datalist = dir([kasprVPTDataDir, 'KASPR_VPT_SWL_MOMENTS_', '*.nc']);

% Empty arrays which will populate with SWL characteristics
layerHeight_km = [];
layerThickness_m = [];
layerMagnitude = [];
layerDuration_s = [];
profileDateTime = [];
stormDate = [];
stormNum = [];
profile_duration_s = [];

% Define the parameters of the convolution
layerThicknesses = [100, 250, 500]; % thicknesses in m of the desired VPT SWLs
SWLThreshold = 0.20; % SW threshold in m/s

% Loop through each KASPR VPT file to apply the SWL detection algorithm
for i=1:10%length(datalist)
    kasprdata = [datalist(i).folder,'/', datalist(i).name];

    % Step 1: Load in KASPR VPT data and extract datetime information
    [timeh, times, ref, spw, snr, rangekm, elev, verticalRes, file_duration_s] =...
        kaspr_variables_VPT(kasprdata);
    profileDate_file = kasprdata(end-17:end-10); % profile date
    profileDateTime_file = kasprdata(end-17:end-3); % profile date and time

    % Step 2: Perform the convolution and identify VPT SWLs
    [layerNumber] = VPT_convolution(layerThicknesses, spw, verticalRes, SWLThreshold);

    % Step 3: Compute VPT SWL properties
    [layerHeight_list, layerThickness_list, layerMagnitude_list, layerDuration_list] =...
        VPT_SWL_properties(layerNumber, spw, rangekm, verticalRes);

    % Step 4: Append SWL properties from current KASPR VPT scan to main
    % list of VPT SWL characteristics
    layerHeight_km = cat(1, layerHeight_km, layerHeight_list); 
    layerThickness_m = cat(1, layerThickness_m, layerThickness_list); 
    layerMagnitude = cat(1, layerMagnitude, layerMagnitude_list); 
    layerDuration_s = cat(1, layerDuration_s, layerDuration_list); 
    dateTime_list = repmat(profileDateTime_file, [length(layerHeight_list),1]);
    profileDateTime = cat(1, profileDateTime, dateTime_list);
    stormDate_list = repmat(profileDate_file, [length(layerHeight_list),1]);
    stormDate = cat(1, stormDate, stormDate_list);
    loc = string(dates.Dates) == profileDate_file;
    stormNum_profile = dates.StormNum(loc, :);
    stormNum_list = repmat(stormNum_profile, [length(layerHeight_list),1]);
    stormNum = cat(1, stormNum, stormNum_list);
    profile_duration_s_list = repmat(file_duration_s, [length(layerHeight_list),1]);
    profile_duration_s = cat(1, profile_duration_s, profile_duration_s_list);
end

% Save VPT SWL dataset as a .csv file
data_save_dir = '/path/to/save/directory/';
T = table(profileDateTime, stormDate, stormNum, layerHeight_km, layerThickness_m, layerMagnitude,...
    layerDuration_s, profile_duration_s);
writetable(T, [data_save_dir, 'VPT_SWL_climatology.csv']);
