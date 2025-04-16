% PPI_SWL_Climatology.m
% Description: Main driver for identifying and analyzing SWLs in KASPR PPI
% data. This script takes ~1h to complete.
% Author: Erin Leghart; erin.leghart@stonybrook.edu
% Last Updated: April 14, 2025

% Directories
kasprPPIDataDir = '/path/to/kasprdata/';
data_save_dir = '/path/to/save/directory/';

% Load Data
filename = 'kaspr_dates.csv';
dates = readtable(filename);
datalist = dir([kasprPPIDataDir, 'KASPR_PPI_SWL_MOMENTS_', '*.nc']);
disp([num2str(length(datalist)), ' KASPR PPI files']);

% Empty arrays which will populate with SWL characteristics
layerHeight_km = [];
layerThickness_m = [];
layerAzimuth_deg = [];
layerMagnitude = [];
profileDateTime = [];
stormDate = [];
stormNum = [];
scanDuration_s = [];

% Define the parameters of the convolution
layerThicknesses = [100, 250, 500]; % thicknesses in m of the desired PPI SWLs
SWLThreshold = 0.25; % SW threshold in m/s

% Loop through each KASPR PPI file to apply the SWL detection algorithm
for i = 1:length(datalist) 
    kasprdata = [datalist(i).folder,'/', datalist(i).name];

    % Step 1: Load in KASPR PPI data and extract datetime information
    [timeh, times, ref, spw, snr, rangekm, xkm, ykm, zkm,...
        elev_deg, az_deg, file_duration_s] = PPI_kaspr_variables(kasprdata);
    profileDate_file = kasprdata(end-17:end-10); % scan date
    profileDateTime_file = kasprdata(end-17:end-3); % scan date and time
    [~,dz] = gradient(zkm); % vertical gradient of altitude in km
    verticalRes = mode(dz, 'all') * 1000; % vertical gradient of altitude m. This is not range resolution.
    verticalRes = repmat(verticalRes, size(spw));

    % Step 2: Perform the convolution and identify PPI SWLs
    [layerNumber] = PPI_convolution(layerThicknesses, spw, verticalRes, SWLThreshold);

    % Step 3: Compute PPI SWL properties
    [layerHeight_list, layerThickness_list, layerAzimuth_list, layerMagnitude_list] = ...
    PPI_SWL_properties(layerNumber, spw, zkm, az_deg, verticalRes);

    % Step 4: Append SWL properties from current KASPR PPI scan to main
    % lists of PPI SWL characyeristics
    layerHeight_km = cat(1, layerHeight_km, layerHeight_list); 
    layerThickness_m = cat(1, layerThickness_m, layerThickness_list); 
    layerAzimuth_deg = cat(1, layerAzimuth_deg, layerAzimuth_list); 
    layerMagnitude = cat(1, layerMagnitude, layerMagnitude_list); 
    file_duration_s_list = repmat(file_duration_s, [length(layerHeight_list),1]);
    scanDuration_s = cat(1, scanDuration_s, file_duration_s_list);
    dateTime_list = repmat(profileDateTime_file, [length(layerHeight_list),1]);
    profileDateTime = cat(1, profileDateTime, dateTime_list);
    stormDate_list = repmat(profileDate_file, [length(layerHeight_list),1]);
    stormDate = cat(1, stormDate, stormDate_list);
    loc = string(dates.Dates) == profileDate_file;
    stormNum_profile = dates.StormNum(loc, :);
    stormNum_list = repmat(stormNum_profile, [length(layerHeight_list),1]);
    stormNum = cat(1, stormNum, stormNum_list);
end

% Save PPI SWL dataset as a .csv file
T = table(profileDateTime, stormDate, stormNum, layerHeight_km, layerThickness_m, layerMagnitude,...
    layerAzimuth_deg, scanDuration_s);
writetable(T, [data_save_dir, 'PPI_SWL_climatology.csv']);

disp('Done');