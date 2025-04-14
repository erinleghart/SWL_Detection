# SWL_Detection
Convolution-based detection algorithm for enhanced spectrum width layers (SWLs) within the KASPR radar in Northeast U.S. coastal winter storms

This repository contains MATLAB scripts used to detect and analyze spectrum width layers (SWLs) from Ka-Band radar observations collected by the Ka-Band Scanning Polarimetric Radar (KASPR) at Stony Brook University. The scripts support the study presented in Leghart et al. (2025, in review), which provides the first climatology of SWLs within cool-season coastal winter storms.
### Repository Contents
- PPI_SWL_Climatology.m: Main script for identifying and extracting SWLs from KASPR plan position indicator (PPI) scans. Outputs SWL height, thickness, azimuth, magnitude, and scan duration.
- VPT_SWL_Climatology.m: Main script for identifying and extracting SWLs from vertically pointing time-height (VPT) radar profiles. Outputs SWL height, thickness, magnitude, and duration.
- Supporting functions (PPI_convolution.m, VPT_convolution.m, kaspr_variables_PPI.m, kaspr_variables_VPT.m, PPI_SWL_properties.m, and VPT_SWL_properties.m) are required to run the full workflow.

### Input Data Format
Each script expects input radar files in NetCDF format named using the convention:
- KASPR_PPI_SWL_MOMENTS_yyyymmdd-HHMMSS.nc for PPI scans
- KASPR_VPT_SWL_MOMENTS_yyyymmdd-HHMMSS.nc for VPT profiles
These files must be stored in subdirectories (kaspr_ppi/ and kaspr_vpt/) relative to the location of the MATLAB scripts. You will also need a CSV file kaspr_dates.csv with storm metadata (e.g., date and storm number).

### Workflow Summary
1. Load radar and metadata: Each script reads in relevant NetCDF files and parses the observation timestamp.
2. Apply SWL detection algorithm: A 2D convolution process is applied to the Doppler spectrum width (SW) field using multiple vertical kernel thicknesses (100, 250, and 500 m). SWL identification is based on local SW enhancement thresholds of 0.25 m/s (PPI) and 0.20 m/s (VPT).
3. Label and characterize SWLs: SW enhancements exceeding the threshold are grouped into contiguous regions and treated as individual SWLs.
4. Export climatology: SWL properties are saved to .csv output files for analysis and visualization.

### Associated Publications
If you use this code or dataset in your work, please cite the associated manuscript(s):
- Leghart, E. C., et al. (2025, in review): Climatology of Spectrum Width Layers in Northeast U.S. Winter Storms Using KASPR Radar Observations. Journal of Applied Meteorology and Climatology.
