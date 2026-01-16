# SWL_Detection
Percentile-based detection algorithm designed to detect locally enhanced, azimuthally continuous Doppler Spectrum Width (SW) features within Northeast U.S. winter precipitation events called SW Layers (SWLs). The SWL identification methodology relies on a percentile-based detection algorithm designed to detect locally enhanced, azimuthally continuous SW features. For each range gate, the 75th percentile of SW was calculated using centered moving windows with vertical (not range) depths of 100, 250, 500, 750, and 1000 m. Multiple vertical thicknesses were used to ensure that  SW enhancements of varying vertical thicknesses were identified. A binary operator is used to conjoin SW enhancements exceeding the 75th percentile into individual SWLs. The algorithm identified 77,955 SWLs in KASPR PPI scans over four winter seasons. A Velocity Azimuth Display (VAD) technique is applied to decompose the flow into resolved and unresolved components, from which proxies for shear and turbulence are derived.

### Repository Contents
This repository contains MATLAB scripts used to detect and analyze spectrum width layers (SWLs) from Ka-Band radar observations collected by the Ka-Band Scanning Polarimetric Radar (KASPR) at Stony Brook University. The scripts support the study presented in Leghart et al. 2026 (Monthly Weather Review, in review), which provides the first climatology of SWLs within Northeast U.S. winter precipitation events.
- PPI_SWL_Climatology.m: Main script for identifying and extracting SWLs from KASPR plan position indicator (PPI) scans. Outputs variables include: profileDateTime, stormDate, stormNum, stormCategory, layerHeight_km, layerThickness_m, layerMagnitude, layerFitVelGrad, layerResidVel, layerAzimuth_deg, layerSquare_km, scanDuration_s, cloudTopHeight, cloudBaseHeight, cloudDepth, within_ml, below_ml
- Supporting functions (PPI_kaspr_variables.m, compute_fourier_coeffs_by_range.m, cloud_depth_PPI.m) are required to run the full workflow.

### Required Input Data
The required data can be found the Dryad dataset cited below.
- KASPR_PPI_SWL_MOMENTS_yyyymmdd-HHMMSS.nc: Files containing the KASPR radar variables needed for the SWL identification algorithm.
- identified_melting_layers_KASPR_PPI.csv: A table containing the dates, times, altitudes, and thicknesses of melting layers within the KASPR 15 degree PPI dataset. This table is used to remove any SWLs which extend into or are contained within the melting layer, as melting influences SW (see Leghart et al. 2026 for more details).
- kaspr_dates.csv: A table which provides required information for event dates, numbers, and types.
- identified_melting_layers_KASPR_PPI.csv: A table which provides required information regarding melting levels in KASPR PPI scans. This is used for the removal of SWLs contaminated by melting.

### Workflow Summary
1. Load KASPR_PPI_SWL_moments_YYYYMMDD-HHmmss.nc file & perform preliminary QC metrics
2. VAD Fitting: This consistists of the function "compute_fourier_coeffs_by_range" which performs a harmonic decomposition of the observed doppler velocity (vel) to be used for the VAD fitting.
3. Cloud Base, Top, and Depth Estimation
4. SWL Detection & QC
- 4.1 Data Prep
- 4.2 Calculate Moving Percentile Windows
- 4.3 Compare Moving Percentile Windows to SW Field
- 4.4 Binary Label Operator & Identification of SWLs
- 4.5 SWL Characteristics & Noise Removal
- 4.6 - Append SWLs to Growing Lists
5. Remove SWLs in the Melting Layer
6. Save SWL dataset as PPI_SWL_Climatology.csv


### Associated Publications
- Leghart, E. C., et al. (2026, in review): Characteristics of Layers of Enhanced Spectrum Width within Northeast United States Winter Precipitation Events. Monthly Weather Review.
- Leghart, Erin; Colle, Brian; Oue, Mariko et al. (2025). Layers of enhanced spectrum width within Northeast United States winter precipitation events [Dataset]. Dryad. https://doi.org/10.5061/dryad.5dv41nshp. 
