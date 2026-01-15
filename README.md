# SWL_Detection
Percentile-based detection algorithm designed to detect locally enhanced, azimuthally continuous SW features within Northeast U.S. winter precipitation events. The SWL identification methodology relies on a percentile-based detection algorithm designed to detect locally enhanced, azimuthally continuous SW features. For each range gate, the 75th percentile of SW was calculated using centered moving windows with vertical (not range) depths of 100, 250, 500, 750, and 1000 m. Multiple vertical thicknesses were used to ensure that  SW enhancements of varying vertical thicknesses were identified. A binary operator is used to conjoin SW enhancements exceeding the 75th percentile into individual SWLs. The algorithm identified 77,955 SWLs in KASPR PPI scans over four winter seasons. A Velocity Azimuth Display (VAD) technique is applied to decompose the flow into resolved and unresolved components, from which proxies for shear and turbulence are derived.

This repository contains MATLAB scripts used to detect and analyze spectrum width layers (SWLs) from Ka-Band radar observations collected by the Ka-Band Scanning Polarimetric Radar (KASPR) at Stony Brook University. The scripts support the study presented in Leghart et al. 2026 (Monthly Weather Review, in review), which provides the first climatology of SWLs within Northeast U.S. winter precipitation events.
### Repository Contents
- PPI_SWL_Climatology.m: Main script for identifying and extracting SWLs from KASPR plan position indicator (PPI) scans. Outputs VARIABLES LIST
- Supporting functions (LIST THEM ALL) are required to run the full workflow.
- kaspr_dates.csv provides required information for event dates, numbers, and types.

### Required Input Data
- KASPR_PPI_SWL_MOMENTS_yyyymmdd-HHMMSS.nc: The KASPR moments files are required to run PPI_SWL_Climatology.m
- identified_melting_layers_KASPR_PPI.csv: A table containing the dates, times, altitudes, and thicknesses of melting layers within the KASPR 15 degree PPI dataset. This table is used to remove any SWLs which extend into or are contained within the melting layer, as melting influences SW (see Leghart et al. 2026 for more details).
- kaspr_dates.csv: A table which provides required information for event dates, numbers, and types.

### Workflow Summary
Add the workflow I have listed within the MATLAB script

### Associated Publications
- Leghart, E. C., et al. (2026, in review): Characteristics of Layers of Enhanced Spectrum Width within Northeast United States Winter Precipitation Events. Monthly Weather Review.
- Leghart, Erin; Colle, Brian; Oue, Mariko et al. (2025). Layers of enhanced spectrum width within Northeast United States winter precipitation events [Dataset]. Dryad. https://doi.org/10.5061/dryad.5dv41nshp. 
