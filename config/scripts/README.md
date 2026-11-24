# Data Processing and Verification Scripts

This directory contains the Python scripts developed for data acquisition, preprocessing, and verification of WRF precipitation simulations.

## Data preprocessing

| Script | Description |
|---|---|
| 01_download_merge.py | Download hourly MERGE precipitation data from CPTEC/INPE. |
| 02_grib_to_netcdf.py | Convert hourly GRIB2 files into a single NetCDF dataset. |
| 03_regrid_merge.py | Interpolate daily and hourly MERGE precipitation onto the WRF grid. |

## Verification

Spatial, station-based, and hourly precipitation verification scripts will be incorporated as the analysis progresses.

## Programming language

Python 3.
