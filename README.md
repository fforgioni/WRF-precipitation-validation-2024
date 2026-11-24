# WRF Precipitation Validation — April–May 2024

![Model](https://img.shields.io/badge/Model-WRF-blue)
![Language](https://img.shields.io/badge/Language-Python-3776AB)
![Status](https://img.shields.io/badge/Status-In%20progress-yellow)

## Overview

This project evaluates the sensitivity of precipitation simulations performed with the Weather Research and Forecasting (WRF) model to different initial and boundary conditions and cumulus parameterizations during the 28 April–3 May 2024 precipitation event.

## Objectives

- Assess the influence of ERA5 and GFS forcing on simulated precipitation.
- Evaluate the performance of Grell-Freitas and Kain-Fritsch cumulus parameterizations.
- Assess the spatial and temporal representation of precipitation.
- Compare model simulations against gridded precipitation estimates and station observations.

## Experimental design

| Experiment | Forcing | Cumulus scheme |
|---|---|---|
| GF-ERA5 | ERA5 | Grell-Freitas |
| GF-GFS | GFS | Grell-Freitas |
| KF-ERA5 | ERA5 | Kain-Fritsch |
| KF-GFS | GFS | Kain-Fritsch |

**Study period:** 28 April – 3 May 2024.

## Validation datasets

| Dataset | Temporal resolution | Source |
|---|---|---|
| MERGE | Daily | CPTEC/INPE |
| MERGE | Hourly | CPTEC/INPE |
| Station observations | To be documented | 60 meteorological stations |

## Verification methodology

### Spatial verification

- Accumulated precipitation patterns.
- Spatial precipitation differences.
- Grid-point statistical metrics.

### Station-based verification

- Comparison against observed precipitation.
- Station-level performance metrics.
- Spatial distribution of model errors.

### Temporal verification

- Hourly precipitation evolution.
- Domain-averaged precipitation.
- Temporal lag analysis.
- Hovmöller diagrams.

### Statistical verification

- Mean Bias (MB).
- Mean Absolute Error (MAE).
- Root Mean Square Error (RMSE).
- Pearson correlation coefficient.
- Fractions Skill Score (FSS).

## Repository structure

The repository will include:

- Model configuration files.
- Data preprocessing scripts.
- Verification and analysis scripts.
- Figures and statistical results.
- Project documentation.

## Installation

The computational environment is managed using Conda.

Create the environment using:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate wrf-validation
```

Run the preprocessing scripts from the repository root:

```bash
python scripts/01_download_merge.py
python scripts/02_grib_to_netcdf.py
python scripts/03_regrid_merge.py
```

Additional input datasets must be placed in the
corresponding directories described in `data/README.md`.

## Data availability

Data sources and access instructions will be documented as the project progresses.

## Academic context

**Course:** Numerical Atmospheric Modelling

**Institution:** Facultad de Ciencias Astronómicas y Geofísicas (FCAG), Universidad Nacional de La Plata (UNLP), Argentina.

## Project status

Work in progress.
