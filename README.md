# WRF Precipitation Validation

### The April–May 2024 Extreme Rainfall Event in Southern Brazil

![Model](https://img.shields.io/badge/Model-WRF-blue)
![Language](https://img.shields.io/badge/Language-Python-3776AB)
![Forcing](https://img.shields.io/badge/Forcing-ERA5%20%26%20GFS-blue)
![Status](https://img.shields.io/badge/Status-Completed-green)

## Overview

This project investigates the sensitivity of precipitation simulations performed with the Weather Research and Forecasting (WRF) model during the extreme rainfall event that affected southern Brazil between 28 April and 3 May 2024.

The experimental design combines two atmospheric forcing datasets (ERA5 and GFS) and two cumulus parameterization schemes (Grell–Freitas and Kain–Fritsch) to evaluate their influence on the spatial and temporal representation of precipitation.

Model simulations are evaluated against gridded precipitation estimates and meteorological station observations using spatial, temporal, and statistical verification techniques.

## Scientific objectives

- Assess the sensitivity of simulated precipitation to initial and lateral boundary conditions.
- Evaluate the influence of cumulus parameterization on precipitation intensity and distribution.
- Examine the ability of WRF to reproduce the spatial and temporal evolution of extreme precipitation.
- Quantify model performance against gridded precipitation estimates and station observations.
- Identify differences and uncertainties associated with the selected model configurations.

## Study area and event

The study focuses on southern Brazil, particularly the state of Rio Grande do Sul, during the extreme precipitation event of late April and early May 2024.

The analysis covers the period from 28 April to 3 May 2024.

## Experimental design

Four numerical experiments were designed to investigate the sensitivity of precipitation simulations to atmospheric forcing and cumulus parameterization.

| Experiment | Atmospheric forcing | Cumulus parameterization |
|:---|:---:|:---|
| GF-ERA5 | ERA5 | Grell–Freitas |
| GF-GFS | GFS | Grell–Freitas |
| KF-ERA5 | ERA5 | Kain–Fritsch |
| KF-GFS | GFS | Kain–Fritsch |

**Analysis period:** 28 April – 3 May 2024.

## Validation datasets

Model precipitation is evaluated using gridded precipitation products and surface meteorological observations.

| Dataset | Temporal resolution | Source |
|:---|:---:|:---|
| MERGE | Daily | CPTEC/INPE |
| MERGE | Hourly | CPTEC/INPE |
| Meteorological station observations | Station-dependent | 60 meteorological stations |

The combination of gridded and station-based observations provides complementary perspectives on model performance.

## Verification methodology

### Spatial verification

- Accumulated precipitation fields.
- Spatial distribution of precipitation.
- Differences between simulated and observed precipitation.
- Grid-point statistical verification.
- Spatial distribution of model errors.

### Station-based verification

- Comparison of simulated and observed precipitation.
- Station-level statistical metrics.
- Spatial variability of model performance.

### Temporal verification

- Hourly precipitation evolution.
- Domain-averaged precipitation time series.
- Temporal lag analysis.
- Hovmöller diagrams.

### Statistical verification

The verification framework considers the following metrics:

| Metric | Abbreviation |
|:---|:---:|
| Mean Bias | MB |
| Mean Absolute Error | MAE |
| Root Mean Square Error | RMSE |
| Pearson correlation coefficient | r |
| Fractions Skill Score | FSS |

These metrics provide complementary information on systematic errors, precipitation magnitude, temporal correspondence, and spatial agreement.

## Computational workflow

The analysis follows a reproducible computational workflow:

1. Acquisition of atmospheric forcing and precipitation datasets.
2. WRF model preprocessing and configuration.
3. Numerical experiments.
4. Processing of model outputs.
5. Spatial and temporal harmonization of datasets.
6. Precipitation verification.
7. Statistical analysis and scientific visualization.

## Repository structure

The repository contains resources associated with the numerical experiments and precipitation verification workflow, including:

- WRF and WPS configuration files.
- Data preprocessing scripts.
- Precipitation verification scripts.
- Statistical analysis tools.
- Scientific visualization scripts.
- Computational environment specifications.
- Project documentation.

## Installation

The Python computational environment is managed using Conda.

Clone the repository:

```bash
git clone https://github.com/fforgioni/WRF-precipitation-validation-2024.git

cd WRF-precipitation-validation-2024
```

Create the computational environment:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate wrf-validation
```

## Data preprocessing

Run the preprocessing scripts from the repository root:

```bash
python scripts/01_download_merge.py
python scripts/02_grib_to_netcdf.py
python scripts/03_regrid_merge.py
```

Additional input datasets must be placed in the corresponding directories described in `data/README.md`.

## Data availability

The study uses atmospheric forcing datasets and precipitation observations obtained from the following institutions and modelling centres:

- ECMWF — ERA5 atmospheric reanalysis.
- NOAA — Global Forecast System (GFS).
- CPTEC/INPE — MERGE precipitation products.
- Surface meteorological station networks.

Due to the size of the atmospheric datasets and WRF output files, users may need to obtain these data separately.

## Reproducibility

This repository provides a framework for reproducing the numerical experiments, precipitation verification procedures, and scientific visualizations.

The computational environment, model configurations, and processing scripts support transparency and reproducibility.

## Academic context

This project was developed within the framework of the Numerical Atmospheric Modelling course.

**Course:** Numerical Atmospheric Modelling  
**Faculty:** Facultad de Ciencias Astronómicas y Geofísicas (FCAG)  
**University:** Universidad Nacional de La Plata (UNLP)  
**Country:** Argentina

## Acknowledgements

The authors acknowledge the institutions responsible for providing atmospheric forcing datasets, precipitation products, meteorological observations, and the numerical modelling tools used in this project.

## Project status

The numerical experiments and precipitation verification study have been completed.
