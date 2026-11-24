# Data Description and Availability

This directory documents the datasets used to evaluate WRF
precipitation simulations during the 28 April – 3 May 2024 event.

## 1. WRF simulations

Four numerical experiments were designed to evaluate the
sensitivity of simulated precipitation to atmospheric forcing
and cumulus parameterizations.

| Experiment | Forcing | Cumulus scheme | File |
|---|---|---|---|
| GF-ERA5 | ERA5 | Grell-Freitas | Grell-Freitas-ERA5.nc |
| GF-GFS | GFS | Grell-Freitas | Grell-Freitas-GFS.nc |
| KF-ERA5 | ERA5 | Kain-Fritsch | Kain-Fritsch-ERA5.nc |
| KF-GFS | GFS | Kain-Fritsch | Kain-Fritsch-GFS.nc |

### Atmospheric forcing datasets

- ERA5: ECMWF atmospheric reanalysis.
- GFS: NOAA Global Forecast System.

Official resources:

- [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/)
- [NOAA Global Forecast System](https://www.ncei.noaa.gov/products/weather-climate-models/global-forecast)

## 2. MERGE precipitation

MERGE precipitation estimates are provided by CPTEC/INPE,
Brazil.

Two products are used for model verification:

| Product | Temporal resolution | Spatial resolution |
|---|---|---|
| MERGE Daily | 24 hours | 0.1° |
| MERGE Hourly | 1 hour | 0.1° |

### Hourly dataset

- Start date: 2024-04-28 00:00 UTC.
- End date: 2024-05-03 23:00 UTC.
- Number of hourly records: 144.
- Original format: GRIB2.
- Processed format: NetCDF.

### Data access

Official CPTEC/INPE repository:

https://ftp.cptec.inpe.br/modelos/tempo/MERGE/GPM/

Hourly data can be downloaded automatically using:

```bash
python scripts/01_download_merge.py
```

GRIB2 files can be converted to NetCDF using:

```bash
python scripts/02_grib_to_netcdf.py
```

Daily and hourly precipitation can be interpolated onto
the WRF grid using:

```bash
python scripts/03_regrid_merge.py
```

## 3. Meteorological station observations

Precipitation observations from 60 meteorological stations
are used for point-based verification.

These observations provide a complementary reference
for evaluating both WRF simulations and MERGE estimates.

Station locations, data providers, temporal resolution,
quality-control procedures, and access conditions will
be documented during the verification stage.

## 4. Spatial harmonization

All MERGE precipitation fields are interpolated onto the
common WRF spatial grid using bilinear interpolation.

The original MERGE timestamps are preserved.

Temporal accumulation intervals must be harmonized before
performing precipitation verification.

## 5. Local data structure

```text
data/
├── README.md
├── raw/
│   ├── wrf/
│   │   ├── Grell-Freitas-ERA5.nc
│   │   ├── Grell-Freitas-GFS.nc
│   │   ├── Kain-Fritsch-ERA5.nc
│   │   └── Kain-Fritsch-GFS.nc
│   ├── merge_hourly/
│   ├── Merge-diario.nc
│   └── stations/
└── processed/
    ├── MERGE_HORARIO_20240428_20240503.nc
    └── regridded/
        ├── MERGE_DIARIO_WRF.nc
        └── MERGE_HORARIO_WRF.nc
```

## 6. Data availability

Large numerical model outputs and raw meteorological
datasets are not directly included in this GitHub repository.

External storage and permanent access links will be
documented upon project completion.

Redistribution of third-party datasets is subject to
their respective licensing and access conditions.
