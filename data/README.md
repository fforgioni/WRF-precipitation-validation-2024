# Data availability and requirements

This directory contains the input datasets required
for the WRF precipitation validation study.

Large NetCDF files are not distributed through this
repository due to their size.

## Required datasets

### WRF simulations

| Filename | Experiment |
|---|---|
| Grell-Freitas-ERA5.nc | GF-ERA5 |
| Grell-Freitas-GFS.nc | GF-GFS |
| Kain-Fritsch-ERA5.nc | KF-ERA5 |
| Kain-Fritsch-GFS.nc | KF-GFS |

Required variables:

- RAINC: Accumulated convective precipitation.
- RAINNC: Accumulated non-convective precipitation.
- XLAT: Latitude.
- XLONG: Longitude.
- Time information.

Nominal horizontal resolution: 12 km.

### Reference precipitation datasets

| Filename | Temporal resolution | Source |
|---|---|---|
| MERGE_HORARIO_WRF.nc | Hourly | CPTEC/INPE |
| Merge-diario.nc | Daily | CPTEC/INPE |
| PERSIANN_WRF.nc | Hourly | CHRS |

The files used in the spatial verification must be
interpolated onto the WRF computational grid.

### Meteorological station observations

Required file:

`Estaciones.xlsx`

Required worksheet:

`Sheet1`

Expected columns:

- nombre
- institucion
- latitud
- longitud
- Daily precipitation observations identified by date.

Station precipitation must be expressed in millimetres.

## Temporal conventions

The main verification period extends from
28 April 2024 at 12 UTC to 3 May 2024 at 12 UTC.

Daily precipitation is evaluated using 24-hour
accumulation windows from 12 UTC to 12 UTC.

The processing scripts use the following
timestamp conventions:

- Hourly MERGE: end of the accumulation interval.
- Daily MERGE: end of the accumulation interval.
- Hourly PERSIANN: beginning of the accumulation interval.
- Station observations: beginning date of the daily interval.

These conventions must be verified against the
original dataset metadata before running the scripts.

## Topographic data

Topographic information is obtained from the
NOAA ETOPO 2022 global relief model.

Geographical boundaries are obtained from
Natural Earth.

The study area script downloads the required
topographic information automatically.

## Data directory structure

```text
data/
├── README.md
├── Estaciones.xlsx
├── Grell-Freitas-ERA5.nc
├── Grell-Freitas-GFS.nc
├── Kain-Fritsch-ERA5.nc
├── Kain-Fritsch-GFS.nc
├── MERGE_HORARIO_WRF.nc
├── Merge-diario.nc
└── PERSIANN_WRF.nc
```

## Data sources

- ERA5: https://cds.climate.copernicus.eu/
- GFS: https://www.ncei.noaa.gov/
- MERGE: https://www.cptec.inpe.br/
- PERSIANN: https://chrsdata.eng.uci.edu/
- NOAA ETOPO: https://www.ncei.noaa.gov/products/etopo-global-relief-model
- Natural Earth: https://www.naturalearthdata.com/

## Data redistribution

Input datasets remain subject to the terms and
conditions established by their respective providers.
