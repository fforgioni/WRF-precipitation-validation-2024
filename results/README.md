# Verification

This directory contains the numerical and statistical
results of the WRF precipitation validation study
for the April–May 2024 extreme rainfall event
in Rio Grande do Sul, Brazil.

## Available results

| File | Description |
|---|---|
| Diagnosticos_horarios_WRF_MERGE_Gauge_CC.xlsx | Hourly precipitation diagnostics and correlations |
| Fractions_Skill_Score_Daily.xlsx | Daily Fractions Skill Score |
| Metricas.xlsx | Station-based verification metrics |
| Metricas_20240430_20mm.xlsx | Categorical verification at 20 mm/day |
| Metricas_20240430_50mm.xlsx | Categorical verification at 50 mm/day |
| Metricas_20240430_100mm.xlsx | Categorical verification at 100 mm/day |
| Precipitacion_acumulada_estaciones.xlsx | Event-total precipitation at meteorological stations |

## Numerical experiments

- Grell–Freitas driven by ERA5 (GF-ERA5).
- Grell–Freitas driven by GFS (GF-GFS).
- Kain–Fritsch driven by ERA5 (KF-ERA5).
- Kain–Fritsch driven by GFS (KF-GFS).

## Reproducibility

The analysis scripts are available in the
[scripts directory](../scripts/).

Intermediate numerical arrays and large NetCDF
files are not included in this repository.
