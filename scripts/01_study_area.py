# -*- coding: utf-8 -*-
"""Plot Rio Grande do Sul topography and meteorological station locations.

Author: Fernando Primo Forgioni (ORCID: 0000-0003-3233-8957).
Input: data/Estaciones.xlsx, Sheet1, with latitud and longitud columns.
External data: NOAA ETOPO 2022 and Natural Earth (internet access required
for downloads not already cached).
Output: figures/study_domain.png.
Run: python scripts/01_study_area.py
"""

import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import requests
import shapely
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from io import BytesIO
from cartopy.io import shapereader
from matplotlib.colors import ListedColormap, BoundaryNorm

# =====================================================
# 1. LEER DATOS DE ESTACIONES
# =====================================================

repository_root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--data-dir", type=Path, default=repository_root / "data")
parser.add_argument("--figures-dir", type=Path, default=repository_root / "figures")
args = parser.parse_args()
carpeta = args.data_dir
args.figures_dir.mkdir(parents=True, exist_ok=True)
archivo = os.path.join(carpeta, "Estaciones.xlsx")

if not Path(archivo).is_file():
    raise FileNotFoundError(f"Station input file not found: {archivo}")

df = pd.read_excel(archivo, sheet_name="Sheet1")

df.columns = df.columns.str.strip().str.lower()

df["latitud"] = pd.to_numeric(df["latitud"], errors="coerce")
df["longitud"] = pd.to_numeric(df["longitud"], errors="coerce")

df = df.dropna(subset=["latitud", "longitud"])

# =====================================================
# 3. CONFIGURACIÓN GRÁFICA
# =====================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.linewidth": 0.8
})

proyeccion = ccrs.PlateCarree()

fig, axs = plt.subplots(
    1, 2,
    figsize=(16, 8),
    subplot_kw={"projection": proyeccion}
)

ax1, ax2 = axs

extent = [-58.5, -48.5, -34.5, -26.5]

for ax in axs:
    ax.set_extent(extent, crs=proyeccion)

# =====================================================
# 4. LÍMITES GEOGRÁFICOS
# =====================================================

estados = cfeature.NaturalEarthFeature(
    category="cultural",
    name="admin_1_states_provinces_lines",
    scale="10m",
    facecolor="none"
)

shp = shapereader.natural_earth(
    resolution="10m",
    category="cultural",
    name="admin_1_states_provinces"
)

registros = shapereader.Reader(shp).records()

rs = [
    registro.geometry
    for registro in registros
    if registro.attributes.get("name") == "Rio Grande do Sul"
]

if not rs:
    raise ValueError(
        "No se encontró el límite de Rio Grande do Sul."
    )

geometria_rs = shapely.union_all(rs)

# =====================================================
# 5. OBTENER TOPOGRAFÍA ETOPO 2022
# =====================================================

print("Descargando datos de elevación de NOAA...")

url = (
    "https://coastwatch.pfeg.noaa.gov/erddap/griddap/"
    "ETOPO_2022_v1_15s.nc"
)

consulta = (
    "z[(-34.5):1:(-26.5)]"
    "[(-58.5):1:(-48.5)]"
)

respuesta = requests.get(
    url + "?" + consulta,
    timeout=300
)

respuesta.raise_for_status()

with xr.open_dataset(
    BytesIO(respuesta.content),
    engine="scipy"
) as ds:

    lon = ds["longitude"].values.copy()
    lat = ds["latitude"].values.copy()

    elevacion = ds["z"].transpose(
        "latitude", "longitude"
    ).values.astype(float)

print("Datos de elevación descargados correctamente.")

# =====================================================
# 6. RECORTAR TOPOGRAFÍA A RIO GRANDE DO SUL
# =====================================================

LON, LAT = np.meshgrid(lon, lat)

shapely.prepare(geometria_rs)

mascara = shapely.contains_xy(
    geometria_rs,
    LON,
    LAT
)

elevacion = np.where(
    mascara,
    np.maximum(elevacion, 0),
    np.nan
)

# =====================================================
# 7. PALETA HIPSOMÉTRICA
# =====================================================

niveles = [
    0, 90, 180, 270, 360, 450,
    540, 630, 720, 810, 900, 990
]

colores = [
    "#185c25",
    "#30853a",
    "#55ab42",
    "#88cc45",
    "#d4e846",
    "#ffff33",
    "#ffdb27",
    "#ffae19",
    "#ff8015",
    "#ed4014",
    "#b71919",
    "#681414"
]

cmap = ListedColormap(colores[:-1])

cmap.set_over(colores[-1])
cmap.set_bad((0, 0, 0, 0))

norm = BoundaryNorm(
    niveles,
    cmap.N,
    clip=False
)

# =====================================================
# 8. PANEL A: REGIÓN DE ESTUDIO
# =====================================================

ax1.add_feature(
    cfeature.LAND.with_scale("10m"),
    facecolor="#eeeeee",
    zorder=0
)

ax1.add_feature(
    cfeature.OCEAN.with_scale("10m"),
    facecolor="#dcecf5",
    zorder=0
)

# Topografía

topografia = ax1.pcolormesh(
    lon,
    lat,
    elevacion,
    cmap=cmap,
    norm=norm,
    shading="auto",
    transform=proyeccion,
    rasterized=True,
    zorder=1
)

# Lagos y lagunas

ax1.add_feature(
    cfeature.LAKES.with_scale("10m"),
    facecolor="#dcecf5",
    edgecolor="#888888",
    linewidth=0.3,
    zorder=2
)

# Costas

ax1.add_feature(
    cfeature.COASTLINE.with_scale("10m"),
    linewidth=0.6,
    zorder=3
)

# Límites internacionales

ax1.add_feature(
    cfeature.BORDERS.with_scale("10m"),
    linewidth=0.6,
    edgecolor="#666666",
    zorder=3
)

# Límites estatales

ax1.add_feature(
    estados,
    edgecolor="#999999",
    linewidth=0.4,
    zorder=3
)

# Límite de Rio Grande do Sul

ax1.add_geometries(
    rs,
    crs=proyeccion,
    facecolor="none",
    edgecolor="#222222",
    linewidth=0.9,
    zorder=4
)

# =====================================================
# 9. COLORBAR HORIZONTAL SOBRE EL OCÉANO
# =====================================================

cax = ax1.inset_axes([
    0.56,   # Posición horizontal
    0.10,   # Posición vertical
    0.39,   # Ancho
    0.025   # Altura
])

cbar = fig.colorbar(
    topografia,
    cax=cax,
    orientation="horizontal",
    ticks=[0, 180, 360, 540, 720, 900],
    extend="max"
)

cbar.ax.tick_params(
    labelsize=9,
    length=2,
    width=0.5,
    pad=2
)

cbar.set_label(
    "Elevation (m)",
    fontsize=10,
    labelpad=4
)

cbar.outline.set_linewidth(0.5)

# Título centrado

ax1.set_title(
    "(a) Study area",
    fontsize=14,
    fontweight="bold",
    loc="center",
    pad=12
)

# =====================================================
# 10. PANEL B: DISTRIBUCIÓN DE ESTACIONES
# =====================================================

ax2.add_feature(
    cfeature.LAND.with_scale("10m"),
    facecolor="#f5f5f5",
    zorder=0
)

ax2.add_feature(
    cfeature.OCEAN.with_scale("10m"),
    facecolor="#dcecf5",
    zorder=0
)

ax2.add_feature(
    cfeature.LAKES.with_scale("10m"),
    facecolor="#dcecf5",
    edgecolor="gray",
    linewidth=0.4,
    zorder=2
)

ax2.add_feature(
    cfeature.COASTLINE.with_scale("10m"),
    linewidth=0.8,
    zorder=3
)

ax2.add_feature(
    cfeature.BORDERS.with_scale("10m"),
    linewidth=0.8,
    zorder=3
)

ax2.add_feature(
    estados,
    edgecolor="gray",
    linewidth=0.6,
    zorder=3
)

# Límite estadual

ax2.add_geometries(
    rs,
    crs=proyeccion,
    facecolor="none",
    edgecolor="#444444",
    linewidth=0.9,
    zorder=4
)

# =====================================================
# 11. ESTACIONES METEOROLÓGICAS
# =====================================================

ax2.scatter(
    df["longitud"],
    df["latitud"],
    s=45,
    marker="o",
    color="#708090",
    edgecolor="black",
    linewidth=0.5,
    transform=proyeccion,
    zorder=5
)

# =====================================================
# 12. TÍTULO DEL PANEL B
# =====================================================

ax2.set_title(
    "(b) Meteorological stations",
    fontsize=14,
    fontweight="bold",
    loc="center",
    pad=12
)

# =====================================================
# 13. AJUSTAR FIGURA
# =====================================================

plt.tight_layout()

fig.subplots_adjust(
    wspace=0.04
)

# =====================================================
# 14. GUARDAR FIGURA
# =====================================================

salida = os.path.join(
    args.figures_dir,
    "study_domain.png"
)

plt.savefig(
    salida,
    dpi=600,
    bbox_inches="tight",
    facecolor="white"
)

print(f"Figura guardada en: {salida}")

plt.show()