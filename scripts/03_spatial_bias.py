# -*- coding: utf-8 -*-
"""Map accumulated precipitation bias (WRF minus MERGE) for four experiments.

Author: Fernando Primo Forgioni (ORCID: 0000-0003-3233-8957).
Inputs: data/Grell-Freitas-ERA5.nc and accumulated arrays produced by
02_accumulated_precipitation.py in results/Mapas_acumulados_evento/.
Outputs: figures/Accumulated_BIAS.png and
results/Mapas_BIAS_espacial/BIAS_<experiment>.npy.
Run: python scripts/03_spatial_bias.py
Natural Earth geographic boundaries may require internet access.
"""

import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.ticker as mticker

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from cartopy.io import shapereader

# =====================================================
# 1. CONFIGURACIÓN
# =====================================================

repository_root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--data-dir", type=Path, default=repository_root / "data")
parser.add_argument("--figures-dir", type=Path, default=repository_root / "figures")
parser.add_argument("--results-dir", type=Path, default=repository_root / "results")
args = parser.parse_args()

carpeta = args.data_dir
carpeta_acumulados = args.results_dir / "Mapas_acumulados_evento"
carpeta_salida = args.results_dir / "Mapas_BIAS_espacial"
archivo_referencia = carpeta / "Grell-Freitas-ERA5.nc"

required_inputs = [archivo_referencia] + [
    carpeta_acumulados / f"Acumulado_{name}.npy"
    for name in ("MERGE", "GF-ERA5", "GF-GFS", "KF-ERA5", "KF-GFS")
]
for input_path in required_inputs:
    if not input_path.is_file():
        raise FileNotFoundError(
            f"Required input file not found: {input_path}. "
            "Run 02_accumulated_precipitation.py first with matching directories."
        )

args.figures_dir.mkdir(parents=True, exist_ok=True)
carpeta_salida.mkdir(parents=True, exist_ok=True)

orden = [
    "GF-ERA5",
    "GF-GFS",
    "KF-ERA5",
    "KF-GFS"
]

titulos = [
    "(a) GF-ERA5",
    "(b) GF-GFS",
    "(c) KF-ERA5",
    "(d) KF-GFS"
]

extent = [-58.5, -48.5, -34.5, -26.5]

proyeccion = ccrs.PlateCarree()

# =====================================================
# 2. CONFIGURACIÓN GRÁFICA
# =====================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 13,
    "axes.linewidth": 0.8
})

fontsize_titulos = 15
fontsize_ticks = 14
fontsize_colorbar = 16

# Escala simétrica del BIAS

vmin = -300
vmax = 300
intervalo = 50

niveles = np.arange(
    vmin,
    vmax + intervalo,
    intervalo
)

ticks = niveles

cmap = plt.get_cmap(
    "RdBu_r"
).copy()

cmap.set_bad("#D9D9D9")

norm = mcolors.BoundaryNorm(
    niveles,
    cmap.N,
    clip=False
)

# =====================================================
# 3. LEER COORDENADAS WRF
# =====================================================

with xr.open_dataset(archivo_referencia) as ds:

    lat = ds["XLAT"]
    lon = ds["XLONG"]

    for dim in ["Time", "time"]:

        if dim in lat.dims:
            lat = lat.isel({dim: 0})

        if dim in lon.dims:
            lon = lon.isel({dim: 0})

    lat, lon = xr.broadcast(lat, lon)

    lon = lon.transpose(*lat.dims)

    lat_ref = lat.values
    lon_ref = lon.values

print(f"Dimensiones de la grilla: {lat_ref.shape}")

# =====================================================
# 4. LEER PRECIPITACIÓN ACUMULADA MERGE
# =====================================================

merge = np.load(
    os.path.join(
        carpeta_acumulados,
        "Acumulado_MERGE.npy"
    )
)

if merge.shape != lat_ref.shape:

    raise ValueError(
        "Las dimensiones de MERGE no coinciden "
        "con la grilla WRF."
    )

# =====================================================
# 5. CALCULAR BIAS ESPACIAL
# =====================================================

campos_bias = {}

for nombre in orden:

    print(f"\nProcesando {nombre}...")

    campo = np.load(
        os.path.join(
            carpeta_acumulados,
            f"Acumulado_{nombre}.npy"
        )
    )

    if campo.shape != merge.shape:

        raise ValueError(
            f"La grilla de {nombre} no coincide con MERGE."
        )

    # Máscara de valores válidos

    validos = (
        np.isfinite(campo) &
        np.isfinite(merge) &
        (campo >= 0) &
        (merge >= 0)
    )

    # BIAS acumulado

    bias = np.where(
        validos,
        campo - merge,
        np.nan
    )

    campos_bias[nombre] = bias

    print(f"BIAS mínimo: {np.nanmin(bias):.2f} mm")
    print(f"BIAS máximo: {np.nanmax(bias):.2f} mm")

# =====================================================
# 6. LÍMITES GEOGRÁFICOS
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

rs = [
    registro.geometry
    for registro in shapereader.Reader(shp).records()
    if registro.attributes.get("name") == "Rio Grande do Sul"
]

if not rs:

    raise ValueError(
        "No se encontró Rio Grande do Sul."
    )

# =====================================================
# 7. FUNCIÓN PARA GENERAR MAPA BASE
# =====================================================

def mapa_base(ax):

    ax.set_extent(
        extent,
        crs=proyeccion
    )

    ax.add_feature(
        cfeature.LAND.with_scale("10m"),
        facecolor="#F5F5F5",
        zorder=0
    )

    ax.add_feature(
        cfeature.OCEAN.with_scale("10m"),
        facecolor="#DCECF5",
        zorder=0
    )

    ax.add_feature(
        cfeature.LAKES.with_scale("10m"),
        facecolor="#DCECF5",
        edgecolor="gray",
        linewidth=0.3,
        zorder=1
    )

    ax.add_feature(
        cfeature.COASTLINE.with_scale("10m"),
        linewidth=0.6,
        zorder=3
    )

    ax.add_feature(
        cfeature.BORDERS.with_scale("10m"),
        linewidth=0.6,
        edgecolor="#555555",
        zorder=3
    )

    ax.add_feature(
        estados,
        edgecolor="#999999",
        linewidth=0.4,
        zorder=3
    )

    # Destacar Rio Grande do Sul

    ax.add_geometries(
        rs,
        crs=proyeccion,
        facecolor="none",
        edgecolor="#333333",
        linewidth=0.9,
        zorder=4
    )

# =====================================================
# 8. GENERAR FIGURA
# =====================================================

fig, axs = plt.subplots(
    2, 2,
    figsize=(13, 10),
    subplot_kw={"projection": proyeccion}
)

axs = axs.flatten()

for i, nombre in enumerate(orden):

    ax = axs[i]

    mapa_base(ax)

    campo = campos_bias[nombre]

    mesh = ax.pcolormesh(
        lon_ref,
        lat_ref,
        campo,
        cmap=cmap,
        norm=norm,
        shading="auto",
        transform=proyeccion,
        zorder=2
    )

    ax.set_title(
        titulos[i],
        fontsize=fontsize_titulos,
        fontweight="bold",
        loc="center",
        pad=10
    )

fig.subplots_adjust(
    left=0.04,
    right=0.86,
    bottom=0.04,
    top=0.96,
    wspace=0.06,
    hspace=0.10
)

# =====================================================
# 9. COLORBAR COMPARTIDA
# =====================================================

sm = cm.ScalarMappable(
    cmap=cmap,
    norm=norm
)

sm.set_array([])

cax = fig.add_axes([
    0.89,
    0.15,
    0.022,
    0.70
])

cbar = fig.colorbar(
    sm,
    cax=cax,
    orientation="vertical",
    boundaries=niveles,
    ticks=ticks,
    spacing="uniform",
    extend="both",
    format=mticker.FormatStrFormatter("%.0f")
)

cbar.set_label(
    "Accumulated precipitation bias (mm)",
    fontsize=fontsize_colorbar,
    labelpad=15
)

cbar.ax.tick_params(
    axis="y",
    which="major",
    labelsize=fontsize_ticks,
    length=4,
    width=0.8,
    pad=5
)

cbar.ax.minorticks_off()

cbar.outline.set_linewidth(0.8)

# =====================================================
# 10. GUARDAR FIGURA
# =====================================================

salida = os.path.join(
    args.figures_dir,
    "Accumulated_BIAS.png"
)

fig.savefig(
    salida,
    dpi=600,
    format="png",
    bbox_inches="tight",
    facecolor="white"
)

print(f"\nFigura guardada: {salida}")

plt.show()
plt.close(fig)

# =====================================================
# 11. GUARDAR CAMPOS DE BIAS
# =====================================================

for nombre in orden:

    np.save(
        os.path.join(
            carpeta_salida,
            f"BIAS_{nombre}.npy"
        ),
        campos_bias[nombre]
    )

print("\nProcesamiento finalizado.")
