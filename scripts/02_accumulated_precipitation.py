# -*- coding: utf-8 -*-
"""Map event-total precipitation for four WRF experiments, PERSIANN and MERGE.

Author: Fernando Primo Forgioni (ORCID: 0000-0003-3233-8957).
Period: 2024-04-28 12:00 to 2024-05-03 12:00 UTC.
Inputs: four WRF NetCDF files, PERSIANN_WRF.nc and MERGE_HORARIO_WRF.nc.
Outputs: figures/accumulated_precipitation.png and
results/Mapas_acumulados_evento/Acumulado_<dataset>.npy.
Run: python scripts/02_accumulated_precipitation.py
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
carpeta_salida = args.results_dir / "Mapas_acumulados_evento"
args.figures_dir.mkdir(parents=True, exist_ok=True)

os.makedirs(carpeta_salida, exist_ok=True)

simulaciones = {
    "GF-ERA5": "Grell-Freitas-ERA5.nc",
    "GF-GFS": "Grell-Freitas-GFS.nc",
    "KF-ERA5": "Kain-Fritsch-ERA5.nc",
    "KF-GFS": "Kain-Fritsch-GFS.nc"
}

archivo_persiann = "PERSIANN_WRF.nc"
archivo_merge = "MERGE_HORARIO_WRF.nc"

# Check all required inputs before loading geographic resources.
for filename in [*simulaciones.values(), archivo_persiann, archivo_merge]:
    input_path = carpeta / filename
    if not input_path.is_file():
        raise FileNotFoundError(f"Required input file not found: {input_path}")

fecha_inicio = pd.Timestamp("2024-04-28 12:00:00")
fecha_final = pd.Timestamp("2024-05-03 12:00:00")

# Convenciones temporales
# True: timestamp representa final del intervalo
# False: timestamp representa inicio del intervalo

persiann_fin_hora = False
merge_fin_hora = True

# Dominio del mapa
extent = [-58.5, -48.5, -34.5, -26.5]

# =====================================================
# 2. COLORES Y CONFIGURACIÓN GRÁFICA
# =====================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.linewidth": 0.8
})

proyeccion = ccrs.PlateCarree()

fontsize_titulos = 13
fontsize_ticks = 14
fontsize_colorbar = 16

# Escala de precipitación acumulada
niveles = np.arange(0, 451, 50)
ticks = np.arange(0, 451, 50)

cmap = plt.get_cmap("YlGnBu", len(niveles) - 1)

norm = mcolors.BoundaryNorm(
    niveles,
    cmap.N,
    clip=False
)

# =====================================================
# 3. TÍTULOS DE LOS PANELES
# =====================================================

orden = [
    "GF-ERA5",
    "GF-GFS",
    "KF-ERA5",
    "KF-GFS",
    "PERSIANN",
    "MERGE"
]

titulos = [
    "(a) Grell-Freitas / ERA5",
    "(b) Grell-Freitas / GFS",
    "(c) Kain-Fritsch / ERA5",
    "(d) Kain-Fritsch / GFS",
    "(e) PERSIANN",
    "(f) MERGE"
]

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

rs = [
    registro.geometry
    for registro in shapereader.Reader(shp).records()
    if registro.attributes.get("name") == "Rio Grande do Sul"
]

if not rs:
    raise ValueError("No se encontró Rio Grande do Sul.")

# =====================================================
# 5. FUNCIONES AUXILIARES
# =====================================================

def leer_fechas(ds):

    if "Times" in ds.variables:

        fechas = []

        for t in ds["Times"].values:

            if isinstance(t, (bytes, np.bytes_)):
                texto = t.decode("utf-8")
            elif isinstance(t, str):
                texto = t
            else:
                texto = "".join(
                    v.decode("utf-8")
                    if isinstance(v, (bytes, np.bytes_))
                    else str(v)
                    for v in t
                )

            texto = texto.replace("_", " ").strip("\x00 ")
            fechas.append(pd.to_datetime(texto))

        return pd.DatetimeIndex(fechas)

    for nombre in ["time", "Time", "XTIME"]:

        if nombre in ds.variables:

            valores = ds[nombre].values

            if np.issubdtype(valores.dtype, np.datetime64):
                return pd.DatetimeIndex(valores)

    raise ValueError("No se pudieron identificar las fechas.")

def leer_coordenadas(ds, nombre_lat, nombre_lon):

    lat = ds[nombre_lat]
    lon = ds[nombre_lon]

    for dim in ["Time", "time"]:
        if dim in lat.dims:
            lat = lat.isel({dim: 0})
        if dim in lon.dims:
            lon = lon.isel({dim: 0})

    lat, lon = xr.broadcast(lat, lon)
    lon = lon.transpose(*lat.dims)

    return lat.values, lon.values, lat.dims

def mapa_base(ax):

    ax.set_extent(extent, crs=proyeccion)

    ax.add_feature(
        cfeature.LAND.with_scale("10m"),
        facecolor="#f5f5f5",
        zorder=0
    )

    ax.add_feature(
        cfeature.OCEAN.with_scale("10m"),
        facecolor="#dcecf5",
        zorder=0
    )

    ax.add_feature(
        cfeature.LAKES.with_scale("10m"),
        facecolor="#dcecf5",
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

    ax.add_geometries(
        rs,
        crs=proyeccion,
        facecolor="none",
        edgecolor="#333333",
        linewidth=0.9,
        zorder=4
    )

# =====================================================
# 6. LEER Y ACUMULAR LAS CORRIDAS WRF
# =====================================================

campos_acumulados = {}

lat_ref = None
lon_ref = None

fechas_wrf_esperadas = pd.date_range(
    start=fecha_inicio,
    end=fecha_final,
    freq="h"
)

for nombre, archivo in simulaciones.items():

    print(f"Procesando {nombre}...")

    ruta = os.path.join(carpeta, archivo)

    with xr.open_dataset(ruta) as ds:

        fechas = leer_fechas(ds)

        lat, lon, dims_esp = leer_coordenadas(
            ds, "XLAT", "XLONG"
        )

        if lat_ref is None:
            lat_ref = lat.copy()
            lon_ref = lon.copy()
        else:
            if (
                lat.shape != lat_ref.shape or
                not np.allclose(lat, lat_ref, atol=1e-4) or
                not np.allclose(lon, lon_ref, atol=1e-4)
            ):
                raise ValueError(
                    f"La grilla de {nombre} no coincide."
                )

        pp_acum = ds["RAINC"] + ds["RAINNC"]

        dim_t = next(
            d for d in pp_acum.dims
            if d not in dims_esp
        )

        idx = np.flatnonzero(
            (fechas >= fecha_inicio) &
            (fechas <= fecha_final)
        )

        pp_acum = pp_acum.isel({dim_t: idx})
        fechas = fechas[idx]

        pp_acum = pp_acum.transpose(
            dim_t,
            *dims_esp
        ).values.astype(float)

    idx_orden = np.argsort(fechas)
    fechas = fechas[idx_orden]
    pp_acum = pp_acum[idx_orden]

    if not fechas.equals(fechas_wrf_esperadas):
        raise ValueError(
            f"{nombre} no tiene la secuencia horaria esperada."
        )

    acumulado_evento = pp_acum[-1] - pp_acum[0]

    acumulado_evento = np.where(
        acumulado_evento >= 0,
        acumulado_evento,
        np.nan
    )

    campos_acumulados[nombre] = acumulado_evento

    print(
        f"Rango {nombre}: "
        f"{np.nanmin(acumulado_evento):.2f} - "
        f"{np.nanmax(acumulado_evento):.2f} mm"
    )

# =====================================================
# 7. LEER Y ACUMULAR PERSIANN
# =====================================================

print("Procesando PERSIANN...")

with xr.open_dataset(
    os.path.join(carpeta, archivo_persiann)
) as ds:

    fechas = leer_fechas(ds)

    if not persiann_fin_hora:
        fechas = fechas + pd.Timedelta(hours=1)

    lat, lon, dims_esp = leer_coordenadas(
        ds, "XLAT", "XLONG"
    )

    if (
        lat.shape != lat_ref.shape or
        not np.allclose(lat, lat_ref, atol=1e-4) or
        not np.allclose(lon, lon_ref, atol=1e-4)
    ):
        raise ValueError("La grilla de PERSIANN no coincide.")

    pp = ds["precipitation"]

    dim_t = next(
        d for d in pp.dims
        if d not in dims_esp
    )

    idx = np.flatnonzero(
        (fechas >= fecha_inicio + pd.Timedelta(hours=1)) &
        (fechas <= fecha_final)
    )

    pp = pp.isel({dim_t: idx})
    fechas = fechas[idx]

    pp = pp.transpose(
        dim_t,
        *dims_esp
    ).values.astype(float)

idx_orden = np.argsort(fechas)
fechas = fechas[idx_orden]
pp = pp[idx_orden]

fechas_horarias = pd.date_range(
    start=fecha_inicio + pd.Timedelta(hours=1),
    end=fecha_final,
    freq="h"
)

if not pd.DatetimeIndex(fechas).equals(fechas_horarias):
    raise ValueError("PERSIANN no tiene la secuencia esperada.")

acumulado_evento = np.nansum(pp, axis=0)

acumulado_evento = np.where(
    acumulado_evento >= 0,
    acumulado_evento,
    np.nan
)

campos_acumulados["PERSIANN"] = acumulado_evento

print(
    f"Rango PERSIANN: "
    f"{np.nanmin(acumulado_evento):.2f} - "
    f"{np.nanmax(acumulado_evento):.2f} mm"
)

# =====================================================
# 8. LEER Y ACUMULAR MERGE HORARIO
# =====================================================

print("Procesando MERGE...")

with xr.open_dataset(
    os.path.join(carpeta, archivo_merge)
) as ds:

    fechas = leer_fechas(ds)

    if not merge_fin_hora:
        fechas = fechas + pd.Timedelta(hours=1)

    lat, lon, dims_esp = leer_coordenadas(
        ds, "latitude", "longitude"
    )

    if (
        lat.shape != lat_ref.shape or
        not np.allclose(lat, lat_ref, atol=1e-4) or
        not np.allclose(lon, lon_ref, atol=1e-4)
    ):
        raise ValueError("La grilla de MERGE no coincide.")

    pp = ds["precipitation"]

    dim_t = next(
        d for d in pp.dims
        if d not in dims_esp
    )

    idx = np.flatnonzero(
        (fechas >= fecha_inicio + pd.Timedelta(hours=1)) &
        (fechas <= fecha_final)
    )

    pp = pp.isel({dim_t: idx})
    fechas = fechas[idx]

    pp = pp.transpose(
        dim_t,
        *dims_esp
    ).values.astype(float)

idx_orden = np.argsort(fechas)
fechas = fechas[idx_orden]
pp = pp[idx_orden]

if not pd.DatetimeIndex(fechas).equals(fechas_horarias):
    raise ValueError("MERGE no tiene la secuencia esperada.")

acumulado_evento = np.nansum(pp, axis=0)

acumulado_evento = np.where(
    acumulado_evento >= 0,
    acumulado_evento,
    np.nan
)

campos_acumulados["MERGE"] = acumulado_evento

print(
    f"Rango MERGE: "
    f"{np.nanmin(acumulado_evento):.2f} - "
    f"{np.nanmax(acumulado_evento):.2f} mm"
)

# =====================================================
# 9. GENERAR FIGURA
# =====================================================

fig, axs = plt.subplots(
    3, 2,
    figsize=(13, 15),
    subplot_kw={"projection": proyeccion}
)

axs = axs.flatten()

for i, nombre in enumerate(orden):

    ax = axs[i]

    mapa_base(ax)

    campo = campos_acumulados[nombre]

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
        pad=9
    )

fig.subplots_adjust(
    left=0.04,
    right=0.86,
    bottom=0.04,
    top=0.96,
    wspace=0.04,
    hspace=0.10
)

# =====================================================
# 10. COLORBAR COMPARTIDA
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
    extend="max",
    format=mticker.FormatStrFormatter("%.0f")
)

cbar.set_label(
    "Accumulated precipitation (mm)",
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
# 11. GUARDAR FIGURA
# =====================================================

salida = os.path.join(
    args.figures_dir,
    "accumulated_precipitation.png"
)

fig.savefig(
    salida,
    dpi=600,
    format="png",
    bbox_inches="tight",
    facecolor="white"
)

print(f"\nFigura guardada en: {salida}")

plt.show()
plt.close(fig)

# =====================================================
# 12. GUARDAR ACUMULADOS EN NETCDF/NPY OPCIONAL
# =====================================================

for nombre in orden:

    np.save(
        os.path.join(
            carpeta_salida,
            f"Acumulado_{nombre}.npy"
        ),
        campos_acumulados[nombre]
    )

print("\nProcesamiento finalizado.")