# -*- coding: utf-8 -*-
"""Evaluate hourly precipitation evolution and spatial correlation over Rio Grande do Sul.

Author: Fernando Primo Forgioni (ORCID: 0000-0003-3233-8957).
Inputs: data/Estaciones.xlsx, four WRF NetCDF files, MERGE_HORARIO_WRF.nc.
Period: 28 April 2024 12 UTC to 3 May 2024 12 UTC.
Outputs: figures/hourly_precipitation.png and
results/Diagnosticos_horarios_WRF_MERGE_Gauge_CC.xlsx.
Run: python scripts/07_hourly_analysis.py
Original calculations and plotting settings are retained.
"""

import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns

from scipy.spatial import cKDTree
from cartopy.io import shapereader
from shapely.geometry import Point
from shapely.prepared import prep

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
args.figures_dir.mkdir(parents=True, exist_ok=True)
args.results_dir.mkdir(parents=True, exist_ok=True)

archivo_estaciones = os.path.join(
    carpeta,
    "Estaciones.xlsx"
)

simulaciones = {
    "GF-ERA5": "Grell-Freitas-ERA5.nc",
    "GF-GFS": "Grell-Freitas-GFS.nc",
    "KF-ERA5": "Kain-Fritsch-ERA5.nc",
    "KF-GFS": "Kain-Fritsch-GFS.nc"
}

archivo_merge = "MERGE_HORARIO_WRF.nc"

# Check required input files before processing.
for input_path in [Path(archivo_estaciones)] + [carpeta / name for name in (*simulaciones.values(), archivo_merge)]:
    if not input_path.is_file():
        raise FileNotFoundError(f"Required input file not found: {input_path}")

# Período del evento

fecha_inicio = pd.Timestamp("2024-04-28 12:00:00")
fecha_final = pd.Timestamp("2024-05-03 12:00:00")

sheet_estaciones = "Sheet1"

# Convenciones temporales
# True: fecha correspondiente al final del acumulado
# False: fecha correspondiente al inicio del acumulado

merge_fin_hora = True
gauge_fin_dia = False

# =====================================================
# 2. COLORES Y ORDEN
# =====================================================

colores = {
    "Gauge": "#708090",
    "GF-ERA5": "#87CEFA",
    "GF-GFS": "#4682B4",
    "KF-ERA5": "#BC8F8F",
    "KF-GFS": "#F5DEB3",
    "PERSIANN": "#5F9EA0",
    "MERGE": "#556B2F"
}

orden = [
    "GF-ERA5",
    "GF-GFS",
    "KF-ERA5",
    "KF-GFS",
    "MERGE"
]

# =====================================================
# 3. FUNCIONES AUXILIARES
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

            texto = texto.replace(
                "_", " "
            ).strip("\x00 ")

            fechas.append(pd.to_datetime(texto))

        return pd.DatetimeIndex(fechas)

    for nombre in ["time", "Time", "XTIME"]:

        if nombre in ds.variables:

            valores = ds[nombre].values

            if np.issubdtype(
                valores.dtype,
                np.datetime64
            ):

                return pd.DatetimeIndex(valores)

    raise ValueError(
        "No se pudieron identificar las fechas."
    )


# -----------------------------------------------------
# Coordenadas de los NetCDF
# -----------------------------------------------------

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


# -----------------------------------------------------
# Coordenadas esféricas
# -----------------------------------------------------

def xyz(lat, lon):

    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    return np.column_stack([
        (np.cos(lat) * np.cos(lon)).ravel(),
        (np.cos(lat) * np.sin(lon)).ravel(),
        np.sin(lat).ravel()
    ])


# -----------------------------------------------------
# Buscar punto de grilla más cercano
# -----------------------------------------------------

def crear_indices_estaciones(
    lat_grid,
    lon_grid,
    lat_est,
    lon_est
):

    xyz_grid = xyz(lat_grid, lon_grid)
    xyz_sta = xyz(lat_est, lon_est)

    validos = np.isfinite(
        xyz_grid
    ).all(axis=1)

    if not validos.any():

        raise ValueError(
            "No hay coordenadas válidas en la grilla."
        )

    arbol = cKDTree(
        xyz_grid[validos]
    )

    distancias, indices = arbol.query(
        xyz_sta,
        k=1
    )

    indices = np.flatnonzero(
        validos
    )[indices]

    ij = np.unravel_index(
        indices,
        lat_grid.shape
    )

    distancia_km = 2 * 6371 * np.arcsin(
        np.clip(distancias / 2, 0, 1)
    )

    print(
        f"Distancia media estación-grilla: "
        f"{distancia_km.mean():.2f} km"
    )

    return ij


# -----------------------------------------------------
# Extraer precipitación en estaciones
# -----------------------------------------------------

def extraer_estaciones_desde_campo(campo, ij):

    iy, ix = ij

    return campo[:, iy, ix]


# -----------------------------------------------------
# Verificar coincidencia de las grillas
# -----------------------------------------------------

def verificar_grilla(lat, lon, lat_ref, lon_ref, nombre):

    if (
        lat.shape != lat_ref.shape or
        lon.shape != lon_ref.shape or
        not np.allclose(
            lat, lat_ref,
            atol=1e-4,
            rtol=0
        ) or
        not np.allclose(
            lon, lon_ref,
            atol=1e-4,
            rtol=0
        )
    ):

        raise ValueError(
            f"La grilla de {nombre} no coincide con WRF."
        )


# =====================================================
# 4. LEER OBSERVACIONES DE ESTACIONES
# =====================================================

df_est = pd.read_excel(
    archivo_estaciones,
    sheet_name=sheet_estaciones
)

df_est.columns = [
    c.strip().lower()
    if isinstance(c, str)
    else c
    for c in df_est.columns
]

df_est["latitud"] = pd.to_numeric(
    df_est["latitud"],
    errors="coerce"
)

df_est["longitud"] = pd.to_numeric(
    df_est["longitud"],
    errors="coerce"
)

df_est = df_est.dropna(
    subset=["latitud", "longitud"]
).reset_index(drop=True)

df_est["id_estacion"] = np.arange(
    1,
    len(df_est) + 1
)

# Identificar columnas con fechas

columnas_fechas = []
fechas_obs = []

for col in df_est.columns:

    fecha = pd.to_datetime(
        str(col),
        errors="coerce"
    )

    if not pd.isna(fecha):

        columnas_fechas.append(col)
        fechas_obs.append(fecha.normalize())

if not columnas_fechas:

    raise ValueError(
        "No se encontraron fechas en el Excel."
    )

# Construir tabla de observaciones

obs = df_est[columnas_fechas].apply(
    pd.to_numeric,
    errors="coerce"
).T

obs.index = pd.DatetimeIndex(fechas_obs)
obs.columns = df_est["id_estacion"].values

obs = obs.sort_index()

# Los ceros son datos válidos

obs = obs.where(
    np.isfinite(obs) & (obs >= 0)
)

if obs.index.duplicated().any():

    raise ValueError(
        "Hay fechas duplicadas en las observaciones."
    )

# Unificar etiquetas al inicio del acumulado diario

if gauge_fin_dia:

    obs.index = obs.index - pd.Timedelta(days=1)

print(f"\nEstaciones disponibles: {len(df_est)}")

# =====================================================
# 5. GEOMETRÍA DE RIO GRANDE DO SUL
# =====================================================

shp = shapereader.natural_earth(
    resolution="10m",
    category="cultural",
    name="admin_1_states_provinces"
)

geometria_rs = None

for registro in shapereader.Reader(shp).records():

    if registro.attributes.get("name") == "Rio Grande do Sul":

        geometria_rs = registro.geometry

        break

if geometria_rs is None:

    raise ValueError(
        "No se encontró Rio Grande do Sul."
    )

geometria_preparada = prep(
    geometria_rs
)


def crear_mascara_rs(lat, lon):

    mascara = (
        np.isfinite(lat) &
        np.isfinite(lon)
    )

    mascara_rs = np.zeros(
        lat.shape,
        dtype=bool
    )

    for i, j in np.ndindex(lat.shape):

        if mascara[i, j]:

            punto = Point(
                float(lon[i, j]),
                float(lat[i, j])
            )

            mascara_rs[i, j] = (
                geometria_preparada.covers(punto)
            )

    return mascara_rs


# =====================================================
# 6. PERÍODO HORARIO COMÚN
# =====================================================

# WRF: fechas de la precipitación acumulada

fechas_wrf_esperadas = pd.date_range(
    start=fecha_inicio,
    end=fecha_final,
    freq="h"
)

# Fechas correspondientes al final de cada hora

fechas_horarias = pd.date_range(
    start=fecha_inicio + pd.Timedelta(hours=1),
    end=fecha_final,
    freq="h"
)

# Finalización de los acumulados diarios

fechas_gauge_plot = pd.date_range(
    start=fecha_inicio + pd.Timedelta(days=1),
    end=fecha_final,
    freq="24h"
)

# Fechas observadas etiquetadas al inicio del día

fechas_gauge = (
    fechas_gauge_plot - pd.Timedelta(days=1)
).normalize()

print("\nPeríodo de análisis")
print(f"Inicio: {fecha_inicio}")
print(f"Final: {fecha_final}")
print(f"Horas: {len(fechas_horarias)}")
print(f"Días: {len(fechas_gauge)}")

# =====================================================
# 7. PROCESAR SIMULACIONES WRF
# =====================================================

campos_horarios_rs = {}
campos_horarios_est = {}

lat_ref = None
lon_ref = None
mascara_rs = None
indices_estaciones = None

for nombre, archivo in simulaciones.items():

    print("\n" + "=" * 60)
    print(f"Procesando {nombre}")
    print("=" * 60)

    ruta = os.path.join(
        carpeta,
        archivo
    )

    with xr.open_dataset(ruta) as ds:

        fechas = leer_fechas(ds)

        lat, lon, dims_esp = leer_coordenadas(
            ds,
            "XLAT",
            "XLONG"
        )

        if lat_ref is None:

            lat_ref = lat.copy()
            lon_ref = lon.copy()

            mascara_rs = crear_mascara_rs(
                lat_ref,
                lon_ref
            )

            indices_estaciones = crear_indices_estaciones(
                lat_ref,
                lon_ref,
                df_est["latitud"].values,
                df_est["longitud"].values
            )

        else:

            verificar_grilla(
                lat,
                lon,
                lat_ref,
                lon_ref,
                nombre
            )

        # Precipitación convectiva y no convectiva

        pp_acum = ds["RAINC"] + ds["RAINNC"]

        dim_t = next(
            d for d in pp_acum.dims
            if d not in dims_esp
        )

        # Seleccionar período de análisis

        idx = np.flatnonzero(
            (fechas >= fecha_inicio) &
            (fechas <= fecha_final)
        )

        pp_acum = pp_acum.isel({
            dim_t: idx
        })

        fechas = fechas[idx]

        pp_acum = pp_acum.transpose(
            dim_t,
            *dims_esp
        ).values.astype(float)

    # Orden cronológico

    idx_orden = np.argsort(fechas)

    fechas = fechas[idx_orden]
    pp_acum = pp_acum[idx_orden]

    if not fechas.equals(fechas_wrf_esperadas):

        raise ValueError(
            f"{nombre} no tiene la secuencia "
            "horaria esperada."
        )

    # -------------------------------------------------
    # Calcular precipitación horaria
    # -------------------------------------------------

    pp_h = np.diff(
        pp_acum,
        axis=0
    )

    # Detectar reinicios del acumulador

    if np.any(pp_h < -0.001):

        raise ValueError(
            f"Disminuciones del acumulador en {nombre}."
        )

    # Corregir pequeños negativos numéricos

    pp_h = np.where(
        (pp_h < 0) & (pp_h >= -0.001),
        0,
        pp_h
    )

    # Guardar campos horarios

    campos_horarios_rs[nombre] = pp_h

    # Extraer precipitación en las estaciones

    campos_horarios_est[nombre] = (
        extraer_estaciones_desde_campo(
            pp_h,
            indices_estaciones
        )
    )

    print(
        f"Horas disponibles: {pp_h.shape[0]}"
    )

    del pp_acum

# =====================================================
# 8. PROCESAR MERGE HORARIO
# =====================================================

print("\n" + "=" * 60)
print("Procesando MERGE")
print("=" * 60)

with xr.open_dataset(
    os.path.join(carpeta, archivo_merge)
) as ds:

    fechas = leer_fechas(ds)

    lat, lon, dims_esp = leer_coordenadas(
        ds,
        "latitude",
        "longitude"
    )

    verificar_grilla(
        lat,
        lon,
        lat_ref,
        lon_ref,
        "MERGE"
    )

    pp = ds["precipitation"]

    dim_t = next(
        d for d in pp.dims
        if d not in dims_esp
    )

    # Unificar las fechas al final de cada hora

    if not merge_fin_hora:

        fechas = fechas + pd.Timedelta(hours=1)

    if fechas.duplicated().any():

        raise ValueError(
            "Hay fechas duplicadas en MERGE."
        )

    idx = fechas.get_indexer(
        fechas_horarias
    )

    if np.any(idx < 0):

        raise ValueError(
            "Faltan registros horarios en MERGE."
        )

    pp = pp.isel({
        dim_t: idx
    })

    pp = pp.transpose(
        dim_t,
        *dims_esp
    ).values.astype(float)

pp = np.where(
    np.isfinite(pp) & (pp >= 0),
    pp,
    np.nan
)

campos_horarios_rs["MERGE"] = pp

campos_horarios_est["MERGE"] = (
    extraer_estaciones_desde_campo(
        pp,
        indices_estaciones
    )
)

print(
    f"Horas disponibles: {pp.shape[0]}"
)

# =====================================================
# 9. MÁSCARAS COMUNES
# =====================================================

# -----------------------------------------------------
# Máscara espacial sobre Rio Grande do Sul
# -----------------------------------------------------

mascara_spatial = mascara_rs.copy()

for nombre in orden:

    campo = campos_horarios_rs[nombre]

    mascara_spatial &= np.all(
        np.isfinite(campo) &
        (campo >= 0),
        axis=0
    )

n_grid = np.sum(mascara_spatial)

print(
    f"\nPuntos de grilla utilizados sobre RS: {n_grid}"
)

if n_grid == 0:

    raise ValueError(
        "No hay puntos de grilla comunes."
    )

# -----------------------------------------------------
# Máscara común de estaciones
# -----------------------------------------------------

obs_evento = obs.reindex(
    fechas_gauge
)

estaciones_validas = (
    obs_evento.notna().all(axis=0).to_numpy(copy=True)
)

# Exigir datos completos durante las 120 horas

for nombre in orden:

    valores = campos_horarios_est[nombre]

    estaciones_validas &= np.all(
        np.isfinite(valores) &
        (valores >= 0),
        axis=0
    )

ids_validos = obs.columns[
    estaciones_validas
]

n_est = len(ids_validos)

print("\n" + "=" * 60)
print(f"Estaciones originales: {len(df_est)}")
print(f"Estaciones utilizadas: {n_est}")
print(f"Estaciones excluidas: {len(df_est) - n_est}")
print("=" * 60)

if n_est == 0:

    raise ValueError(
        "No hay estaciones comunes para "
        "la comparación acumulada."
    )

# =====================================================
# 10. PRECIPITACIÓN MEDIA HORARIA SOBRE RS
# =====================================================

media_horaria = pd.DataFrame(
    index=fechas_horarias
)

for nombre in orden:

    valores = campos_horarios_rs[nombre][
        :, mascara_spatial
    ]

    media_horaria[nombre] = np.mean(
        valores,
        axis=1
    )

# =====================================================
# 11. PRECIPITACIÓN ACUMULADA EN ESTACIONES
# =====================================================

acumulada_est = pd.DataFrame(
    index=fechas_horarias
)

for nombre in orden:

    valores = campos_horarios_est[nombre][
        :, estaciones_validas
    ]

    # Precipitación horaria media

    media_est = np.mean(
        valores,
        axis=1
    )

    # Precipitación media acumulada

    acumulada_est[nombre] = np.cumsum(
        media_est
    )

# -----------------------------------------------------
# Gauge: precipitación media diaria
# -----------------------------------------------------

gauge_media_diaria = obs_evento[
    ids_validos
].mean(axis=1)

gauge_acumulado = gauge_media_diaria.cumsum()

# Asignar fechas correspondientes al final
# de los intervalos diarios

gauge_media_diaria.index = fechas_gauge_plot
gauge_acumulado.index = fechas_gauge_plot

print("\nPrecipitación acumulada final:")

print(
    acumulada_est.iloc[-1].round(2)
)

print(
    f"Gauge: {gauge_acumulado.iloc[-1]:.2f} mm"
)

# =====================================================
# 12. CORRELACIÓN ESPACIAL HORARIA RESPECTO A MERGE
# =====================================================

correlacion_espacial = pd.DataFrame(
    index=fechas_horarias
)

referencia = campos_horarios_rs["MERGE"][
    :, mascara_spatial
]

for nombre in simulaciones:

    valores = campos_horarios_rs[nombre][
        :, mascara_spatial
    ]

    correlaciones = []

    for i in range(len(fechas_horarias)):

        o = referencia[i, :]
        s = valores[i, :]

        # Seleccionar pares válidos
        # Los ceros se conservan

        validos = (
            np.isfinite(o) &
            np.isfinite(s)
        )

        o = o[validos]
        s = s[validos]

        # Correlación espacial de Pearson

        if (
            len(o) >= 3 and
            np.std(o) > 0 and
            np.std(s) > 0
        ):

            cc = np.corrcoef(
                o,
                s
            )[0, 1]

        else:

            cc = np.nan

        correlaciones.append(cc)

    correlacion_espacial[nombre] = correlaciones

print("\nCorrelación espacial media respecto a MERGE:")

print(
    correlacion_espacial.mean().round(3)
)

# =====================================================
# 13. CONFIGURACIÓN GRÁFICA
# =====================================================

sns.set_theme(
    style="ticks",
    context="paper",
    font="DejaVu Sans"
)

plt.rcParams.update({
    "font.size": 13,
    "axes.linewidth": 0.9
})

fontsize_tit = 15
fontsize_ejes = 14
fontsize_ticks = 13
fontsize_leyenda = 13

grosor_wrf = 1.9
grosor_merge = 2.4

tam_gauge = 65

# =====================================================
# 14. FIGURA DE TRES PANELES
# =====================================================

fig, axs = plt.subplots(
    3, 1,
    figsize=(13, 11),
    sharex=True
)

configuracion = [
    {
        "datos": media_horaria,
        "titulo": "(a) Mean hourly precipitation",
        "ylabel": r"Precipitation (mm h$^{-1}$)"
    },
    {
        "datos": acumulada_est,
        "titulo": "(b) Accumulated precipitation at station locations",
        "ylabel": "Precipitation (mm)"
    },
    {
        "datos": correlacion_espacial,
        "titulo": "(c) Spatial correlation against MERGE",
        "ylabel": "Correlation coefficient"
    }
]

for i, config in enumerate(configuracion):

    ax = axs[i]
    datos = config["datos"]

    # -------------------------------------------------
    # Seleccionar conjuntos de datos
    # -------------------------------------------------

    if i == 2:

        nombres = list(
            simulaciones.keys()
        )

    else:

        nombres = orden

    # -------------------------------------------------
    # Graficar las series
    # -------------------------------------------------

    for nombre in nombres:

        x = datos.index
        y = datos[nombre].values

        # Iniciar la precipitación acumulada en cero

        if i == 1:

            x = x.insert(
                0,
                fecha_inicio
            )

            y = np.insert(
                y,
                0,
                0
            )

        # Estilos de línea

        estilo = (
            "--" if nombre == "MERGE" else "-"
        )

        grosor = (
            grosor_merge
            if nombre == "MERGE"
            else grosor_wrf
        )

        ax.plot(
            x,
            y,
            color=colores[nombre],
            linestyle=estilo,
            linewidth=grosor,
            label=nombre,
            zorder=3
        )

    # -------------------------------------------------
    # Observaciones: Gauge
    # -------------------------------------------------

    if i == 1:

        ax.scatter(
            fechas_gauge_plot,
            gauge_acumulado.values,
            s=tam_gauge,
            marker="o",
            facecolor=colores["Gauge"],
            edgecolor="white",
            linewidth=0.8,
            zorder=10,
            label="Gauge"
        )

    # -------------------------------------------------
    # Títulos
    # -------------------------------------------------

    ax.set_title(
        config["titulo"],
        fontsize=fontsize_tit,
        fontweight="bold",
        loc="left",
        pad=12
    )

    ax.set_ylabel(
        config["ylabel"],
        fontsize=fontsize_ejes,
        labelpad=10
    )

    # -------------------------------------------------
    # Límites horizontales
    # -------------------------------------------------

    ax.set_xlim(
        fecha_inicio,
        fecha_final
    )

    # -------------------------------------------------
    # Límites verticales
    # -------------------------------------------------

    if i == 2:

        ax.set_ylim(-1, 1)

        ax.set_yticks(
            np.linspace(-1, 1, 11)
        )

        # Línea de correlación nula

        ax.axhline(
            0,
            color="#555555",
            linestyle="--",
            linewidth=1.0,
            zorder=1
        )

    else:

        ax.set_ylim(bottom=0)

    # -------------------------------------------------
    # Configuración temporal
    # -------------------------------------------------

    ax.xaxis.set_major_locator(
        mdates.HourLocator(byhour=12)
    )

    ax.xaxis.set_major_formatter(
        mdates.DateFormatter("%d %b")
    )

    ax.tick_params(
        axis="both",
        labelsize=fontsize_ticks,
        direction="out",
        length=4,
        width=0.8
    )

    # -------------------------------------------------
    # Grilla horizontal
    # -------------------------------------------------

    ax.grid(
        axis="y",
        linestyle="--",
        linewidth=0.6,
        alpha=0.30
    )

    ax.grid(
        axis="x",
        visible=False
    )

    ax.set_axisbelow(True)

    sns.despine(ax=ax)

# Etiqueta temporal únicamente en el último panel

axs[2].set_xlabel(
    "Date (UTC)",
    fontsize=fontsize_ejes,
    labelpad=10
)

# =====================================================
# 15. LEYENDA COMPARTIDA
# =====================================================

handles, labels = axs[1].get_legend_handles_labels()

fig.legend(
    handles,
    labels,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.02),
    ncol=6,
    fontsize=fontsize_leyenda,
    frameon=False,
    handlelength=3
)

fig.subplots_adjust(
    left=0.10,
    right=0.97,
    bottom=0.12,
    top=0.95,
    hspace=0.25
)

# =====================================================
# 16. GUARDAR FIGURA
# =====================================================

salida = os.path.join(
    args.figures_dir,
    "hourly_precipitation.png"
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
# 17. GUARDAR DIAGNÓSTICOS EN EXCEL
# =====================================================

archivo_excel = os.path.join(
    args.results_dir,
    "Diagnosticos_horarios_WRF_MERGE_Gauge_CC.xlsx"
)

with pd.ExcelWriter(
    archivo_excel,
    engine="openpyxl"
) as writer:

    media_horaria.to_excel(
        writer,
        sheet_name="Mean_hourly_RS"
    )

    acumulada_est.to_excel(
        writer,
        sheet_name="Accumulated_station_locs"
    )

    correlacion_espacial.to_excel(
        writer,
        sheet_name="Spatial_correlation"
    )

    pd.DataFrame({
        "Gauge_daily_mean": gauge_media_diaria,
        "Gauge_accumulated": gauge_acumulado
    }).to_excel(
        writer,
        sheet_name="Gauge_points"
    )

    correlacion_espacial.describe().to_excel(
        writer,
        sheet_name="Spatial_CC_summary"
    )

    pd.DataFrame({
        "Parameter": [
            "Start UTC",
            "End UTC",
            "Number of hours",
            "Number of stations",
            "Number of grid points",
            "MERGE end-hour timestamp",
            "Gauge end-day timestamp"
        ],
        "Value": [
            str(fecha_inicio),
            str(fecha_final),
            len(fechas_horarias),
            n_est,
            n_grid,
            merge_fin_hora,
            gauge_fin_dia
        ]
    }).to_excel(
        writer,
        sheet_name="Configuration",
        index=False
    )

print("\nProcesamiento finalizado.")