# -*- coding: utf-8 -*-
"""Compare event-total and daily precipitation estimates against station observations.

Author: Fernando Primo Forgioni (ORCID: 0000-0003-3233-8957).
Inputs: data/Estaciones.xlsx, four WRF NetCDF files, PERSIANN_WRF.nc,
and Merge-diario.nc. Daily intervals start at 12 UTC.
Outputs: figures/boxplot_precipitation_event.png and
figures/density_scatter_daily_precipitation.png.
Run: python scripts/06_statistical_comparison.py
The original common-date and complete-station selection is retained.
"""

import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree

# =====================================================
# 1. CONFIGURACIÓN
# =====================================================

repository_root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--data-dir", type=Path, default=repository_root / "data")
parser.add_argument("--figures-dir", type=Path, default=repository_root / "figures")
args = parser.parse_args()
carpeta = args.data_dir

archivo_estaciones = os.path.join(
    carpeta, "Estaciones.xlsx"
)

carpeta_salida = args.figures_dir

os.makedirs(carpeta_salida, exist_ok=True)

simulaciones = {
    "GF-ERA5": "Grell-Freitas-ERA5.nc",
    "GF-GFS": "Grell-Freitas-GFS.nc",
    "KF-ERA5": "Kain-Fritsch-ERA5.nc",
    "KF-GFS": "Kain-Fritsch-GFS.nc"
}

productos = {
    "PERSIANN": "PERSIANN_WRF.nc",
    "MERGE": "Merge-diario.nc"
}

# Check required input files before processing.
for input_path in [Path(archivo_estaciones)] + [carpeta / name for name in (*simulaciones.values(), *productos.values())]:
    if not input_path.is_file():
        raise FileNotFoundError(f"Required input file not found: {input_path}")

# Hora UTC de inicio del acumulado diario

hora_inicio = 12

# Convención temporal de los productos

persiann_fin_hora = False
merge_fin_dia = True

# =====================================================
# 2. LEER ESTACIONES METEOROLÓGICAS
# =====================================================

df = pd.read_excel(
    archivo_estaciones,
    sheet_name="Sheet1"
)

df.columns = [
    c.strip().lower() if isinstance(c, str) else c
    for c in df.columns
]

df["latitud"] = pd.to_numeric(
    df["latitud"], errors="coerce"
)

df["longitud"] = pd.to_numeric(
    df["longitud"], errors="coerce"
)

df = df.dropna(
    subset=["latitud", "longitud"]
).reset_index(drop=True)

df["id_estacion"] = np.arange(1, len(df) + 1)

# Identificar columnas de fechas

columnas_fechas = []
fechas = []

for columna in df.columns:

    fecha = pd.to_datetime(
        str(columna),
        errors="coerce"
    )

    if not pd.isna(fecha):

        columnas_fechas.append(columna)
        fechas.append(fecha.normalize())

if not columnas_fechas:
    raise ValueError(
        "No se encontraron fechas en el Excel."
    )

# Construir tabla de observaciones

obs = df[columnas_fechas].apply(
    pd.to_numeric,
    errors="coerce"
).T

obs.index = pd.DatetimeIndex(fechas)
obs.columns = df["id_estacion"].values

obs = obs.sort_index()

obs = obs.where(obs >= 0)

if obs.index.duplicated().any():
    raise ValueError("Hay fechas duplicadas en el Excel.")

print(f"Estaciones: {len(df)}")
print(f"Días del evento: {len(obs)}")
print(f"Fecha inicial: {obs.index.min()}")
print(f"Fecha final: {obs.index.max()}")

# =====================================================
# 3. FUNCIONES AUXILIARES
# =====================================================

def xyz(lat, lon):

    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    return np.column_stack([
        (np.cos(lat) * np.cos(lon)).ravel(),
        (np.cos(lat) * np.sin(lon)).ravel(),
        np.sin(lat).ravel()
    ])


xyz_estaciones = xyz(
    df["latitud"].values,
    df["longitud"].values
)


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

            if np.issubdtype(
                valores.dtype, np.datetime64
            ):

                return pd.DatetimeIndex(valores)

    raise ValueError(
        "No se pudieron identificar las fechas."
    )


# =====================================================
# 4. EXTRAER DATOS EN LAS ESTACIONES
# =====================================================

def extraer_estaciones(ds, variable, nombre_lat, nombre_lon):

    lat = ds[nombre_lat]
    lon = ds[nombre_lon]

    # Eliminar dimensión temporal de las coordenadas

    for dim in ["Time", "time"]:

        if dim in lat.dims:
            lat = lat.isel({dim: 0})

        if dim in lon.dims:
            lon = lon.isel({dim: 0})

    # Coordenadas unidimensionales o bidimensionales

    lat, lon = xr.broadcast(lat, lon)

    lon = lon.transpose(*lat.dims)

    dimensiones = lat.dims

    lat_grid = lat.values
    lon_grid = lon.values

    # Buscar punto más cercano

    coordenadas = xyz(
        lat_grid,
        lon_grid
    )

    validos = np.isfinite(coordenadas).all(axis=1)

    if not validos.any():
        raise ValueError("No hay coordenadas válidas.")

    arbol = cKDTree(coordenadas[validos])

    distancias, indices = arbol.query(
        xyz_estaciones,
        k=1
    )

    indices = np.flatnonzero(validos)[indices]

    indices_grid = np.unravel_index(
        indices,
        lat_grid.shape
    )

    seleccion = {
        dim: xr.DataArray(indice, dims="estacion")
        for dim, indice in zip(
            dimensiones,
            indices_grid
        )
    }

    # Extraer precipitación

    pp = variable.isel(seleccion)

    dimensiones_t = [
        d for d in pp.dims if d != "estacion"
    ]

    if len(dimensiones_t) != 1:
        raise ValueError(
            f"Dimensiones inesperadas: {pp.dims}"
        )

    dim_t = dimensiones_t[0]

    valores = pp.transpose(
        dim_t,
        "estacion"
    ).values

    fechas = leer_fechas(ds)

    tabla = pd.DataFrame(
        valores,
        index=fechas,
        columns=df["id_estacion"].values
    )

    tabla = tabla.sort_index()

    if tabla.index.duplicated().any():
        raise ValueError(
            "Hay fechas duplicadas en el NetCDF."
        )

    distancia_km = 2 * 6371 * np.arcsin(
        np.clip(distancias / 2, 0, 1)
    )

    print(
        f"Distancia media estación-grilla: "
        f"{distancia_km.mean():.2f} km"
    )

    return tabla


# =====================================================
# 5. PROCESAR SIMULACIONES WRF
# =====================================================

series_diarias = {}

for nombre, archivo in simulaciones.items():

    print(f"\nProcesando {nombre}...")

    ruta = os.path.join(carpeta, archivo)

    with xr.open_dataset(ruta) as ds:

        tabla = extraer_estaciones(
            ds,
            ds["RAINC"] + ds["RAINNC"],
            "XLAT",
            "XLONG"
        )

    # Detectar reinicios del acumulador

    if np.any(tabla.diff().values < -0.001):

        raise ValueError(
            f"Disminuciones del acumulador en {nombre}."
        )

    # Inicio y final de cada día

    inicio = obs.index + pd.Timedelta(
        hours=hora_inicio
    )

    final = inicio + pd.Timedelta(days=1)

    pp_inicio = tabla.reindex(inicio).values
    pp_final = tabla.reindex(final).values

    diferencia = pp_final - pp_inicio

    # Corregir pequeños errores numéricos

    diferencia = np.where(
        (diferencia < 0) &
        (diferencia >= -0.001),
        0,
        diferencia
    )

    pp_diaria = pd.DataFrame(
        diferencia,
        index=obs.index,
        columns=obs.columns
    )

    pp_diaria = pp_diaria.where(pp_diaria >= 0)

    series_diarias[nombre] = pp_diaria

    print(f"Datos válidos: {pp_diaria.count().sum()}")

# =====================================================
# 6. PROCESAR PERSIANN
# =====================================================

print("\nProcesando PERSIANN...")

with xr.open_dataset(
    os.path.join(carpeta, productos["PERSIANN"])
) as ds:

    tabla = extraer_estaciones(
        ds,
        ds["precipitation"],
        "XLAT",
        "XLONG"
    )

tabla = tabla.where(tabla >= 0)

# Si los registros representan la hora anterior

if persiann_fin_hora:

    tabla.index = tabla.index - pd.Timedelta(hours=1)

# Acumulados diarios de 24 horas

pp_diaria = []

for fecha in obs.index:

    inicio = fecha + pd.Timedelta(hours=hora_inicio)

    horas = pd.date_range(
        start=inicio,
        periods=24,
        freq="h"
    )

    acumulado = tabla.reindex(horas).sum(
        axis=0,
        min_count=24
    )

    pp_diaria.append(acumulado)

pp_diaria = pd.DataFrame(
    pp_diaria,
    index=obs.index,
    columns=obs.columns
)

series_diarias["PERSIANN"] = pp_diaria

print(f"Datos válidos: {pp_diaria.count().sum()}")

# =====================================================
# 7. PROCESAR MERGE
# =====================================================

print("\nProcesando MERGE...")

with xr.open_dataset(
    os.path.join(carpeta, productos["MERGE"])
) as ds:

    tabla = extraer_estaciones(
        ds,
        ds["precipitation"],
        "latitude",
        "longitude"
    )

tabla = tabla.where(tabla >= 0)

# Verificar horario de los registros diarios

if not np.all(
    tabla.index.hour == hora_inicio
):

    raise ValueError(
        "El horario de MERGE no coincide con hora_inicio."
    )

# MERGE: fecha de finalización del acumulado

if merge_fin_dia:

    tabla.index = tabla.index - pd.Timedelta(days=1)

tabla.index = tabla.index.normalize()

if tabla.index.duplicated().any():
    raise ValueError("Hay fechas diarias duplicadas en MERGE.")

pp_diaria = tabla.reindex(obs.index)

series_diarias["MERGE"] = pp_diaria

print(f"Datos válidos: {pp_diaria.count().sum()}")

# =====================================================
# 8. SELECCIONAR FECHAS Y ESTACIONES COMUNES
# =====================================================

# Primero seleccionar las fechas comunes

fechas_comunes = obs.index[obs.notna().any(axis=1)]

for nombre, tabla in series_diarias.items():

    fechas_validas = tabla.index[tabla.notna().any(axis=1)]

    fechas_comunes = fechas_comunes.intersection(
        fechas_validas
    )

if len(fechas_comunes) == 0:
    raise ValueError("No hay fechas comunes.")

obs = obs.loc[fechas_comunes]

for nombre in series_diarias:

    series_diarias[nombre] = series_diarias[nombre].loc[
        fechas_comunes
    ]

print("\nFechas comunes utilizadas:")
print(fechas_comunes)
print(f"Días analizados: {len(fechas_comunes)}")

# Ahora seleccionar las estaciones con datos completos

mascara_comun = obs.notna()

for tabla in series_diarias.values():

    mascara_comun &= tabla.notna()

estaciones_validas = mascara_comun.all(axis=0)

ids_validos = estaciones_validas[
    estaciones_validas
].index

n_estaciones = len(ids_validos)

print("\n" + "=" * 60)
print(f"Estaciones originales: {len(df)}")
print(f"Estaciones utilizadas: {n_estaciones}")
print(f"Estaciones excluidas: {len(df) - n_estaciones}")
print("=" * 60)

if n_estaciones == 0:

    raise ValueError(
        "No hay estaciones con datos completos "
        "en todos los conjuntos."
    )
# =====================================================
# 9. PRECIPITACIÓN ACUMULADA DEL EVENTO
# =====================================================

datos_acumulados = {}

datos_acumulados["Gauge"] = obs[
    ids_validos
].sum(
    axis=0,
    min_count=len(obs)
)

for nombre, tabla in series_diarias.items():

    datos_acumulados[nombre] = tabla[
        ids_validos
    ].sum(
        axis=0,
        min_count=len(obs)
    )

acumulados = pd.DataFrame(datos_acumulados)

orden = [
    "Gauge",
    "GF-ERA5",
    "GF-GFS",
    "KF-ERA5",
    "KF-GFS",
    "PERSIANN",
    "MERGE"
]

acumulados = acumulados[orden]

print("\nResumen de la precipitación acumulada:")
print(acumulados.describe().round(2))

# =====================================================
# 10. CONFIGURACIÓN GRÁFICA
# =====================================================

sns.set_theme(
    style="whitegrid",
    context="paper",
    font="DejaVu Sans"
)

colores = {
    "Gauge": "#708090",
    "GF-ERA5": "#87CEFA",
    "GF-GFS": "#4682B4",
    "KF-ERA5": "#BC8F8F",
    "KF-GFS": "#F5DEB3",
    "PERSIANN": "#5F9EA0",
    "MERGE": "#556B2F"
}

# Transformar los datos a formato largo

df_plot = acumulados.melt(
    var_name="Dataset",
    value_name="Precipitation"
)

# =====================================================
# 11. GENERAR BOXPLOT HORIZONTAL
# =====================================================

fig, ax = plt.subplots(figsize=(10, 7))

sns.boxplot(
    data=df_plot,
    x="Precipitation",
    y="Dataset",
    hue="Dataset",
    order=orden,
    hue_order=orden,
    palette=colores,
    orient="h",
    dodge=False,
    legend=False,
    width=0.55,
    saturation=1,
    linewidth=1.2,
    showfliers=False,
    showmeans=True,
    meanprops={
        "marker": "D",
        "markerfacecolor": "white",
        "markeredgecolor": "black",
        "markersize": 5
    },
    medianprops={
        "color": "black",
        "linewidth": 1.8
    },
    ax=ax
)


# =====================================================
# 12. EJES Y ESTÉTICA
# =====================================================

ax.set_ylabel("")

ax.set_xlabel(
    "Accumulated precipitation (mm)",
    fontsize=16,
    labelpad=12
)

ax.tick_params(
    axis="x",
    labelsize=14
)

ax.tick_params(
    axis="y",
    labelsize=14,
    length=0
)

ax.set_xlim(left=0)

ax.grid(
    axis="x",
    linestyle="--",
    linewidth=0.6,
    alpha=0.3
)

ax.grid(axis="y", visible=False)

ax.set_axisbelow(True)

sns.despine(ax=ax)

fig.tight_layout()
# =====================================================
# 13. GUARDAR FIGURA
# =====================================================

salida = os.path.join(
    carpeta_salida,
    "boxplot_precipitation_event.png"
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
# DIAGRAMAS DE DISPERSIÓN POR DENSIDAD
# =====================================================

import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.ticker as mticker

from scipy.stats import gaussian_kde

# =====================================================
# 1. CONFIGURACIÓN GRÁFICA
# =====================================================

sns.set_theme(
    style="ticks",
    context="paper",
    font="DejaVu Sans"
)

orden_scatter = [
    "GF-ERA5",
    "GF-GFS",
    "KF-ERA5",
    "KF-GFS",
    "PERSIANN",
    "MERGE"
]

titulos_scatter = [
    "(a) Grell-Freitas / ERA5",
    "(b) Grell-Freitas / GFS",
    "(c) Kain-Fritsch / ERA5",
    "(d) Kain-Fritsch / GFS",
    "(e) PERSIANN",
    "(f) MERGE"
]

# Tamaños de fuente y puntos

fontsize_titulos = 15
fontsize_ejes = 15
fontsize_ticks = 13
fontsize_colorbar = 15

tam_puntos = 12

# =====================================================
# 2. PREPARAR LOS DATOS Y CALCULAR DENSIDAD
# =====================================================

datos_scatter = {}

observado = obs[ids_validos].values.ravel()

for nombre in orden_scatter:

    estimado = series_diarias[nombre][
        ids_validos
    ].values.ravel()

    # Seleccionar pares válidos
    # Los ceros se conservan

    validos = (
        np.isfinite(observado) &
        np.isfinite(estimado) &
        (observado >= 0) &
        (estimado >= 0)
    )

    x = observado[validos]
    y = estimado[validos]

    # Calcular densidad bidimensional

    xy = np.vstack([x, y])

    densidad = gaussian_kde(xy)(xy)

    datos_scatter[nombre] = {
        "x": x,
        "y": y,
        "densidad": densidad
    }

    print(f"{nombre}: N = {len(x)}")

# =====================================================
# 3. NORMALIZACIÓN COMÚN DE DENSIDADES
# =====================================================

densidad_max = max(
    np.max(datos_scatter[nombre]["densidad"])
    for nombre in orden_scatter
)

for nombre in orden_scatter:

    datos_scatter[nombre]["densidad"] /= densidad_max

# =====================================================
# 4. LÍMITES COMUNES DE LOS EJES
# =====================================================

pp_max = max(
    max(
        np.max(datos_scatter[nombre]["x"]),
        np.max(datos_scatter[nombre]["y"])
    )
    for nombre in orden_scatter
)

intervalo_ticks = 50

limite = max(
    intervalo_ticks,
    np.ceil(pp_max / intervalo_ticks) * intervalo_ticks
)

# =====================================================
# 5. PALETA DE COLORES
# =====================================================

niveles = np.linspace(0, 1, 11)

cmap = plt.get_cmap(
    "rainbow",
    len(niveles) - 1
)

norm = mcolors.BoundaryNorm(
    niveles,
    cmap.N,
    clip=True
)

# =====================================================
# 6. GENERAR FIGURA DE SEIS PANELES
# =====================================================

fig, axs = plt.subplots(
    3, 2,
    figsize=(12, 16),
    sharex=True,
    sharey=True
)

axs = axs.flatten()

for i, nombre in enumerate(orden_scatter):

    ax = axs[i]

    x = datos_scatter[nombre]["x"]
    y = datos_scatter[nombre]["y"]
    densidad = datos_scatter[nombre]["densidad"]

    # Ordenar los puntos para que los de mayor
    # densidad queden por encima de los demás

    indices = np.argsort(densidad)

    x = x[indices]
    y = y[indices]
    densidad = densidad[indices]

    # Línea de referencia 1:1

    ax.plot(
        [0, limite],
        [0, limite],
        color="#555555",
        linewidth=1.0,
        linestyle="-",
        zorder=1
    )

    # Diagrama de dispersión

    ax.scatter(
        x,
        y,
        c=densidad,
        cmap=cmap,
        norm=norm,
        s=tam_puntos,
        edgecolors="none",
        alpha=1,
        zorder=3
    )

    # Títulos

    ax.set_title(
        titulos_scatter[i],
        fontsize=fontsize_titulos,
        fontweight="bold",
        pad=10
    )

    # Límites idénticos

    ax.set_xlim(0, limite)
    ax.set_ylim(0, limite)

    ax.set_aspect(
        "equal",
        adjustable="box"
    )

    # Ticks

    ax.xaxis.set_major_locator(
        mticker.MultipleLocator(intervalo_ticks)
    )

    ax.yaxis.set_major_locator(
        mticker.MultipleLocator(intervalo_ticks)
    )

    ax.tick_params(
        axis="both",
        direction="in",
        labelsize=fontsize_ticks,
        length=4,
        width=0.8,
        top=True,
        right=True
    )

    # Bordes de los paneles

    for spine in ax.spines.values():

        spine.set_visible(True)
        spine.set_linewidth(0.8)

    ax.grid(False)

    # Etiquetas solamente en la columna izquierda

    if i in [0, 2, 4]:

        ax.set_ylabel(
            r"Estimated precipitation (mm day$^{-1}$)",
            fontsize=fontsize_ejes,
            labelpad=10
        )

    # Etiquetas solamente en la última fila

    if i in [4, 5]:

        ax.set_xlabel(
            r"Gauge (mm day$^{-1}$)",
            fontsize=fontsize_ejes,
            labelpad=10
        )

# =====================================================
# 7. AJUSTAR ESPACIADO
# =====================================================

fig.subplots_adjust(
    left=0.12,
    right=0.85,
    bottom=0.06,
    top=0.96,
    wspace=0.12,
    hspace=0.15
)

# =====================================================
# 8. COLORBAR COMPARTIDA
# =====================================================

sm = cm.ScalarMappable(
    cmap=cmap,
    norm=norm
)

sm.set_array([])

cax = fig.add_axes([
    0.88,
    0.20,
    0.025,
    0.60
])

cbar = fig.colorbar(
    sm,
    cax=cax,
    orientation="vertical",
    boundaries=niveles,
    ticks=niveles,
    spacing="uniform",
    format=mticker.FormatStrFormatter("%.1f")
)

cbar.set_label(
    "Relative density",
    fontsize=fontsize_colorbar,
    labelpad=12
)

cbar.ax.tick_params(
    labelsize=fontsize_ticks,
    length=4,
    width=0.8
)

cbar.ax.minorticks_off()

cbar.outline.set_linewidth(0.8)

# =====================================================
# 9. GUARDAR FIGURA
# =====================================================

salida_scatter = os.path.join(
    carpeta_salida,
    "density_scatter_daily_precipitation.png"
)

fig.savefig(
    salida_scatter,
    dpi=600,
    format="png",
    bbox_inches="tight",
    facecolor="white"
)

print(f"\nFigura guardada en: {salida_scatter}")

plt.show()
plt.close(fig)