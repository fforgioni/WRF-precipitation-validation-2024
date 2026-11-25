# -*- coding: utf-8 -*-
"""Categorical precipitation verification against MERGE over Rio Grande do Sul.

Author: Fernando Primo Forgioni (ORCID: 0000-0003-3233-8957).
Period: 30 April 2024 12 UTC to 1 May 2024 12 UTC.
Thresholds: 20, 50 and 100 mm/day by default.
Inputs: four WRF NetCDF files and MERGE_HORARIO_WRF.nc in data/.
Outputs: categorical maps in figures/ and metrics/arrays in results/.
Run: python scripts/08_categorical_verification.py
Original calculations and plotting settings are retained.
"""

import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

import cartopy.crs as ccrs
import cartopy.feature as cfeature

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
parser.add_argument("--threshold", type=int, choices=[20, 50, 100], nargs="+", default=[20, 50, 100])
args = parser.parse_args()
carpeta = args.data_dir
carpeta_salida = args.results_dir / "Mapas_CSI_Categoricos"
carpeta_salida.mkdir(parents=True, exist_ok=True)
args.figures_dir.mkdir(parents=True, exist_ok=True)

simulaciones = {
    "GF-ERA5": "Grell-Freitas-ERA5.nc",
    "GF-GFS": "Grell-Freitas-GFS.nc",
    "KF-ERA5": "Kain-Fritsch-ERA5.nc",
    "KF-GFS": "Kain-Fritsch-GFS.nc"
}

archivo_merge = "MERGE_HORARIO_WRF.nc"

for name in (*simulaciones.values(), archivo_merge):
    input_path = carpeta / name
    if not input_path.is_file():
        raise FileNotFoundError(f"Required input file not found: {input_path}")

# =====================================================
# FECHA Y UMBRAL DE PRECIPITACIÓN
# =====================================================

fecha_inicio = pd.Timestamp("2024-04-30 12:00:00")
fecha_final = fecha_inicio + pd.Timedelta(days=1)


# True: timestamp al final de la hora
# False: timestamp al inicio de la hora

merge_fin_hora = True

# Extensión geográfica

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
fontsize_metricas = 13
fontsize_leyenda = 14

orden = [
    "GF-ERA5",
    "GF-GFS",
    "KF-ERA5",
    "KF-GFS"
]

titulos = {
    "GF-ERA5": "(a) GF-ERA5",
    "GF-GFS": "(b) GF-GFS",
    "KF-ERA5": "(c) KF-ERA5",
    "KF-GFS": "(d) KF-GFS"
}

# =====================================================
# COLORES DE LAS CATEGORÍAS
# =====================================================

colores_categorias = {
    "No event": "#D9D9D9",
    "Hit": "#4DAF4A",
    "Miss": "#E41A1C",
    "False alarm": "#FDB462"
}

cmap = mcolors.ListedColormap([
    colores_categorias["No event"],
    colores_categorias["Hit"],
    colores_categorias["Miss"],
    colores_categorias["False alarm"]
])

cmap.set_bad(
    color="white",
    alpha=0
)

norm = mcolors.BoundaryNorm(
    boundaries=[-0.5, 0.5, 1.5, 2.5, 3.5],
    ncolors=cmap.N
)

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

            fechas.append(
                pd.to_datetime(texto)
            )

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
# Leer coordenadas
# -----------------------------------------------------

def leer_coordenadas(
    ds,
    nombre_lat,
    nombre_lon
):

    lat = ds[nombre_lat]
    lon = ds[nombre_lon]

    for dim in ["Time", "time"]:

        if dim in lat.dims:

            lat = lat.isel({dim: 0})

        if dim in lon.dims:

            lon = lon.isel({dim: 0})

    lat, lon = xr.broadcast(
        lat,
        lon
    )

    lon = lon.transpose(
        *lat.dims
    )

    return (
        lat.values,
        lon.values,
        lat.dims
    )


# -----------------------------------------------------
# Verificar coincidencia de grillas
# -----------------------------------------------------

def verificar_grilla(
    lat,
    lon,
    lat_ref,
    lon_ref,
    nombre
):

    if (
        lat.shape != lat_ref.shape or
        lon.shape != lon_ref.shape or
        not np.allclose(
            lat,
            lat_ref,
            atol=1e-4,
            rtol=0
        ) or
        not np.allclose(
            lon,
            lon_ref,
            atol=1e-4,
            rtol=0
        )
    ):

        raise ValueError(
            f"La grilla de {nombre} "
            "no coincide con la referencia."
        )


# -----------------------------------------------------
# Clasificación categórica
# -----------------------------------------------------

def categorizar(
    sim,
    obs,
    mascara,
    umbral
):

    cat = np.full(
        sim.shape,
        np.nan
    )

    validos = (
        mascara &
        np.isfinite(sim) &
        np.isfinite(obs) &
        (sim >= 0) &
        (obs >= 0)
    )

    # 0: No event
    # Ninguno alcanza el umbral

    cat[
        validos &
        (sim < umbral) &
        (obs < umbral)
    ] = 0

    # 1: Hit
    # Ambos alcanzan el umbral

    cat[
        validos &
        (sim >= umbral) &
        (obs >= umbral)
    ] = 1

    # 2: Miss
    # MERGE alcanza el umbral, WRF no

    cat[
        validos &
        (sim < umbral) &
        (obs >= umbral)
    ] = 2

    # 3: False alarm
    # WRF alcanza el umbral, MERGE no

    cat[
        validos &
        (sim >= umbral) &
        (obs < umbral)
    ] = 3

    return cat


# -----------------------------------------------------
# Calcular métricas categóricas
# -----------------------------------------------------

def calcular_metricas(cat):

    H = np.sum(cat == 1)

    M = np.sum(cat == 2)

    FA = np.sum(cat == 3)

    CN = np.sum(cat == 0)

    # Critical Success Index

    if (H + M + FA) > 0:

        csi = H / (H + M + FA)

    else:

        csi = np.nan

    # Probability of Detection

    if (H + M) > 0:

        pod = H / (H + M)

    else:

        pod = np.nan

    # False Alarm Ratio

    if (H + FA) > 0:

        far = FA / (H + FA)

    else:

        far = np.nan

    return H, M, FA, CN, csi, pod, far


# =====================================================
# 4. PERÍODO TEMPORAL
# =====================================================

# WRF: 25 registros acumulados para obtener 24 horas

fechas_wrf = pd.date_range(
    start=fecha_inicio,
    end=fecha_final,
    freq="h"
)

# MERGE: 24 registros horarios

fechas_horarias = pd.date_range(
    start=fecha_inicio + pd.Timedelta(hours=1),
    end=fecha_final,
    freq="h"
)

print("\n" + "=" * 60)
print("PERÍODO DE ANÁLISIS")
print("=" * 60)

print(f"Inicio: {fecha_inicio}")
print(f"Final: {fecha_final}")
print(f"Horas: {len(fechas_horarias)}")
print(f"Thresholds: {args.threshold} mm/day")

# =====================================================
# 5. LEER SIMULACIONES WRF
# =====================================================

campos = {}

lat_ref = None
lon_ref = None

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

        else:

            verificar_grilla(
                lat,
                lon,
                lat_ref,
                lon_ref,
                nombre
            )

        # Precipitación acumulada total

        pp_acum = (
            ds["RAINC"].astype(np.float64) +
            ds["RAINNC"].astype(np.float64)
        )

        dim_t = next(
            d for d in pp_acum.dims
            if d not in dims_esp
        )

        if fechas.duplicated().any():

            raise ValueError(
                f"Fechas duplicadas en {nombre}."
            )

        indices = fechas.get_indexer(
            fechas_wrf
        )

        if np.any(indices < 0):

            raise ValueError(
                f"Faltan horas en {nombre}."
            )

        pp_acum = pp_acum.isel({
            dim_t: indices
        })

        pp_acum = pp_acum.transpose(
            dim_t,
            *dims_esp
        ).values

    # Calcular precipitación horaria

    pp_h = np.diff(
        pp_acum,
        axis=0
    )

    # Verificar reinicios del acumulador

    if np.any(pp_h < -0.001):

        raise ValueError(
            f"Disminuciones del acumulador en {nombre}."
        )

    # Corregir pequeños negativos numéricos

    pp_h = np.where(
        (pp_h < 0) &
        (pp_h >= -0.001),
        0,
        pp_h
    )

    # Acumulado diario de las 24 horas

    pp_diaria = np.sum(
        pp_h,
        axis=0,
        dtype=np.float64
    )

    campos[nombre] = pp_diaria

    print(
        f"Horas disponibles: {pp_h.shape[0]}"
    )

    print(
        f"Máximo diario: {np.nanmax(pp_diaria):.2f} mm"
    )

    del pp_acum, pp_h

# =====================================================
# 6. LEER MERGE HORARIO
# =====================================================

print("\n" + "=" * 60)
print("Procesando MERGE")
print("=" * 60)

with xr.open_dataset(
    os.path.join(
        carpeta,
        archivo_merge
    )
) as ds:

    fechas = leer_fechas(ds)

    if not merge_fin_hora:

        fechas = fechas + pd.Timedelta(hours=1)

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

    print(
        "Unidades MERGE:",
        pp.attrs.get(
            "units",
            "No especificadas"
        )
    )

    dim_t = next(
        d for d in pp.dims
        if d not in dims_esp
    )

    if fechas.duplicated().any():

        raise ValueError(
            "Fechas duplicadas en MERGE."
        )

    indices = fechas.get_indexer(
        fechas_horarias
    )

    if np.any(indices < 0):

        raise ValueError(
            "Faltan horas en MERGE."
        )

    pp = pp.isel({
        dim_t: indices
    })

    pp = pp.transpose(
        dim_t,
        *dims_esp
    ).values.astype(np.float64)

# Eliminar valores negativos y no válidos

pp = np.where(
    np.isfinite(pp) & (pp >= 0),
    pp,
    np.nan
)

# Acumular las 24 horas
# np.sum conserva los NaN

campos["MERGE"] = np.sum(
    pp,
    axis=0,
    dtype=np.float64
)

print(f"Horas disponibles: {pp.shape[0]}")

print(
    f"Máximo diario MERGE: "
    f"{np.nanmax(campos['MERGE']):.2f} mm"
)

del pp

# =====================================================
# 7. GEOMETRÍA DE RIO GRANDE DO SUL
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

# =====================================================
# 8. MÁSCARA ESPACIAL
# =====================================================

mascara_rs = np.zeros(
    lat_ref.shape,
    dtype=bool
)

for i, j in np.ndindex(lat_ref.shape):

    if (
        np.isfinite(lat_ref[i, j]) and
        np.isfinite(lon_ref[i, j])
    ):

        punto = Point(
            float(lon_ref[i, j]),
            float(lat_ref[i, j])
        )

        mascara_rs[i, j] = (
            geometria_preparada.covers(punto)
        )

# Máscara común para todos los conjuntos

mascara_valida = mascara_rs.copy()

for nombre, campo in campos.items():

    mascara_valida &= (
        np.isfinite(campo) &
        (campo >= 0)
    )

n_grid = mascara_valida.sum()

print(
    f"\nPuntos de grilla válidos sobre RS: {n_grid}"
)

if n_grid == 0:

    raise ValueError(
        "No hay puntos válidos sobre RS."
    )

# =====================================================
# 9. LÍMITES GEOGRÁFICOS
# =====================================================

estados = cfeature.NaturalEarthFeature(
    category="cultural",
    name="admin_1_states_provinces_lines",
    scale="10m",
    facecolor="none"
)

rs_geom = [geometria_rs]

# =====================================================
# 10. FUNCIÓN PARA GENERAR MAPA BASE
# =====================================================

def mapa_base(ax):

    ax.set_extent(
        extent,
        crs=proyeccion
    )

    # Continente

    ax.add_feature(
        cfeature.LAND.with_scale("10m"),
        facecolor="#F5F5F5",
        zorder=0
    )

    # Océano

    ax.add_feature(
        cfeature.OCEAN.with_scale("10m"),
        facecolor="#DCECF5",
        zorder=0
    )

    # Lagos

    ax.add_feature(
        cfeature.LAKES.with_scale("10m"),
        facecolor="#DCECF5",
        edgecolor="gray",
        linewidth=0.3,
        zorder=1
    )

    # Costas

    ax.add_feature(
        cfeature.COASTLINE.with_scale("10m"),
        linewidth=0.6,
        zorder=3
    )

    # Límites internacionales

    ax.add_feature(
        cfeature.BORDERS.with_scale("10m"),
        linewidth=0.6,
        edgecolor="#555555",
        zorder=3
    )

    # Límites administrativos

    ax.add_feature(
        estados,
        edgecolor="#999999",
        linewidth=0.4,
        zorder=3
    )

    # Destacar Rio Grande do Sul

    ax.add_geometries(
        rs_geom,
        crs=proyeccion,
        facecolor="none",
        edgecolor="#333333",
        linewidth=1.0,
        zorder=4
    )

# =====================================================
# 11. GENERAR FIGURA
# =====================================================

for umbral in args.threshold:
    print(f"\nThreshold: {umbral} mm/day")
    fig, axs = plt.subplots(
        2, 2,
        figsize=(13, 10),
        subplot_kw={
            "projection": proyeccion
        }
    )

    axs = axs.flatten()

    resumen = []

    campos_categoricos = {}

    obs = campos["MERGE"]

    for i, nombre in enumerate(orden):

        ax = axs[i]

        mapa_base(ax)

        sim = campos[nombre]

        # -------------------------------------------------
        # Clasificación categórica
        # -------------------------------------------------

        cat = categorizar(
            sim=sim,
            obs=obs,
            mascara=mascara_valida,
            umbral=umbral
        )

        campos_categoricos[nombre] = cat

        # -------------------------------------------------
        # Graficar categorías
        # -------------------------------------------------

        ax.pcolormesh(
            lon_ref,
            lat_ref,
            cat,
            cmap=cmap,
            norm=norm,
            shading="auto",
            transform=proyeccion,
            zorder=2
        )

        # -------------------------------------------------
        # Calcular métricas
        # -------------------------------------------------

        H, M, FA, CN, csi, pod, far = (
            calcular_metricas(cat)
        )

        resumen.append({
            "Start_UTC": fecha_inicio,
            "End_UTC": fecha_final,
            "Threshold_mm_day": umbral,
            "Model": nombre,
            "Valid_gridpoints": n_grid,
            "Hits": H,
            "Misses": M,
            "False_alarms": FA,
            "Correct_negatives": CN,
            "CSI": csi,
            "POD": pod,
            "FAR": far
        })

        # -------------------------------------------------
        # Título de cada panel
        # -------------------------------------------------

        ax.set_title(
            titulos[nombre],
            fontsize=fontsize_titulos,
            fontweight="bold",
            pad=25
        )

        # -------------------------------------------------
        # Métricas debajo del nombre de la simulación
        # -------------------------------------------------

        texto_metricas = (
            f"CSI = {csi:.2f}    "
            f"POD = {pod:.2f}    "
            f"FAR = {far:.2f}"
        )

        ax.text(
            0.5,
            1.015,
            texto_metricas,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=fontsize_metricas,
            fontweight="normal"
        )

        print("\n" + nombre)

        print(f"Hits: {H}")
        print(f"Misses: {M}")
        print(f"False alarms: {FA}")
        print(f"CSI: {csi:.3f}")
        print(f"POD: {pod:.3f}")
        print(f"FAR: {far:.3f}")

    # =====================================================
    # 12. LEYENDA COMPARTIDA
    # =====================================================

    handles = [

        mpatches.Patch(
            facecolor=colores_categorias["No event"],
            edgecolor="#999999",
            linewidth=0.4,
            label="No event"
        ),

        mpatches.Patch(
            facecolor=colores_categorias["Hit"],
            edgecolor="none",
            label="Hit"
        ),

        mpatches.Patch(
            facecolor=colores_categorias["Miss"],
            edgecolor="none",
            label="Miss"
        ),

        mpatches.Patch(
            facecolor=colores_categorias["False alarm"],
            edgecolor="none",
            label="False alarm"
        )
    ]

    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.025),
        ncol=4,
        frameon=False,
        fontsize=fontsize_leyenda,
        handlelength=2,
        handleheight=1.2,
        columnspacing=2.5
    )

    # =====================================================
    # 13. AJUSTAR ESPACIADO
    # =====================================================

    fig.subplots_adjust(
        left=0.04,
        right=0.98,
        bottom=0.10,
        top=0.92,
        wspace=0.06,
        hspace=0.20
    )

    # =====================================================
    # 14. GUARDAR FIGURA
    # =====================================================

    salida = os.path.join(
        args.figures_dir,
        f"Categorical_skill_{umbral}mm.png"
    )

    fig.savefig(
        salida,
        dpi=600,
        format="png",
        bbox_inches="tight",
        facecolor="white"
    )

    print(
        f"\nFigura guardada: {salida}"
    )

    plt.show()
    plt.close(fig)

    # =====================================================
    # 15. GUARDAR MÉTRICAS EN EXCEL
    # =====================================================

    df_resumen = pd.DataFrame(
        resumen
    )

    archivo_excel = os.path.join(
        carpeta_salida,
        f"Metricas_20240430_{umbral}mm.xlsx"
    )

    with pd.ExcelWriter(
        archivo_excel,
        engine="openpyxl"
    ) as writer:

        df_resumen.to_excel(
            writer,
            sheet_name="Metrics",
            index=False
        )

        pd.DataFrame({
            "Parameter": [
                "Start UTC",
                "End UTC",
                "Threshold (mm/day)",
                "Valid grid points",
                "MERGE end-hour timestamp"
            ],
            "Value": [
                str(fecha_inicio),
                str(fecha_final),
                umbral,
                n_grid,
                merge_fin_hora
            ]
        }).to_excel(
            writer,
            sheet_name="Configuration",
            index=False
        )

    print(
        f"\nExcel guardado: {archivo_excel}"
    )

    # =====================================================
    # 16. GUARDAR CAMPOS CATEGÓRICOS
    # =====================================================

    for nombre in orden:

        np.save(
            os.path.join(
                carpeta_salida,
                f"Categorical_{nombre}_20240430_{umbral}mm.npy"
            ),
            campos_categoricos[nombre]
        )

    print("\nProcesamiento finalizado.")
