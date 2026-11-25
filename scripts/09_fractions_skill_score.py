# -*- coding: utf-8 -*-
"""Evaluate daily precipitation Fractions Skill Score against MERGE.

Author: Fernando Primo Forgioni (ORCID: 0000-0003-3233-8957).
Inputs: four WRF NetCDF files and MERGE_HORARIO_WRF.nc in data/.
Period: 28 April 2024 12 UTC to 3 May 2024 12 UTC.
Outputs: figures/FSS_Daily_WRF_MERGE.png and
results/Fractions_Skill_Score_Daily.xlsx.
Run: python scripts/09_fractions_skill_score.py
Original calculations and plotting settings are retained.
"""

import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.ndimage import uniform_filter
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

simulaciones = {
    "GF-ERA5": "Grell-Freitas-ERA5.nc",
    "GF-GFS": "Grell-Freitas-GFS.nc",
    "KF-ERA5": "Kain-Fritsch-ERA5.nc",
    "KF-GFS": "Kain-Fritsch-GFS.nc"
}

archivo_merge = "MERGE_HORARIO_WRF.nc"

# Check required input files before processing.
for name in (*simulaciones.values(), archivo_merge):
    input_path = carpeta / name
    if not input_path.is_file():
        raise FileNotFoundError(f"Required input file not found: {input_path}")

# Período del evento

fecha_inicio = pd.Timestamp("2024-04-28 12:00:00")
fecha_final = pd.Timestamp("2024-05-03 12:00:00")

# Resolución horizontal nominal (km)

dx = 12

# Umbrales de precipitación diaria (mm/día)

umbrales = [10, 25, 50, 100]

# Tamaño de las ventanas espaciales
# Número de puntos de grilla por lado

ventanas = [1, 3, 5, 7, 9, 11, 15, 21, 25]

# Convención temporal de MERGE
# True: timestamp al final del intervalo horario
# False: timestamp al inicio del intervalo horario

merge_fin_hora = True

# =====================================================
# 2. PALETA DE COLORES
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

orden = list(simulaciones.keys())

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
# Leer coordenadas espaciales
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
# Verificar coincidencia de grillas
# -----------------------------------------------------

def verificar_grilla(
    lat, lon, lat_ref, lon_ref, nombre
):

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
            f"La grilla de {nombre} no coincide."
        )


# =====================================================
# 4. PERÍODO TEMPORAL COMÚN
# =====================================================

fechas_wrf = pd.date_range(
    fecha_inicio,
    fecha_final,
    freq="h"
)

fechas_horarias = pd.date_range(
    fecha_inicio + pd.Timedelta(hours=1),
    fecha_final,
    freq="h"
)

fechas_diarias = pd.date_range(
    fecha_inicio + pd.Timedelta(days=1),
    fecha_final,
    freq="24h"
)

n_horas = len(fechas_horarias)

if n_horas % 24 != 0:

    raise ValueError(
        "El período no contiene días completos."
    )

n_dias = n_horas // 24

print("\n" + "=" * 60)
print("PERÍODO DEL EVENTO")
print("=" * 60)

print(f"Inicio: {fecha_inicio}")
print(f"Final: {fecha_final}")
print(f"Horas: {n_horas}")
print(f"Días: {n_dias}")

print("\nFinal de los acumulados diarios:")

for fecha in fechas_diarias:

    print(fecha)

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

        pp_acum = ds["RAINC"] + ds["RAINNC"]

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
        ).values.astype(np.float64)

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

    campos[nombre] = pp_h.astype(
        np.float32
    )

    print(
        f"Horas disponibles: {pp_h.shape[0]}"
    )

    del pp_acum, pp_h

# =====================================================
# 6. LEER MERGE HORARIO
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

    # Unificar fechas al final del intervalo horario

    if not merge_fin_hora:

        fechas = fechas + pd.Timedelta(hours=1)

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

    pp = ds["precipitation"]

    print(
        "Unidades de MERGE:",
        pp.attrs.get("units", "No especificadas")
    )

    dim_t = next(
        d for d in pp.dims
        if d not in dims_esp
    )

    pp = pp.isel({
        dim_t: indices
    })

    pp = pp.transpose(
        dim_t,
        *dims_esp
    ).values.astype(np.float32)

campos["MERGE"] = np.where(
    np.isfinite(pp) & (pp >= 0),
    pp,
    np.nan
)

del pp

print(
    f"Horas disponibles: "
    f"{campos['MERGE'].shape[0]}"
)

# =====================================================
# 7. CONVERTIR PRECIPITACIÓN HORARIA A DIARIA
# =====================================================

print("\n" + "=" * 60)
print("CALCULANDO ACUMULADOS DIARIOS")
print("=" * 60)

for nombre in campos:

    campo = campos[nombre]

    nt, ny, nx = campo.shape

    if nt != n_dias * 24:

        raise ValueError(
            f"Número incorrecto de horas en {nombre}."
        )

    # Convertir:
    # (120, 180, 200)
    # a:
    # (5, 24, 180, 200)
    #
    # Luego sumar las 24 horas de cada día.
    #
    # np.sum conserva los NaN, evitando acumulados
    # diarios incompletos.

    campo_diario = campo.reshape(
        n_dias,
        24,
        ny,
        nx
    ).sum(
        axis=1,
        dtype=np.float64
    )

    campos[nombre] = campo_diario.astype(
        np.float32
    )

    print(
        f"{nombre}: "
        f"{campos[nombre].shape[0]} días"
    )

    print(
        f"Máximo diario: "
        f"{np.nanmax(campos[nombre]):.2f} mm"
    )

    del campo, campo_diario

# =====================================================
# 8. MÁSCARA DE RIO GRANDE DO SUL
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

# =====================================================
# 9. MÁSCARA COMÚN DE DATOS VÁLIDOS
# =====================================================

mascara_valida = (
    np.isfinite(lat_ref) &
    np.isfinite(lon_ref)
)

for nombre, campo in campos.items():

    mascara_valida &= np.all(
        np.isfinite(campo) &
        (campo >= 0),
        axis=0
    )

# Centros de las ventanas dentro de Rio Grande do Sul

mascara_evaluacion = (
    mascara_rs &
    mascara_valida
)

n_grid = np.sum(
    mascara_evaluacion
)

print(
    f"\nPuntos de grilla válidos sobre RS: {n_grid}"
)

if n_grid == 0:

    raise ValueError(
        "No hay puntos comunes para evaluar."
    )

# =====================================================
# 10. FUNCIÓN PARA CALCULAR FRACCIONES
# =====================================================

def calcular_fracciones(
    campo_binario,
    ventana,
    centros_validos
):

    # El filtro se aplica solamente en el espacio.
    # No se realiza ningún suavizado temporal.

    fracciones = uniform_filter(
        campo_binario.astype(np.float64),
        size=(1, ventana, ventana),
        mode="constant",
        cval=0
    )

    return np.clip(
        fracciones[:, centros_validos],
        0,
        1
    )

# =====================================================
# 11. CENTROS VÁLIDOS POR TAMAÑO DE VENTANA
# =====================================================

centros_por_ventana = {}

print("\n" + "=" * 60)
print("VENTANAS ESPACIALES")
print("=" * 60)

for ventana in ventanas:

    cobertura = uniform_filter(
        mascara_valida.astype(np.float64),
        size=ventana,
        mode="constant",
        cval=0
    )

    # Exigir ventanas completas con datos válidos.
    # El centro debe estar dentro de RS.
    # La ventana puede extenderse fuera del estado.

    centros = (
        mascara_evaluacion &
        (cobertura >= 1 - 1e-10)
    )

    centros_por_ventana[ventana] = centros

    print(
        f"Ventana {ventana:2d} x {ventana:2d}: "
        f"{centros.sum()} centros válidos"
    )

# =====================================================
# 12. CALCULAR FRACTIONS SKILL SCORE
# =====================================================

resultados = []
referencias = []

observado = campos["MERGE"]

for umbral in umbrales:

    print("\n" + "=" * 60)
    print(f"Umbral diario: {umbral} mm/día")
    print("=" * 60)

    # -------------------------------------------------
    # Campo binario de MERGE
    # -------------------------------------------------

    obs_binario = (
        (observado >= umbral) &
        mascara_valida[None, :, :]
    )

    # Frecuencia observada de excedencia
    # en los centros de evaluación de RS

    frecuencia_obs = np.mean(
        obs_binario[:, mascara_evaluacion]
    )

    # -------------------------------------------------
    # Criterio convencional de habilidad útil
    # -------------------------------------------------

    fss_util = 0.5 + 0.5 * frecuencia_obs

    referencias.append({
        "Threshold_mm_day": umbral,
        "Observed_frequency": frecuencia_obs,
        "Observed_frequency_pct": frecuencia_obs * 100,
        "Useful_FSS": fss_util
    })

    print(
        f"Frecuencia observada: "
        f"{frecuencia_obs * 100:.3f} %"
    )

    print(
        f"FSS útil: {fss_util:.3f}"
    )

    # -------------------------------------------------
    # Campos binarios de las simulaciones
    # -------------------------------------------------

    binarios_wrf = {}

    for nombre in orden:

        binarios_wrf[nombre] = (
            (campos[nombre] >= umbral) &
            mascara_valida[None, :, :]
        )

    # -------------------------------------------------
    # Evaluar diferentes ventanas espaciales
    # -------------------------------------------------

    for ventana in ventanas:

        centros = centros_por_ventana[ventana]

        n_centros = centros.sum()

        print(
            f"\nVentana: {ventana} x {ventana}"
        )

        if n_centros == 0:

            for nombre in orden:

                resultados.append({
                    "Model": nombre,
                    "Threshold_mm_day": umbral,
                    "Window_gridpoints": ventana,
                    "Total_gridpoints": ventana ** 2,
                    "Scale_km": ventana * dx,
                    "FSS": np.nan,
                    "Useful_FSS": fss_util,
                    "Valid_centers": 0
                })

            continue

        # ---------------------------------------------
        # Fracciones observadas
        # ---------------------------------------------

        fracciones_obs = calcular_fracciones(
            obs_binario,
            ventana,
            centros
        )

        suma_obs2 = np.sum(
            fracciones_obs ** 2,
            dtype=np.float64
        )

        # ---------------------------------------------
        # Evaluar cada simulación
        # ---------------------------------------------

        for nombre in orden:

            fracciones_sim = calcular_fracciones(
                binarios_wrf[nombre],
                ventana,
                centros
            )

            # Numerador de la ecuación del FSS

            numerador = np.sum(
                (fracciones_sim - fracciones_obs) ** 2,
                dtype=np.float64
            )

            # Denominador de la ecuación del FSS

            denominador = (
                np.sum(
                    fracciones_sim ** 2,
                    dtype=np.float64
                ) +
                suma_obs2
            )

            # FSS agregado sobre los cinco días

            if denominador > 0:

                fss = 1 - (
                    numerador / denominador
                )

                fss = np.clip(
                    fss,
                    0,
                    1
                )

            else:

                # FSS indefinido cuando ninguno
                # de los campos supera el umbral

                fss = np.nan

            resultados.append({
                "Model": nombre,
                "Threshold_mm_day": umbral,
                "Window_gridpoints": ventana,
                "Total_gridpoints": ventana ** 2,
                "Scale_km": ventana * dx,
                "FSS": fss,
                "Useful_FSS": fss_util,
                "Valid_centers": n_centros
            })

            print(
                f"{nombre}: FSS = {fss:.3f}"
            )

            del fracciones_sim

        del fracciones_obs

    del binarios_wrf, obs_binario

# =====================================================
# 13. ORGANIZAR RESULTADOS
# =====================================================

df_fss = pd.DataFrame(
    resultados
)

df_referencias = pd.DataFrame(
    referencias
)

# =====================================================
# 14. VENTANA MÍNIMA CON HABILIDAD ÚTIL
# =====================================================

escalas_utiles = []

for umbral in umbrales:

    for nombre in orden:

        datos = df_fss[
            (df_fss["Model"] == nombre) &
            (df_fss["Threshold_mm_day"] == umbral)
        ].sort_values(
            "Window_gridpoints"
        )

        seleccion = datos[
            datos["FSS"] >= datos["Useful_FSS"]
        ]

        if not seleccion.empty:

            ventana_minima = seleccion[
                "Window_gridpoints"
            ].iloc[0]

            escala_minima_km = ventana_minima * dx

            puntos_totales = ventana_minima ** 2

        else:

            ventana_minima = np.nan
            escala_minima_km = np.nan
            puntos_totales = np.nan

        escalas_utiles.append({
            "Model": nombre,
            "Threshold_mm_day": umbral,
            "Minimum_window_gridpoints": ventana_minima,
            "Total_gridpoints": puntos_totales,
            "Minimum_useful_scale_km": escala_minima_km
        })

df_escalas = pd.DataFrame(
    escalas_utiles
)

print("\n" + "=" * 60)
print("VENTANA MÍNIMA CON HABILIDAD ÚTIL")
print("=" * 60)

print(
    df_escalas.pivot(
        index="Model",
        columns="Threshold_mm_day",
        values="Minimum_window_gridpoints"
    )
)

# =====================================================
# 15. CONFIGURACIÓN GRÁFICA
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

fontsize_titulos = 16
fontsize_ejes = 15
fontsize_ticks = 13
fontsize_leyenda = 14

marcadores = {
    "GF-ERA5": "o",
    "GF-GFS": "s",
    "KF-ERA5": "^",
    "KF-GFS": "D"
}

letras = [
    "a",
    "b",
    "c",
    "d"
]

# =====================================================
# 16. FIGURA DE CUATRO PANELES
# =====================================================

fig, axs = plt.subplots(
    2, 2,
    figsize=(13, 10),
    sharex=True,
    sharey=True
)

axs = axs.flatten()

for i, umbral in enumerate(umbrales):

    ax = axs[i]

    datos_umbral = df_fss[
        df_fss["Threshold_mm_day"] == umbral
    ]

    # -------------------------------------------------
    # Curvas de FSS
    # -------------------------------------------------

    for nombre in orden:

        datos = datos_umbral[
            datos_umbral["Model"] == nombre
        ].sort_values(
            "Window_gridpoints"
        )

        ax.plot(
            datos["Window_gridpoints"],
            datos["FSS"],
            color=colores[nombre],
            linewidth=2.2,
            marker=marcadores[nombre],
            markersize=6.5,
            markeredgecolor="#555555",
            markeredgewidth=0.4,
            label=nombre,
            zorder=4
        )

    # -------------------------------------------------
    # Línea de habilidad útil
    # -------------------------------------------------

    fss_util = df_referencias.loc[
        df_referencias["Threshold_mm_day"] == umbral,
        "Useful_FSS"
    ].iloc[0]

    ax.axhline(
        fss_util,
        color="#555555",
        linestyle="--",
        linewidth=1.4,
        label="Useful FSS",
        zorder=2
    )

    # -------------------------------------------------
    # Títulos
    # -------------------------------------------------

    ax.set_title(
        rf"({letras[i]}) $\geq$ {umbral} mm day$^{{-1}}$",
        fontsize=fontsize_titulos,
        fontweight="bold",
        loc="left",
        pad=12
    )

    # -------------------------------------------------
    # Límites y ticks
    # -------------------------------------------------

    ax.set_xlim(0, 26)

    ax.set_ylim(0, 1.02)

    ax.set_xticks(
        ventanas
    )

    ax.set_yticks(
        np.arange(0, 1.01, 0.2)
    )

    ax.tick_params(
        axis="both",
        labelsize=fontsize_ticks,
        direction="out",
        length=4,
        width=0.8
    )

    # -------------------------------------------------
    # Grilla
    # -------------------------------------------------

    ax.grid(
        axis="both",
        linestyle="--",
        linewidth=0.6,
        alpha=0.25
    )

    ax.set_axisbelow(True)

    sns.despine(ax=ax)

# =====================================================
# 17. ETIQUETAS DE LOS EJES
# =====================================================

axs[0].set_ylabel(
    "Fractions Skill Score",
    fontsize=fontsize_ejes
)

axs[2].set_ylabel(
    "Fractions Skill Score",
    fontsize=fontsize_ejes
)

axs[2].set_xlabel(
    "Grid points per side",
    fontsize=fontsize_ejes
)

axs[3].set_xlabel(
    "Grid points per side",
    fontsize=fontsize_ejes
)

# =====================================================
# 18. LEYENDA COMPARTIDA
# =====================================================

handles, labels = axs[0].get_legend_handles_labels()

fig.legend(
    handles,
    labels,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.02),
    ncol=5,
    fontsize=fontsize_leyenda,
    frameon=False,
    handlelength=2.8
)

fig.subplots_adjust(
    left=0.10,
    right=0.97,
    bottom=0.13,
    top=0.94,
    wspace=0.15,
    hspace=0.22
)

# =====================================================
# 19. GUARDAR FIGURA
# =====================================================

salida = os.path.join(
    args.figures_dir,
    "FSS_Daily_WRF_MERGE.png"
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
# 20. GUARDAR RESULTADOS EN EXCEL
# =====================================================

archivo_excel = os.path.join(
    args.results_dir,
    "Fractions_Skill_Score_Daily.xlsx"
)

with pd.ExcelWriter(
    archivo_excel,
    engine="openpyxl"
) as writer:

    df_fss.to_excel(
        writer,
        sheet_name="FSS",
        index=False
    )

    df_referencias.to_excel(
        writer,
        sheet_name="Useful_FSS",
        index=False
    )

    df_escalas.to_excel(
        writer,
        sheet_name="Minimum_useful_window",
        index=False
    )

    # Guardar una hoja por umbral

    for umbral in umbrales:

        tabla = df_fss[
            df_fss["Threshold_mm_day"] == umbral
        ].pivot(
            index="Window_gridpoints",
            columns="Model",
            values="FSS"
        )

        tabla = tabla.reindex(
            columns=orden
        )

        tabla.to_excel(
            writer,
            sheet_name=f"Threshold_{umbral}mm"
        )

    # Configuración del análisis

    pd.DataFrame({
        "Parameter": [
            "Start UTC",
            "End UTC",
            "Number of days",
            "Number of hours",
            "Grid spacing (km)",
            "Valid grid points over RS",
            "MERGE end-hour timestamp"
        ],
        "Value": [
            str(fecha_inicio),
            str(fecha_final),
            n_dias,
            n_horas,
            dx,
            n_grid,
            merge_fin_hora
        ]
    }).to_excel(
        writer,
        sheet_name="Configuration",
        index=False
    )

print(f"\nExcel guardado: {archivo_excel}")

print("\nProcesamiento finalizado.")