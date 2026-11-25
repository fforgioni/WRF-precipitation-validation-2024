# -*- coding: utf-8 -*-
"""Calculate daily station-based precipitation verification metrics.

Author: Fernando Primo Forgioni (ORCID: 0000-0003-3233-8957).
Inputs: data/Estaciones.xlsx, four WRF NetCDF files, PERSIANN_WRF.nc,
and Merge-diario.nc. The original fallback station workbook is retained.
Outputs: results/Metricas.xlsx (metrics, daily series, grid matches, settings).
Run: python scripts/04_station_metrics.py
Daily intervals start at 12 UTC; POD and FAR use a 1 mm/day threshold.
"""

import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

from scipy.spatial import cKDTree

# =====================================================
# 1. CONFIGURACIÓN
# =====================================================

repository_root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--data-dir", type=Path, default=repository_root / "data")
parser.add_argument("--results-dir", type=Path, default=repository_root / "results")
args = parser.parse_args()
carpeta = args.data_dir
args.results_dir.mkdir(parents=True, exist_ok=True)

archivo_estaciones = os.path.join(
    carpeta, "Estaciones.xlsx"
)

archivo_salida = os.path.join(
    args.results_dir, "Metricas.xlsx"
)

simulaciones = {
    "GF-ERA5": "Grell-Freitas-ERA5.nc",
    "GF-GFS": "Grell-Freitas-GFS.nc",
    "KF-ERA5": "Kain-Fritsch-ERA5.nc",
    "KF-GFS": "Kain-Fritsch-GFS.nc"
}

# Umbral de precipitación para POD y FAR (mm/día)

umbral = 1.0

productos = {
    "PERSIANN": "PERSIANN_WRF.nc",
    "MERGE": "Merge-diario.nc"
}

# Check required inputs before processing.
for input_path in [Path(archivo_estaciones)] + [carpeta / name for name in (*simulaciones.values(), *productos.values())]:
    if not input_path.is_file():
        raise FileNotFoundError(f"Required input file not found: {input_path}")

# PERSIANN: supuesto de acumulado horario en mm.
# False: etiqueta al inicio de la hora; True: al final.
persiann_fin_hora = False

# MERGE: supuesto de etiqueta al final del acumulado de 24 h.
merge_fin_dia = True

# Hora UTC de inicio del acumulado diario
# Ajustar según las observaciones

hora_inicio = 12  # Verificar que las estaciones correspondan a 12–12 UTC

# Si los NetCDF ya contienen precipitación diaria,
# indicar aquí el nombre de la variable.
# Dejar None para utilizar RAINC + RAINNC.

variable_diaria = None

# =====================================================
# 2. LEER ESTACIONES
# =====================================================

df = pd.read_excel(
    archivo_estaciones,
    sheet_name="Sheet1"
)

# Alternativa si Estaciones.xlsx no contiene
# las precipitaciones observadas

if df.empty or len(df.columns) <= 5:

    df = pd.read_excel(
        os.path.join(
            carpeta,
            "Datos-estaciones_unificadas.xlsx"
        ),
        sheet_name="Estaciones_unificadas"
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

# Identificar columnas con fechas

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

# Observaciones: filas = fechas, columnas = estaciones

obs = df[columnas_fechas].apply(
    pd.to_numeric,
    errors="coerce"
)

obs = obs.T

obs.index = pd.DatetimeIndex(fechas)
obs.columns = df["id_estacion"].values

obs = obs.sort_index()

# Eliminar valores negativos inválidos

obs = obs.where(obs >= 0)

print(f"Estaciones: {len(df)}")
print(f"Fecha inicial: {obs.index.min()}")
print(f"Fecha final: {obs.index.max()}")

# Información de las estaciones

metadata = df[[
    "id_estacion",
    "nombre",
    "institucion",
    "latitud",
    "longitud"
]].copy()

# =====================================================
# 3. FUNCIÓN PARA BUSCAR EL PUNTO MÁS CERCANO
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

# =====================================================
# 4. FUNCIÓN PARA LEER FECHAS WRF
# =====================================================

def leer_fechas(ds):

    if "Times" in ds.variables:

        fechas_wrf = []

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

            fechas_wrf.append(
                pd.to_datetime(texto)
            )

        return pd.DatetimeIndex(fechas_wrf)

    for nombre in ["time", "Time", "XTIME"]:

        if nombre in ds.coords:

            valores = ds[nombre].values

            if np.issubdtype(
                valores.dtype, np.datetime64
            ):

                return pd.DatetimeIndex(valores)

    raise ValueError(
        "No se pudieron identificar las fechas del NetCDF."
    )

# =====================================================
# 5. FUNCIÓN PARA CALCULAR MÉTRICAS
# =====================================================

def calcular_metricas(observado, simulado):

    observado = np.asarray(observado, dtype=float)
    simulado = np.asarray(simulado, dtype=float)

    validos = (
        np.isfinite(observado) &
        np.isfinite(simulado)
    )

    o = observado[validos]
    s = simulado[validos]

    n = len(o)

    if n == 0:

        return {
            "CC": np.nan,
            "RMSE": np.nan,
            "BIAS": np.nan,
            "POD": np.nan,
            "FAR": np.nan,
            "N": 0
        }

    # Correlación de Pearson

    if n >= 3 and np.std(o) > 0 and np.std(s) > 0:

        cc = np.corrcoef(o, s)[0, 1]

    else:

        cc = np.nan

    # RMSE

    rmse = np.sqrt(
        np.mean((s - o)**2)
    )

    # BIAS

    bias = np.mean(s - o)

    # Tabla de contingencia

    evento_obs = o >= umbral
    evento_sim = s >= umbral

    H = np.sum(evento_obs & evento_sim)
    M = np.sum(evento_obs & ~evento_sim)
    F = np.sum(~evento_obs & evento_sim)

    # Probability of Detection

    pod = H / (H + M) if (H + M) > 0 else np.nan

    # False Alarm Ratio

    far = F / (H + F) if (H + F) > 0 else np.nan

    return {
        "CC": cc,
        "RMSE": rmse,
        "BIAS": bias,
        "POD": pod,
        "FAR": far,
        "N": n
    }

# =====================================================
# 6. PROCESAR SIMULACIONES
# =====================================================

series_simuladas = {}
puntos_grilla = []

for nombre, archivo in simulaciones.items():

    print("\n" + "=" * 60)
    print(f"Procesando: {nombre}")
    print("=" * 60)

    ruta = os.path.join(carpeta, archivo)

    with xr.open_dataset(ruta) as ds:

        fechas_wrf = leer_fechas(ds)

        # ---------------------------------------------
        # Coordenadas WRF
        # ---------------------------------------------

        lat = ds["XLAT"]
        lon = ds["XLONG"]

        # Eliminar dimensión temporal de las coordenadas

        for dim in ["Time", "time"]:

            if dim in lat.dims:
                lat = lat.isel({dim: 0})

            if dim in lon.dims:
                lon = lon.isel({dim: 0})

        dim_y, dim_x = lat.dims

        lat_wrf = lat.values
        lon_wrf = lon.values

        # ---------------------------------------------
        # Vecino más cercano
        # ---------------------------------------------

        xyz_wrf = xyz(lat_wrf, lon_wrf)

        arbol = cKDTree(xyz_wrf)

        distancias, indices = arbol.query(
            xyz_estaciones,
            k=1
        )

        iy, ix = np.unravel_index(
            indices,
            lat_wrf.shape
        )

        # Distancia estación-grilla en km

        distancia_km = 2 * 6371 * np.arcsin(
            np.clip(distancias / 2, 0, 1)
        )

        tabla_grilla = metadata.copy()

        tabla_grilla["Simulacion"] = nombre
        tabla_grilla["Lat_WRF"] = lat_wrf[iy, ix]
        tabla_grilla["Lon_WRF"] = lon_wrf[iy, ix]
        tabla_grilla["Distancia_km"] = distancia_km

        puntos_grilla.append(tabla_grilla)

        # ---------------------------------------------
        # Extraer precipitación de las estaciones
        # ---------------------------------------------

        seleccion = {
            dim_y: xr.DataArray(iy, dims="estacion"),
            dim_x: xr.DataArray(ix, dims="estacion")
        }

        if variable_diaria is not None:

            pp = ds[variable_diaria].isel(seleccion)

        else:

            # Precipitación convectiva y no convectiva

            rainc = ds["RAINC"].isel(seleccion)
            rainnc = ds["RAINNC"].isel(seleccion)

            pp = rainc + rainnc

        dim_t = next(
            d for d in pp.dims if d != "estacion"
        )

        valores = pp.transpose(
            dim_t, "estacion"
        ).values

    tabla = pd.DataFrame(
        valores,
        index=fechas_wrf,
        columns=df["id_estacion"].values
    )

    tabla = tabla.sort_index()

    if tabla.index.duplicated().any():

        raise ValueError(
            f"Hay fechas duplicadas en {nombre}."
        )

    # ---------------------------------------------
    # Calcular precipitación diaria
    # ---------------------------------------------

    if variable_diaria is None:

        # Restar el acumulado al inicio del día
        # del acumulado al final del día.

        inicio = obs.index + pd.Timedelta(
            hours=hora_inicio
        )

        final = inicio + pd.Timedelta(days=1)

        pp_inicio = tabla.reindex(inicio).values
        pp_final = tabla.reindex(final).values

        diferencia = pp_final - pp_inicio

        # Detectar posibles reinicios del acumulador

        if np.any(tabla.diff().values < -0.001):

            raise ValueError(
                f"Se detectaron reinicios o disminuciones "
                f"del acumulador en {nombre}."
            )

        # Corregir pequeños negativos numéricos

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

    else:

        # Los datos ya representan acumulados diarios

        tabla.index = tabla.index.normalize()

        if tabla.index.duplicated().any():

            raise ValueError(
                f"{nombre} tiene más de un registro diario."
            )

        pp_diaria = tabla.reindex(obs.index)

    pp_diaria = pp_diaria.where(pp_diaria >= 0)

    series_simuladas[nombre] = pp_diaria

    print(f"Datos válidos: {pp_diaria.count().sum()}")
    print(f"Distancia media: {distancia_km.mean():.2f} km")

# =====================================================
# 6B. PROCESAR PERSIANN Y MERGE
# =====================================================

for nombre, archivo in productos.items():
    print(f"\nProcesando: {nombre}")
    with xr.open_dataset(os.path.join(carpeta, archivo)) as ds:
        lat = ds["XLAT" if nombre == "PERSIANN" else "latitude"]
        lon = ds["XLONG" if nombre == "PERSIANN" else "longitude"]
        dim_y, dim_x = lat.dims
        lat_g, lon_g = lat.values, lon.values
        arbol = cKDTree(xyz(lat_g, lon_g))
        distancias, indices = arbol.query(xyz_estaciones, k=1)
        iy, ix = np.unravel_index(indices, lat_g.shape)
        distancia_km = 2 * 6371 * np.arcsin(np.clip(distancias / 2, 0, 1))
        tabla_grilla = metadata.copy()
        tabla_grilla["Simulacion"] = nombre
        tabla_grilla["Lat_WRF"] = lat_g[iy, ix]
        tabla_grilla["Lon_WRF"] = lon_g[iy, ix]
        tabla_grilla["Distancia_km"] = distancia_km
        puntos_grilla.append(tabla_grilla)

        pp = ds["precipitation"].isel({
            dim_y: xr.DataArray(iy, dims="estacion"),
            dim_x: xr.DataArray(ix, dims="estacion")
        }).transpose("time", "estacion")
        tabla = pd.DataFrame(pp.values, index=pd.DatetimeIndex(ds.time.values),
                             columns=obs.columns).sort_index()
    tabla = tabla.where(tabla >= 0)
    if tabla.index.duplicated().any():
        raise ValueError(f"Fechas duplicadas en {nombre}")

    if nombre == "PERSIANN":
        # Sumar exactamente 24 registros horarios por ventana diaria.
        if persiann_fin_hora:
            tabla.index -= pd.Timedelta(hours=1)
        diarios = []
        for fecha in obs.index:
            inicio = fecha + pd.Timedelta(hours=hora_inicio)
            horas = pd.date_range(inicio, periods=24, freq="h")
            diarios.append(tabla.reindex(horas).sum(min_count=24).values)
        pp_diaria = pd.DataFrame(diarios, index=obs.index, columns=obs.columns)
    else:
        # MERGE conserva sus horas originales (12 UTC).
        # Asociar al inicio o al final del intervalo, según configuración.
        horas_merge = tabla.index.hour.unique()
        if len(horas_merge) != 1 or horas_merge[0] != hora_inicio:
            raise ValueError("MERGE tiene etiquetas a 12 UTC. Configure hora_inicio=12 "
                             "si las estaciones representan ventanas 12–12 UTC; "
                             "no mezcle acumulados diarios con ventanas distintas.")
        if merge_fin_dia:
            tabla.index -= pd.Timedelta(days=1)
        tabla.index = tabla.index.normalize()
        pp_diaria = tabla.reindex(obs.index)
    series_simuladas[nombre] = pp_diaria
    print(f"Datos válidos: {pp_diaria.count().sum()}")
    print(f"Distancia media: {distancia_km.mean():.2f} km")

# =====================================================
# 7. CALCULAR MÉTRICAS POR ESTACIÓN
# =====================================================

metricas = ["CC", "RMSE", "BIAS", "POD", "FAR", "N"]

resultados = {
    metrica: metadata.copy()
    for metrica in metricas
}

# Utilizar las mismas fechas válidas para los seis conjuntos

mascara_comun = obs.notna()

for tabla in series_simuladas.values():

    mascara_comun &= tabla.notna()

detalle = []

for nombre, pp_diaria in series_simuladas.items():

    valores_metricas = {
        metrica: []
        for metrica in metricas
    }

    for estacion in df["id_estacion"]:

        validos = mascara_comun[estacion]

        observado = obs[estacion][validos]

        simulado = pp_diaria[estacion][validos]

        resultado = calcular_metricas(
            observado.values,
            simulado.values
        )

        for metrica in metricas:

            valores_metricas[metrica].append(
                resultado[metrica]
            )

        detalle.append({
            "id_estacion": estacion,
            "Simulacion": nombre,
            **resultado
        })

    for metrica in metricas:

        resultados[metrica][nombre] = (
            valores_metricas[metrica]
        )

# =====================================================
# 8. CONSTRUIR TABLA DE PRECIPITACIONES DIARIAS
# =====================================================

series = []

for estacion in df["id_estacion"]:

    temporal = pd.DataFrame({
        "id_estacion": estacion,
        "Fecha": obs.index,
        "OBS": obs[estacion].values,
        "Valido_comun": mascara_comun[estacion].values
    })

    for nombre, tabla in series_simuladas.items():

        temporal[nombre] = tabla[estacion].values

    series.append(temporal)

series_diarias = pd.concat(
    series,
    ignore_index=True
)

series_diarias = series_diarias.merge(
    metadata,
    on="id_estacion",
    how="left"
)

detalle = pd.DataFrame(detalle).merge(
    metadata,
    on="id_estacion",
    how="left"
)

# =====================================================
# 9. GUARDAR EXCEL
# =====================================================

with pd.ExcelWriter(
    archivo_salida,
    engine="openpyxl"
) as writer:

    for metrica, tabla in resultados.items():

        tabla.to_excel(
            writer,
            sheet_name=metrica,
            index=False
        )

    detalle.to_excel(
        writer,
        sheet_name="Todas_metricas",
        index=False
    )

    series_diarias.to_excel(
        writer,
        sheet_name="Series_diarias",
        index=False
    )

    pd.concat(
        puntos_grilla,
        ignore_index=True
    ).to_excel(
        writer,
        sheet_name="Puntos_grilla",
        index=False
    )

    pd.DataFrame({
        "Parametro": [
            "Umbral POD y FAR (mm/día)",
            "Hora de inicio UTC",
            "Fechas válidas",
            "BIAS",
            "PERSIANN etiqueta fin de hora",
            "MERGE etiqueta fin de día"
        ],
        "Valor": [
            umbral,
            hora_inicio,
            "Comunes a los seis conjuntos",
            "Simulado menos observado",
            persiann_fin_hora,
            merge_fin_dia
        ]
    }).to_excel(
        writer,
        sheet_name="Configuracion",
        index=False
    )

    # Formato de las hojas

    for hoja in writer.sheets.values():

        hoja.freeze_panes = "A2"
        hoja.auto_filter.ref = hoja.dimensions

        for columna in hoja.columns:

            ancho = max(
                len(str(celda.value))
                if celda.value is not None else 0
                for celda in columna
            )

            hoja.column_dimensions[
                columna[0].column_letter
            ].width = min(ancho + 3, 40)

        for fila in hoja.iter_rows(min_row=2):

            for celda in fila:

                if isinstance(celda.value, float):

                    celda.number_format = "0.000"

print("\n" + "=" * 60)
print("PROCESAMIENTO FINALIZADO")
print(f"Excel guardado en: {archivo_salida}")
print("=" * 60)