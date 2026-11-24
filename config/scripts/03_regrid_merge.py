# -*- coding: utf-8 -*-

"""
Regrid daily and hourly MERGE precipitation onto the WRF grid.

Study period: 28 April - 3 May 2024

Experiments:
    GF-ERA5: Grell-Freitas with ERA5 forcing
    GF-GFS:  Grell-Freitas with GFS forcing
    KF-ERA5: Kain-Fritsch with ERA5 forcing
    KF-GFS:  Kain-Fritsch with GFS forcing

Regridding method: Bilinear interpolation

The script verifies spatial consistency across the four
WRF experiments and preserves the original MERGE timestamps.
"""

import xarray as xr
import numpy as np

from pathlib import Path


# =====================================================
# 1. CONFIGURATION
# =====================================================

ROOT_DIR = Path(__file__).resolve().parents[1]

RAW_DIR = ROOT_DIR / "data" / "raw"

PROCESSED_DIR = ROOT_DIR / "data" / "processed"

WRF_DIR = RAW_DIR / "wrf"

OUTPUT_DIR = PROCESSED_DIR / "regridded"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# WRF experiments

WRF_FILES = [
    WRF_DIR / "Grell-Freitas-ERA5.nc",
    WRF_DIR / "Grell-Freitas-GFS.nc",
    WRF_DIR / "Kain-Fritsch-ERA5.nc",
    WRF_DIR / "Kain-Fritsch-GFS.nc"
]

REFERENCE_FILE = WRF_FILES[0]

# MERGE datasets

MERGE_DAILY = (
    RAW_DIR / "Merge-diario.nc"
)

MERGE_HOURLY = (
    PROCESSED_DIR /
    "MERGE_HORARIO_20240428_20240503.nc"
)

# Spatial interpolation settings

INTERPOLATION_METHOD = "linear"

SPATIAL_MARGIN = 0.3

GRID_TOLERANCE = 1e-5


# =====================================================
# 2. READ WRF GRID
# =====================================================

def read_wrf_grid(filepath):
    """Extract geographical coordinates from a WRF file."""

    with xr.open_dataset(
        filepath,
        decode_times=False
    ) as ds:

        lat = ds["XLAT"]
        lon = ds["XLONG"]

        for dim in ["Time", "time"]:

            if dim in lat.dims:
                lat = lat.isel({dim: 0})

            if dim in lon.dims:
                lon = lon.isel({dim: 0})

        lat = lat.transpose(
            "south_north",
            "west_east"
        ).values.copy()

        lon = lon.transpose(
            "south_north",
            "west_east"
        ).values.copy()

        dx = float(ds.attrs.get("DX", np.nan))
        dy = float(ds.attrs.get("DY", np.nan))

    # Standardize longitude convention

    lon = (lon + 180) % 360 - 180

    if not (
        np.isfinite(lat).all()
        and np.isfinite(lon).all()
    ):

        raise ValueError(
            f"Invalid grid coordinates: {filepath.name}"
        )

    return lat, lon, dx, dy


# =====================================================
# 3. VERIFY WRF GRID CONSISTENCY
# =====================================================

def verify_wrf_grids(lat_ref, lon_ref):
    """Verify that all experiments share the same grid."""

    print("\nVerifying WRF grids...")

    for filepath in WRF_FILES:

        lat, lon, _, _ = read_wrf_grid(filepath)

        same_grid = (
            lat.shape == lat_ref.shape
            and lon.shape == lon_ref.shape
            and np.allclose(
                lat,
                lat_ref,
                atol=GRID_TOLERANCE,
                rtol=0
            )
            and np.allclose(
                lon,
                lon_ref,
                atol=GRID_TOLERANCE,
                rtol=0
            )
        )

        if not same_grid:

            raise ValueError(
                f"Grid mismatch: {filepath.name}"
            )

        print(f"Grid verified: {filepath.name}")

    print("All WRF experiments share the same grid.")


# =====================================================
# 4. PRECIPITATION UNIT STANDARDIZATION
# =====================================================

def standardize_units(pp):
    """Standardize precipitation amounts to millimeters."""

    original_units = pp.attrs.get("units", "")

    units = (
        original_units.lower()
        .replace(" ", "")
        .replace("**", "^")
    )

    if units in [
        "mm",
        "kgm^-2",
        "kgm-2",
        "kg/m^2",
        "kg/m2",
        "kgm^(-2)"
    ]:

        # 1 kg/m2 of liquid water equals 1 mm

        pp.attrs["units"] = "mm"

    else:

        raise ValueError(
            f"Unsupported precipitation units: "
            f"{original_units}"
        )

    return pp, original_units


# =====================================================
# 5. REGRID MERGE DATASETS
# =====================================================

def regrid_merge(
    input_file,
    output_name,
    frequency,
    lat_wrf,
    lon_wrf,
    dx,
    dy
):
    """Interpolate MERGE precipitation onto the WRF grid."""

    print("\n" + "=" * 60)
    print(f"Processing: {input_file.name}")
    print("=" * 60)

    ny, nx = lat_wrf.shape

    target_lat = xr.DataArray(
        lat_wrf,
        dims=("south_north", "west_east")
    )

    target_lon = xr.DataArray(
        lon_wrf,
        dims=("south_north", "west_east")
    )

    with xr.open_dataset(input_file) as ds:

        # Standardize geographical coordinates

        ds = ds.assign_coords(
            longitude=(ds.longitude + 180) % 360 - 180
        )

        ds = ds.sortby("longitude")
        ds = ds.sortby("latitude")

        pp = ds["precipitation"]

        print("\nOriginal dimensions:")
        print(pp.sizes)

        # Standardize precipitation units

        pp, original_units = standardize_units(pp)

        print(f"\nOriginal units: {original_units}")
        print("Output units: mm")

        # Verify spatial coverage

        inside = (
            (lat_wrf >= float(ds.latitude.min()))
            & (lat_wrf <= float(ds.latitude.max()))
            & (lon_wrf >= float(ds.longitude.min()))
            & (lon_wrf <= float(ds.longitude.max()))
        )

        coverage = inside.mean() * 100

        print(f"\nSpatial coverage: {coverage:.2f}%")

        if not inside.all():

            raise ValueError(
                "MERGE does not fully cover the WRF domain."
            )

        # Subset MERGE to the WRF geographical extent

        pp = pp.sel(
            latitude=slice(
                lat_wrf.min() - SPATIAL_MARGIN,
                lat_wrf.max() + SPATIAL_MARGIN
            ),
            longitude=slice(
                lon_wrf.min() - SPATIAL_MARGIN,
                lon_wrf.max() + SPATIAL_MARGIN
            )
        )

        # Mask physically invalid precipitation values

        pp = pp.where(
            np.isfinite(pp) & (pp >= 0)
        )

        # Bilinear interpolation

        print("\nInterpolating onto the WRF grid...")

        pp_interp = pp.interp(
            latitude=target_lat,
            longitude=target_lon,
            method=INTERPOLATION_METHOD
        )

        pp_interp = pp_interp.transpose(
            "time",
            "south_north",
            "west_east"
        ).load()

        # Build output dataset

        ds_final = xr.Dataset(

            data_vars={

                "precipitation": (
                    ("time", "south_north", "west_east"),
                    pp_interp.values.astype("float32")
                )

            },

            coords={

                "time": pp_interp.time.values,

                "south_north": np.arange(ny),

                "west_east": np.arange(nx),

                "latitude": (
                    ("south_north", "west_east"),
                    lat_wrf
                ),

                "longitude": (
                    ("south_north", "west_east"),
                    lon_wrf
                )

            }

        )

    # =================================================
    # 6. METADATA
    # =================================================

    ds_final.precipitation.attrs = {

        "units": "mm",

        "long_name": (
            "MERGE precipitation on the WRF grid"
        ),

        "original_units": original_units,

        "temporal_resolution": frequency

    }

    ds_final.latitude.attrs = {
        "units": "degrees_north",
        "standard_name": "latitude"
    }

    ds_final.longitude.attrs = {
        "units": "degrees_east",
        "standard_name": "longitude"
    }

    ds_final.attrs = {

        "title": "MERGE precipitation regridded to WRF",

        "source": "CPTEC/INPE MERGE",

        "target_grid": "WRF",

        "reference_file": REFERENCE_FILE.name,

        "interpolation_method": "bilinear",

        "DX": float(dx),

        "DY": float(dy),

        "temporal_resolution": frequency,

        "time_zone": "UTC",

        "time_reference": "Original MERGE timestamps",

        "processing": (
            "Spatial interpolation without temporal aggregation"
        )

    }

    # =================================================
    # 7. SAVE NETCDF
    # =================================================

    output_file = OUTPUT_DIR / output_name

    encoding = {

        "precipitation": {

            "zlib": True,

            "complevel": 4,

            "dtype": "float32",

            "_FillValue": np.float32(np.nan)

        }

    }

    print("\nSaving NetCDF file...")

    ds_final.to_netcdf(
        output_file,
        engine="netcdf4",
        encoding=encoding
    )

    # =================================================
    # 8. QUALITY CONTROL
    # =================================================

    print("\nOutput dimensions:")
    print(ds_final.sizes)

    print(
        "\nFirst timestamp:",
        ds_final.time.values[0]
    )

    print(
        "Last timestamp:",
        ds_final.time.values[-1]
    )

    print(
        "Minimum precipitation:",
        float(ds_final.precipitation.min())
    )

    print(
        "Maximum precipitation:",
        float(ds_final.precipitation.max())
    )

    missing = float(
        ds_final.precipitation.isnull().mean()
    ) * 100

    print(f"Missing values: {missing:.2f}%")

    print(f"\nOutput file: {output_file}")

    ds_final.close()


# =====================================================
# 9. MAIN EXECUTION
# =====================================================

def main():

    print("=" * 60)
    print("MERGE REGRIDDING TO WRF")
    print("=" * 60)

    lat_wrf, lon_wrf, dx, dy = read_wrf_grid(
        REFERENCE_FILE
    )

    ny, nx = lat_wrf.shape

    print(f"\nReference file: {REFERENCE_FILE.name}")
    print(f"Grid dimensions: {ny} x {nx}")

    print(f"DX: {dx / 1000:.2f} km")
    print(f"DY: {dy / 1000:.2f} km")

    verify_wrf_grids(lat_wrf, lon_wrf)

    # Daily precipitation

    regrid_merge(
        input_file=MERGE_DAILY,
        output_name="MERGE_DIARIO_WRF.nc",
        frequency="daily",
        lat_wrf=lat_wrf,
        lon_wrf=lon_wrf,
        dx=dx,
        dy=dy
    )

    # Hourly precipitation

    regrid_merge(
        input_file=MERGE_HOURLY,
        output_name="MERGE_HORARIO_WRF.nc",
        frequency="hourly",
        lat_wrf=lat_wrf,
        lon_wrf=lon_wrf,
        dx=dx,
        dy=dy
    )

    print("\n" + "=" * 60)
    print("REGRIDDING COMPLETED")
    print("=" * 60)

    print(f"Output directory: {OUTPUT_DIR}")


if __name__ == "__main__":

    main()
