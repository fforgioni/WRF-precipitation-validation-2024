# -*- coding: utf-8 -*-

"""
Convert hourly MERGE GRIB2 files into a single NetCDF dataset.

Study period: 28 April - 3 May 2024
Source: CPTEC/INPE MERGE GPM HOURLY
Temporal resolution: 1 hour
Spatial resolution: Original MERGE grid

The script preserves the original precipitation units
and spatial coordinates. No interpolation is applied.
"""

import xarray as xr
import pandas as pd
import numpy as np

from pathlib import Path


# =====================================================
# 1. CONFIGURATION
# =====================================================

START_DATE = "2024-04-28 00:00"
END_DATE = "2024-05-03 23:00"

ROOT_DIR = Path(__file__).resolve().parents[1]

INPUT_DIR = ROOT_DIR / "data" / "raw" / "merge_hourly"

OUTPUT_DIR = ROOT_DIR / "data" / "processed"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FILE = OUTPUT_DIR / "MERGE_HORARIO_20240428_20240503.nc"

VARIABLE = "rdp"


# =====================================================
# 2. IDENTIFY INPUT FILES
# =====================================================

def get_input_files():
    """Identify and validate the expected hourly files."""

    dates = pd.date_range(
        start=START_DATE,
        end=END_DATE,
        freq="h"
    )

    files = []
    missing = []

    for date in dates:

        timestamp = date.strftime("%Y%m%d%H")

        filename = (
            f"MERGE_CPTEC_{timestamp}.grib2"
        )

        filepath = INPUT_DIR / filename

        if filepath.is_file():

            files.append((date, filepath))

        else:

            missing.append(filename)

    if missing:

        print("\nMissing files:")

        for filename in missing:
            print(filename)

        raise FileNotFoundError(
            f"{len(missing)} hourly files are missing."
        )

    return files


# =====================================================
# 3. READ INDIVIDUAL GRIB2 FILES
# =====================================================

def read_grib(filepath, date):
    """Extract precipitation from an individual GRIB2 file."""

    with xr.open_dataset(
        filepath,
        engine="cfgrib",
        backend_kwargs={"indexpath": ""}
    ) as ds:

        if VARIABLE not in ds.data_vars:

            raise ValueError(
                f"Variable {VARIABLE} not found in "
                f"{filepath.name}. "
                f"Available variables: {list(ds.data_vars)}"
            )

        pp = ds[[VARIABLE]].load()

        # Remove original temporal coordinates

        pp = pp.drop_vars(
            ["time", "step", "valid_time"],
            errors="ignore"
        )

        # Assign timestamp from the original filename

        pp = pp.expand_dims(time=[date])

        # Standardize precipitation variable name

        pp = pp.rename(
            {VARIABLE: "precipitation"}
        )

    return pp


# =====================================================
# 4. CONVERSION FUNCTION
# =====================================================

def convert_merge():

    files = get_input_files()

    datasets = []

    print("=" * 60)
    print("MERGE GRIB2 TO NETCDF CONVERSION")
    print("=" * 60)

    print(f"Input directory: {INPUT_DIR}")
    print(f"Expected files: {len(files)}")

    for i, (date, filepath) in enumerate(
        files,
        start=1
    ):

        print(
            f"[{i}/{len(files)}] "
            f"Processing: {filepath.name}"
        )

        ds = read_grib(filepath, date)

        if not datasets:

            print("\nOriginal precipitation attributes:")
            print(ds.precipitation.attrs)

            print("\nSpatial dimensions:")
            print(ds.sizes)

        datasets.append(ds)

    # =================================================
    # 5. CONCATENATE HOURLY DATASETS
    # =================================================

    print("\nConcatenating hourly datasets...")

    ds_final = xr.concat(
        datasets,
        dim="time",
        join="exact",
        compat="equals",
        coords="minimal"
    )

    ds_final = ds_final.sortby("time")

    # =================================================
    # 6. VALIDATE TEMPORAL COVERAGE
    # =================================================

    expected_dates = pd.date_range(
        START_DATE,
        END_DATE,
        freq="h"
    )

    available_dates = pd.DatetimeIndex(
        ds_final.time.values
    )

    if available_dates.has_duplicates:

        raise ValueError(
            "Duplicate timestamps were detected."
        )

    if not available_dates.equals(expected_dates):

        raise ValueError(
            "Temporal coverage does not match "
            "the expected study period."
        )

    # =================================================
    # 7. DATASET METADATA
    # =================================================

    ds_final.attrs = {
        "title": "MERGE hourly precipitation",
        "source": "CPTEC/INPE MERGE GPM HOURLY",
        "study_period": "2024-04-28 to 2024-05-03",
        "temporal_resolution": "1 hour",
        "spatial_resolution": "Original MERGE grid",
        "time_zone": "UTC",
        "time_reference": (
            "Timestamp from original GRIB2 filename"
        ),
        "processing": (
            "GRIB2 to NetCDF conversion without "
            "spatial interpolation or unit conversion"
        )
    }

    ds_final.precipitation.attrs["long_name"] = (
        "Hourly precipitation"
    )

    # =================================================
    # 8. SAVE COMPRESSED NETCDF
    # =================================================

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
        OUTPUT_FILE,
        engine="netcdf4",
        encoding=encoding
    )

    # =================================================
    # 9. OUTPUT SUMMARY
    # =================================================

    print("\n" + "=" * 60)
    print("CONVERSION COMPLETED")
    print("=" * 60)

    print(ds_final)

    print(f"\nExpected hours: {len(expected_dates)}")
    print(f"Available hours: {len(available_dates)}")

    print(
        "\nPrecipitation units:",
        ds_final.precipitation.attrs.get(
            "units",
            "Not specified"
        )
    )

    print(
        "Minimum precipitation:",
        float(ds_final.precipitation.min())
    )

    print(
        "Maximum precipitation:",
        float(ds_final.precipitation.max())
    )

    print(
        "Missing values (%):",
        float(
            ds_final.precipitation.isnull().mean()
        ) * 100
    )

    print(f"\nOutput file: {OUTPUT_FILE}")

    ds_final.close()


# =====================================================
# 10. MAIN EXECUTION
# =====================================================

if __name__ == "__main__":

    convert_merge()
