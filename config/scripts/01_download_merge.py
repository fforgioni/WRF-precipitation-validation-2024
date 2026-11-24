# -*- coding: utf-8 -*-

"""
Download hourly MERGE precipitation data from CPTEC/INPE.

Study period: 28 April - 3 May 2024
Product: MERGE GPM HOURLY
Format: GRIB2
Time reference: UTC

The script supports automatic retries and skips previously
downloaded valid files.
"""

import requests
import pandas as pd

from pathlib import Path
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# =====================================================
# 1. CONFIGURATION
# =====================================================

START_DATE = "2024-04-28 00:00"
END_DATE = "2024-05-03 23:00"

BASE_URL = (
    "https://ftp.cptec.inpe.br/"
    "modelos/tempo/MERGE/GPM/HOURLY"
)

# Repository root directory

ROOT_DIR = Path(__file__).resolve().parents[1]

# Output directory

OUTPUT_DIR = ROOT_DIR / "data" / "raw" / "merge_hourly"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# =====================================================
# 2. GRIB FILE VALIDATION
# =====================================================

def is_valid_grib(filepath):
    """Check GRIB header and end marker."""

    try:

        if filepath.stat().st_size < 20:
            return False

        with filepath.open("rb") as file:

            header = file.read(4)

            file.seek(-4, 2)
            footer = file.read(4)

        return header == b"GRIB" and footer == b"7777"

    except OSError:

        return False


# =====================================================
# 3. HTTP SESSION
# =====================================================

def create_session():
    """Create an HTTP session with automatic retries."""

    session = requests.Session()

    retry = Retry(
        total=3,
        backoff_factor=2,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"]
    )

    adapter = HTTPAdapter(max_retries=retry)

    session.mount("https://", adapter)

    return session


# =====================================================
# 4. DOWNLOAD FUNCTION
# =====================================================

def download_merge():

    dates = pd.date_range(
        start=START_DATE,
        end=END_DATE,
        freq="h"
    )

    downloaded = 0
    existing = 0
    failed = []

    print("=" * 60)
    print("MERGE HOURLY PRECIPITATION DOWNLOAD")
    print("=" * 60)

    print(f"Start date: {START_DATE} UTC")
    print(f"End date: {END_DATE} UTC")
    print(f"Expected files: {len(dates)}")
    print(f"Output directory: {OUTPUT_DIR}")

    with create_session() as session:

        for i, date in enumerate(dates, start=1):

            year = date.strftime("%Y")
            month = date.strftime("%m")
            day = date.strftime("%d")
            timestamp = date.strftime("%Y%m%d%H")

            filename = f"MERGE_CPTEC_{timestamp}.grib2"

            url = (
                f"{BASE_URL}/{year}/{month}/{day}/"
                f"{filename}"
            )

            destination = OUTPUT_DIR / filename

            temporary = OUTPUT_DIR / f"{filename}.part"

            # Skip previously downloaded valid files

            if is_valid_grib(destination):

                existing += 1

                print(
                    f"[{i}/{len(dates)}] "
                    f"Already exists: {filename}"
                )

                continue

            # Download file

            try:

                with session.get(
                    url,
                    stream=True,
                    timeout=(15, 120)
                ) as response:

                    response.raise_for_status()

                    with temporary.open("wb") as file:

                        for chunk in response.iter_content(
                            chunk_size=1024 * 1024
                        ):

                            if chunk:
                                file.write(chunk)

                # Validate downloaded GRIB file

                if not is_valid_grib(temporary):

                    raise ValueError("Invalid GRIB file")

                temporary.replace(destination)

                downloaded += 1

                print(
                    f"[{i}/{len(dates)}] "
                    f"Downloaded: {filename}"
                )

            except (
                requests.RequestException,
                ValueError,
                OSError
            ) as error:

                failed.append(filename)

                temporary.unlink(missing_ok=True)

                print(
                    f"[{i}/{len(dates)}] "
                    f"Failed: {filename} - {error}"
                )

    # =================================================
    # 5. DOWNLOAD SUMMARY
    # =================================================

    print("\n" + "=" * 60)
    print("DOWNLOAD SUMMARY")
    print("=" * 60)

    print(f"Expected files: {len(dates)}")
    print(f"Downloaded files: {downloaded}")
    print(f"Existing files: {existing}")
    print(f"Failed downloads: {len(failed)}")

    log_file = OUTPUT_DIR / "failed_downloads.txt"

    if failed:

        log_file.write_text(
            "\n".join(failed),
            encoding="utf-8"
        )

        print(f"\nFailed downloads saved to: {log_file}")

    else:

        log_file.unlink(missing_ok=True)

        print("\nAll requested files are available.")


# =====================================================
# 6. MAIN EXECUTION
# =====================================================

if __name__ == "__main__":

    download_merge()
