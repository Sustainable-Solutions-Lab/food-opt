# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Download SoilGrids organic carbon stock data via WCS at custom resolution.

This script uses the ISRIC WCS (Web Coverage Service) to download global soil
organic carbon stock data at a specified resolution. The native resolution is
250m, but this script resamples to a coarser resolution (e.g., 10km) to reduce
file size while maintaining global coverage.

Snakemake passes the ``snakemake`` object into this module; no standalone CLI
usage is supported.
"""

from pathlib import Path

import requests


def download_via_wcs(
    coverage_id: str,
    target_resolution_m: float,
    output_path: Path,
) -> None:
    """Download SoilGrids coverage via WCS and resample to target resolution.

    Parameters
    ----------
    coverage_id : str
        WCS coverage ID (e.g., "ocs_0-30cm_mean").
    target_resolution_m : float
        Target resolution in meters (e.g., 10000 for 10km).
    output_path : Path
        Output GeoTIFF file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # WCS endpoint
    # Coverage bounds from DescribeCoverage: x: -19949750 to 19861750, y: -6147500 to 8361000
    # Native grid: 159246 x 58034 pixels at 250m resolution
    # in EPSG:152160 (Interrupted Goode Homolosine)
    wcs_url = "https://maps.isric.org/mapserv"

    # Native resolution is 250m, calculate scaling factor for target resolution
    native_res_m = 250.0
    scale_factor = target_resolution_m / native_res_m

    # Calculate output dimensions
    native_width = 159246
    native_height = 58034
    output_width = int(native_width / scale_factor)
    output_height = int(native_height / scale_factor)

    print(f"Downloading {coverage_id} from ISRIC WCS...")
    print(f"Native resolution: {native_res_m}m, target: {target_resolution_m}m")
    print(f"Scale factor: {scale_factor:.2f}x")
    print(f"Output dimensions: {output_width} x {output_height} pixels")

    # WCS 2.0 ScaleSize parameter to request at specific dimensions
    params = {
        "map": "/map/ocs.map",
        "service": "WCS",
        "version": "2.0.1",
        "request": "GetCoverage",
        "coverageid": coverage_id,
        "format": "image/tiff",
        # Request full global extent
        "subset": [
            "x(-19949750,19861750)",
            "y(-6147500,8361000)",
        ],
        # Scale to target resolution using ScaleSize extension
        "scalesize": f"x({output_width}),y({output_height})",
    }

    # Download the coverage
    print("Requesting data from WCS...")
    response = requests.get(wcs_url, params=params, stream=True, timeout=1200)
    response.raise_for_status()

    # Save directly to output
    print(f"Saving to {output_path}...")
    with open(output_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

    print(f"Downloaded and saved {output_path}")


if __name__ == "__main__":
    download_via_wcs(
        coverage_id=snakemake.params.coverage_id,
        target_resolution_m=snakemake.params.target_resolution_m,
        output_path=Path(snakemake.output[0]),
    )
