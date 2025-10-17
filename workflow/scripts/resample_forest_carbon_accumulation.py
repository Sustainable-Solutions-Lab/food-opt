"""
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: GPL-3.0-or-later
"""

from pathlib import Path

import numpy as np
import xarray as xr
from affine import Affine
from rasterio.crs import CRS
from rasterio.warp import Resampling, reproject
import rasterio


NO_DATA = -9999.0


def _load_target_grid(
    grid_path: str,
) -> tuple[Affine, CRS, tuple[int, int], dict[str, np.ndarray]]:
    ds = xr.load_dataset(grid_path)
    try:
        transform = Affine(*ds.attrs["transform"])
    except KeyError as exc:
        raise ValueError("grid definition missing affine transform metadata") from exc
    try:
        crs = CRS.from_wkt(ds.attrs["crs_wkt"])
    except KeyError as exc:
        raise ValueError("grid definition missing CRS metadata") from exc

    height = int(ds.attrs.get("height", ds.sizes["y"]))
    width = int(ds.attrs.get("width", ds.sizes["x"]))

    cols = np.arange(width, dtype=np.float64)
    rows = np.arange(height, dtype=np.float64)
    lon = transform.c + (cols + 0.5) * transform.a
    lat = transform.f + (rows + 0.5) * transform.e

    coords = {
        "y": lat.astype(np.float32),
        "x": lon.astype(np.float32),
    }

    return transform, crs, (height, width), coords


def main() -> None:
    grid_path: str = snakemake.input.grid  # type: ignore[name-defined]
    regrowth_path: str = snakemake.input.regrowth_raw  # type: ignore[name-defined]
    output_path: str = snakemake.output.regrowth  # type: ignore[name-defined]
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    (
        target_transform,
        target_crs,
        target_shape,
        coords,
    ) = _load_target_grid(grid_path)

    with rasterio.open(regrowth_path) as src:
        src_crs = src.crs
        if src_crs is None:
            raise ValueError("regrowth raster missing CRS information")
        src_transform = src.transform

        dst = np.full(target_shape, NO_DATA, dtype=np.float32)
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src_transform,
            src_crs=src_crs,
            dst_transform=target_transform,
            dst_crs=target_crs,
            src_nodata=src.nodata,
            dst_nodata=NO_DATA,
            resampling=Resampling.average,
        )

    dst[dst == NO_DATA] = np.nan

    ds_out = xr.Dataset(
        {
            "regrowth_tc_per_ha_yr": (("y", "x"), dst.astype(np.float32)),
        },
        coords=coords,
        attrs={
            "transform": target_transform.to_gdal(),
            "crs_wkt": target_crs.to_wkt(),
        },
    )

    encoding = {
        "regrowth_tc_per_ha_yr": {"zlib": True, "complevel": 4, "dtype": "float32"}
    }
    ds_out.to_netcdf(output_path, encoding=encoding)


if __name__ == "__main__":
    main()
