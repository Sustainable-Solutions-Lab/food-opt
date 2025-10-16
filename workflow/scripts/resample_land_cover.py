"""
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: GPL-3.0-or-later
"""

from pathlib import Path

import numpy as np
import xarray as xr


COARSEN_FACTOR = 30
CHUNKS = {"lat": 900, "lon": 1800}

FOREST_CLASSES: tuple[int, ...] = (50, 60, 70, 80, 90, 160, 170, 100)
GRASSLAND_CLASSES: tuple[int, ...] = (110, 130, 140, 150)
CROPLAND_CLASSES: tuple[int, ...] = (10, 20, 30, 40)


def _target_coords(resource_classes_path: str) -> tuple[xr.DataArray, xr.DataArray]:
    ds = xr.load_dataset(resource_classes_path)
    try:
        y = ds["y"].astype(np.float32)
        x = ds["x"].astype(np.float32)
    except KeyError as exc:
        raise ValueError(
            "resource_classes.nc must expose 'y' and 'x' coordinates"
        ) from exc
    return y, x


def _open_land_cover(path: str) -> xr.DataArray:
    ds = xr.open_dataset(path, chunks="auto")
    data = ds["lccs_class"]
    if "time" in data.dims:
        data = data.isel(time=0, drop=True)
    nodata = data.attrs.get("_FillValue")
    if nodata is None:
        nodata = data.encoding.get("_FillValue")
    if nodata is None:
        nodata = 0
    return data.chunk(CHUNKS).where(data != nodata)


def _fraction(
    lc: xr.DataArray,
    classes: tuple[int, ...],
) -> xr.DataArray:
    valid = lc.notnull()
    hits = xr.where(lc.isin(classes), 1.0, 0.0)
    mask = xr.where(valid, hits, np.nan)
    frac = mask.coarsen(
        lat=COARSEN_FACTOR,
        lon=COARSEN_FACTOR,
        boundary="trim",
    ).mean(skipna=True)
    return frac.clip(0.0, 1.0).astype(np.float32)


def main() -> None:
    classes_path: str = snakemake.input.classes  # type: ignore[name-defined]
    land_cover_path: str = snakemake.input.land_cover  # type: ignore[name-defined]
    output_path: str = snakemake.output.fractions  # type: ignore[name-defined]

    y_coords, x_coords = _target_coords(classes_path)
    land_cover = _open_land_cover(land_cover_path)

    lat_size = land_cover.sizes.get("lat")
    lon_size = land_cover.sizes.get("lon")
    if lat_size is None or lon_size is None:
        raise ValueError("land cover input must have 'lat' and 'lon' dimensions")
    if lat_size % COARSEN_FACTOR != 0 or lon_size % COARSEN_FACTOR != 0:
        raise ValueError("land cover grid size must be divisible by the coarsen factor")

    forest_fraction = _fraction(land_cover, FOREST_CLASSES)
    cropland_fraction = _fraction(land_cover, CROPLAND_CLASSES)
    grassland_fraction = _fraction(land_cover, GRASSLAND_CLASSES)

    ds_out = xr.Dataset(
        {
            "forest_fraction": forest_fraction.rename(lat="y", lon="x"),
            "cropland_fraction": cropland_fraction.rename(lat="y", lon="x"),
            "grassland_fraction": grassland_fraction.rename(lat="y", lon="x"),
        }
    ).assign_coords(y=y_coords, x=x_coords)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    encoding = {
        "forest_fraction": {"zlib": True, "complevel": 4, "dtype": "float32"},
        "cropland_fraction": {"zlib": True, "complevel": 4, "dtype": "float32"},
        "grassland_fraction": {"zlib": True, "complevel": 4, "dtype": "float32"},
    }
    ds_out.to_netcdf(output_path, encoding=encoding)


if __name__ == "__main__":
    main()
