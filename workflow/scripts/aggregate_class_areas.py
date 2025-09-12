"""
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: GPL-3.0-or-later
"""

from pathlib import Path
import numpy as np
import geopandas as gpd
import rasterio
from exactextract import exact_extract
from exactextract.raster import NumPyRasterSource
import xarray as xr

from raster_utils import calculate_all_cell_areas, scale_fraction


def read_raster_float(path: str):
    src = rasterio.open(path)
    arr = src.read(1).astype(float)
    if src.nodata is not None:
        arr = np.where(arr == src.nodata, np.nan, arr)
    return arr, src


def raster_bounds(transform, width: int, height: int):
    xmin = transform.c
    ymax = transform.f
    xmax = xmin + width * transform.a
    ymin = ymax + height * transform.e
    return xmin, ymin, xmax, ymax


if __name__ == "__main__":
    # Inputs
    regions_path: str = snakemake.input.regions  # type: ignore[name-defined]
    classes_nc: str = snakemake.input.classes  # type: ignore[name-defined]
    # Suitability/area inputs as lists of file paths
    sr_files: list[str] = list(snakemake.input.sr)  # type: ignore[attr-defined]
    si_files: list[str] = list(snakemake.input.si)  # type: ignore[attr-defined]
    irrigated_share_path: str | None = getattr(snakemake.input, "irrigated_share", None)  # type: ignore[attr-defined]

    land_limit_mode: str = snakemake.params.land_limit_dataset  # type: ignore[name-defined]

    # Load classes
    ds = xr.load_dataset(classes_nc)
    classes = ds["resource_class"].values.astype(np.int16)

    # Reference grid parameters from a suitability raster (rainfed)
    # Use first rainfed suitability file as reference
    if not sr_files:
        raise ValueError("No rainfed suitability files provided")
    sr0, src0 = read_raster_float(sr_files[0])
    height, width = sr0.shape
    transform = src0.transform
    crs = src0.crs
    xmin, ymin, xmax, ymax = raster_bounds(transform, width, height)
    crs_wkt = crs.to_wkt() if crs else None

    # Regions
    regions_gdf = gpd.read_file(regions_path)
    if regions_gdf.crs and crs and regions_gdf.crs != crs:
        regions_gdf = regions_gdf.to_crs(crs)
    regions_for_extract = regions_gdf.reset_index()

    # Cell areas
    cell_area_ha = calculate_all_cell_areas(src0)

    # Build max suitability per pixel across crops for each ws
    def max_suitability(files: list[str]) -> np.ndarray:
        it = iter(files)
        base = scale_fraction(read_raster_float(next(it))[0])
        for p in it:
            base = np.fmax(base, scale_fraction(read_raster_float(p)[0]))
        return base

    sr_max = max_suitability(sr_files) if sr_files else np.zeros((height, width))
    if land_limit_mode == "suitability":
        si_max = max_suitability(si_files) if si_files else np.zeros((height, width))
        area_i = si_max * cell_area_ha
    else:
        if not irrigated_share_path:
            raise ValueError(
                "irrigated_share input required when land_limit_dataset='irrigated'"
            )
        irr_share = scale_fraction(read_raster_float(irrigated_share_path)[0])
        area_i = irr_share * cell_area_ha
    area_r = sr_max * cell_area_ha

    # Aggregate sums per region and resource class via exactextract, for each ws
    import pandas as pd

    def aggregate_area(area: np.ndarray, ws: str) -> pd.DataFrame:
        out = []
        n_classes = int(np.nanmax(classes)) + 1 if np.isfinite(classes).any() else 0
        for cls in range(n_classes):
            mask = classes == cls
            if not np.any(mask):
                continue
            arr = np.where(mask, area, np.nan)
            a_src = NumPyRasterSource(
                arr,
                xmin=xmin,
                ymin=ymin,
                xmax=xmax,
                ymax=ymax,
                nodata=np.nan,
                srs_wkt=crs_wkt,
            )
            a_stats = exact_extract(
                a_src,
                regions_for_extract,
                ["sum"],
                include_cols=["region"],
                output="pandas",
            )
            if a_stats.empty:
                continue
            a_stats = a_stats.rename(columns={"sum": "area_ha"})
            a_stats["resource_class"] = cls
            a_stats["water_supply"] = ws
            out.append(a_stats)
        if not out:
            return pd.DataFrame(
                columns=["region", "resource_class", "water_supply", "area_ha"]
            )
        return pd.concat(out, ignore_index=True)

    df_r = aggregate_area(area_r, "r")
    df_i = aggregate_area(area_i, "i")
    out_df = pd.concat([df_r, df_i], ignore_index=True)
    out_df = out_df.set_index(["region", "water_supply", "resource_class"]).sort_index()

    out_path = Path(snakemake.output[0])  # type: ignore[name-defined]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path)
