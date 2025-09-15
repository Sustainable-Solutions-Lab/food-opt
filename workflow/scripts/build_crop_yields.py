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
import pandas as pd


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
    classes_nc: str = snakemake.input.classes  # type: ignore[name-defined]
    yield_path: str = snakemake.input.yield_raster  # type: ignore[name-defined]
    suit_path: str = snakemake.input.suitability_raster  # type: ignore[name-defined]
    regions_path: str = snakemake.input.regions  # type: ignore[name-defined]
    crop_code: str = snakemake.wildcards.crop  # type: ignore[name-defined]
    conv_csv: str | None = getattr(snakemake.input, "yield_unit_conversions", None)  # type: ignore[attr-defined]

    # Load classes
    ds = xr.load_dataset(classes_nc)
    class_labels = ds["resource_class"].values.astype(np.int16)

    # Load rasters
    y_raw, y_src = read_raster_float(yield_path)
    # Determine conversion factor from CSV (default 0.001 t per kg)
    factor = 0.001
    if conv_csv:
        try:
            df_conv = pd.read_csv(conv_csv, comment="#")
            df_conv = df_conv.set_index("code")
            if crop_code in df_conv.index and pd.notna(
                df_conv.at[crop_code, "factor_to_t_per_ha"]
            ):
                factor = float(df_conv.at[crop_code, "factor_to_t_per_ha"])
        except Exception:
            # Fall back to default if the table cannot be read
            factor = 0.001
    y_tpha = y_raw * factor
    s_raw, _ = read_raster_float(suit_path)
    s_frac = scale_fraction(s_raw)

    height, width = y_tpha.shape
    transform = y_src.transform
    crs = y_src.crs
    crs_wkt = crs.to_wkt() if crs else None
    xmin, ymin, xmax, ymax = raster_bounds(transform, width, height)
    cell_area_ha = calculate_all_cell_areas(y_src)

    area_ha = s_frac * cell_area_ha

    # Regions
    regions_gdf = gpd.read_file(regions_path)
    if regions_gdf.crs and crs and regions_gdf.crs != crs:
        regions_gdf = regions_gdf.to_crs(crs)
    regions_for_extract = regions_gdf.reset_index()

    # Aggregate mean yield and sum area per class
    import pandas as pd

    out = []
    n_classes = (
        int(np.nanmax(class_labels)) + 1 if np.isfinite(class_labels).any() else 0
    )
    for cls in range(n_classes):
        mask = class_labels == cls
        if not np.any(mask):
            continue
        y = np.where(mask, y_tpha, np.nan)
        a = np.where(mask, area_ha, np.nan)

        y_src_np = NumPyRasterSource(
            y,
            xmin=xmin,
            ymin=ymin,
            xmax=xmax,
            ymax=ymax,
            nodata=np.nan,
            srs_wkt=crs_wkt,
        )
        a_src_np = NumPyRasterSource(
            a,
            xmin=xmin,
            ymin=ymin,
            xmax=xmax,
            ymax=ymax,
            nodata=np.nan,
            srs_wkt=crs_wkt,
        )

        y_stats = exact_extract(
            y_src_np,
            regions_for_extract,
            ["mean"],
            include_cols=["region"],
            output="pandas",
        )
        a_stats = exact_extract(
            a_src_np,
            regions_for_extract,
            ["sum"],
            include_cols=["region"],
            output="pandas",
        )
        if y_stats.empty or a_stats.empty:
            continue
        merged = y_stats.rename(columns={"mean": "yield"}).merge(
            a_stats.rename(columns={"sum": "suitable_area"}), on="region", how="inner"
        )
        merged["resource_class"] = cls
        out.append(merged)

    if out:
        df = (
            pd.concat(out, ignore_index=True)
            .set_index(["region", "resource_class"])
            .sort_index()
        )
    else:
        df = pd.DataFrame(
            columns=["region", "resource_class", "yield", "suitable_area"]
        ).set_index(["region", "resource_class"])  # type: ignore[name-defined]

    Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)  # type: ignore[name-defined]
    df.to_csv(snakemake.output[0])  # type: ignore[name-defined]
