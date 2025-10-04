#! /usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Plot cropland fraction maps with resource class detail at pixel resolution.

Inputs (via snakemake):
- input.network: solved PyPSA network (NetCDF)
- input.regions: regions GeoJSON with a 'region' column
- input.land_area_by_class: CSV with columns [region, water_supply, resource_class, area_ha]
- input.resource_classes: NetCDF with variables 'resource_class' and 'region_id'

Outputs:
- PDF map under results/{name}/plots/
- CSV summary alongside the PDF (same stem, `_by_region_class.csv`)

Optional parameters (snakemake.params):
- water_supply: "i"/"irrigated", "r"/"rainfed", or omitted for all.
- title: figure title; defaults derived from `water_supply`.
- colorbar_label: colorbar caption; defaults derived from `water_supply`.
- cmap: Matplotlib colormap name (default: "YlOrRd").
- csv_prefix: prefix for CSV column names (default determined from `water_supply`).

Notes:
- Cropland use is computed from actual land flows supplied by land generators
  (carrier 'land'), i.e. n.generators_t.p for generators named like
  'land_{region}_class{k}_{ws}'. Water supply filtering (irrigated vs rainfed)
  happens via the optional parameter above. This reflects realized land use,
  not capacity.
- Total land area per (region, resource class) pair comes from
  land_area_by_class.csv, filtered to the requested water supply when
  specified, matching the model’s land availability basis.
- Each pixel inherits the cropland fraction of its (region, resource class)
  combination, so within-region spatial patterns remain visible.
"""

from pathlib import Path
import logging
import re

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
import geopandas as gpd
import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from affine import Affine
from rasterio.transform import array_bounds


_LAND_GEN_RE = re.compile(
    r"^land_(?P<region>.+?)_class(?P<resource_class>\d+)_?(?P<water_supply>[a-z]*)$"
)

logger = logging.getLogger(__name__)


def _used_cropland_area_by_region_class(
    n: pypsa.Network, water_supply: str | None
) -> pd.Series:
    """Return used cropland area by region and resource class.

    Extracts positive output from land generators (carrier 'land') at snapshot
    'now'. Generator names follow the pattern
    'land_{region}_class{resource_class}_{water_supply}'. Water supply letters
    (e.g. r/i) are ignored for aggregation.

    Returns area in hectares.
    """

    if n.generators.empty or n.generators_t.p.empty:
        return pd.Series(dtype=float)

    land_gen = n.generators[n.generators["carrier"] == "land"]
    if land_gen.empty:
        return pd.Series(dtype=float)

    names = land_gen.index
    if "now" in n.snapshots:
        p_now = n.generators_t.p.loc["now", names]
    else:
        p_now = n.generators_t.p.iloc[0][names]
    p_now = p_now.fillna(0.0)

    rows: list[tuple[str, int, float]] = []
    for name, value in p_now.items():
        match = _LAND_GEN_RE.match(str(name))
        if not match:
            continue
        region = match.group("region")
        resource_class = int(match.group("resource_class"))
        ws = (match.group("water_supply") or "").lower()
        if water_supply is not None and ws != water_supply:
            continue
        rows.append((region, resource_class, max(float(value), 0.0) * 1e6))

    if not rows:
        return pd.Series(dtype=float)

    df = pd.DataFrame(rows, columns=["region", "resource_class", "used_ha"])
    used = (
        df.groupby(["region", "resource_class"], sort=False)["used_ha"]
        .sum()
        .astype(float)
    )
    used.index = used.index.set_levels(
        used.index.levels[1].astype(int), level="resource_class"
    )
    return used


def _normalize_water_supply(value: object) -> str | None:
    if value is None:
        return None

    text = str(value).strip().lower()
    if not text or text in {"all", "total", "both", "any"}:
        return None

    mapping = {
        "i": "i",
        "irr": "i",
        "irrigated": "i",
        "r": "r",
        "rainfed": "r",
    }

    try:
        return mapping[text]
    except KeyError as exc:
        raise ValueError(
            "Unsupported water_supply parameter (expected irrigated/rainfed)"
        ) from exc


def main() -> None:
    n = pypsa.Network(snakemake.input.network)  # type: ignore[name-defined]
    regions_path: str = snakemake.input.regions  # type: ignore[name-defined]
    land_area_csv: str = snakemake.input.land_area_by_class  # type: ignore[name-defined]
    classes_path: str = snakemake.input.resource_classes  # type: ignore[name-defined]
    out_pdf = Path(snakemake.output.pdf)  # type: ignore[name-defined]

    params = getattr(snakemake, "params", None)  # type: ignore[name-defined]
    water_param = getattr(params, "water_supply", None) if params is not None else None
    water_supply = _normalize_water_supply(water_param)

    base_prefix = {
        None: "cropland",
        "i": "irrigated_cropland",
        "r": "rainfed_cropland",
    }[water_supply]

    csv_prefix = getattr(params, "csv_prefix", None) if params is not None else None
    if not csv_prefix:
        csv_prefix = base_prefix

    default_titles = {
        None: "Cropland Fraction by Region and Resource Class",
        "i": "Irrigated Cropland Fraction by Region and Resource Class",
        "r": "Rainfed Cropland Fraction by Region and Resource Class",
    }
    default_colorbars = {
        None: "Cropland / total model land area",
        "i": "Irrigated cropland / irrigable land area",
        "r": "Rainfed cropland / rainfed land area",
    }

    title = (
        getattr(params, "title", None) if params is not None else None
    ) or default_titles[water_supply]
    colorbar_label = (
        getattr(params, "colorbar_label", None) if params is not None else None
    ) or default_colorbars[water_supply]
    cmap_name = (
        getattr(params, "cmap", None) if params is not None else None
    ) or "YlOrRd"

    gdf = gpd.read_file(regions_path)
    if "region" not in gdf.columns:
        raise ValueError("regions input must contain a 'region' column")
    gdf = gdf.reset_index(drop=True).to_crs(4326)
    region_name_to_id = {region: idx for idx, region in enumerate(gdf["region"])}
    gdf = gdf.set_index("region", drop=False)

    used_ha = _used_cropland_area_by_region_class(n, water_supply)

    df_land = pd.read_csv(land_area_csv)
    required_cols = {"region", "resource_class", "area_ha"}
    if not required_cols.issubset(df_land.columns):
        raise ValueError("land_area_by_class.csv must contain required columns")

    df_land = df_land.dropna(subset=list(required_cols))
    df_land["resource_class"] = df_land["resource_class"].astype(int)
    if "water_supply" in df_land.columns:
        df_land["water_supply"] = (
            df_land["water_supply"].astype(str).str.strip().str.lower()
        )
        if water_supply is not None:
            df_land = df_land[df_land["water_supply"] == water_supply]

    total_ha = (
        df_land.groupby(["region", "resource_class"])["area_ha"].sum().astype(float)
    )

    classes = sorted(
        set(total_ha.index.get_level_values("resource_class"))
        | set(used_ha.index.get_level_values("resource_class"))
    )
    if not classes:
        raise ValueError("No resource classes found")

    full_index = pd.MultiIndex.from_product(
        [gdf.index, classes], names=["region", "resource_class"]
    )
    used_ha = used_ha.reindex(full_index).fillna(0.0)
    total_ha = total_ha.reindex(full_index).fillna(0.0)

    with np.errstate(divide="ignore", invalid="ignore"):
        fractions = (
            (used_ha / total_ha).replace([np.inf, -np.inf], np.nan).clip(0.0, 1.0)
        )

    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    classes_ds = xr.open_dataset(classes_path, mode="r")
    if "resource_class" not in classes_ds or "region_id" not in classes_ds:
        raise ValueError("resource_classes input must contain required variables")

    class_grid = classes_ds["resource_class"].values.astype(np.int16)
    region_grid = classes_ds["region_id"].values.astype(np.int32)
    transform = Affine(*classes_ds.attrs["transform"])
    height, width = class_grid.shape
    lon_min, lat_min, lon_max, lat_max = array_bounds(height, width, transform)
    extent = (lon_min, lon_max, lat_min, lat_max)  # Fixed orientation!
    classes_ds.close()

    fraction_grid = np.full(class_grid.shape, np.nan, dtype=float)
    for (region, cls), frac in fractions.items():
        ridx = region_name_to_id.get(region)
        if ridx is not None:
            mask = (region_grid == ridx) & (class_grid == cls)
            fraction_grid[mask] = frac

    vmax = (
        max(0.5, min(1.0, np.nanmax(fraction_grid)))
        if not np.all(np.isnan(fraction_grid))
        else 0.5
    )

    cmap = plt.colormaps[cmap_name]
    fig, ax = plt.subplots(
        figsize=(13, 6.5), dpi=150, subplot_kw={"projection": ccrs.EqualEarth()}
    )
    ax.set_facecolor("#f7f9fb")
    ax.set_global()

    plate = ccrs.PlateCarree()
    im = ax.imshow(
        fraction_grid,
        origin="upper",
        extent=extent,
        transform=plate,
        cmap=cmap,
        vmin=0,
        vmax=vmax,
        alpha=0.8,
        zorder=1,
    )

    ax.add_geometries(
        gdf.geometry,
        crs=plate,
        facecolor="none",
        edgecolor="#666666",
        linewidth=0.3,
        zorder=2,
    )

    gl = ax.gridlines(
        draw_labels=True,
        crs=plate,
        linewidth=0.35,
        color="#888888",
        alpha=0.45,
        linestyle="--",
    )
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
    gl.ylocator = mticker.FixedLocator(np.arange(-60, 61, 15))
    gl.xformatter = LongitudeFormatter(number_format=".0f")
    gl.yformatter = LatitudeFormatter(number_format=".0f")
    gl.xlabel_style = {"size": 8, "color": "#555555"}
    gl.ylabel_style = {"size": 8, "color": "#555555"}
    gl.top_labels = gl.right_labels = False

    cbar = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.02)
    cbar.set_label(colorbar_label)

    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    plt.close(fig)

    data = pd.DataFrame(
        {
            f"{csv_prefix}_used_ha": used_ha,
            (
                "land_total_ha"
                if water_supply is None
                else (
                    "irrigable_land_total_ha"
                    if water_supply == "i"
                    else "rainfed_land_total_ha"
                )
            ): total_ha,
            f"{csv_prefix}_fraction": fractions,
        }
    )
    csv_out = out_pdf.with_suffix("").parent / f"{out_pdf.stem}_by_region_class.csv"
    data.reset_index().to_csv(csv_out, index=False)
    logger.info("Saved cropland fraction map to %s and CSV to %s", out_pdf, csv_out)


if __name__ == "__main__":
    main()
