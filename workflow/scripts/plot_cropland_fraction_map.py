#! /usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Plot cropland fraction map with resource class detail at pixel resolution.

Inputs (via snakemake):
- input.network: solved PyPSA network (NetCDF)
- input.regions: regions GeoJSON with a 'region' column
- input.land_area_by_class: CSV with columns [region, water_supply, resource_class, area_ha]
- input.resource_classes: NetCDF with variables 'resource_class' and 'region_id'

Output:
- results/{name}/plots/cropland_fraction_map.pdf (single map at pixel resolution)

Notes:
- Cropland use is computed from actual land flows supplied by land generators
  (carrier 'land'), i.e. n.generators_t.p for generators named like
  'land_{region}_class{k}_{ws}', aggregated per resource class across water
  supply options. This reflects realized land use, not capacity.
- Total land area per (region, resource class) pair is the sum of area_ha over
  water supplies from land_area_by_class.csv, matching the model’s land
  availability basis.
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


def _used_cropland_area_by_region_class(n: pypsa.Network) -> pd.Series:
    """Return used cropland area by region and resource class.

    Extracts positive output from land generators (carrier 'land') at snapshot
    'now'. Generator names follow the pattern
    'land_{region}_class{resource_class}_{water_supply}'. Water supply letters
    (e.g. r/i) are ignored for aggregation.
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
        rows.append((region, resource_class, max(float(value), 0.0)))

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


def main() -> None:
    n = pypsa.Network(snakemake.input.network)  # type: ignore[name-defined]
    regions_path: str = snakemake.input.regions  # type: ignore[name-defined]
    land_area_csv: str = snakemake.input.land_area_by_class  # type: ignore[name-defined]
    classes_path: str = snakemake.input.resource_classes  # type: ignore[name-defined]
    out_pdf = Path(snakemake.output.pdf)  # type: ignore[name-defined]

    gdf = gpd.read_file(regions_path)
    if "region" not in gdf.columns:
        raise ValueError("regions input must contain a 'region' column")
    gdf = gdf.reset_index(drop=True).to_crs(4326)
    region_name_to_id = {region: idx for idx, region in enumerate(gdf["region"])}
    gdf = gdf.set_index("region", drop=False)

    used_ha = _used_cropland_area_by_region_class(n)

    df_land = pd.read_csv(land_area_csv)
    required_cols = {"region", "resource_class", "area_ha"}
    if not required_cols.issubset(df_land.columns):
        raise ValueError("land_area_by_class.csv must contain required columns")

    df_land = df_land.dropna(subset=list(required_cols))
    df_land["resource_class"] = df_land["resource_class"].astype(int)
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

    classes_ds = xr.load_dataset(classes_path)
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

    cmap = plt.get_cmap("YlOrRd")
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
    cbar.set_label("Cropland / total model land area")

    ax.set_title("Cropland Fraction by Region and Resource Class")
    plt.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    plt.close(fig)

    data = pd.DataFrame(
        {
            "cropland_used_ha": used_ha,
            "land_total_ha": total_ha,
            "cropland_fraction": fractions,
        }
    )
    csv_out = out_pdf.with_suffix("").parent / f"{out_pdf.stem}_by_region_class.csv"
    data.reset_index().to_csv(csv_out, index=False)
    logger.info("Saved cropland fraction map to %s and CSV to %s", out_pdf, csv_out)


if __name__ == "__main__":
    main()
