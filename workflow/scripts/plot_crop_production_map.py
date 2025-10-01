#! /usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
from pathlib import Path
from typing import Dict, Iterable, Tuple

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
import geopandas as gpd
import matplotlib

matplotlib.use("pdf")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import pypsa

logger = logging.getLogger(__name__)


def _region_from_bus0(bus0: str) -> str:
    parts = bus0.split("_")
    return parts[1] if len(parts) >= 2 else "unknown"


def _dict_to_df(data: Dict[Tuple[str, str], float]) -> pd.DataFrame:
    if not data:
        return pd.DataFrame()
    series = pd.Series(data).sort_index()
    df = series.unstack(fill_value=0.0).sort_index(axis=0).sort_index(axis=1)
    df.index.name = "region"
    return df


def _aggregate_production_by_region(n: pypsa.Network, snapshot: str) -> pd.DataFrame:
    data: Dict[Tuple[str, str], float] = {}

    def add(region: str, crop: str, value: float) -> None:
        if not np.isfinite(value) or value <= 0:
            return
        key = (region, crop)
        data[key] = data.get(key, 0.0) + float(value)

    produce_links = [name for name in n.links.index if str(name).startswith("produce_")]
    if produce_links:
        p1 = n.links_t.p1.loc[snapshot, produce_links]
        bus0 = n.links.loc[produce_links, "bus0"]
        for name, value in p1.items():
            crop = str(name).split("_")[1] if "_" in str(name) else "unknown"
            region = _region_from_bus0(str(bus0.loc[name]))
            add(region, crop, abs(float(value)))

    graze_links = [name for name in n.links.index if str(name).startswith("graze_")]
    if graze_links:
        p1 = n.links_t.p1.loc[snapshot, graze_links]
        bus0 = n.links.loc[graze_links, "bus0"]
        for name, value in p1.items():
            region = _region_from_bus0(str(bus0.loc[name]))
            add(region, "grassland", abs(float(value)))

    df = _dict_to_df(data)
    if "grassland" not in df.columns:
        df["grassland"] = 0.0
    return df


def _aggregate_land_use_by_region(n: pypsa.Network, snapshot: str) -> pd.DataFrame:
    """Aggregate land use by region and crop.

    Returns land area in hectares.
    """
    data: Dict[Tuple[str, str], float] = {}

    def add(region: str, crop: str, value: float) -> None:
        if not np.isfinite(value) or value <= 0:
            return
        key = (region, crop)
        data[key] = data.get(key, 0.0) + float(value)

    produce_links = [name for name in n.links.index if str(name).startswith("produce_")]
    if produce_links:
        p0 = n.links_t.p0.loc[snapshot, produce_links]
        bus0 = n.links.loc[produce_links, "bus0"]
        for name, value in p0.items():
            crop = str(name).split("_")[1] if "_" in str(name) else "unknown"
            region = _region_from_bus0(str(bus0.loc[name]))
            add(region, crop, max(float(value), 0.0) * 1e6)

    graze_links = [name for name in n.links.index if str(name).startswith("graze_")]
    if graze_links:
        p0 = n.links_t.p0.loc[snapshot, graze_links]
        bus0 = n.links.loc[graze_links, "bus0"]
        for name, value in p0.items():
            region = _region_from_bus0(str(bus0.loc[name]))
            add(region, "grassland", max(float(value), 0.0) * 1e6)

    df = _dict_to_df(data)
    if "grassland" not in df.columns:
        df["grassland"] = 0.0
    return df


def _setup_regions(regions_path: str) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    gdf = gpd.read_file(regions_path)
    if gdf.crs is None:
        logger.warning("Regions GeoDataFrame missing CRS; assuming EPSG:4326")
        gdf = gdf.set_crs(4326, allow_override=True)
    else:
        gdf = gdf.to_crs(4326)

    if "region" not in gdf.columns:
        raise ValueError("Regions GeoDataFrame must contain a 'region' column")

    gdf = gdf.set_index("region", drop=False)
    gdf_eq = gdf.to_crs("+proj=eqearth")
    gdf_eq = gdf_eq.set_index("region", drop=False)
    return gdf, gdf_eq


def _plot_pie_map(
    by_region: pd.DataFrame,
    gdf: gpd.GeoDataFrame,
    gdf_eq: gpd.GeoDataFrame,
    colors: Dict[str, str],
    output_path: str,
    title: str,
    legend_title: str,
    pie_scale_title: str,
    pie_unit: str,
    min_total: float,
) -> None:
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(
        figsize=(13, 6.5),
        dpi=150,
        subplot_kw={"projection": ccrs.EqualEarth()},
    )
    ax.set_facecolor("#f7f9fb")
    ax.set_global()
    plate = ccrs.PlateCarree()

    ax.add_geometries(
        gdf.geometry,
        crs=plate,
        facecolor="#e6eef2",
        edgecolor="#666666",
        linewidth=0.3,
        zorder=1,
    )

    model_regions = gdf.index
    present_regions = by_region.index if not by_region.empty else pd.Index([])
    missing_regions = model_regions.difference(present_regions)

    if len(missing_regions) > 0:
        ax.add_geometries(
            gdf.loc[missing_regions].geometry,
            crs=plate,
            facecolor="#f0f0f0",
            edgecolor="#666666",
            linewidth=0.3,
            hatch="..",
            zorder=1.5,
        )

    filtered = by_region.copy()
    if not filtered.empty:
        filtered = filtered.reindex(model_regions.intersection(filtered.index)).fillna(
            0.0
        )
        crop_totals = filtered.sum(axis=0).sort_values(ascending=False)
        crop_totals = crop_totals[crop_totals >= min_total]
        if crop_totals.empty:
            filtered = filtered.iloc[:, 0:0]
        else:
            filtered = filtered.loc[:, crop_totals.index]

    if not filtered.empty:
        crops = list(filtered.columns)
        color_list = [colors[c] for c in crops]
        totals = filtered.sum(axis=1)
        if totals.max() > 0:
            xmin, ymin, xmax, ymax = gdf_eq.total_bounds
            width = xmax - xmin
            height = ymax - ymin
            r_max = 0.024 * max(width, height)
            radii = (np.sqrt(totals / totals.max()) * r_max).fillna(0.0)

            centroids = gdf_eq.geometry.representative_point()
            for region in filtered.index:
                point = centroids.loc[region]
                x, y = point.x, point.y
                values = filtered.loc[region].values.tolist()
                _draw_pie(ax, x, y, values, color_list, float(radii.get(region, 0.0)))

            handles = [mpatches.Patch(facecolor=colors[c], label=c) for c in crops]
            legend1 = ax.legend(
                handles=handles,
                title=legend_title,
                loc="lower left",
                bbox_to_anchor=(0.15, 0.03),
                fontsize=8,
                title_fontsize=9,
                frameon=True,
                borderpad=0.8,
                labelspacing=0.6,
                handletextpad=0.6,
            )
            legend1._legend_box.align = "left"
            ax.add_artist(legend1)

            ref_fracs = np.array([0.25, 0.5, 1.0])
            ref_vals = np.unique(totals.max() * ref_fracs)
            handle_scale = 900.0
            size_handles = [
                ax.scatter(
                    [],
                    [],
                    s=float((val / totals.max()) * handle_scale),
                    facecolors="#d0d7de",
                    edgecolors="#555555",
                    linewidths=0.6,
                    alpha=0.7,
                )
                for val in ref_vals
            ]
            size_labels = [f"{val:,.0f} {pie_unit}" for val in ref_vals]
            legend2 = ax.legend(
                size_handles,
                size_labels,
                title=pie_scale_title,
                loc="lower left",
                bbox_to_anchor=(0.6, 0.03),
                fontsize=8,
                title_fontsize=9,
                frameon=True,
                scatterpoints=1,
                handlelength=1.5,
                borderpad=0.8,
                labelspacing=2,
            )
            legend2._legend_box.align = "left"
            ax.add_artist(legend2)

    if len(missing_regions) > 0:
        hatch_handle = mpatches.Patch(
            facecolor="#f0f0f0",
            edgecolor="#666666",
            hatch="..",
            label="No activity",
        )
        ax.legend(
            handles=[hatch_handle],
            loc="lower right",
            bbox_to_anchor=(0.99, 0.02),
            fontsize=8,
            frameon=True,
        )

    for name, spine in ax.spines.items():
        if name == "geo":
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_edgecolor("#555555")
            spine.set_alpha(0.7)
        else:
            spine.set_visible(False)

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
    gl.top_labels = False
    gl.right_labels = False

    ax.set_xlabel("Longitude", fontsize=8, color="#555555")
    ax.set_ylabel("Latitude", fontsize=8, color="#555555")
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)

    if not by_region.empty:
        csv_path = out.with_suffix("")
        csv_out = csv_path.parent / f"{csv_path.name}_by_region.csv"
        by_region.sort_index(axis=1).to_csv(csv_out, index=True)
        logger.info("Saved regional totals to %s", csv_out)

    logger.info("Saved map to %s", out)


def _draw_pie(
    ax: plt.Axes,
    x: float,
    y: float,
    sizes: Iterable[float],
    colors: Iterable[str],
    radius: float,
) -> None:
    total = float(sum(sizes))
    if total <= 0 or radius <= 0:
        return
    sizes = list(sizes)
    colors = list(colors)
    angles = np.cumsum([0.0] + [s / total * 360.0 for s in sizes])
    for i, size in enumerate(sizes):
        if size <= 0:
            continue
        wedge = mpatches.Wedge(
            center=(x, y),
            r=radius,
            theta1=angles[i],
            theta2=angles[i + 1],
            facecolor=colors[i],
            edgecolor="white",
            linewidth=0.4,
            alpha=0.85,
            zorder=10,
        )
        ax.add_patch(wedge)
    circ = mpatches.Circle(
        (x, y),
        radius=radius,
        facecolor="none",
        edgecolor="#444444",
        linewidth=0.3,
        alpha=0.7,
        zorder=11,
    )
    ax.add_patch(circ)


def main() -> None:
    n = pypsa.Network(snakemake.input.network)  # type: ignore[name-defined]
    regions_path: str = snakemake.input.regions  # type: ignore[name-defined]
    prod_pdf: str = snakemake.output.production_pdf  # type: ignore[name-defined]
    land_pdf: str = snakemake.output.land_pdf  # type: ignore[name-defined]

    snapshot = "now" if "now" in n.snapshots else n.snapshots[0]

    gdf, gdf_eq = _setup_regions(regions_path)

    production = _aggregate_production_by_region(n, snapshot)
    land_use = _aggregate_land_use_by_region(n, snapshot)

    all_regions = gdf.index
    all_columns = sorted(set(production.columns) | set(land_use.columns))
    if production.empty:
        production = pd.DataFrame(
            index=pd.Index([], name="region"), columns=all_columns
        )
    else:
        production = production.reindex(
            all_regions.intersection(production.index)
        ).fillna(0.0)
        production = production.reindex(columns=all_columns, fill_value=0.0)

    if land_use.empty:
        land_use = pd.DataFrame(index=pd.Index([], name="region"), columns=all_columns)
    else:
        land_use = land_use.reindex(all_regions.intersection(land_use.index)).fillna(
            0.0
        )
        land_use = land_use.reindex(columns=all_columns, fill_value=0.0)

    # Define logical color scheme by crop category
    crop_colors = {
        # Grains - shades of gold/tan
        "wheat": "#DAA520",
        "rice": "#F4A460",
        "maize": "#FFD700",
        "barley": "#D2B48C",
        "sorghum": "#CD853F",
        "millet": "#DEB887",
        "oats": "#F5DEB3",
        "rye": "#BC8F8F",
        # Legumes - shades of brown/maroon
        "soybean": "#8B4513",
        "groundnut": "#A0522D",
        "chickpea": "#D2691E",
        "cowpea": "#8B5A2B",
        "pigeonpea": "#9B6B47",
        "bean": "#6B4423",
        "lentil": "#704214",
        "pea": "#8B7355",
        # Oil crops - shades of yellow/orange
        "rapeseed": "#FFB347",
        "sunflower": "#FFA500",
        "sesame": "#FF8C00",
        "palm": "#FF7F50",
        # Vegetables - shades of green
        "cassava": "#6B8E23",
        "potato": "#8FBC8F",
        "sweetpotato": "#9ACD32",
        "yam": "#556B2F",
        "tomato": "#DC143C",
        # Sugar crops - purple/pink
        "sugarcane": "#BA55D3",
        "sugarbeet": "#DA70D6",
        # Fruits/other
        "banana": "#FFE135",
        "plantain": "#F0E68C",
        "coconut": "#8B7D6B",
        # Grassland - light green
        "grassland": "#90EE90",
    }

    # Fallback colors for any crops not in the predefined list
    fallback_cmap = plt.colormaps["Set3"]
    colors = {}
    fallback_idx = 0
    for crop in all_columns:
        if crop in crop_colors:
            colors[crop] = crop_colors[crop]
        else:
            colors[crop] = mcolors.to_hex(fallback_cmap(fallback_idx % 12))
            fallback_idx += 1

    _plot_pie_map(
        production,
        gdf,
        gdf_eq,
        colors,
        prod_pdf,
        title="Crop and Grassland Output by Region",
        legend_title="Crops / grassland",
        pie_scale_title="Pie size ∝ total production",
        pie_unit="t",
        min_total=10_000.0,
    )

    _plot_pie_map(
        land_use,
        gdf,
        gdf_eq,
        colors,
        land_pdf,
        title="Land Use by Crop and Grassland",
        legend_title="Crops / grassland",
        pie_scale_title="Pie size ∝ total land area",
        pie_unit="ha",
        min_total=1_000.0,
    )


if __name__ == "__main__":
    main()
