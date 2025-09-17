#! /usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
from pathlib import Path
from typing import Dict, List, Tuple

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
import geopandas as gpd
import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import pypsa

logger = logging.getLogger(__name__)


def _extract_crop_by_region(n: pypsa.Network) -> pd.DataFrame:
    """Return dataframe with index=region, columns=crop, values=tonnes.

    Compatible with link names like:
    - produce_{crop}_{region}_class{resource_class}  (legacy)
    - produce_{crop}_{irrigated|rainfed}_{region}_class{resource_class}  (current)
    Region is derived robustly from bus0: land_{region}_class{k}_{i|r}.
    """
    data: Dict[Tuple[str, str], float] = {}
    links = [link for link in n.links.index if str(link).startswith("produce_")]
    if not links:
        return pd.DataFrame()

    for link in links:
        parts = str(link).split("_")
        crop = parts[1] if len(parts) >= 2 else "unknown"
        # Derive region from bus0 to avoid ambiguity with ws tokens
        try:
            bus0 = n.links.at[link, "bus0"]
            bparts = str(bus0).split("_")
            # Expect: land_{region}_class{k}_{ws}
            region = bparts[1] if len(bparts) >= 2 else "unknown"
        except Exception:
            # Fallback to legacy parsing: third token as region
            region = parts[2] if len(parts) >= 3 else "unknown"

        flow = float(n.links_t.p1.loc["now", link])  # t at crop bus
        data[(region, crop)] = data.get((region, crop), 0.0) + abs(flow)

    if not data:
        return pd.DataFrame()

    s = pd.Series(data)
    df = s.unstack(fill_value=0.0)
    df.index.name = "region"
    return df.sort_index(axis=0).sort_index(axis=1)


def _draw_pie(
    ax: plt.Axes,
    x: float,
    y: float,
    sizes: List[float],
    colors: List[str],
    radius: float,
) -> None:
    total = float(sum(sizes))
    if total <= 0 or radius <= 0:
        return
    angles = np.cumsum([0.0] + [s / total * 360.0 for s in sizes])
    for i in range(len(sizes)):
        if sizes[i] <= 0:
            continue
        wedge = mpatches.Wedge(
            center=(x, y),
            r=radius,
            theta1=angles[i],
            theta2=angles[i + 1],
            facecolor=colors[i],
            edgecolor="white",
            linewidth=0.4,
        )
        ax.add_patch(wedge)
    # subtle outline
    circ = mpatches.Circle(
        (x, y),
        radius=radius,
        facecolor="none",
        edgecolor="#444444",
        linewidth=0.3,
        alpha=0.7,
    )
    ax.add_patch(circ)


def plot_crop_production_pies(
    n: pypsa.Network, regions_path: str, output_path: str
) -> None:
    logger.info("Loading regions and network")
    gdf = gpd.read_file(regions_path)
    if gdf.crs is None:
        logger.warning("Input CRS missing; assuming EPSG:4326 (WGS84)")
        gdf = gdf.set_crs(4326, allow_override=True)
    else:
        gdf = gdf.to_crs(4326)

    if "region" not in gdf.columns:
        raise ValueError("Regions GeoDataFrame must contain a 'region' column")

    gdf = gdf.set_index("region", drop=False)
    gdf_eq = gdf.to_crs("+proj=eqearth")

    # Compute crop production by region (only regions that have links)
    by_region = _extract_crop_by_region(n)

    # Prepare indices
    gdf_eq = gdf_eq.set_index("region", drop=False)
    model_regions = gdf.index
    present_regions = by_region.index if not by_region.empty else pd.Index([])
    missing_regions = model_regions.difference(present_regions)
    # Filter by_region to modeled regions only
    if not by_region.empty:
        by_region = by_region.reindex(
            model_regions.intersection(present_regions)
        ).fillna(0.0)

        # Keep only crops with meaningful production and order by total output
        crop_totals = by_region.sum(axis=0).sort_values(ascending=False)
        crop_totals = crop_totals[crop_totals >= 10_000.0]
        if crop_totals.empty:
            by_region = by_region.iloc[:, 0:0]
        else:
            by_region = by_region.loc[:, crop_totals.index]

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

    # Hatch overlay for regions without production links
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

    # Draw pies for regions that have production links
    if not by_region.empty:
        crops = list(by_region.columns)
        cmap = plt.get_cmap("tab20")
        colors = {
            crop: matplotlib.colors.to_hex(cmap(i % 20)) for i, crop in enumerate(crops)
        }

        totals = by_region.sum(axis=1)
        if totals.max() > 0:
            xmin, ymin, xmax, ymax = gdf_eq.total_bounds
            width = xmax - xmin
            height = ymax - ymin
            r_max = 0.03 * max(width, height)
            radii = (np.sqrt(totals / totals.max()) * r_max).fillna(0.0)

            centroids = gdf_eq.geometry.representative_point()
            for region in by_region.index:
                point = centroids.loc[region]
                x, y = point.x, point.y
                values = by_region.loc[region].values.tolist()
                cols = [colors[c] for c in crops]
                _draw_pie(ax, x, y, values, cols, float(radii.get(region, 0.0)))

            # Legend: crops (colors)
            handles = [mpatches.Patch(facecolor=colors[c], label=c) for c in crops]
            legend1 = ax.legend(
                handles=handles,
                title="Crops",
                loc="lower left",
                bbox_to_anchor=(0.01, 0.02),
                fontsize=8,
                title_fontsize=9,
                frameon=True,
            )
            ax.add_artist(legend1)

            # Legend: pie size scale
            ref_vals = np.linspace(totals.max() / 5.0, totals.max(), 3)
            ref_r = np.sqrt(ref_vals / totals.max()) * r_max
            x0 = xmax - 0.06 * width
            y0 = ymax - 0.1 * height
            for i, (rv, rr) in enumerate(zip(ref_vals, ref_r)):
                circ = mpatches.Circle(
                    (x0, y0 - i * 1.5 * r_max),
                    radius=float(rr),
                    facecolor="#d0d7de",
                    edgecolor="#555555",
                    linewidth=0.5,
                    alpha=0.7,
                )
                ax.add_patch(circ)
                ax.text(
                    x0 + 1.2 * r_max,
                    y0 - i * 1.5 * r_max,
                    f"{rv:,.0f} t",
                    va="center",
                    fontsize=8,
                )
            ax.text(
                x0,
                y0 - 1.5 * r_max * len(ref_vals) + 0.01 * height,
                "Pie size ∝ total crop production",
                fontsize=8,
            )

    # Legend: hatched regions indicating no production links
    if len(missing_regions) > 0:
        hatch_handle = mpatches.Patch(
            facecolor="#f0f0f0",
            edgecolor="#666666",
            hatch="..",
            label="No production links",
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
    ax.set_title("Crop Production by Region")
    plt.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)

    # Also write CSV for reference
    if not by_region.empty:
        csv_path = out.with_suffix("")  # remove .pdf
        csv_out = csv_path.parent / f"{csv_path.name}_by_region.csv"
        by_region.to_csv(csv_out, index=True)
        logger.info("Saved crop production by region to %s", csv_out)

    logger.info("Saved crop production pie map to %s", output_path)


if __name__ == "__main__":
    # snakemake inputs
    n = pypsa.Network(snakemake.input.network)
    plot_crop_production_pies(n, snakemake.input.regions, snakemake.output.pdf)
