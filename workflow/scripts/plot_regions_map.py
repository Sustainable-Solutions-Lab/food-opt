#! /usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def plot_regions_equal_earth(regions_path: str, output_path: str) -> None:
    logger.info("Loading regions from %s", regions_path)
    gdf = gpd.read_file(regions_path)

    if gdf.crs is None:
        logger.warning("Input CRS missing; assuming EPSG:4326 (WGS84)")
        gdf = gdf.set_crs(4326, allow_override=True)

    # Equal Earth projection via PROJ
    ee_crs = "+proj=eqearth"
    gdf_ee = gdf.to_crs(ee_crs)

    # Prepare output directory
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6), dpi=150)
    gdf_ee.plot(ax=ax, linewidth=0.3, edgecolor="#444444", facecolor="#cfd8dc")
    ax.set_axis_off()
    ax.set_title("Regions (Equal Earth)", fontsize=12)
    plt.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)

    logger.info("Saved regions map to %s", output_path)


if __name__ == "__main__":
    plot_regions_equal_earth(snakemake.input.regions, snakemake.output.pdf)
