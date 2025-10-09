#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Generate global coverage map for introduction documentation.

Shows all modeled regions with country boundaries to illustrate
the global scope and spatial resolution of the model.
"""

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))
from doc_figures_config import apply_doc_style, COLORS, FIGURE_SIZES, save_doc_figure


def main(regions_path: str, output_path: str):
    """Generate global coverage map.

    Args:
        regions_path: Path to regions GeoJSON file
        output_path: Path for output SVG file
    """
    # Apply documentation styling
    apply_doc_style()

    # Load regions
    regions = gpd.read_file(regions_path)

    # Ensure CRS is WGS84
    if regions.crs is None:
        regions = regions.set_crs(4326, allow_override=True)
    else:
        regions = regions.to_crs(4326)

    # Create figure with EqualEarth projection
    fig, ax = plt.subplots(
        figsize=FIGURE_SIZES["map_wide"],
        subplot_kw={"projection": ccrs.EqualEarth()},
    )

    ax.set_global()
    ax.set_facecolor("#f7f9fb")

    # Plot regions using cartopy
    ax.add_geometries(
        regions.geometry,
        crs=ccrs.PlateCarree(),
        facecolor=COLORS["primary"],
        edgecolor="white",
        linewidth=0.3,
        alpha=0.7,
    )

    # Add coastlines for context (very subtle)
    ax.coastlines(linewidth=0.3, color="#888888", alpha=0.3)

    # Style the map frame
    for name, spine in ax.spines.items():
        if name == "geo":
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_edgecolor("#555555")
            spine.set_alpha(0.7)
        else:
            spine.set_visible(False)

    ax.set_title("Global Model Coverage and Regional Aggregation", fontsize=12, pad=10)

    # Add text annotation with region count
    n_regions = len(regions)
    ax.text(
        0.02,
        0.98,
        f"{n_regions} optimization regions",
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(
            boxstyle="round,pad=0.5", facecolor="white", alpha=0.8, edgecolor="none"
        ),
    )

    # Save figure
    save_doc_figure(fig, output_path, format="svg")
    plt.close(fig)


if __name__ == "__main__":
    # Snakemake integration
    main(
        regions_path=snakemake.input.regions,
        output_path=snakemake.output.svg,
    )
