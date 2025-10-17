#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Visualise annualised land-use change emission factors (LEFs)."""

import sys
from pathlib import Path

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).parent.parent))
from doc_figures_config import (
    COLORMAPS,
    FIGURE_SIZES,
    apply_doc_style,
    save_doc_figure,
)


def _load_lef(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    ds = xr.load_dataset(path)
    lef = ds["LEF_tCO2_per_ha_yr"].astype(np.float32)
    lat = ds["y"].astype(np.float32).values
    lon = ds["x"].astype(np.float32).values
    uses = [str(u) for u in lef.coords["use"].values]

    data = lef.values
    data = np.where(np.isfinite(data), data, np.nan)

    if lat[0] > lat[-1]:
        lat = lat[::-1]
        data = data[:, ::-1, :]
    if lon[0] > lon[-1]:
        lon = lon[::-1]
        data = data[:, :, ::-1]

    return data, lat, lon, uses


def _symmetric_limits(arrays: list[np.ndarray], percentile: float = 99.0) -> float:
    finite_vals = np.concatenate(
        [a[np.isfinite(a)] for a in arrays if np.any(np.isfinite(a))]
    )
    if finite_vals.size == 0:
        return 1.0
    limit = float(np.nanpercentile(np.abs(finite_vals), percentile))
    return max(limit, 0.1)


def main(
    annualized_path: str,
    regions_path: str,
    svg_output_path: str,
    png_output_path: str,
) -> None:
    apply_doc_style()

    lef_cube, lat, lon, uses = _load_lef(annualized_path)

    regions = gpd.read_file(regions_path)
    if regions.crs is None:
        regions = regions.set_crs(4326, allow_override=True)
    else:
        regions = regions.to_crs(4326)

    use_to_label = {
        "cropland": "Cropland land-use emission factor",
        "pasture": "Pasture land-use emission factor",
        "spared": "Spared land-use emission factor",
    }

    panels = []
    for use, data in zip(uses, lef_cube):
        panels.append((use, data))

    vmax = _symmetric_limits([arr for _, arr in panels])

    ncols = 2
    nrows = 2
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(FIGURE_SIZES["map_wide"][0], FIGURE_SIZES["map_wide"][1] * 1.4),
        subplot_kw={"projection": ccrs.EqualEarth()},
    )
    axes = axes.flatten()

    extent = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]

    for ax, (use, data) in zip(axes, panels):
        ax.set_global()
        ax.set_facecolor("white")
        im = ax.imshow(
            data,
            extent=extent,
            transform=ccrs.PlateCarree(),
            origin="lower",
            cmap=COLORMAPS["diverging"],
            vmin=-vmax,
            vmax=vmax,
            interpolation="nearest",
        )
        ax.coastlines(linewidth=0.3, color="#666666", alpha=0.4)
        ax.add_geometries(
            regions.geometry,
            crs=ccrs.PlateCarree(),
            facecolor="none",
            edgecolor="black",
            linewidth=0.2,
            alpha=0.3,
        )
        ax.set_title(use_to_label.get(use, use.title()), fontsize=11, pad=8)

    # If there are fewer panels than axes (e.g. last panel missing), hide extras
    for ax in axes[len(panels) :]:
        ax.set_visible(False)

    fig.subplots_adjust(
        left=0.05, right=0.97, top=0.91, bottom=0.16, wspace=0.18, hspace=0.25
    )

    cbar = fig.colorbar(
        im,
        ax=[a for a in axes if a.get_visible()],
        orientation="horizontal",
        fraction=0.035,
        pad=0.08,
    )
    cbar.set_label("Land-use emission factor (tCO₂ per ha per year)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    fig.text(
        0.05,
        0.09,
        "Positive values = emissions cost    Negative values = sequestration credit",
        fontsize=8,
    )

    save_doc_figure(fig, svg_output_path, format="svg")
    save_doc_figure(fig, png_output_path, format="png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main(
        annualized_path=snakemake.input.annualized,  # type: ignore[name-defined]
        regions_path=snakemake.input.regions,  # type: ignore[name-defined]
        svg_output_path=snakemake.output.svg,  # type: ignore[name-defined]
        png_output_path=snakemake.output.png,  # type: ignore[name-defined]
    )
