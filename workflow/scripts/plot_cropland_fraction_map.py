#! /usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Plot a world map of cropland fraction per region.

Inputs (via snakemake):
- input.network: solved PyPSA network (NetCDF)
- input.regions: regions GeoJSON with a 'region' column
- input.land_area_by_class: CSV with columns [region, water_supply, resource_class, area_ha]

Output:
- results/{name}/plots/cropland_fraction_map.pdf

Notes:
- Cropland use is computed from actual land flows supplied by the land
  generators (carrier 'land'), i.e. n.generators_t.p for generators named
  like 'land_{region}_class{k}_{ws}', summed over classes and water supplies.
  This reflects realized land use, not capacity.
- Total land area per region is the sum of area_ha over all resource classes
  and water supplies from land_area_by_class.csv, matching the model’s land
  availability basis.
"""

from pathlib import Path
import logging
import geopandas as gpd
import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa

logger = logging.getLogger(__name__)


def _used_cropland_area_by_region(n: pypsa.Network) -> pd.Series:
    """Return a Series mapping region -> used cropland area [ha].

    Sums positive output from land generators (carrier 'land') at snapshot
    'now'. Generator names equal their land bus names: land_{region}_class{k}_{ws}.
    """
    if n.generators.empty or n.generators_t.p.empty:
        return pd.Series(dtype=float)
    # Filter land generators
    land_gen = n.generators[n.generators["carrier"] == "land"]
    if land_gen.empty:
        return pd.Series(dtype=float)
    names = land_gen.index
    # Extract flows at 'now' (use 0 if snapshot missing)
    if "now" in n.snapshots:
        p_now = n.generators_t.p.loc["now", names]
    else:
        # fallback: take first snapshot
        p_now = n.generators_t.p.iloc[0][names]
    p_now = p_now.fillna(0.0)

    # Group by region parsed from name
    def region_from_name(s: str) -> str:
        parts = s.split("_")
        return parts[1] if len(parts) >= 2 else "unknown"

    regions = [region_from_name(str(nm)) for nm in names]
    df = pd.DataFrame({"region": regions, "p": p_now.values})
    df["p"] = df["p"].clip(lower=0.0)
    used = df.groupby("region")["p"].sum()
    return used.astype(float)


def main() -> None:
    n = pypsa.Network(snakemake.input.network)  # type: ignore[name-defined]
    regions_path: str = snakemake.input.regions  # type: ignore[name-defined]
    land_area_csv: str = snakemake.input.land_area_by_class  # type: ignore[name-defined]
    out_pdf = Path(snakemake.output.pdf)  # type: ignore[name-defined]

    # Load regions and set projection
    gdf = gpd.read_file(regions_path)
    if gdf.crs is None:
        gdf = gdf.set_crs(4326, allow_override=True)
    gdf = gdf.set_index("region", drop=False)
    gdf_ee = gdf.to_crs("+proj=eqearth")

    # Used cropland area from solved network
    used_ha = _used_cropland_area_by_region(n)

    # Total available land area from preprocessing
    df_land = pd.read_csv(land_area_csv)
    if not {"region", "area_ha"}.issubset(df_land.columns):
        raise ValueError("land_area_by_class.csv must contain 'region' and 'area_ha'")
    total_ha = df_land.groupby("region")["area_ha"].sum().astype(float)

    # Align to regions
    idx = gdf.index
    used_ha = used_ha.reindex(idx).fillna(0.0)
    total_ha = total_ha.reindex(idx).fillna(0.0)

    # Fraction and mask for undefined denominators
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = (used_ha / total_ha).replace([np.inf, -np.inf], np.nan)
    frac = frac.clip(lower=0.0, upper=1.0)

    gdf_ee = gdf_ee.copy()
    gdf_ee["cropland_used_ha"] = used_ha.values
    gdf_ee["land_total_ha"] = total_ha.values
    gdf_ee["cropland_fraction"] = frac.values

    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    # Plot
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.set_axis_off()

    # Base map
    gdf_ee.plot(ax=ax, linewidth=0.3, edgecolor="#666666", facecolor="#f5f7f9")

    # Valid fraction regions
    valid = gdf_ee[~gdf_ee["cropland_fraction"].isna()]
    if not valid.empty:
        vmin, vmax = 0.0, 0.5
        cmap = plt.get_cmap("YlGn")
        coll = valid.plot(
            ax=ax,
            column="cropland_fraction",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            linewidth=0.3,
            edgecolor="#666666",
            legend=False,
        )
        # Add colorbar
        sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
        cbar.set_label("Cropland / total model land area")

    # Hatch regions with zero total area (undefined fraction)
    zero_total = gdf_ee[gdf_ee["land_total_ha"] <= 0]
    if not zero_total.empty:
        zero_total.plot(
            ax=ax,
            facecolor="#f0f0f0",
            edgecolor="#666666",
            linewidth=0.3,
            hatch="//",
        )

    ax.set_title("Cropland Fraction by Region")
    plt.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    plt.close(fig)

    # Optional CSV sidecar for reference
    csv_out = out_pdf.with_suffix("")
    csv_out = csv_out.parent / f"{csv_out.name}_by_region.csv"
    gdf_ee[["region", "cropland_used_ha", "land_total_ha", "cropland_fraction"]].to_csv(
        csv_out, index=False
    )
    logger.info("Saved cropland fraction map to %s and CSV to %s", out_pdf, csv_out)


if __name__ == "__main__":
    main()
