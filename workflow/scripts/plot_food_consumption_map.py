#! /usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Plot food consumption by health cluster using pie charts on a map."""

import logging
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple

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

GRAMS_PER_MEGATONNE = 1e12
DAYS_PER_YEAR = 365
PLATE_CARREE = ccrs.PlateCarree()


def _select_snapshot(network: pypsa.Network) -> pd.Index | str:
    if "now" in network.snapshots:
        return "now"
    if len(network.snapshots) == 1:
        return network.snapshots[0]
    raise ValueError("Expected snapshot 'now' or single snapshot in solved network")


def _parse_group_and_iso(bus: str) -> Tuple[str, str] | None:
    if not bus.startswith("group_"):
        return None
    remainder = bus[len("group_") :]
    if "_" not in remainder:
        return None
    group, iso = remainder.rsplit("_", 1)
    return group, iso.upper()


def _bus_column_to_leg(column: str) -> int | None:
    if not column.startswith("bus"):
        return None
    suffix = column[len("bus") :]
    if not suffix:
        return 0
    if suffix.isdigit():
        return int(suffix)
    return None


def _link_dispatch_at_snapshot(
    network: pypsa.Network, snapshot
) -> dict[int, pd.Series]:
    dispatch: dict[int, pd.Series] = {}
    if not hasattr(network, "links_t"):
        return dispatch
    for attr in dir(network.links_t):
        if not attr.startswith("p"):
            continue
        suffix = attr[1:]
        if not suffix.isdigit():
            continue
        series = getattr(network.links_t, attr)
        if snapshot not in series.index:
            continue
        dispatch[int(suffix)] = series.loc[snapshot]
    return dispatch


def _aggregate_consume_link_totals(
    network: pypsa.Network, snapshot
) -> Dict[Tuple[str, str], float]:
    links = network.links
    if links.empty:
        return {}
    consume_links = links[links.index.str.startswith("consume_")]
    if consume_links.empty:
        return {}

    bus_columns = [col for col in consume_links.columns if col.startswith("bus")]
    if not bus_columns:
        return {}

    dispatch_lookup = _link_dispatch_at_snapshot(network, snapshot)
    if not dispatch_lookup:
        return {}

    totals: Dict[Tuple[str, str], float] = {}
    bus_values = consume_links[bus_columns]
    for link_name, row in bus_values.iterrows():
        for column, bus_value in row.items():
            if not isinstance(bus_value, str):
                continue
            parsed = _parse_group_and_iso(bus_value)
            if parsed is None:
                continue
            leg = _bus_column_to_leg(column)
            if leg is None:
                continue
            dispatch = dispatch_lookup.get(leg)
            if dispatch is None:
                continue
            value = float(dispatch.get(link_name, 0.0))
            if value == 0.0 or not np.isfinite(value):
                continue
            group, iso = parsed
            key = (iso, group)
            totals[key] = totals.get(key, 0.0) + abs(value)
    return totals


def _aggregate_country_group_mass(
    network: pypsa.Network,
    snapshot,
) -> pd.DataFrame:
    totals = _aggregate_consume_link_totals(network, snapshot)

    if not totals:
        return pd.DataFrame()

    series = pd.Series(totals, dtype=float)
    df = series.unstack(fill_value=0.0).sort_index(axis=0).sort_index(axis=1)
    df.index.name = "iso3"
    return df


def _aggregate_cluster_group_mass(
    country_group: pd.DataFrame,
    iso_to_cluster: Mapping[str, int],
) -> pd.DataFrame:
    if country_group.empty:
        return pd.DataFrame()

    data: Dict[Tuple[int, str], float] = {}
    for iso, row in country_group.iterrows():
        cluster = iso_to_cluster.get(str(iso).upper())
        if cluster is None:
            continue
        for group, value in row.items():
            if value <= 0.0 or not np.isfinite(value):
                continue
            key = (int(cluster), str(group))
            data[key] = data.get(key, 0.0) + float(value)

    if not data:
        return pd.DataFrame()

    series = pd.Series(data, dtype=float)
    df = series.unstack(fill_value=0.0).sort_index(axis=0).sort_index(axis=1)
    df.index.name = "health_cluster"
    return df


def _cluster_population(
    population_df: pd.DataFrame,
    iso_to_cluster: Mapping[str, int],
) -> Dict[int, float]:
    if population_df.empty:
        return {}
    pop_df = population_df.copy()
    if "iso3" not in pop_df.columns:
        raise ValueError("Population table must contain an 'iso3' column")
    pop_df = pop_df.assign(iso3=lambda df: df["iso3"].str.upper())
    pop_df = pop_df[pop_df["iso3"].isin(iso_to_cluster.keys())]
    pop_df = pop_df.dropna(subset=["population"])
    pop_df["population"] = pop_df["population"].astype(float)

    result: Dict[int, float] = {}
    for iso, group_df in pop_df.groupby("iso3"):
        cluster = iso_to_cluster.get(iso)
        if cluster is None:
            continue
        value = float(group_df["population"].sum())
        if value <= 0.0:
            continue
        result[int(cluster)] = result.get(int(cluster), 0.0) + value
    return result


def _colors_for_groups(groups: Sequence[str]) -> Dict[str, str]:
    cmap = plt.get_cmap("tab20")
    colors: Dict[str, str] = {}
    for idx, group in enumerate(groups):
        colors[group] = mcolors.to_hex(cmap(idx % cmap.N))
    return colors


def _prepare_cluster_geodata(
    regions_path: str,
    iso_to_cluster: Mapping[str, int],
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    regions = gpd.read_file(regions_path)
    if regions.crs is None:
        logger.warning("Regions GeoDataFrame missing CRS; assuming EPSG:4326")
        regions = regions.set_crs(4326, allow_override=True)
    else:
        regions = regions.to_crs(4326)

    if "country" not in regions.columns:
        raise ValueError("Regions GeoDataFrame must contain a 'country' column")

    regions = regions.assign(country=lambda df: df["country"].str.upper())
    regions["health_cluster"] = regions["country"].map(iso_to_cluster)
    cluster_regions = regions.dropna(subset=["health_cluster"])
    if cluster_regions.empty:
        return regions, gpd.GeoDataFrame(columns=["health_cluster"], crs=regions.crs)

    cluster_regions = cluster_regions.assign(
        health_cluster=lambda df: df["health_cluster"].astype(int)
    )
    dissolved = (
        cluster_regions.dissolve(by="health_cluster", as_index=False)
        .set_index("health_cluster", drop=False)
        .sort_index()
    )
    dissolved_eq = dissolved.to_crs("+proj=eqearth")
    return dissolved, dissolved_eq


def _draw_pie(
    ax: plt.Axes,
    x: float,
    y: float,
    sizes: Iterable[float],
    colors: Iterable[str],
    radius: float,
) -> None:
    sizes_list = list(sizes)
    total = float(sum(sizes_list))
    if radius <= 0 or not np.isfinite(total) or total <= 0:
        return
    colors_list = list(colors)
    cumulative = np.cumsum([0.0] + [max(s, 0.0) for s in sizes_list])
    norm = cumulative[-1]
    if norm <= 0:
        return
    angles = cumulative / norm * 360.0
    for idx, size in enumerate(sizes_list):
        if size <= 0:
            continue
        wedge = mpatches.Wedge(
            center=(x, y),
            r=radius,
            theta1=angles[idx],
            theta2=angles[idx + 1],
            facecolor=colors_list[idx],
            edgecolor="white",
            linewidth=0.4,
            alpha=0.9,
            zorder=5,
        )
        ax.add_patch(wedge)
    outline = mpatches.Circle(
        (x, y),
        radius,
        facecolor="none",
        edgecolor="#444444",
        linewidth=0.35,
        alpha=0.7,
        zorder=6,
    )
    ax.add_patch(outline)


def _plot_cluster_pies(
    cluster_data: pd.DataFrame,
    cluster_gdf: gpd.GeoDataFrame,
    cluster_gdf_eq: gpd.GeoDataFrame,
    colors: Mapping[str, str],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(
        figsize=(13, 6.5),
        dpi=150,
        subplot_kw={"projection": ccrs.EqualEarth()},
    )
    ax.set_facecolor("#f7f9fb")
    ax.set_global()

    has_geometry = cluster_gdf is not None and not cluster_gdf.empty
    has_geometry_eq = cluster_gdf_eq is not None and not cluster_gdf_eq.empty

    if has_geometry:
        ax.add_geometries(
            cluster_gdf.geometry,
            crs=PLATE_CARREE,
            facecolor="#e4edf4",
            edgecolor="#666666",
            linewidth=0.2,
            zorder=1,
        )

    filtered = cluster_data.copy()
    if not filtered.empty:
        if has_geometry:
            missing_clusters = sorted(set(filtered.index).difference(cluster_gdf.index))
            if missing_clusters:
                logger.warning(
                    "Dropping %d clusters without geometry: %s",
                    len(missing_clusters),
                    ", ".join(map(str, missing_clusters[:10])),
                )
            filtered = filtered.loc[
                cluster_gdf.index.intersection(filtered.index)
            ].fillna(0.0)
        group_totals = filtered.sum(axis=0).sort_values(ascending=False)
        filtered = filtered.loc[:, group_totals[group_totals > 0].index]

    if filtered.empty:
        ax.text(
            0.5,
            0.5,
            "No food consumption data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
            color="#555555",
        )
        ax.set_title("Food Consumption by Health Cluster")
        fig.savefig(output_path, bbox_inches="tight", dpi=300)
        plt.close(fig)
        return

    totals = filtered.sum(axis=1)
    if totals.max() <= 0:
        ax.text(
            0.5,
            0.5,
            "No positive totals to plot",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
            color="#555555",
        )
        ax.set_title("Food Consumption by Health Cluster")
        fig.savefig(output_path, bbox_inches="tight", dpi=300)
        plt.close(fig)
        return

    if not has_geometry_eq:
        raise ValueError("Cluster geometries required for pie placement")

    centroid_points = cluster_gdf_eq.geometry.representative_point()
    xmin, ymin, xmax, ymax = cluster_gdf_eq.total_bounds
    width = xmax - xmin
    height = ymax - ymin
    scale_extent = max(width, height)
    if not np.isfinite(scale_extent) or scale_extent <= 0:
        scale_extent = 1.0
    radius_scale = 0.025 * scale_extent
    radii = np.sqrt(totals / totals.max()) * radius_scale

    ordered_groups = filtered.columns.tolist()
    color_list = [colors[group] for group in ordered_groups]

    for cluster in filtered.index:
        point = centroid_points.loc[cluster]
        radius = float(radii.get(cluster, 0.0))
        values = filtered.loc[cluster, ordered_groups].tolist()
        _draw_pie(ax, point.x, point.y, values, color_list, radius)

    handles = [
        mpatches.Patch(facecolor=colors[group], label=group) for group in ordered_groups
    ]
    legend1 = ax.legend(
        handles,
        [f"{g}" for g in ordered_groups],
        title="Food groups",
        loc="lower left",
        bbox_to_anchor=(0.03, 0.02),
        fontsize=8,
        title_fontsize=9,
        frameon=True,
        borderpad=0.8,
        labelspacing=0.6,
        handlelength=1.2,
    )
    legend1._legend_box.align = "left"
    ax.add_artist(legend1)

    ref_fracs = np.array([0.25, 0.5, 1.0])
    ref_vals = np.unique(totals.max() * ref_fracs)
    marker_scale = 900.0
    size_handles = [
        ax.scatter(
            [],
            [],
            s=float((val / totals.max()) * marker_scale),
            facecolors="#d0d7de",
            edgecolors="#555555",
            linewidths=0.6,
            alpha=0.8,
        )
        for val in ref_vals
    ]
    size_labels = [f"{val:,.1f} g/person/day" for val in ref_vals]
    legend2 = ax.legend(
        size_handles,
        size_labels,
        title="Pie size ∝ total",
        loc="lower left",
        bbox_to_anchor=(0.38, 0.02),
        fontsize=8,
        title_fontsize=9,
        frameon=True,
        scatterpoints=1,
        borderpad=0.8,
        labelspacing=1.2,
    )
    legend2._legend_box.align = "left"
    ax.add_artist(legend2)

    gridlines = ax.gridlines(
        draw_labels=True,
        crs=PLATE_CARREE,
        linewidth=0.35,
        color="#888888",
        alpha=0.4,
        linestyle="--",
    )
    gridlines.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
    gridlines.ylocator = mticker.FixedLocator(np.arange(-60, 61, 15))
    gridlines.xformatter = LongitudeFormatter(number_format=".0f")
    gridlines.yformatter = LatitudeFormatter(number_format=".0f")
    gridlines.xlabel_style = {"size": 8, "color": "#555555"}
    gridlines.ylabel_style = {"size": 8, "color": "#555555"}
    gridlines.top_labels = False
    gridlines.right_labels = False

    ax.set_title("Food Consumption by Health Cluster (g/person/day)")
    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close(fig)


def main() -> None:
    try:
        snakemake  # type: ignore[name-defined]
    except NameError as exc:  # pragma: no cover
        raise RuntimeError("This script must be run via Snakemake") from exc

    network_path = snakemake.input.network  # type: ignore[attr-defined]
    population_path = snakemake.input.population  # type: ignore[attr-defined]
    clusters_path = snakemake.input.clusters  # type: ignore[attr-defined]
    regions_path = snakemake.input.regions  # type: ignore[attr-defined]
    output_pdf = Path(snakemake.output.pdf)  # type: ignore[attr-defined]
    output_csv = Path(snakemake.output.csv)

    logger.info("Loading network from %s", network_path)
    network = pypsa.Network(network_path)
    snapshot = _select_snapshot(network)
    logger.info("Using snapshot '%s'", snapshot)

    logger.info("Aggregating country-level consumption")
    country_group = _aggregate_country_group_mass(network, snapshot)
    logger.info("Found %d countries with food group loads", country_group.shape[0])

    cluster_df = pd.read_csv(clusters_path)
    if (
        "country_iso3" not in cluster_df.columns
        or "health_cluster" not in cluster_df.columns
    ):
        raise ValueError(
            "Cluster table must contain 'country_iso3' and 'health_cluster' columns"
        )
    iso_to_cluster = (
        cluster_df.assign(country_iso3=lambda df: df["country_iso3"].str.upper())
        .set_index("country_iso3")["health_cluster"]
        .astype(int)
        .to_dict()
    )

    cluster_mass = _aggregate_cluster_group_mass(country_group, iso_to_cluster)
    logger.info("Aggregated to %d clusters", cluster_mass.shape[0])

    population_df = pd.read_csv(population_path)
    pop_map = _cluster_population(population_df, iso_to_cluster)
    population_series = pd.Series(pop_map, dtype=float)
    missing_population = sorted(
        cluster
        for cluster in cluster_mass.index
        if population_series.get(int(cluster), 0.0) <= 0.0
    )
    if missing_population:
        raise ValueError(
            f"Missing population for clusters: {', '.join(map(str, missing_population))}"
        )

    consumption = cluster_mass.mul(GRAMS_PER_MEGATONNE).div(DAYS_PER_YEAR)
    consumption = consumption.div(population_series, axis=0)

    consumption = consumption.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    consumption = consumption.loc[(consumption.sum(axis=1) > 0)]

    if consumption.empty:
        logger.warning("No positive consumption found after conversion")

    groups = consumption.sum(axis=0).sort_values(ascending=False).index.tolist()
    colors = _colors_for_groups(groups)

    cluster_gdf, cluster_gdf_eq = _prepare_cluster_geodata(regions_path, iso_to_cluster)
    _plot_cluster_pies(consumption, cluster_gdf, cluster_gdf_eq, colors, output_pdf)

    consumption.sort_index(axis=0).sort_index(axis=1).to_csv(output_csv, index=True)
    logger.info("Saved cluster consumption data to %s", output_csv)
    logger.info("Saved map to %s", output_pdf)


if __name__ == "__main__":
    main()
