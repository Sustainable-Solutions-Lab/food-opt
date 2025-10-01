#! /usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Plot objective breakdown and visualize health risk factors by region."""

import logging
from collections import defaultdict
from dataclasses import dataclass
from math import exp
from pathlib import Path
from typing import Dict, Iterable, Mapping

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
import geopandas as gpd
import matplotlib

matplotlib.use("pdf")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import pypsa

logger = logging.getLogger(__name__)


@dataclass
class HealthInputs:
    risk_breakpoints: pd.DataFrame
    cluster_cause: pd.DataFrame
    cause_log_breakpoints: pd.DataFrame
    cluster_summary: pd.DataFrame
    clusters: pd.DataFrame
    food_map: pd.DataFrame
    population: pd.DataFrame
    cluster_risk_baseline: pd.DataFrame


@dataclass
class HealthResults:
    cause_costs: pd.DataFrame
    risk_costs: pd.DataFrame
    intake: pd.DataFrame
    cluster_population: Mapping[int, float]


def sanitize_identifier(value: str) -> str:
    return (
        value.replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("/", "_")
        .replace("-", "_")
    )


def sanitize_food_name(food: str) -> str:
    return sanitize_identifier(food)


def objective_category(n: pypsa.Network, component: str, **_: object) -> pd.Series:
    """Group assets into high-level categories for system cost aggregation."""

    static = n.components[component].static
    if static.empty:
        return pd.Series(dtype="object")

    index = static.index
    if component == "Link":
        mapping = {
            "produce": "Crop production",
            "trade": "Trade",
            "convert": "Processing",
            "consume": "Consumption",
        }
        categories = [
            mapping.get(str(name).split("_", 1)[0], "Other") for name in index
        ]
        return pd.Series(categories, index=index, name="category")

    return pd.Series(component, index=index, name="category")


def _build_food_lookup(
    food_map: pd.DataFrame,
) -> Dict[str, list[dict[str, float | str]]]:
    lookup: Dict[str, list[dict[str, float | str]]] = {}
    for sanitized, group in food_map.groupby("sanitized"):
        lookup[sanitized] = group[["risk_factor", "share"]].to_dict("records")
    return lookup


def _cluster_population(
    cluster_summary: pd.DataFrame,
    clusters: pd.DataFrame,
    population: pd.DataFrame,
) -> Dict[int, float]:
    clusters = clusters.assign(country_iso3=lambda df: df["country_iso3"].str.upper())
    cluster_lookup = (
        clusters.set_index("country_iso3")["health_cluster"].astype(int).to_dict()
    )

    cluster_summary = cluster_summary.assign(
        health_cluster=lambda df: df["health_cluster"].astype(int)
    )
    baseline = (
        cluster_summary.set_index("health_cluster")["population_persons"]
        .astype(float)
        .to_dict()
    )
    population = population.assign(iso3=lambda df: df["iso3"].str.upper())
    population_map = population.set_index("iso3")["population"].astype(float).to_dict()

    result: Dict[int, float] = {}
    for cluster, base_value in baseline.items():
        members = [iso for iso, c in cluster_lookup.items() if c == cluster]
        planning = sum(population_map.get(iso, 0.0) for iso in members)
        result[int(cluster)] = planning if planning > 0 else float(base_value)

    return result


def _prepare_health_inputs(
    inputs: HealthInputs,
) -> tuple[
    Dict[str, pd.DataFrame],
    Dict[str, pd.DataFrame],
    Dict[str, list[dict[str, float | str]]],
    Dict[str, int],
    Dict[int, float],
    Dict[int, float],
]:
    risk_tables = {}
    for risk, group in inputs.risk_breakpoints.groupby("risk_factor"):
        pivot = (
            group.sort_values(["intake_g_per_day", "cause"])
            .pivot_table(
                index="intake_g_per_day",
                columns="cause",
                values="log_rr",
                aggfunc="first",
            )
            .sort_index()
        )
        risk_tables[str(risk)] = pivot

    cause_tables = {
        str(cause): df.sort_values("log_rr_total")
        for cause, df in inputs.cause_log_breakpoints.groupby("cause")
    }

    food_map = inputs.food_map.copy()
    food_map["sanitized"] = food_map["food"].apply(sanitize_food_name)
    food_lookup = _build_food_lookup(food_map)

    clusters = inputs.clusters.assign(
        country_iso3=lambda df: df["country_iso3"].str.upper()
    )
    cluster_lookup = (
        clusters.set_index("country_iso3")["health_cluster"].astype(int).to_dict()
    )

    cluster_summary = inputs.cluster_summary.assign(
        health_cluster=lambda df: df["health_cluster"].astype(int)
    )
    value_per_yll = (
        cluster_summary.set_index("health_cluster")["value_per_yll_usd_per_yll"]
        .astype(float)
        .to_dict()
    )
    cluster_population = _cluster_population(
        cluster_summary,
        clusters,
        inputs.population,
    )

    return (
        risk_tables,
        cause_tables,
        food_lookup,
        cluster_lookup,
        value_per_yll,
        cluster_population,
    )


def compute_health_results(n: pypsa.Network, inputs: HealthInputs) -> HealthResults:
    (
        risk_tables,
        cause_tables,
        food_lookup,
        cluster_lookup,
        value_per_yll,
        cluster_population,
    ) = _prepare_health_inputs(inputs)

    cluster_cause = inputs.cluster_cause.assign(
        health_cluster=lambda df: df["health_cluster"].astype(int)
    )

    intake_totals: dict[tuple[int, str], float] = defaultdict(float)
    p_now = n.links_t.p0.loc["now"]

    for link_name in n.links.index:
        name = str(link_name)
        if not name.startswith("consume_"):
            continue
        base, _, country = name.rpartition("_")
        if len(country) != 3:
            continue
        sanitized_food = base[len("consume_") :]
        entries = food_lookup.get(sanitized_food)
        if not entries:
            continue
        cluster = cluster_lookup.get(country.upper())
        if cluster is None:
            continue
        population = cluster_population.get(int(cluster), 0.0)
        if population <= 0:
            continue
        flow = float(p_now.get(link_name, 0.0))
        if flow == 0:
            continue
        scale = 1_000_000.0 / (365.0 * population)
        for entry in entries:
            share = float(entry["share"])
            if share <= 0:
                continue
            risk = str(entry["risk_factor"])
            intake_totals[(int(cluster), risk)] += flow * share * scale

    intake_series = (
        pd.Series(intake_totals, dtype=float).rename("intake_g_per_day").sort_index()
    )
    intake_df = (
        intake_series.reset_index()
        if not intake_series.empty
        else pd.DataFrame(columns=["cluster", "risk_factor", "intake_g_per_day"])
    )

    cause_records: list[dict[str, float | int]] = []
    risk_costs: dict[tuple[int, str], float] = defaultdict(float)

    for (cluster, cause), row in cluster_cause.set_index(
        ["health_cluster", "cause"]
    ).iterrows():
        cluster = int(cluster)
        cause = str(cause)
        value = float(value_per_yll.get(cluster, 0.0))
        yll_base = float(row.get("yll_base", 0.0))
        if value <= 0 or yll_base <= 0:
            continue
        coeff = value * yll_base
        rr_ref = exp(float(row.get("log_rr_total_ref", 0.0)))

        risk_contribs: dict[str, float] = {}
        total_log = 0.0

        for risk, table in risk_tables.items():
            if cause not in table.columns:
                continue
            xs = table.index.to_numpy(dtype=float)
            if xs.size == 0:
                continue
            ys = table[cause].to_numpy(dtype=float)
            intake_value = float(intake_totals.get((cluster, risk), 0.0))
            contribution = float(np.interp(intake_value, xs, ys))
            risk_contribs[risk] = contribution
            total_log += contribution

        cause_bp = cause_tables.get(cause)
        if cause_bp is None or cause_bp.empty:
            continue

        log_points = cause_bp["log_rr_total"].to_numpy(dtype=float)
        rr_points = cause_bp["rr_total"].to_numpy(dtype=float)
        rr_total = float(np.interp(total_log, log_points, rr_points))
        cost = coeff / rr_ref * rr_total - coeff

        cause_records.append(
            {
                "cluster": cluster,
                "cause": cause,
                "cost": cost,
                "log_total": total_log,
                "rr_total": rr_total,
                "coeff": coeff,
            }
        )

        if abs(total_log) > 1e-12 and cost != 0.0:
            for risk, contribution in risk_contribs.items():
                share = contribution / total_log if total_log != 0 else 0.0
                risk_costs[(cluster, risk)] += cost * share

    cause_df = pd.DataFrame(cause_records)
    risk_df = pd.DataFrame(
        (
            {
                "cluster": cluster,
                "risk_factor": risk,
                "cost": value,
            }
            for (cluster, risk), value in risk_costs.items()
        )
    )

    return HealthResults(
        cause_costs=cause_df,
        risk_costs=risk_df,
        intake=intake_df,
        cluster_population=cluster_population,
    )


def compute_baseline_risk_costs(inputs: HealthInputs) -> pd.DataFrame:
    (
        risk_tables,
        cause_tables,
        _food_lookup,
        _cluster_lookup,
        value_per_yll,
        _cluster_population,
    ) = _prepare_health_inputs(inputs)

    baseline = inputs.cluster_risk_baseline.assign(
        health_cluster=lambda df: df["health_cluster"].astype(int),
        risk_factor=lambda df: df["risk_factor"].astype(str),
    )
    baseline_intake = {
        (int(row.health_cluster), str(row.risk_factor)): float(
            row.baseline_intake_g_per_day
        )
        for row in baseline.itertuples(index=False)
    }

    cluster_cause = inputs.cluster_cause.assign(
        health_cluster=lambda df: df["health_cluster"].astype(int)
    )

    records: list[dict[str, float | int | str]] = []

    for (cluster, cause), row in cluster_cause.set_index(
        [
            "health_cluster",
            "cause",
        ]
    ).iterrows():
        cluster = int(cluster)
        value = float(value_per_yll.get(cluster, 0.0))
        yll_base = float(row.get("yll_base", 0.0))
        if value <= 0 or yll_base <= 0:
            continue
        coeff = value * yll_base

        contributions: dict[str, float] = {}
        total_log = 0.0

        for risk, table in risk_tables.items():
            if cause not in table.columns:
                continue
            xs = table.index.to_numpy(dtype=float)
            if xs.size == 0:
                continue
            ys = table[cause].to_numpy(dtype=float)
            intake_value = float(baseline_intake.get((cluster, risk), 0.0))
            contribution = float(np.interp(intake_value, xs, ys))
            contributions[risk] = contribution
            total_log += contribution

        if not contributions or abs(total_log) <= 1e-12:
            continue

        for risk, contribution in contributions.items():
            share = contribution / total_log if total_log != 0 else 0.0
            records.append(
                {
                    "cluster": cluster,
                    "risk_factor": risk,
                    "cost": coeff * share,
                }
            )

    result = pd.DataFrame(records)
    if not result.empty:
        result["cluster"] = result["cluster"].astype(int)
    return result.reset_index(drop=True)


def compute_system_costs(n: pypsa.Network) -> pd.Series:
    costs = n.statistics.system_cost(groupby=objective_category)
    if isinstance(costs, pd.DataFrame):
        costs = costs.iloc[:, 0]
    if costs.empty:
        return pd.Series(dtype=float)
    idx = costs.index
    if "category" not in idx.names:
        idx = idx.set_names(list(idx.names[:-1]) + ["category"])
        costs.index = idx
    return costs.groupby("category").sum().sort_values(ascending=False)


def ghg_cost_from_network(n: pypsa.Network, ghg_price: float) -> float:
    co2 = (
        float(n.stores_t.e.loc["now", "co2"]) if "co2" in n.stores_t.e.columns else 0.0
    )
    ch4 = (
        float(n.stores_t.e.loc["now", "ch4"]) if "ch4" in n.stores_t.e.columns else 0.0
    )
    return ghg_price * (co2 + 25.0 * ch4)


def choose_scale(values: Iterable[float]) -> tuple[float, str]:
    max_val = max((abs(v) for v in values), default=1.0)
    if max_val >= 1e12:
        return 1e12, "trillion USD"
    if max_val >= 1e9:
        return 1e9, "billion USD"
    if max_val >= 1e6:
        return 1e6, "million USD"
    return 1.0, "USD"


def plot_cost_breakdown(series: pd.Series, output_path: Path) -> None:
    if series.empty:
        logger.warning("No cost data available for plotting")
        return

    scale, label = choose_scale(series.values)
    values = series / scale

    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)
    bars = ax.bar(series.index, values, color="#4e79a7")

    ax.set_ylabel(f"Cost ({label})")
    ax.set_title("Objective Breakdown by Category")
    ax.grid(axis="y", linestyle="--", alpha=0.3)

    for bar, raw in zip(bars, series.values):
        height = bar.get_height()
        ax.annotate(
            f"{raw / scale:,.2f}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close(fig)


def build_cluster_risk_tables(
    risk_costs_df: pd.DataFrame,
    cluster_population: Mapping[int, float],
) -> tuple[dict[str, dict[int, float]], dict[str, dict[int, float]]]:
    if risk_costs_df.empty:
        return {}, {}

    risk_costs = risk_costs_df.copy()
    risk_costs["cluster"] = risk_costs["cluster"].astype(int)

    populations = pd.Series(cluster_population, name="population")
    risk_costs = risk_costs.merge(
        populations.rename_axis("cluster").reset_index(), on="cluster", how="left"
    )
    risk_costs["cost_per_capita"] = risk_costs["cost"] / risk_costs["population"]

    cost_map: dict[str, dict[int, float]] = defaultdict(dict)
    per_capita_map: dict[str, dict[int, float]] = defaultdict(dict)

    for row in risk_costs.itertuples(index=False):
        risk = str(row.risk_factor)
        cluster = int(row.cluster)
        cost_map[risk][cluster] = float(row.cost)
        per_capita_map[risk][cluster] = float(row.cost_per_capita)

    return cost_map, per_capita_map


def plot_health_map(
    gdf: gpd.GeoDataFrame,
    cluster_lookup: Mapping[str, int],
    per_capita_by_risk: Mapping[str, Mapping[int, float]],
    output_path: Path,
    top_risks: Iterable[str],
    *,
    diverging: bool = True,
    value_label: str = "Health cost per capita (USD)",
) -> None:
    risks = list(top_risks)
    if not risks:
        logger.warning(
            "No risk factors available for mapping; creating placeholder figure"
        )
        fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No health risk data available",
            ha="center",
            va="center",
            fontsize=12,
            color="#555555",
        )
        fig.savefig(output_path, bbox_inches="tight", dpi=300)
        plt.close(fig)
        return

    plate = ccrs.PlateCarree()
    fig, axes = plt.subplots(
        1,
        len(risks),
        figsize=(5.2 * len(risks), 5.4),
        dpi=150,
        subplot_kw={"projection": ccrs.EqualEarth()},
    )

    if len(risks) == 1:
        axes = [axes]  # type: ignore[list-item]

    for ax, risk in zip(axes, risks):
        data = gdf.copy()
        cluster_map = per_capita_by_risk.get(risk, {})
        data["value"] = data["country"].map(
            lambda iso: cluster_map.get(cluster_lookup.get(iso, -1))
        )

        values = data["value"].dropna()
        if diverging:
            if values.empty:
                vmin, vmax = -1.0, 1.0
            else:
                vmin, vmax = values.min(), values.max()
            bound = max(abs(vmin), abs(vmax), 1e-6)
            norm = mcolors.TwoSlopeNorm(vmin=-bound, vcenter=0, vmax=bound)
            cmap = matplotlib.colormaps["RdBu_r"]
        else:
            vmax = float(values.max()) if not values.empty else 1.0
            if not np.isfinite(vmax) or vmax <= 0:
                vmax = 1.0
            norm = mcolors.Normalize(vmin=0.0, vmax=vmax)
            cmap = matplotlib.colormaps["Blues"]

        ax.set_facecolor("#f7f9fb")
        ax.set_global()

        data.plot(
            ax=ax,
            transform=plate,
            column="value",
            cmap=cmap,
            norm=norm,
            linewidth=0.2,
            edgecolor="#666666",
            missing_kwds={
                "color": "#eeeeee",
                "edgecolor": "#999999",
                "hatch": "///",
                "label": "No data",
            },
        )

        gl = ax.gridlines(
            draw_labels=True,
            crs=plate,
            linewidth=0.3,
            color="#aaaaaa",
            alpha=0.5,
            linestyle="--",
        )
        gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 60))
        gl.ylocator = mticker.FixedLocator(np.arange(-60, 61, 30))
        gl.xformatter = LongitudeFormatter(number_format=".0f")
        gl.yformatter = LatitudeFormatter(number_format=".0f")
        gl.xlabel_style = {"size": 7, "color": "#555555"}
        gl.ylabel_style = {"size": 7, "color": "#555555"}
        gl.top_labels = False
        gl.right_labels = False

        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", pad=0.05, shrink=0.75)
        cbar.ax.set_xlabel(value_label, fontsize=8)
        cbar.ax.tick_params(labelsize=7)

        ax.set_title(risk.replace("_", " ").title(), fontsize=11)

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close(fig)


def build_health_region_table(
    gdf: gpd.GeoDataFrame,
    cluster_lookup: Mapping[str, int],
    cost_by_risk: Mapping[str, Mapping[int, float]],
    per_capita_by_risk: Mapping[str, Mapping[int, float]],
) -> pd.DataFrame:
    records: list[dict[str, object]] = []

    for _, row in gdf[["region", "country"]].iterrows():
        iso = str(row.country)
        cluster = cluster_lookup.get(iso)
        if cluster is None:
            continue
        for risk, cluster_costs in cost_by_risk.items():
            records.append(
                {
                    "region": row.region,
                    "country": iso,
                    "health_cluster": cluster,
                    "risk_factor": risk,
                    "cost_usd": cluster_costs.get(cluster, float("nan")),
                    "cost_per_capita_usd": per_capita_by_risk.get(risk, {}).get(
                        cluster, float("nan")
                    ),
                }
            )

    return pd.DataFrame(records)


def main() -> None:
    n = pypsa.Network(snakemake.input.network)  # type: ignore[name-defined]
    logger.info("Loaded network with objective %.3e", n.objective)

    health_inputs = HealthInputs(
        risk_breakpoints=pd.read_csv(snakemake.input.risk_breakpoints),  # type: ignore[attr-defined]
        cluster_cause=pd.read_csv(snakemake.input.health_cluster_cause),
        cause_log_breakpoints=pd.read_csv(snakemake.input.health_cause_log),
        cluster_summary=pd.read_csv(snakemake.input.health_cluster_summary),
        clusters=pd.read_csv(snakemake.input.health_clusters),
        food_map=pd.read_csv(
            snakemake.input.food_risk_map,  # type: ignore[attr-defined]
            comment="#",
            names=["food", "risk_factor", "share"],
        ),
        population=pd.read_csv(snakemake.input.population),
        cluster_risk_baseline=pd.read_csv(snakemake.input.health_cluster_risk_baseline),
    )

    health_results = compute_health_results(n, health_inputs)

    system_costs = compute_system_costs(n)
    ghg_price = float(snakemake.params.ghg_price)
    ghg_cost = ghg_cost_from_network(n, ghg_price)

    total_series = system_costs.copy()
    if not health_results.cause_costs.empty:
        health_total = health_results.cause_costs["cost"].sum()
        total_series.loc["Health burden"] = health_total
        logger.info("Computed total health contribution %.3e USD", health_total)
    else:
        logger.warning("Health results are empty; skipping health contribution")

    if ghg_cost != 0.0:
        total_series.loc["GHG pricing"] = ghg_cost
        logger.info("Computed GHG contribution %.3e USD", ghg_cost)

    total_series = total_series.sort_values(key=np.abs, ascending=False)

    breakdown_csv = Path(snakemake.output.breakdown_csv)  # type: ignore[attr-defined]
    breakdown_pdf = Path(snakemake.output.breakdown_pdf)
    breakdown_csv.parent.mkdir(parents=True, exist_ok=True)

    total_series.rename("cost_usd").to_csv(breakdown_csv, header=True)
    plot_cost_breakdown(total_series, breakdown_pdf)
    logger.info("Wrote objective breakdown plot to %s", breakdown_pdf)

    (
        cost_by_risk,
        per_capita_by_risk,
    ) = build_cluster_risk_tables(
        health_results.risk_costs, health_results.cluster_population
    )

    regions_gdf = gpd.read_file(snakemake.input.regions)  # type: ignore[attr-defined]
    if regions_gdf.crs is None:
        regions_gdf = regions_gdf.set_crs(4326, allow_override=True)
    else:
        regions_gdf = regions_gdf.to_crs(4326)

    cluster_lookup = (
        health_inputs.clusters.assign(
            country_iso3=lambda df: df["country_iso3"].str.upper()
        )
        .set_index("country_iso3")["health_cluster"]
        .astype(int)
        .to_dict()
    )
    regions_gdf = regions_gdf.assign(country=lambda df: df["country"].str.upper())

    risk_totals = (
        health_results.risk_costs.groupby("risk_factor")["cost"].sum()
        if not health_results.risk_costs.empty
        else pd.Series(dtype=float)
    )
    if not risk_totals.empty:
        top_risks = risk_totals.reindex(
            risk_totals.abs().sort_values(ascending=False).index
        ).index[: min(3, len(risk_totals))]
    else:
        top_risks = list(per_capita_by_risk.keys())[:3]

    map_pdf = Path(snakemake.output.health_map_pdf)
    plot_health_map(
        regions_gdf,
        cluster_lookup,
        per_capita_by_risk,
        map_pdf,
        top_risks,
        diverging=True,
        value_label="Health cost per capita (USD)",
    )
    logger.info("Saved health risk map to %s", map_pdf)

    region_table = build_health_region_table(
        regions_gdf,
        cluster_lookup,
        cost_by_risk,
        per_capita_by_risk,
    )
    region_table.to_csv(Path(snakemake.output.health_map_csv), index=False)
    logger.info("Wrote regional health table to %s", snakemake.output.health_map_csv)

    # Baseline health burden maps
    baseline_risk_costs = compute_baseline_risk_costs(health_inputs)
    (
        baseline_cost_by_risk,
        baseline_per_capita_by_risk,
    ) = build_cluster_risk_tables(
        baseline_risk_costs,
        health_results.cluster_population,
    )

    baseline_totals = (
        baseline_risk_costs.groupby("risk_factor")["cost"].sum()
        if not baseline_risk_costs.empty
        else pd.Series(dtype=float)
    )
    if not baseline_totals.empty:
        baseline_top = baseline_totals.reindex(
            baseline_totals.abs().sort_values(ascending=False).index
        ).index[: min(3, len(baseline_totals))]
    else:
        baseline_top = list(baseline_per_capita_by_risk.keys())[:3]

    baseline_map_pdf = Path(snakemake.output.health_baseline_map_pdf)
    plot_health_map(
        regions_gdf,
        cluster_lookup,
        baseline_per_capita_by_risk,
        baseline_map_pdf,
        baseline_top,
        diverging=True,
        value_label="Baseline health cost per capita (USD)",
    )
    logger.info("Saved baseline health risk map to %s", baseline_map_pdf)

    baseline_region_table = build_health_region_table(
        regions_gdf,
        cluster_lookup,
        baseline_cost_by_risk,
        baseline_per_capita_by_risk,
    )
    baseline_region_table.to_csv(
        Path(snakemake.output.health_baseline_map_csv), index=False
    )
    logger.info(
        "Wrote baseline regional health table to %s",
        snakemake.output.health_baseline_map_csv,
    )


if __name__ == "__main__":
    main()
