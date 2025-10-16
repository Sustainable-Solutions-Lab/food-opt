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
    risk_factors: list[str],
    value_per_yll: float,
    food_groups_df: pd.DataFrame,
) -> tuple[
    Dict[str, pd.DataFrame],
    Dict[str, pd.DataFrame],
    Dict[str, list[dict[str, float | str]]],
    Dict[str, int],
    Dict[int, float],
    float,
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

    # Load food→risk factor mapping for configured risk factors
    food_map = food_groups_df[food_groups_df["group"].isin(risk_factors)].copy()
    food_map = food_map.rename(columns={"group": "risk_factor"})
    food_map["share"] = 1.0
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
        cluster_population,
        float(value_per_yll),
    )


def compute_health_results(
    n: pypsa.Network,
    inputs: HealthInputs,
    risk_factors: list[str],
    value_per_yll: float,
    tmrel_g_per_day: dict[str, float],
    food_groups_df: pd.DataFrame,
) -> HealthResults:
    """Compute health costs from optimized network, relative to TMREL intake levels."""
    (
        risk_tables,
        cause_tables,
        food_lookup,
        cluster_lookup,
        cluster_population,
        value_per_yll_const,
    ) = _prepare_health_inputs(inputs, risk_factors, value_per_yll, food_groups_df)

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
        value = float(value_per_yll_const)
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

        # Compute cost for each risk factor in isolation (comparing intake vs TMREL)
        # This measures the cost of deviation from optimal (TMREL) intake
        for risk, log_rr_at_intake in risk_contribs.items():
            # Get log(RR) at TMREL for this risk factor
            risk_table = risk_tables[risk]
            if cause not in risk_table.columns:
                continue
            xs = risk_table.index.to_numpy(dtype=float)
            ys = risk_table[cause].to_numpy(dtype=float)
            tmrel_intake = float(tmrel_g_per_day.get(risk, 0.0))
            log_rr_at_tmrel = float(np.interp(tmrel_intake, xs, ys))

            # Cost if we change this risk factor from current intake to TMREL,
            # holding all other risk factors at their current levels
            log_total_at_tmrel = total_log - log_rr_at_intake + log_rr_at_tmrel
            rr_at_tmrel = float(np.interp(log_total_at_tmrel, log_points, rr_points))
            cost_at_tmrel = coeff / rr_ref * rr_at_tmrel - coeff
            risk_cost_contribution = cost - cost_at_tmrel
            # Cap at zero: intake beyond TMREL for protective foods provides no additional benefit
            risk_cost_contribution = max(0.0, risk_cost_contribution)
            risk_costs[(cluster, risk)] += risk_cost_contribution

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


def compute_baseline_risk_costs(
    inputs: HealthInputs,
    risk_factors: list[str],
    value_per_yll: float,
    tmrel_g_per_day: dict[str, float],
    food_groups_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute baseline health costs by risk factor and by cause, relative to TMREL intake levels.

    Health costs represent the monetized burden from deviations from optimal (TMREL) intake.
    Both total costs and individual risk factor contributions are measured relative to TMREL.

    Returns:
        (risk_costs_df, cause_costs_df) where risk_costs has columns
        (cluster, risk_factor, cost) and cause_costs has columns
        (cluster, cause, cost, log_total, rr_total, coeff).
    """
    (
        risk_tables,
        cause_tables,
        _food_lookup,
        _cluster_lookup,
        _cluster_population,
        value_per_yll_const,
    ) = _prepare_health_inputs(inputs, risk_factors, value_per_yll, food_groups_df)

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

    risk_records: list[dict[str, float | int | str]] = []
    cause_records: list[dict[str, float | int]] = []

    for (cluster, cause), row in cluster_cause.set_index(
        [
            "health_cluster",
            "cause",
        ]
    ).iterrows():
        cluster = int(cluster)
        value = float(value_per_yll_const)
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

        if not contributions:
            continue

        # Get cause breakpoints for RR interpolation
        cause_bp = cause_tables.get(cause)
        if cause_bp is None or cause_bp.empty:
            continue

        log_points = cause_bp["log_rr_total"].to_numpy(dtype=float)
        rr_points = cause_bp["rr_total"].to_numpy(dtype=float)

        # Compute total baseline cost for this cause
        rr_ref = exp(float(row.get("log_rr_total_ref", 0.0)))
        rr_total = float(np.interp(total_log, log_points, rr_points))
        total_cost = coeff / rr_ref * rr_total - coeff

        # Record cause-level cost
        cause_records.append(
            {
                "cluster": cluster,
                "cause": cause,
                "cost": total_cost,
                "log_total": total_log,
                "rr_total": rr_total,
                "coeff": coeff,
            }
        )

        # Compute cost for each risk factor in isolation (comparing baseline vs TMREL)
        # This measures the cost of deviation from optimal (TMREL) intake
        for risk, log_rr_at_baseline in contributions.items():
            # Get log(RR) at TMREL for this risk factor
            risk_table = risk_tables[risk]
            if cause not in risk_table.columns:
                continue
            xs = risk_table.index.to_numpy(dtype=float)
            ys = risk_table[cause].to_numpy(dtype=float)
            tmrel_intake = float(tmrel_g_per_day.get(risk, 0.0))
            log_rr_at_tmrel = float(np.interp(tmrel_intake, xs, ys))

            # Cost if we change this risk factor from baseline to TMREL,
            # holding all other risk factors at their baseline levels
            log_total_at_tmrel = total_log - log_rr_at_baseline + log_rr_at_tmrel
            rr_at_tmrel = float(np.interp(log_total_at_tmrel, log_points, rr_points))
            cost_at_tmrel = coeff / rr_ref * rr_at_tmrel - coeff
            risk_cost_contribution = total_cost - cost_at_tmrel
            # Cap at zero: intake beyond TMREL for protective foods provides no additional benefit
            risk_cost_contribution = max(0.0, risk_cost_contribution)

            risk_records.append(
                {
                    "cluster": cluster,
                    "risk_factor": risk,
                    "cost": risk_cost_contribution,
                }
            )

    risk_df = pd.DataFrame(risk_records)
    if not risk_df.empty:
        risk_df["cluster"] = risk_df["cluster"].astype(int)

    cause_df = pd.DataFrame(cause_records)
    if not cause_df.empty:
        cause_df["cluster"] = cause_df["cluster"].astype(int)

    return risk_df.reset_index(drop=True), cause_df.reset_index(drop=True)


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


def compute_total_health_costs_per_capita(
    cause_costs: pd.DataFrame,
    cluster_population: Mapping[int, float],
) -> dict[int, float]:
    """
    Compute total health costs per capita for each cluster.

    The total cost per cluster is the sum across causes (health costs are
    additive across different health outcomes).
    """
    if cause_costs.empty:
        return {}

    total_per_cluster = cause_costs.groupby("cluster")["cost"].sum().to_dict()

    result: dict[int, float] = {}
    for cluster, total_cost in total_per_cluster.items():
        pop = cluster_population.get(int(cluster), 0.0)
        if pop > 0:
            result[int(cluster)] = float(total_cost) / pop

    return result


def plot_health_map(
    gdf: gpd.GeoDataFrame,
    cluster_lookup: Mapping[str, int],
    per_capita_by_risk: Mapping[str, Mapping[int, float]],
    output_path: Path,
    top_risks: Iterable[str],
    *,
    diverging: bool = True,
    value_label: str = "Health cost per capita (USD)",
    total_per_capita: Mapping[int, float] | None = None,
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
    # Add one extra panel for total if provided
    num_panels = len(risks) + (1 if total_per_capita is not None else 0)

    # Calculate grid dimensions (prefer 3 columns)
    ncols = min(3, num_panels)
    nrows = (num_panels + ncols - 1) // ncols  # ceiling division

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5.2 * ncols, 5.4 * nrows),
        dpi=150,
        subplot_kw={"projection": ccrs.EqualEarth()},
    )

    # Flatten axes array for easy iteration
    if num_panels == 1:
        axes = [axes]  # type: ignore[list-item]
    else:
        axes = axes.flatten()  # type: ignore[union-attr]

    # Plot individual risk factors
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

    # Plot total health costs in the final panel if provided
    if total_per_capita is not None:
        ax = axes[-1]
        data = gdf.copy()
        data["value"] = data["country"].map(
            lambda iso: total_per_capita.get(cluster_lookup.get(iso, -1))
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

        ax.set_title("Total Health Cost", fontsize=11)

    # Hide any unused subplots
    for idx in range(num_panels, len(axes)):
        axes[idx].set_visible(False)

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
        population=pd.read_csv(snakemake.input.population),
        cluster_risk_baseline=pd.read_csv(snakemake.input.health_cluster_risk_baseline),
    )

    food_groups_df = pd.read_csv(snakemake.input.food_groups)

    value_per_yll = float(snakemake.params.health_value_per_yll)
    tmrel_g_per_day: dict[str, float] = dict(snakemake.params.health_tmrel_g_per_day)

    health_results = compute_health_results(
        n,
        health_inputs,
        snakemake.params.health_risk_factors,
        value_per_yll,
        tmrel_g_per_day,
        food_groups_df,
    )

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

    # Use all risk factors from config instead of just top 3
    all_risks = snakemake.params.health_risk_factors

    # Compute total health costs per capita
    total_per_capita = compute_total_health_costs_per_capita(
        health_results.cause_costs, health_results.cluster_population
    )

    map_pdf = Path(snakemake.output.health_map_pdf)
    plot_health_map(
        regions_gdf,
        cluster_lookup,
        per_capita_by_risk,
        map_pdf,
        all_risks,
        diverging=True,
        value_label="Health cost per capita (USD)",
        total_per_capita=total_per_capita,
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
    baseline_risk_costs, baseline_cause_costs = compute_baseline_risk_costs(
        health_inputs,
        snakemake.params.health_risk_factors,
        value_per_yll,
        tmrel_g_per_day,
        food_groups_df,
    )
    (
        baseline_cost_by_risk,
        baseline_per_capita_by_risk,
    ) = build_cluster_risk_tables(
        baseline_risk_costs,
        health_results.cluster_population,
    )

    # Compute baseline total costs from cause-level costs
    baseline_total_per_capita = compute_total_health_costs_per_capita(
        baseline_cause_costs, health_results.cluster_population
    )

    baseline_map_pdf = Path(snakemake.output.health_baseline_map_pdf)
    plot_health_map(
        regions_gdf,
        cluster_lookup,
        baseline_per_capita_by_risk,
        baseline_map_pdf,
        all_risks,
        diverging=True,
        value_label="Baseline health cost per capita (USD)",
        total_per_capita=baseline_total_per_capita,
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
