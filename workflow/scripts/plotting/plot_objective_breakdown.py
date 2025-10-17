#! /usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Generate an objective breakdown plot with updated emissions accounting."""

import logging
import sys
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_health_impacts import (
    HealthInputs,
    compute_health_results,
)

logger = logging.getLogger(__name__)

TONNE_TO_MEGATONNE = 1e-6


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


def compute_system_costs(n: pypsa.Network) -> pd.Series:
    """Aggregate system costs by the objective categories defined above."""

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


def compute_ghg_cost_breakdown(n: pypsa.Network, ghg_price: float) -> dict[str, float]:
    """Return the objective contribution from the priced GHG store."""

    if len(n.snapshots) == 0:
        return {}

    snapshot = n.snapshots[-1]
    if snapshot not in n.stores_t.e.index:
        return {}

    store_levels = n.stores_t.e.loc[snapshot]
    if "ghg" not in store_levels.index:
        return {}

    level_mt = float(store_levels["ghg"])
    if level_mt == 0.0:
        return {}

    contribution = ghg_price * level_mt / TONNE_TO_MEGATONNE
    label = "GHG pricing (CO₂-eq)"
    logger.info(
        "Computed %s contribution %.3e USD (level %.3e MtCO2-eq, price %.2f USD/tCO2-eq)",
        label,
        contribution,
        level_mt,
        ghg_price,
    )
    return {label: contribution}


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


def load_health_inputs() -> HealthInputs:
    """Load the health-related input tables referenced by Snakemake."""

    return HealthInputs(
        risk_breakpoints=pd.read_csv(snakemake.input.risk_breakpoints),  # type: ignore[name-defined]
        cluster_cause=pd.read_csv(snakemake.input.health_cluster_cause),
        cause_log_breakpoints=pd.read_csv(snakemake.input.health_cause_log),
        cluster_summary=pd.read_csv(snakemake.input.health_cluster_summary),
        clusters=pd.read_csv(snakemake.input.health_clusters),
        population=pd.read_csv(snakemake.input.population),
        cluster_risk_baseline=pd.read_csv(snakemake.input.health_cluster_risk_baseline),
    )


def compute_health_total(
    n: pypsa.Network,
    health_inputs: HealthInputs,
    risk_factors: list[str],
    value_per_yll: float,
    tmrel_g_per_day: dict[str, float],
    food_groups_df: pd.DataFrame,
) -> float:
    """Return aggregate health burden contribution to the objective."""

    health_results = compute_health_results(
        n,
        health_inputs,
        risk_factors,
        value_per_yll,
        tmrel_g_per_day,
        food_groups_df,
    )
    if health_results.cause_costs.empty:
        logger.warning("Health results are empty; skipping health objective term")
        return 0.0

    total = float(health_results.cause_costs["cost"].sum())
    logger.info("Computed total health contribution %.3e USD", total)
    return total


def main() -> None:
    n = pypsa.Network(snakemake.input.network)  # type: ignore[name-defined]
    logger.info("Loaded network with objective %.3e", n.objective)

    system_costs = compute_system_costs(n)
    store_cost = 0.0
    if "Store" in system_costs.index:
        store_cost = float(system_costs.loc["Store"])
        system_costs = system_costs.drop(index="Store")
        if store_cost != 0.0:
            logger.info(
                "Removed raw 'Store' category contribution %.3e USD to avoid double counting GHG costs",
                store_cost,
            )

    total_series = system_costs.copy()

    # Add health-related objective contribution
    food_groups_df = pd.read_csv(snakemake.input.food_groups)
    health_total = compute_health_total(
        n,
        load_health_inputs(),
        snakemake.params.health_risk_factors,  # type: ignore[attr-defined]
        float(snakemake.params.health_value_per_yll),
        dict(snakemake.params.health_tmrel_g_per_day),
        food_groups_df,
    )
    if health_total != 0.0:
        total_series.loc["Health burden"] = health_total

    # Add priced greenhouse gas emissions
    ghg_price = float(snakemake.params.ghg_price)
    ghg_terms = compute_ghg_cost_breakdown(n, ghg_price)
    if ghg_terms:
        total_series = total_series.combine_first(pd.Series(dtype=float))
        for label, value in ghg_terms.items():
            total_series.loc[label] = value

    total_series = total_series.sort_values(key=np.abs, ascending=False)

    breakdown_csv = Path(snakemake.output.breakdown_csv)  # type: ignore[attr-defined]
    breakdown_pdf = Path(snakemake.output.breakdown_pdf)
    breakdown_csv.parent.mkdir(parents=True, exist_ok=True)

    total_series.rename("cost_usd").to_csv(breakdown_csv, header=True)
    logger.info("Wrote objective breakdown table to %s", breakdown_csv)
    plot_cost_breakdown(total_series, breakdown_pdf)
    logger.info("Wrote objective breakdown plot to %s", breakdown_pdf)


if __name__ == "__main__":
    main()
