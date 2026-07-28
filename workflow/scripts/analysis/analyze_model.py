# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Run all model analysis steps with the network loaded once.

Combines extract_statistics, extract_net_emissions, extract_objective_breakdown,
extract_ghg_attribution, and extract_health_impacts into a single pass over the
solved network. Intermediate DataFrames (food_consumption, food_group_consumption)
are passed in memory rather than written and re-read from disk.
"""

import logging
from pathlib import Path

import pandas as pd
import pypsa

from workflow.scripts.analysis.extract_baseline_deviation import (
    extract_baseline_deviation,
)
from workflow.scripts.analysis.extract_food_prices import extract_food_prices
from workflow.scripts.analysis.extract_ghg_attribution import (
    add_monetary_value as add_ghg_monetary_value,
)
from workflow.scripts.analysis.extract_ghg_attribution import (
    compute_bus_intensities,
    compute_ghg_totals,
    join_intensities_to_consumption,
)
from workflow.scripts.analysis.extract_health_impacts import (
    add_monetary_value as add_health_monetary_value,
)
from workflow.scripts.analysis.extract_health_impacts import (
    compute_health_attribution,
    compute_health_marginals,
    extract_yll_totals,
    load_health_data,
)
from workflow.scripts.analysis.extract_net_emissions import extract_net_emissions
from workflow.scripts.analysis.extract_objective_breakdown import (
    extract_objective_breakdown,
)
from workflow.scripts.analysis.extract_statistics import (
    extract_animal_production,
    extract_consumption_tables,
    extract_crop_production,
    extract_feed_by_animal,
    extract_feed_by_category,
    extract_feed_by_source,
    extract_land_use,
    extract_luc_breakdown,
)
from workflow.scripts.analysis.extract_water_metrics import extract_water_by_region
from workflow.scripts.logging_config import setup_script_logging
from workflow.scripts.snakemake_utils import load_solved_network
from workflow.scripts.solve_namespace import ANALYSIS_OUTPUT_NAMES


def analysis_output_paths(output) -> dict[str, str]:
    """Return the canonical analysis output paths from a namespace."""
    return {name: str(getattr(output, name)) for name in ANALYSIS_OUTPUT_NAMES}


def write_empty_outputs(output) -> None:
    """Write empty Parquet files for all canonical analysis outputs."""
    for path in analysis_output_paths(output).values():
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame().to_parquet(p)


def run_analysis(
    n: pypsa.Network,
    *,
    output_paths: dict[str, str],
    food_groups_path: str,
    m49_codes_path: str,
    population_path: str,
    ghg_price: float,
    ch4_gwp: float,
    n2o_gwp: float,
    value_per_yll: float,
    health_risk_factors: list[str],
    logger: logging.Logger,
    health_enabled: bool,
    risk_breakpoints_path: str | None = None,
    health_cluster_cause_path: str | None = None,
    health_cause_log_path: str | None = None,
    health_clusters_path: str | None = None,
    tmrel_path: str | None = None,
) -> None:
    """Run all analysis extractions on a solved network.

    Parameters
    ----------
    n
        Solved PyPSA network with solution assigned.
    output_paths
        Mapping of output name to file path (e.g. ``{"crop_production": "/.../crop_production.parquet"}``).
    food_groups_path, m49_codes_path, ...
        Paths to input data files needed by analysis.
    ghg_price, ch4_gwp, n2o_gwp, value_per_yll
        Scalar parameters for emissions and health valuation.
    health_risk_factors
        List of risk factor names for health impact computation.
    logger
        Logger instance.
    """
    # --- Statistics ---
    logger.info("Extracting statistics...")
    crop_production = extract_crop_production(n)
    land_use = extract_land_use(n)
    animal_production = extract_animal_production(n)
    food_consumption, food_group_consumption = extract_consumption_tables(n)
    feed_by_category = extract_feed_by_category(n)
    feed_by_animal = extract_feed_by_animal(n)
    feed_by_source = extract_feed_by_source(n)

    # --- Food prices ---
    logger.info("Extracting food prices...")
    food_prices = extract_food_prices(n)

    # --- Baseline deviation ---
    logger.info("Extracting baseline deviation...")
    baseline_deviation = extract_baseline_deviation(n)

    # --- LUC breakdown ---
    logger.info("Extracting LUC breakdown...")
    m49 = pd.read_csv(m49_codes_path, sep=";", comment="#")
    country_to_continent = {}
    for _, row in m49.iterrows():
        iso3 = row.get("ISO-alpha3 Code")
        region = row.get("Region Name")
        if pd.notna(iso3) and pd.notna(region) and str(iso3).strip():
            country_to_continent[str(iso3).strip()] = str(region).strip()
    luc_breakdown = extract_luc_breakdown(n, country_to_continent)

    # --- Net emissions ---
    logger.info("Extracting net emissions...")
    net_emissions = extract_net_emissions(n, ch4_gwp, n2o_gwp)

    # --- Objective breakdown ---
    logger.info("Extracting objective breakdown...")
    objective_breakdown = extract_objective_breakdown(n)

    # --- Water metrics (scarcity / withdrawal / groundwater per region) ---
    logger.info("Extracting water metrics...")
    water_metrics = extract_water_by_region(n)

    # --- GHG attribution ---
    logger.info("Computing GHG attribution...")
    food_groups = pd.read_csv(food_groups_path)
    bus_intensities = compute_bus_intensities(n, ch4_gwp, n2o_gwp)
    ghg_attribution = join_intensities_to_consumption(
        food_consumption, food_groups, bus_intensities
    )
    ghg_attribution = add_ghg_monetary_value(ghg_attribution, ghg_price)
    ghg_attribution = ghg_attribution.sort_values(["country", "food"]).reset_index(
        drop=True
    )
    ghg_attribution_totals = compute_ghg_totals(ghg_attribution)

    # --- Health impacts ---
    # Only computed when health is enabled for this scenario; otherwise the
    # network carries no YLL stores and the GBD-derived health processing
    # inputs are absent, so the health outputs are written empty.
    if health_enabled:
        logger.info("Extracting health impacts...")
        risk_factors = list(health_risk_factors)
        health_data = load_health_data(
            {
                "risk_breakpoints": risk_breakpoints_path,
                "health_cluster_cause": health_cluster_cause_path,
                "health_cause_log": health_cause_log_path,
                "health_clusters": health_clusters_path,
                "population": population_path,
                "tmrel": tmrel_path,
            }
        )
        health_marginals = compute_health_marginals(n, health_data, risk_factors)
        health_marginals = add_health_monetary_value(health_marginals, value_per_yll)
        health_marginals = health_marginals.sort_values(
            ["country", "food_group"]
        ).reset_index(drop=True)
        health_totals = extract_yll_totals(n)
        health_attribution = compute_health_attribution(n, health_data, risk_factors)
    else:
        logger.info("Health disabled for this scenario; writing empty health outputs.")
        health_marginals = pd.DataFrame()
        health_totals = pd.DataFrame()
        health_attribution = pd.DataFrame()

    # Write all outputs
    output_dir = Path(output_paths["crop_production"]).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {
        "crop_production": crop_production,
        "land_use": land_use,
        "animal_production": animal_production,
        "food_consumption": food_consumption,
        "food_group_consumption": food_group_consumption,
        "net_emissions": net_emissions,
        "objective_breakdown": objective_breakdown,
        "ghg_attribution": ghg_attribution,
        "ghg_attribution_totals": ghg_attribution_totals,
        "health_marginals": health_marginals,
        "health_totals": health_totals,
        "health_attribution": health_attribution,
        "feed_by_category": feed_by_category,
        "feed_by_animal": feed_by_animal,
        "feed_by_source": feed_by_source,
        "luc_breakdown": luc_breakdown,
        "baseline_deviation": baseline_deviation,
        "food_prices": food_prices,
        "water_metrics": water_metrics,
    }
    # Producer-side drift guard: every output declared in the
    # canonical list must be produced here, and we must not produce
    # extras the downstream consumers don't know about.
    if set(results) != set(ANALYSIS_OUTPUT_NAMES):
        missing = sorted(set(ANALYSIS_OUTPUT_NAMES) - set(results))
        extra = sorted(set(results) - set(ANALYSIS_OUTPUT_NAMES))
        raise RuntimeError(
            "analyze_model results mismatch ANALYSIS_OUTPUT_NAMES "
            f"(missing={missing}, extra={extra}). Update the canonical "
            "list in workflow/scripts/solve_namespace.py."
        )
    for name, df in results.items():
        df.to_parquet(output_paths[name])

    logger.info("Wrote all analysis outputs to %s", output_dir)


def run_analysis_from_namespace(smk, n: pypsa.Network, logger: logging.Logger) -> None:
    """Run analysis using the shared Snakemake/manifest namespace contract."""
    health_enabled = bool(smk.params.health_enabled)
    run_analysis(
        n,
        output_paths=analysis_output_paths(smk.output),
        food_groups_path=smk.input.food_groups,
        m49_codes_path=smk.input.m49,
        population_path=smk.input.population,
        ghg_price=float(smk.params.ghg_price),
        ch4_gwp=float(smk.params.ch4_gwp),
        n2o_gwp=float(smk.params.n2o_gwp),
        value_per_yll=float(smk.params.health_value_per_yll),
        health_risk_factors=list(smk.params.health_risk_factors),
        logger=logger,
        health_enabled=health_enabled,
        risk_breakpoints_path=(
            smk.input.health_risk_breakpoints if health_enabled else None
        ),
        health_cluster_cause_path=(
            smk.input.health_cluster_cause if health_enabled else None
        ),
        health_cause_log_path=(smk.input.health_cause_log if health_enabled else None),
        health_clusters_path=smk.input.health_clusters if health_enabled else None,
        tmrel_path=smk.input.health_tmrel if health_enabled else None,
    )


def main() -> None:
    """Snakemake entry point: load solved network and run analysis."""
    logger = setup_script_logging(snakemake.log[0])

    network_path = Path(snakemake.input.network)
    if network_path.stat().st_size == 0:
        logger.warning(
            "Network file is empty (solve failed or timed out). Writing empty outputs."
        )
        write_empty_outputs(snakemake.output)
        return

    try:
        n = load_solved_network(snakemake.input.network)
    except Exception as e:
        logger.warning(
            "Failed to load network (%s) — likely an unsolved model. "
            "Writing empty outputs.",
            e,
        )
        write_empty_outputs(snakemake.output)
        return

    if n.links.empty:
        logger.warning(
            "Network has no links — likely an unsolved model. Writing empty outputs."
        )
        write_empty_outputs(snakemake.output)
        return

    logger.info("Loaded network with %d links", len(n.links))

    run_analysis_from_namespace(snakemake, n, logger)


if __name__ == "__main__":
    main()
