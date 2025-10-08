# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Pre-compute health data for SOS2 linearisation in the solver."""

import logging
import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import geopandas as gpd
import numpy as np
import numpy.typing as npt
import pandas as pd
from sklearn.cluster import KMeans

AGE_BUCKETS = [
    "<1",
    "1-4",
    "5-9",
    "10-14",
    "15-19",
    "20-24",
    "25-29",
    "30-34",
    "35-39",
    "40-44",
    "45-49",
    "50-54",
    "55-59",
    "60-64",
    "65-69",
    "70-74",
    "75-79",
    "80-84",
    "85-89",
    "90-94",
    "95+",
]


logger = logging.getLogger(__name__)


def _load_life_expectancy(path: str) -> pd.Series:
    """Load processed life expectancy data from prepare_life_table.py output."""
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError("Life table file is empty")

    required_cols = {"age", "life_exp"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Life table missing required columns: {required_cols}")

    # Validate all expected age buckets are present
    missing = [bucket for bucket in AGE_BUCKETS if bucket not in df["age"].values]
    if missing:
        raise ValueError(
            "Life table missing life expectancy entries for age buckets: "
            + ", ".join(missing)
        )

    series = df.set_index("age")["life_exp"]
    series.name = "life_exp"
    return series


def _build_country_clusters(
    regions_path: str,
    countries: Iterable[str],
    n_clusters: int,
) -> Tuple[pd.Series, Dict[int, List[str]]]:
    regions = gpd.read_file(regions_path)

    regions_equal_area = regions.to_crs(6933)
    dissolved = regions_equal_area.dissolve(by="country", as_index=True)
    centroids = dissolved.geometry.centroid

    coords = np.column_stack([centroids.x.values, centroids.y.values])
    k = max(1, min(int(n_clusters), len(coords)))
    if k < int(n_clusters):
        logger.info(
            f"Requested {n_clusters} clusters but only {len(coords)} countries available; using {k}."
        )

    if len(coords) == 1:
        labels = np.array([0])
    else:
        km = KMeans(n_clusters=k, n_init=20, random_state=0)
        labels = km.fit_predict(coords)

    dissolved["health_cluster"] = labels
    cluster_series = dissolved["health_cluster"].astype(int)
    grouped = cluster_series.groupby(cluster_series).groups
    cluster_to_countries = {
        int(cluster): sorted(list(indexes)) for cluster, indexes in grouped.items()
    }
    return cluster_series, cluster_to_countries


class RelativeRiskTable(Dict[Tuple[str, str], Dict[str, np.ndarray]]):
    """Container mapping (risk, cause) to exposure grids and log RR values."""


def _build_rr_tables(
    rr_df: pd.DataFrame, risk_factors: Iterable[str]
) -> Tuple[RelativeRiskTable, Dict[str, float]]:
    """Build lookup tables for relative risk curves by (risk, cause) pairs.

    Returns:
        table: Dict mapping (risk, cause) to exposure arrays and log(RR) values
        max_exposure_g_per_day: Dict mapping risk factor to maximum exposure level in data
    """
    table: RelativeRiskTable = RelativeRiskTable()
    max_exposure_g_per_day: Dict[str, float] = {}
    allowed = set(risk_factors)

    for (risk, cause), grp in rr_df.groupby(["risk_factor", "cause"], sort=True):
        if risk not in allowed:
            continue

        grp = grp.sort_values("exposure_g_per_day")
        exposures = grp["exposure_g_per_day"].astype(float).values
        if len(exposures) == 0:
            continue
        log_rr_mean = np.log(grp["rr_mean"].astype(float).values)

        table[(risk, cause)] = {
            "exposures": exposures,
            "log_rr_mean": log_rr_mean,
        }
        max_exposure_g_per_day[risk] = max(
            max_exposure_g_per_day.get(risk, 0.0), float(exposures.max())
        )

    return table, max_exposure_g_per_day


def _evaluate_rr(
    table: RelativeRiskTable, risk: str, cause: str, intake: float
) -> float:
    """Interpolate relative risk for given intake using log-linear interpolation."""
    data = table[(risk, cause)]
    exposures: npt.NDArray[np.floating] = data["exposures"]
    log_rr: npt.NDArray[np.floating] = data["log_rr_mean"]

    if intake <= exposures[0]:
        return float(math.exp(log_rr[0]))
    if intake >= exposures[-1]:
        return float(math.exp(log_rr[-1]))

    log_val = float(np.interp(intake, exposures, log_rr))
    return float(math.exp(log_val))


def _load_input_data(
    snakemake,
    cfg_countries: List[str],
    reference_year: int,
    use_constant_vsl: bool,
) -> tuple:
    """Load and perform initial processing of all input datasets."""
    cluster_series, cluster_to_countries = _build_country_clusters(
        snakemake.input["regions"],
        cfg_countries,
        int(snakemake.params["health"]["region_clusters"]),
    )

    cluster_map = cluster_series.rename("health_cluster").reset_index()
    cluster_map.columns = ["country_iso3", "health_cluster"]
    cluster_map = cluster_map.sort_values("country_iso3")

    diet = pd.read_csv(snakemake.input["diet"])
    rr_df = pd.read_csv(snakemake.input["relative_risks"])
    dr = pd.read_csv(
        snakemake.input["dr"],
        header=None,
        names=["age", "cause", "country", "year", "value"],
    )
    pop = pd.read_csv(snakemake.input["population"])
    pop["value"] = pd.to_numeric(pop["value"], errors="coerce") / 1_000.0
    life_exp = _load_life_expectancy(snakemake.input["life_table"])

    vsl_input = snakemake.input.get("vsl") if hasattr(snakemake.input, "get") else None
    vsl: pd.DataFrame | None = None
    if not use_constant_vsl and vsl_input:
        vsl = pd.read_csv(
            vsl_input,
            header=None,
            names=["income_group", "country", "year", "gdp_scenario", "stat", "value"],
        )

    return (
        cluster_series,
        cluster_to_countries,
        cluster_map,
        diet,
        rr_df,
        dr,
        pop,
        life_exp,
        vsl,
    )


def _filter_and_prepare_data(
    diet: pd.DataFrame,
    dr: pd.DataFrame,
    pop: pd.DataFrame,
    rr_df: pd.DataFrame,
    cfg_countries: List[str],
    reference_year: int,
    life_exp: pd.Series,
    risk_factors: List[str],
) -> tuple:
    """Filter datasets to reference year and compute derived quantities."""
    # Filter dietary intake data
    diet = diet[
        (diet["scenario"] == "BMK")
        & (diet["unit"] == "g/d_w")
        & (diet["year"] == reference_year)
        & (diet["country"].isin(cfg_countries))
    ]

    # Build relative risk lookup tables
    rr_lookup, max_exposure_g_per_day = _build_rr_tables(rr_df, risk_factors)

    # Filter mortality and population data
    dr = dr[(dr["year"] == reference_year) & (dr["country"].isin(cfg_countries))].copy()
    pop = pop[
        (pop["year"] == reference_year) & (pop["country"].isin(cfg_countries))
    ].copy()

    valid_ages = life_exp.index
    dr = dr[dr["age"].isin(valid_ages)].copy()
    pop_age = pop[pop["age"].isin(valid_ages)].copy()

    pop_total = (
        pop[pop["age"] == "all-a"]
        .groupby("country")["value"]
        .sum()
        .astype(float)
        .reindex(cfg_countries)
    )

    # Determine relevant risk-cause pairs
    relevant_pairs = {
        (risk, cause) for (risk, cause) in rr_lookup.keys() if risk in risk_factors
    }
    relevant_causes = sorted({cause for _, cause in relevant_pairs})
    risk_to_causes: Dict[str, List[str]] = {}
    for risk, cause in relevant_pairs:
        risk_to_causes.setdefault(risk, set()).add(cause)
    risk_to_causes = {
        risk: sorted(list(causes)) for risk, causes in risk_to_causes.items()
    }

    dr = dr[dr["cause"].isin(relevant_causes)].copy()

    # Map diet items to risk factors
    item_to_risk = {
        "whole_grains": "whole_grains",
        "legumes": "legumes",
        "soybeans": "legumes",
        "nuts_seeds": "nuts_seeds",
        "vegetables": "vegetables",
        "fruits_trop": "fruits",
        "fruits_temp": "fruits",
        "fruits_starch": "fruits",
        "fruits": "fruits",
        "beef": "red_meat",
        "lamb": "red_meat",
        "pork": "red_meat",
        "red_meat": "red_meat",
        "prc_meat": "prc_meat",
        "shellfish": "fish",
        "fish_freshw": "fish",
        "fish_pelag": "fish",
        "fish_demrs": "fish",
        "fish": "fish",
    }
    diet["risk_factor"] = diet["item"].map(item_to_risk)
    diet = diet.dropna(subset=["risk_factor"])
    intake_by_country = (
        diet.groupby(["country", "risk_factor"])["value"].sum().unstack(fill_value=0.0)
    )

    return (
        dr,
        pop_age,
        pop_total,
        rr_lookup,
        max_exposure_g_per_day,
        relevant_causes,
        risk_to_causes,
        intake_by_country,
    )


def _compute_baseline_health_metrics(
    dr: pd.DataFrame,
    pop_age: pd.DataFrame,
    life_exp: pd.Series,
    vsl: pd.DataFrame | None,
    cfg_countries: List[str],
    reference_year: int,
    use_constant_vsl: bool,
    constant_vsl_value: float | None,
) -> pd.DataFrame:
    """Compute YLL and value-of-statistical-life statistics by country."""
    pop_age = pop_age.rename(columns={"value": "population"})
    dr = dr.rename(columns={"value": "death_rate"})
    combo = dr.merge(pop_age, on=["age", "country", "year"], how="left").merge(
        life_exp.rename("life_exp"), left_on="age", right_index=True, how="left"
    )
    combo["population"] = combo["population"].fillna(0.0)
    combo["death_rate"] = combo["death_rate"].fillna(0.0)
    combo["death_count"] = combo["death_rate"] * combo["population"]
    combo["yll"] = combo["death_count"] * combo["life_exp"]

    deaths_by_country = combo.groupby("country")["death_count"].sum()
    yll_by_country = combo.groupby("country")["yll"].sum()
    avg_yll_per_death = (
        yll_by_country / deaths_by_country.replace({0.0: np.nan})
    ).replace([np.inf, -np.inf], np.nan)

    if not use_constant_vsl and vsl is None:
        raise ValueError(
            "Regional VSL dataset is not available; set health.value_of_statistical_life"
            " to a numeric constant or reintroduce a regional dataset."
        )

    if use_constant_vsl or vsl is None:
        vsl_per_country = pd.Series(constant_vsl_value, index=cfg_countries)
    else:
        vsl_filtered = vsl[
            (vsl["year"].astype(int) == reference_year)
            & (vsl["gdp_scenario"] == "mean")
            & (vsl["stat"] == "mean")
            & (vsl["country"].isin(cfg_countries))
        ].copy()
        vsl_filtered["value"] = vsl_filtered["value"].astype(float) * 1e6
        vsl_per_country = (
            vsl_filtered.groupby("country")["value"].mean().reindex(cfg_countries)
        )

    pop_total = combo.groupby("country")["population"].sum()
    index = pd.Index(cfg_countries, name="country")
    vsl_stats = pd.DataFrame(index=index)
    vsl_stats["vsl"] = vsl_per_country
    vsl_stats["avg_yll_per_death"] = avg_yll_per_death.reindex(index)
    vsl_stats["population_total"] = pop_total.reindex(index)
    vsl_stats = vsl_stats.replace({0.0: np.nan})
    vsl_stats["value_per_yll"] = vsl_stats["vsl"] / vsl_stats["avg_yll_per_death"]

    return vsl_stats, combo


def _process_health_clusters(
    cluster_to_countries: Dict[int, List[str]],
    pop_total: pd.Series,
    combo: pd.DataFrame,
    vsl_stats: pd.DataFrame,
    risk_factors: List[str],
    intake_by_country: pd.DataFrame,
    max_exposure_g_per_day: Dict[str, float],
    rr_lookup: RelativeRiskTable,
    risk_to_causes: Dict[str, List[str]],
    relevant_causes: List[str],
) -> tuple:
    """Process each health cluster to compute baseline metrics and intakes."""
    cluster_summary_rows = []
    cluster_cause_rows = []
    cluster_risk_baseline_rows = []
    baseline_intake_registry: Dict[str, set] = {risk: set() for risk in risk_factors}

    for cluster_id, members in cluster_to_countries.items():
        pop_weights = pop_total.reindex(members).fillna(0.0)
        total_pop_thousand = float(pop_weights.sum())
        if total_pop_thousand <= 0:
            continue

        total_population_persons = total_pop_thousand * 1_000.0
        cluster_combo = combo[combo["country"].isin(members)]
        yll_by_cause_cluster = cluster_combo.groupby("cause")["yll"].sum()

        vsl_subset = vsl_stats.reindex(members)
        weight_sum = vsl_subset["population_total"].sum()
        if weight_sum > 0:
            value_per_yll = (
                vsl_subset["value_per_yll"] * vsl_subset["population_total"]
            ).sum() / weight_sum
        else:
            value_per_yll = float("nan")

        cluster_summary_rows.append(
            {
                "health_cluster": int(cluster_id),
                "population_persons": total_population_persons,
                "value_per_yll_usd_per_yll": value_per_yll,
            }
        )

        log_rr_ref_totals: Dict[str, float] = {cause: 0.0 for cause in relevant_causes}

        for risk in risk_factors:
            if risk not in intake_by_country.columns:
                baseline_intake = 0.0
            else:
                baseline_intake = (
                    intake_by_country[risk].reindex(members).fillna(0.0) * pop_weights
                ).sum() / total_pop_thousand
            baseline_intake = float(baseline_intake)
            if not math.isfinite(baseline_intake):
                baseline_intake = 0.0
            max_exposure = float(max_exposure_g_per_day.get(risk, baseline_intake))
            baseline_intake = max(0.0, min(baseline_intake, max_exposure))
            baseline_intake_registry.setdefault(risk, set()).add(baseline_intake)

            cluster_risk_baseline_rows.append(
                {
                    "health_cluster": int(cluster_id),
                    "risk_factor": risk,
                    "baseline_intake_g_per_day": baseline_intake,
                }
            )

            causes = risk_to_causes.get(risk, [])
            for cause in causes:
                if (risk, cause) not in rr_lookup:
                    continue
                rr_val = _evaluate_rr(rr_lookup, risk, cause, baseline_intake)
                log_rr = math.log(rr_val)
                log_rr_ref_totals[cause] = log_rr_ref_totals.get(cause, 0.0) + log_rr

        for cause in relevant_causes:
            cluster_cause_rows.append(
                {
                    "health_cluster": int(cluster_id),
                    "cause": cause,
                    "log_rr_total_ref": log_rr_ref_totals.get(cause, 0.0),
                    "yll_base": yll_by_cause_cluster.get(cause, 0.0),
                }
            )

    cluster_summary = pd.DataFrame(
        cluster_summary_rows,
        columns=["health_cluster", "population_persons", "value_per_yll_usd_per_yll"],
    )
    cluster_cause_baseline = pd.DataFrame(
        cluster_cause_rows,
        columns=["health_cluster", "cause", "log_rr_total_ref", "yll_base"],
    )
    cluster_risk_baseline = pd.DataFrame(
        cluster_risk_baseline_rows,
        columns=["health_cluster", "risk_factor", "baseline_intake_g_per_day"],
    )

    return (
        cluster_summary,
        cluster_cause_baseline,
        cluster_risk_baseline,
        baseline_intake_registry,
    )


def _generate_breakpoint_tables(
    risk_factors: List[str],
    max_exposure_g_per_day: Dict[str, float],
    baseline_intake_registry: Dict[str, set],
    intake_step: float,
    rr_lookup: RelativeRiskTable,
    risk_to_causes: Dict[str, List[str]],
    relevant_causes: List[str],
    log_rr_points: int,
) -> tuple:
    """Generate SOS2 linearization breakpoint tables for risks and causes."""
    risk_breakpoint_rows = []
    cause_log_min: Dict[str, float] = {cause: 0.0 for cause in relevant_causes}
    cause_log_max: Dict[str, float] = {cause: 0.0 for cause in relevant_causes}

    for risk in risk_factors:
        max_exposure = float(max_exposure_g_per_day.get(risk, 0.0))
        if max_exposure <= 0:
            continue
        grid_points = set(np.arange(0.0, max_exposure + intake_step, intake_step))
        causes = risk_to_causes.get(risk, [])
        for cause in causes:
            exposures = rr_lookup.get((risk, cause), {}).get("exposures")
            if exposures is not None:
                grid_points.update(float(x) for x in exposures)
        grid_points.add(0.0)
        grid_points.add(max_exposure)
        for val in baseline_intake_registry.get(risk, set()):
            grid_points.add(float(val))
        grid = sorted(grid_points)

        for cause in causes:
            key = (risk, cause)
            if key not in rr_lookup:
                continue
            log_values: List[float] = []
            for intake in grid:
                rr_val = _evaluate_rr(rr_lookup, risk, cause, intake)
                log_rr = math.log(rr_val)
                log_values.append(log_rr)
                risk_breakpoint_rows.append(
                    {
                        "risk_factor": risk,
                        "cause": cause,
                        "intake_g_per_day": float(intake),
                        "log_rr": log_rr,
                    }
                )
            if log_values:
                cause_log_min[cause] += min(log_values)
                cause_log_max[cause] += max(log_values)

    risk_breakpoints = pd.DataFrame(risk_breakpoint_rows)

    cause_breakpoint_rows = []
    for cause in relevant_causes:
        min_total = cause_log_min.get(cause)
        max_total = cause_log_max.get(cause)
        if min_total is None or max_total is None:
            continue
        if not math.isfinite(min_total):
            min_total = 0.0
        if not math.isfinite(max_total):
            max_total = 0.0
        if max_total < min_total:
            min_total, max_total = max_total, min_total
        if abs(max_total - min_total) < 1e-6:
            log_vals = np.array([min_total])
        else:
            log_vals = np.linspace(min_total, max_total, max(log_rr_points, 2))
        for log_val in log_vals:
            cause_breakpoint_rows.append(
                {
                    "cause": cause,
                    "log_rr_total": float(log_val),
                    "rr_total": math.exp(float(log_val)),
                }
            )

    cause_log_breakpoints = pd.DataFrame(cause_breakpoint_rows)

    return risk_breakpoints, cause_log_breakpoints


def main() -> None:
    """Main entry point for health cost preparation."""
    cfg_countries: List[str] = list(snakemake.params["countries"])
    health_cfg = snakemake.params["health"]
    risk_factors: List[str] = list(health_cfg["risk_factors"])
    reference_year = int(health_cfg["reference_year"])
    intake_step = float(health_cfg["intake_grid_step"])
    log_rr_points = int(health_cfg["log_rr_points"])
    vosl_setting = health_cfg["value_of_statistical_life"]

    # Parse value-of-statistical-life setting
    use_constant_vsl = False
    constant_vsl_value = None
    if isinstance(vosl_setting, str):
        if vosl_setting.lower() != "regional":
            try:
                constant_vsl_value = float(vosl_setting)
            except ValueError as exc:
                raise ValueError(
                    "health.value_of_statistical_life must be 'regional' or a number"
                ) from exc
            else:
                use_constant_vsl = True
    elif vosl_setting is not None:
        constant_vsl_value = float(vosl_setting)
        use_constant_vsl = True

    # Load input data
    (
        cluster_series,
        cluster_to_countries,
        cluster_map,
        diet,
        rr_df,
        dr,
        pop,
        life_exp,
        vsl,
    ) = _load_input_data(snakemake, cfg_countries, reference_year, use_constant_vsl)

    # Filter and prepare datasets
    (
        dr,
        pop_age,
        pop_total,
        rr_lookup,
        max_exposure_g_per_day,
        relevant_causes,
        risk_to_causes,
        intake_by_country,
    ) = _filter_and_prepare_data(
        diet, dr, pop, rr_df, cfg_countries, reference_year, life_exp, risk_factors
    )

    # Compute baseline health metrics
    vsl_stats, combo = _compute_baseline_health_metrics(
        dr,
        pop_age,
        life_exp,
        vsl,
        cfg_countries,
        reference_year,
        use_constant_vsl,
        constant_vsl_value,
    )

    # Process health clusters
    (
        cluster_summary,
        cluster_cause_baseline,
        cluster_risk_baseline,
        baseline_intake_registry,
    ) = _process_health_clusters(
        cluster_to_countries,
        pop_total,
        combo,
        vsl_stats,
        risk_factors,
        intake_by_country,
        max_exposure_g_per_day,
        rr_lookup,
        risk_to_causes,
        relevant_causes,
    )

    # Generate breakpoint tables for SOS2 linearization
    risk_breakpoints, cause_log_breakpoints = _generate_breakpoint_tables(
        risk_factors,
        max_exposure_g_per_day,
        baseline_intake_registry,
        intake_step,
        rr_lookup,
        risk_to_causes,
        relevant_causes,
        log_rr_points,
    )

    # Write outputs
    output_dir = Path(snakemake.output["risk_breakpoints"]).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    risk_breakpoints.sort_values(["risk_factor", "cause", "intake_g_per_day"]).to_csv(
        snakemake.output["risk_breakpoints"], index=False
    )
    cluster_cause_baseline.sort_values(["health_cluster", "cause"]).to_csv(
        snakemake.output["cluster_cause"], index=False
    )
    cause_log_breakpoints.sort_values(["cause", "log_rr_total"]).to_csv(
        snakemake.output["cause_log"], index=False
    )
    cluster_summary.sort_values("health_cluster").to_csv(
        snakemake.output["cluster_summary"], index=False
    )
    cluster_map.to_csv(snakemake.output["clusters"], index=False)
    cluster_risk_baseline.sort_values(["health_cluster", "risk_factor"]).to_csv(
        snakemake.output["cluster_risk_baseline"], index=False
    )


if __name__ == "__main__":
    main()
