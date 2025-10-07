# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Pre-compute health data for SOS2 linearisation in the solver."""

import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import geopandas as gpd
import numpy as np
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

_AGE_BUCKET_PATTERN = re.compile(r"^(?P<start>\d+)\s*-\s*(?P<end>\d+)$")


def _load_csv(path: str, names: List[str]) -> pd.DataFrame:
    return pd.read_csv(path, header=None, names=names)


def _normalize_age_bucket(label: object) -> str | None:
    text = str(label).strip().lower()
    if text in {"<1", "under age 1", "under 1", "0", "0-0", "0 – 0"}:
        return "<1"
    if text in {"1-4", "1 – 4", "01-04", "01–04", "1 to 4"}:
        return "1-4"

    match = _AGE_BUCKET_PATTERN.match(text.replace("–", "-").replace("to", "-"))
    if match:
        start = int(match.group("start"))
        end = int(match.group("end"))
        if start == 0:
            return "<1"
        if start == 1 and end in {4, 5}:
            return "1-4"
        if start >= 5 and end == start + 4 and start <= 90:
            return f"{start}-{end}"
        if start >= 95:
            return "95+"

    if text in {"95-99", "95 – 99", "100+", "95+", "95 plus", "100 plus", "100+ years"}:
        return "95+"

    return None


def _load_wpp_life_expectancy(path: str, reference_year: int) -> pd.Series:
    df = pd.read_csv(path, low_memory=False)
    if df.empty:
        raise ValueError("WPP life table file is empty")

    variant_col = df["Variant"].astype(str).str.lower()
    df = df[variant_col == "medium"]
    if df.empty:
        raise ValueError("WPP life table missing 'Medium' variant entries")

    sex_col = df["Sex"].astype(str).str.lower()
    df = df[sex_col == "total"]
    if df.empty:
        raise ValueError("WPP life table missing 'Total' sex entries")

    df["Time"] = pd.to_numeric(df["Time"], errors="coerce")
    df = df.dropna(subset=["Time"])
    if df.empty:
        raise ValueError("WPP life table missing valid 'Time' values")

    available_years = sorted({int(value) for value in df["Time"].unique()})
    if not available_years:
        raise ValueError("WPP life table has no available years")

    if reference_year in available_years:
        target_year = reference_year
    else:
        target_year = min(
            available_years, key=lambda year: (abs(year - reference_year), year)
        )
        print(
            "[prepare_health_costs] WPP life table missing year"
            f" {reference_year}; using {target_year} instead."
        )

    df = df[df["Time"].astype(int) == int(target_year)]
    if df.empty:
        raise ValueError(f"WPP life table has no data for year {target_year}")

    df = df[df["Location"].astype(str) == "World"]
    if df.empty:
        raise ValueError("WPP life table missing 'World' aggregate records")

    age_to_life_exp: Dict[str, float] = {}
    for _, row in df.iterrows():
        bucket = _normalize_age_bucket(row.get("AgeGrp"))
        if bucket is None or bucket in age_to_life_exp:
            continue
        try:
            age_to_life_exp[bucket] = float(row["ex"])
        except (TypeError, ValueError):
            continue

    if "95+" not in age_to_life_exp:
        candidates = df[df["AgeGrp"].astype(str).isin(["95-99", "95+", "100+"])]
        if not candidates.empty:
            first = candidates.iloc[0]
            try:
                age_to_life_exp["95+"] = float(first["ex"])
            except (TypeError, ValueError):
                pass

    missing = [bucket for bucket in AGE_BUCKETS if bucket not in age_to_life_exp]
    if missing:
        raise ValueError(
            "WPP life table missing life expectancy entries for age buckets: "
            + ", ".join(missing)
        )

    ordered = {bucket: age_to_life_exp[bucket] for bucket in AGE_BUCKETS}
    series = pd.Series(ordered, name="life_exp")
    series.index.name = "age"
    return series


def _build_country_clusters(
    regions_path: str,
    countries: Iterable[str],
    n_clusters: int,
) -> Tuple[pd.Series, Dict[int, List[str]]]:
    regions = gpd.read_file(regions_path)
    regions = regions[regions["country"].isin(countries)]
    if regions.empty:
        raise ValueError("No regions match configured countries for health clustering")

    regions_equal_area = regions.to_crs(6933)
    dissolved = regions_equal_area.dissolve(by="country", as_index=True)
    centroids = dissolved.geometry.centroid

    coords = np.column_stack([centroids.x.values, centroids.y.values])
    k = max(1, min(int(n_clusters), len(coords)))
    if k < int(n_clusters):
        print(
            f"[prepare_health_costs] Requested {n_clusters} clusters but only {len(coords)} countries available; using {k}."
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


def _build_rr_lookup(rr_int: pd.DataFrame) -> Dict[Tuple[str, str], pd.Series]:
    lookup: Dict[Tuple[str, str], pd.Series] = {}
    for (risk, cause), grp in rr_int.groupby(["risk_factor", "cause"]):
        grp = grp.sort_values("intake")
        series = grp.set_index("intake")["rr_mean"].astype(float)
        full_index = pd.Index(
            range(int(series.index.min()), int(series.index.max()) + 1)
        )
        series = series.reindex(full_index).interpolate(limit_direction="both")
        lookup[(risk, cause)] = series
    return lookup


def _evaluate_rr(
    lookup: Dict[Tuple[str, str], pd.Series], risk: str, cause: str, intake: float
) -> float:
    series = lookup[(risk, cause)]
    level = int(round(float(intake)))
    level = max(series.index.min(), min(series.index.max(), level))
    value = float(series.loc[level])
    return max(value, 1e-9)


def main() -> None:
    snakemake = globals().get("snakemake")  # type: ignore
    if snakemake is None:
        raise RuntimeError("This script must run via Snakemake")

    cfg_countries: List[str] = list(snakemake.params["countries"])
    health_cfg = snakemake.params["health"]
    risk_factors: List[str] = list(health_cfg["risk_factors"])
    region_clusters = int(health_cfg["region_clusters"])
    reference_year = int(health_cfg["reference_year"])
    intake_step = float(health_cfg["intake_grid_step"])
    log_rr_points = int(health_cfg.get("log_rr_points", 40))
    vosl_setting = health_cfg.get("value_of_statistical_life", "regional")

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

    cluster_series, cluster_to_countries = _build_country_clusters(
        snakemake.input["regions"], cfg_countries, region_clusters
    )

    cluster_map = cluster_series.rename("health_cluster").reset_index()
    cluster_map.columns = ["country_iso3", "health_cluster"]
    cluster_map = cluster_map.sort_values("country_iso3")

    diet = _load_csv(
        snakemake.input["diet"],
        ["scenario", "unit", "item", "country", "year", "value"],
    )
    rr_int = _load_csv(
        snakemake.input["rr_int"],
        ["risk_factor", "cause", "intake", "stat", "value"],
    )
    rr_max = _load_csv(snakemake.input["rr_max"], ["risk_factor", "max_intake"])
    dr = _load_csv(snakemake.input["dr"], ["age", "cause", "country", "year", "value"])
    pop = pd.read_csv(snakemake.input["population"])
    pop["value"] = pd.to_numeric(pop["value"], errors="coerce") / 1_000.0
    life_exp = _load_wpp_life_expectancy(snakemake.input["life_table"], reference_year)
    vsl = _load_csv(
        snakemake.input["vsl"],
        ["income_group", "country", "year", "gdp_scenario", "stat", "value"],
    )

    diet = diet[
        (diet["scenario"] == "BMK")
        & (diet["unit"] == "g/d_w")
        & (diet["year"] == reference_year)
        & (diet["country"].isin(cfg_countries))
    ]

    rr_int = rr_int[rr_int["stat"] == "mean"].copy()
    rr_int["intake"] = rr_int["intake"].astype(int)
    rr_int["rr_mean"] = rr_int["value"].astype(float)
    rr_lookup = _build_rr_lookup(rr_int)

    rr_max = rr_max.set_index("risk_factor")["max_intake"].astype(float)

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

    if use_constant_vsl:
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

    index = pd.Index(cfg_countries, name="country")
    vsl_stats = pd.DataFrame(index=index)
    vsl_stats["vsl"] = vsl_per_country
    vsl_stats["avg_yll_per_death"] = avg_yll_per_death.reindex(index)
    vsl_stats["population_total"] = pop_total
    vsl_stats = vsl_stats.replace({0.0: np.nan})
    vsl_stats["value_per_yll"] = vsl_stats["vsl"] / vsl_stats["avg_yll_per_death"]

    item_to_risk = {
        "whole_grains": "whole_grains",
        "legumes": "legumes",
        "soybeans": "legumes",
        "nuts_seeds": "nuts_seeds",
        "vegetables": "vegetables",
        "fruits_trop": "fruits",
        "fruits_temp": "fruits",
        "fruits_starch": "fruits",
        "beef": "red_meat",
        "lamb": "red_meat",
        "pork": "red_meat",
        "prc_meat": "prc_meat",
        "shellfish": "fish",
        "fish_freshw": "fish",
        "fish_pelag": "fish",
        "fish_demrs": "fish",
    }

    diet["risk_factor"] = diet["item"].map(item_to_risk)
    diet = diet.dropna(subset=["risk_factor"])
    intake_by_country = (
        diet.groupby(["country", "risk_factor"])["value"].sum().unstack(fill_value=0.0)
    )

    cluster_summary_rows = []
    cluster_cause_rows = []
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
            baseline_intake = max(
                0.0, min(baseline_intake, float(rr_max.get(risk, baseline_intake)))
            )
            baseline_intake_registry.setdefault(risk, set()).add(baseline_intake)

            causes = risk_to_causes.get(risk, [])
            for cause in causes:
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

    risk_breakpoint_rows = []
    cause_log_min: Dict[str, float] = {cause: 0.0 for cause in relevant_causes}
    cause_log_max: Dict[str, float] = {cause: 0.0 for cause in relevant_causes}
    for cause in relevant_causes:
        cause_log_min[cause] = 0.0
        cause_log_max[cause] = 0.0

    for risk in risk_factors:
        max_intake = float(rr_max.get(risk, 0.0))
        if max_intake <= 0:
            continue
        grid_points = set(np.arange(0.0, max_intake + intake_step, intake_step))
        grid_points.add(0.0)
        grid_points.add(max_intake)
        for val in baseline_intake_registry.get(risk, set()):
            grid_points.add(float(val))
        grid = sorted(grid_points)

        causes = risk_to_causes.get(risk, [])
        for cause in causes:
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


if __name__ == "__main__":
    main()
