# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
import math
from collections import defaultdict
from typing import Dict

import numpy as np
import pandas as pd
import pypsa
import xarray as xr

try:  # Used for type annotations / documentation; fallback when unavailable
    import linopy  # type: ignore
except Exception:  # pragma: no cover - documentation build without linopy

    class linopy:  # type: ignore
        class Model:  # Minimal stub to satisfy type checkers and autodoc
            pass


from pypsa._options import options
from pypsa.optimization.optimize import _optimize_guard


HEALTH_AUX_MAP: dict[int, set[str]] = {}


def _register_health_variable(m: "linopy.Model", name: str) -> None:
    aux = HEALTH_AUX_MAP.setdefault(id(m), set())
    aux.add(name)


_SOS2_COUNTER = [0]  # Use list for mutable counter


def _add_sos2_with_fallback(m, variable, sos_dim: str, solver_name: str) -> list[str]:
    """Add SOS2 or binary fallback depending on solver support."""

    if solver_name.lower() != "highs":
        m.add_sos_constraints(variable, sos_type=2, sos_dim=sos_dim)
        return []

    coords = variable.coords[sos_dim]
    n_points = len(coords)
    if n_points <= 1:
        return []

    other_dims = [dim for dim in variable.dims if dim != sos_dim]

    interval_dim = f"{sos_dim}_interval"
    suffix = 1
    while interval_dim in variable.dims:
        interval_dim = f"{sos_dim}_interval{suffix}"
        suffix += 1

    interval_index = pd.Index(range(n_points - 1), name=interval_dim)
    binary_coords = [variable.coords[d] for d in other_dims] + [interval_index]

    # Use counter instead of checking existing variables
    _SOS2_COUNTER[0] += 1
    base_name = f"{variable.name}_segment" if variable.name else "health_segment"
    binary_name = f"{base_name}_{_SOS2_COUNTER[0]}"

    binaries = m.add_variables(coords=binary_coords, binary=True, name=binary_name)

    m.add_constraints(binaries.sum(interval_dim) == 1)

    # Vectorize SOS2 constraints: variable[i] <= binary[i-1] + binary[i]
    if n_points >= 2:
        # Build adjacency matrix using numpy indexing
        adjacency_data = np.zeros((n_points, n_points - 1))

        # Fill adjacency matrix using vectorized operations
        # binary[i] affects variable[i] and variable[i+1]
        indices = np.arange(n_points - 1)
        adjacency_data[indices, indices] = 1  # binary[i] -> variable[i]
        adjacency_data[indices + 1, indices] = 1  # binary[i] -> variable[i+1]

        # Convert to DataArray with proper coordinates
        adjacency = xr.DataArray(
            adjacency_data,
            coords={sos_dim: coords, interval_dim: range(n_points - 1)},
            dims=[sos_dim, interval_dim],
        )

        # Create constraint: variables <= adjacency @ binaries
        rhs = (adjacency * binaries).sum(interval_dim)
        m.add_constraints(variable <= rhs)

    return [binary_name]


OBJECTIVE_COEFF_TARGET = 1e8

logger = logging.getLogger(__name__)


def rescale_objective(
    m: "linopy.Model", target_max_coeff: float = OBJECTIVE_COEFF_TARGET
) -> float:
    """Scale objective to keep coefficients within numerical comfort zone."""

    dataset = m.objective.data
    if "coeffs" not in dataset:
        return 1.0

    coeffs = dataset["coeffs"].values
    if coeffs.size == 0:
        return 1.0

    max_coeff = float(np.nanmax(np.abs(coeffs)))
    if not math.isfinite(max_coeff) or max_coeff == 0.0:
        return 1.0

    if max_coeff <= target_max_coeff:
        return 1.0

    exponent = math.ceil(math.log10(max_coeff / target_max_coeff))
    if exponent <= 0:
        return 1.0

    scale_factor = 10.0**exponent
    m.objective = m.objective / scale_factor
    logger.warning(
        "Scaled objective by 1/%s to reduce max coefficient from %.3g to %.3g",
        scale_factor,
        max_coeff,
        max_coeff / scale_factor,
    )
    return scale_factor


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


def add_health_objective(
    n: pypsa.Network,
    risk_breakpoints_path: str,
    cluster_cause_path: str,
    cause_log_path: str,
    cluster_summary_path: str,
    clusters_path: str,
    population_totals_path: str,
    risk_factors: list[str],
    solver_name: str,
    value_per_yll: float,
) -> None:
    """Add SOS2-based health costs with log-linear aggregation."""

    m = n.model

    risk_breakpoints = pd.read_csv(risk_breakpoints_path)
    cluster_cause = pd.read_csv(cluster_cause_path)
    cause_log_breakpoints = pd.read_csv(cause_log_path)
    cluster_summary = pd.read_csv(cluster_summary_path)
    if "health_cluster" in cluster_summary.columns:
        cluster_summary["health_cluster"] = cluster_summary["health_cluster"].astype(
            int
        )
    cluster_map = pd.read_csv(clusters_path)
    population_totals = pd.read_csv(population_totals_path)

    # Load food→risk factor mapping from food_groups.csv (only GBD risk factors)
    food_groups_df = pd.read_csv(snakemake.input.food_groups)
    food_map = food_groups_df[food_groups_df["group"].isin(risk_factors)].copy()
    food_map = food_map.rename(columns={"group": "risk_factor"})
    food_map["share"] = 1.0
    food_map["sanitized"] = food_map["food"].apply(sanitize_food_name)
    food_map = food_map.set_index("sanitized")[["risk_factor", "share"]]

    cluster_lookup = cluster_map.set_index("country_iso3")["health_cluster"].to_dict()
    cluster_population_baseline = cluster_summary.set_index("health_cluster")[
        "population_persons"
    ].to_dict()

    cluster_cause_metadata = cluster_cause.set_index(["health_cluster", "cause"])

    population_totals = population_totals.dropna(subset=["iso3", "population"]).copy()
    population_totals["iso3"] = population_totals["iso3"].astype(str).str.upper()
    population_map = population_totals.set_index("iso3")["population"].to_dict()

    cluster_population_planning: dict[int, float] = defaultdict(float)
    for iso3, cluster in cluster_lookup.items():
        value = float(population_map.get(iso3, 0.0))
        if value <= 0:
            continue
        cluster_population_planning[int(cluster)] += value

    cluster_population: dict[int, float] = {}
    for cluster, baseline_pop in cluster_population_baseline.items():
        planning_pop = cluster_population_planning.get(int(cluster), 0.0)
        if planning_pop > 0:
            cluster_population[int(cluster)] = planning_pop
        else:
            cluster_population[int(cluster)] = float(baseline_pop)

    risk_breakpoints = risk_breakpoints.sort_values(
        ["risk_factor", "intake_g_per_day", "cause"]
    )
    cause_log_breakpoints = cause_log_breakpoints.sort_values(["cause", "log_rr_total"])

    p = m.variables["Link-p"].sel(snapshot="now")

    terms_by_key: dict[tuple[int, str], list[tuple[float, object]]] = defaultdict(list)

    for link_name in p.coords["name"].values:
        if not isinstance(link_name, str) or not link_name.startswith("consume_"):
            continue
        base, _, country = link_name.rpartition("_")
        if not base.startswith("consume_") or len(country) != 3:
            continue
        sanitized_food = base[len("consume_") :]
        if sanitized_food not in food_map.index:
            continue
        risk_factor = food_map.at[sanitized_food, "risk_factor"]
        share = float(food_map.at[sanitized_food, "share"])
        if share <= 0:
            continue
        cluster = cluster_lookup.get(country)
        if cluster is None:
            continue
        population = float(cluster_population.get(cluster, 0.0))
        if population <= 0:
            continue
        coeff = share * 1_000_000.0 / (365.0 * population)
        var = p.sel(name=link_name)
        terms_by_key[(int(cluster), str(risk_factor))].append((coeff, var))

    if not terms_by_key:
        logger.info("No consumption links map to health risk factors; skipping")
        return

    risk_data = {}
    for risk, grp in risk_breakpoints.groupby("risk_factor"):
        intakes = pd.Index(sorted(grp["intake_g_per_day"].unique()), name="intake")
        if intakes.empty:
            continue
        pivot = (
            grp.pivot_table(
                index="intake_g_per_day",
                columns="cause",
                values="log_rr",
                aggfunc="first",
            )
            .reindex(intakes, axis=0)
            .sort_index()
        )
        risk_data[risk] = {"intakes": intakes, "log_rr": pivot}

    cause_breakpoint_data = {
        cause: df.sort_values("log_rr_total")
        for cause, df in cause_log_breakpoints.groupby("cause")
    }

    log_rr_totals: Dict[tuple[int, str], object] = {}

    # Group (cluster, risk) pairs by their intake coordinate patterns
    intake_groups: dict[tuple[float, ...], list[tuple[int, str]]] = defaultdict(list)
    for (cluster, risk), terms in terms_by_key.items():
        risk_table = risk_data.get(risk)
        if risk_table is None or not terms:
            continue
        coords = risk_table["intakes"]
        if len(coords) == 0:
            continue
        coords_key = tuple(coords.values)
        intake_groups[coords_key].append((cluster, risk))

    # Process each group with vectorized operations
    for coords_key, cluster_risk_pairs in intake_groups.items():
        coords = pd.Index(coords_key, name="intake")

        # Create flattened index for this group
        cluster_risk_labels = [
            f"c{cluster}_r{sanitize_identifier(risk)}"
            for cluster, risk in cluster_risk_pairs
        ]
        cluster_risk_index = pd.Index(cluster_risk_labels, name="cluster_risk")

        # Single vectorized variable creation
        lambdas_group = m.add_variables(
            lower=0,
            upper=1,
            coords=[cluster_risk_index, coords],
            name=f"health_lambda_group_{hash(coords_key)}",
        )

        # Register all variables
        _register_health_variable(m, lambdas_group.name)

        # Single SOS2 constraint call for entire group
        aux_names = _add_sos2_with_fallback(
            m, lambdas_group, sos_dim="intake", solver_name=solver_name
        )
        for aux_name in aux_names:
            _register_health_variable(m, aux_name)

        # Vectorized convexity constraints
        m.add_constraints(lambdas_group.sum("intake") == 1)

        # Process each (cluster, risk) for balance constraints and log_rr contributions
        coeff_intake = xr.DataArray(
            coords.values, coords={"intake": coords.values}, dims=["intake"]
        )

        for (cluster, risk), label in zip(cluster_risk_pairs, cluster_risk_labels):
            terms = terms_by_key[(cluster, risk)]
            lambdas = lambdas_group.sel(cluster_risk=label)

            intake_expr = m.linexpr((coeff_intake, lambdas)).sum("intake")
            flow_expr = m.linexpr(*terms)
            m.add_constraints(
                flow_expr == intake_expr,
                name=f"health_intake_balance_c{cluster}_r{sanitize_identifier(risk)}",
            )

            risk_table = risk_data[risk]
            log_rr_matrix = risk_table["log_rr"]
            for cause in log_rr_matrix.columns:
                values = log_rr_matrix[cause].to_numpy()
                coeff_log = xr.DataArray(
                    values, coords={"intake": coords.values}, dims=["intake"]
                )
                contrib = m.linexpr((coeff_log, lambdas)).sum("intake")
                key = (cluster, cause)
                if key in log_rr_totals:
                    log_rr_totals[key] = log_rr_totals[key] + contrib
                else:
                    log_rr_totals[key] = contrib

    constant_objective = 0.0
    objective_expr = None

    # Group (cluster, cause) pairs by their log_total coordinate patterns
    log_total_groups: dict[tuple[float, ...], list[tuple[int, str]]] = defaultdict(list)
    cluster_cause_data: dict[tuple[int, str], dict] = {}

    for (cluster, cause), row in cluster_cause_metadata.iterrows():
        cluster = int(cluster)
        cause = str(cause)
        yll_base = float(row.get("yll_base", 0.0))
        if not math.isfinite(value_per_yll) or value_per_yll <= 0:
            continue
        if yll_base == 0 or not math.isfinite(yll_base):
            continue

        cause_bp = cause_breakpoint_data[cause]
        coords_key = tuple(cause_bp["log_rr_total"].values)
        if len(coords_key) == 1:
            raise ValueError(
                "Need at least two breakpoints for piecewise linear approximation"
            )

        log_total_groups[coords_key].append((cluster, cause))

        # Store data for later use
        log_rr_total_ref = float(row.get("log_rr_total_ref", 0.0))
        cluster_cause_data[(cluster, cause)] = {
            "value_per_yll": value_per_yll,
            "yll_base": yll_base,
            "log_rr_total_ref": log_rr_total_ref,
            "rr_ref": math.exp(log_rr_total_ref),
            "cause_bp": cause_bp,
        }

    # Process each group with vectorized operations
    for coords_key, cluster_cause_pairs in log_total_groups.items():
        coords = pd.Index(coords_key, name="log_total")

        # Create flattened index for this group
        cluster_cause_labels = [
            f"c{cluster}_cause{sanitize_identifier(cause)}"
            for cluster, cause in cluster_cause_pairs
        ]
        cluster_cause_index = pd.Index(cluster_cause_labels, name="cluster_cause")

        # Single vectorized variable creation
        lambda_total_group = m.add_variables(
            lower=0,
            upper=1,
            coords=[cluster_cause_index, coords],
            name=f"health_lambda_total_group_{hash(coords_key)}",
        )

        # Register all variables
        _register_health_variable(m, lambda_total_group.name)

        # Single SOS2 constraint call for entire group
        aux_names = _add_sos2_with_fallback(
            m, lambda_total_group, sos_dim="log_total", solver_name=solver_name
        )
        for aux_name in aux_names:
            _register_health_variable(m, aux_name)

        # Vectorized convexity constraints
        m.add_constraints(lambda_total_group.sum("log_total") == 1)

        # Process each (cluster, cause) for balance constraints and objective
        coeff_log_total = xr.DataArray(
            coords.values,
            coords={"log_total": coords.values},
            dims=["log_total"],
        )

        for (cluster, cause), label in zip(cluster_cause_pairs, cluster_cause_labels):
            data = cluster_cause_data[(cluster, cause)]
            lambda_total = lambda_total_group.sel(cluster_cause=label)

            total_expr = log_rr_totals.get((cluster, cause), m.linexpr(0.0))
            cause_bp = data["cause_bp"]

            log_interp = m.linexpr((coeff_log_total, lambda_total)).sum("log_total")
            coeff_rr = xr.DataArray(
                cause_bp["rr_total"].values,
                coords={"log_total": coords.values},
                dims=["log_total"],
            )
            rr_interp = m.linexpr((coeff_rr, lambda_total)).sum("log_total")

            sanitized_cause = sanitize_identifier(cause)
            m.add_constraints(
                log_interp == total_expr,
                name=f"health_total_balance_c{cluster}_cause{sanitized_cause}",
            )

            coeff = data["value_per_yll"] * data["yll_base"]
            scaled_expr = rr_interp * (coeff / data["rr_ref"])
            objective_expr = (
                scaled_expr if objective_expr is None else objective_expr + scaled_expr
            )
            constant_objective -= coeff

    if objective_expr is not None:
        m.objective = m.objective + objective_expr
        logger.info("Added health cost objective")

    if constant_objective != 0.0:
        adjustments = n.meta.setdefault("objective_constant_terms", {})
        adjustments["health"] = adjustments.get("health", 0.0) + constant_objective
        logger.debug(
            "Recorded health objective constant %.3g in network metadata",
            constant_objective,
        )

    if objective_expr is None and constant_objective == 0.0:
        logger.info("No health objective terms added (missing overlaps)")


if __name__ == "__main__":
    n = pypsa.Network(snakemake.input.network)

    # Create the linopy model
    logger.info("Creating linopy model...")
    n.optimize.create_model()
    logger.info("Linopy model created.")

    solver_name = snakemake.params.solver
    solver_options = snakemake.params.solver_options or {}

    # Add health impacts to the objective if data is available
    add_health_objective(
        n,
        snakemake.input.health_risk_breakpoints,
        snakemake.input.health_cluster_cause,
        snakemake.input.health_cause_log,
        snakemake.input.health_cluster_summary,
        snakemake.input.health_clusters,
        snakemake.input.population,
        snakemake.params.health_risk_factors,
        solver_name,
        float(snakemake.params.health_value_per_yll),
    )

    scaling_factor = rescale_objective(n.model)
    if scaling_factor != 1.0:
        previous = n.meta.get("objective_scaling_factor", 1.0)
        n.meta["objective_scaling_factor"] = previous * scaling_factor
        n.meta["objective_scaling_target_coeff"] = OBJECTIVE_COEFF_TARGET

    status, condition = n.model.solve(
        solver_name=solver_name,
        **solver_options,
    )
    result = (status, condition)

    if status == "ok":
        aux_names = HEALTH_AUX_MAP.pop(id(n.model), set())
        variables_container = n.model.variables
        removed = {}
        for name in aux_names:
            if name in variables_container.data:
                removed[name] = variables_container.data.pop(name)

        try:
            n.optimize.assign_solution()
            n.optimize.assign_duals(False)
            n.optimize.post_processing()
        finally:
            if removed:
                variables_container.data.update(removed)

        if options.debug.runtime_verification:
            _optimize_guard(n)

        n.export_to_netcdf(snakemake.output.network)
    elif condition in {"infeasible", "infeasible_or_unbounded"}:
        logger.error("Model is infeasible or unbounded!")
        if solver_name == "gurobi":
            try:
                logger.error("Infeasible constraints:")
                n.model.print_infeasibilities()
            except Exception as exc:
                logger.error("Could not compute infeasibilities: %s", exc)
        else:
            logger.error("Infeasibility diagnosis only available with Gurobi solver")
    else:
        logger.error("Optimization unsuccessful: %s", result)
