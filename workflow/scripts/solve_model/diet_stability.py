# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Diet-side deviation penalty constraints for the food systems model.

Anchors per-(food, country) ``food_consumption`` link dispatch toward the
observed baseline-year diet (the same matched_baseline that
:func:`fix_food_consumption_to_baseline` consumes when
``validation.enforce_baseline_diet`` is true).

Land-side and feed-side penalties only anchor what is *produced*; they leave
the model free to reroute the same hectares toward a different *diet* (e.g.
less feed grain -> more legumes for direct human consumption). The diet
component fills that gap by penalising consumption deviations.

The shared ``penalty_mode`` and ``deviation_type`` from the parent
``deviation_penalty`` block apply:

- ``l1``: linear absolute-value penalty per link, scaled by
  ``deviation_penalty.diet.l1_cost`` (bn USD/Mt).
- ``quadratic``: ``0.5 * deviation_penalty.quadratic_cost * sum((p - baseline)^2)``.
- ``hard``: per-country diet churn budget (see
  :func:`_add_diet_churn_constraints`) mirroring the land/feed regional
  churn budgets, so all deviation components share one flexibility
  geometry.

The baseline ``target_mt`` per link is always the observed-diet anchor for
the configured ``baseline_year``, independent of the scenario's GHG/YLL
pricing, so diet penalties compose cleanly with piecewise consumer-values
utility.
"""

import logging

import numpy as np
import pandas as pd
import pypsa
import xarray as xr

from workflow.scripts.solve_model.production_stability import add_churn_abs_split

logger = logging.getLogger(__name__)


def add_diet_stability_constraints(
    n: pypsa.Network,
    matched_baseline: pd.DataFrame,
    dp_cfg: dict,
    slack_marginal_cost: float,
) -> None:
    """Add diet-component deviation penalty constraints.

    Parameters
    ----------
    n : pypsa.Network
        Network whose ``n.model`` already exists.
    matched_baseline : pd.DataFrame
        Output of ``_match_baseline_to_consume_links`` with columns
        ``name`` and ``target_mt``.
    dp_cfg : dict
        The resolved ``deviation_penalty`` block. ``l1_cost`` must already
        be numeric (see ``resolve_calibrated_l1_costs``).
    slack_marginal_cost : float
        Cost (bn USD/Mt) of the per-country churn-budget slack in hard
        mode when ``diet.enable_slack`` is true; unused otherwise.
    """
    if not dp_cfg["enabled"]:
        return
    diet_cfg = dp_cfg["diet"]
    if not diet_cfg["enabled"]:
        return

    if matched_baseline is None or matched_baseline.empty:
        logger.warning(
            "deviation_penalty.diet is enabled but no baseline diet matched "
            "any food_consumption links; skipping."
        )
        return

    consume_links = n.links.static[n.links.static["carrier"] == "food_consumption"]
    if consume_links.empty:
        logger.info("No food_consumption links present; skipping diet penalty.")
        return

    targets = (
        matched_baseline.set_index("name")["target_mt"]
        .reindex(consume_links.index)
        .fillna(0.0)
        .astype(float)
    )

    min_baseline = float(diet_cfg["min_baseline"])
    if min_baseline <= 0:
        raise ValueError(
            f"deviation_penalty.diet.min_baseline must be > 0; got {min_baseline}"
        )
    deviation_type = dp_cfg["deviation_type"]
    penalty_mode = dp_cfg["penalty_mode"]

    link_p = n.model.variables["Link-p"].sel(snapshot="now", name=consume_links.index)
    baselines = xr.DataArray(
        targets.to_numpy(),
        coords={"name": consume_links.index},
        dims="name",
    )

    if penalty_mode == "hard":
        _add_diet_churn_constraints(
            n, consume_links, link_p, baselines, targets, diet_cfg, slack_marginal_cost
        )
        return

    if deviation_type == "relative":
        denominator = xr.where(baselines > min_baseline, baselines, min_baseline)
        deviation = (link_p - baselines) / denominator
    elif deviation_type == "absolute":
        deviation = link_p - baselines
    else:
        raise ValueError(
            f"deviation_penalty.deviation_type must be 'absolute' or "
            f"'relative', got {deviation_type!r}"
        )

    if penalty_mode == "l1":
        cost = float(diet_cfg["l1_cost"])
        dev_pos = n.model.add_variables(
            lower=0,
            coords=[consume_links.index],
            dims=["name"],
            name="diet_stability_dev_pos",
        )
        dev_neg = n.model.add_variables(
            lower=0,
            coords=[consume_links.index],
            dims=["name"],
            name="diet_stability_dev_neg",
        )
        n.model.add_constraints(
            dev_pos - dev_neg == deviation,
            name="GlobalConstraint-diet_stability_dev_split",
        )
        n.model.objective += cost * (dev_pos.sum() + dev_neg.sum())
        logger.info(
            "Added %d per-(food, country) diet L1 penalties "
            "(cost=%.4f bn USD/Mt, mode=%s)",
            len(consume_links),
            cost,
            deviation_type,
        )
    elif penalty_mode == "quadratic":
        cost = float(dp_cfg["quadratic_cost"])
        dev = n.model.add_variables(
            coords=[consume_links.index],
            dims=["name"],
            name="diet_stability_dev",
        )
        n.model.add_constraints(
            dev == deviation,
            name="GlobalConstraint-diet_stability_dev",
        )
        n.model.objective += 0.5 * cost * (dev * dev).sum()
        logger.info(
            "Added %d per-(food, country) diet quadratic penalties "
            "(cost=%.4f bn USD per (Mt)^2, mode=%s)",
            len(consume_links),
            cost,
            deviation_type,
        )
    else:
        raise ValueError(
            f"deviation_penalty.penalty_mode must be 'hard', 'l1' or "
            f"'quadratic', got {penalty_mode!r}"
        )


def _add_diet_churn_constraints(
    n: pypsa.Network,
    consume_links: pd.DataFrame,
    link_p,
    baselines: xr.DataArray,
    targets: pd.Series,
    diet_cfg: dict,
    slack_marginal_cost: float,
) -> None:
    """Add a per-country diet churn budget (hard mode).

    For each country, the total absolute consumption churn summed over all
    ``food_consumption`` links is bounded::

        sum_{foods} |consumption - target_mt|  <=  2 * delta * C_country

    where ``delta = max_relative_deviation`` and ``C_country`` is the
    country's total baseline consumption (Mt). Replacing baseline mass with
    an equal mass of other foods costs twice that mass in churn (a removal
    plus an addition), so the factor 2 makes ``delta`` a 0-1 dial:
    ``delta=0`` freezes the baseline diet and ``delta=1`` allows a complete
    reorganization. This mirrors the land/feed churn budgets in
    :mod:`production_stability`; the churn is in native Mt
    (``deviation_type`` does not apply). When ``enable_slack`` is true a
    per-country slack variable allows exceeding the budget at
    ``slack_marginal_cost`` per Mt.
    """
    m = n.model
    delta = float(diet_cfg["max_relative_deviation"])

    if "country" not in consume_links.columns:
        raise ValueError(
            "food_consumption links lack a 'country' column; cannot form a "
            "diet churn budget."
        )
    countries = consume_links["country"].astype(str)

    # Per-link absolute churn via a nonneg equality split (Mt).
    churn_abs = add_churn_abs_split(
        m, consume_links.index, link_p - baselines, "diet_churn"
    )

    # Per-country budget: twice the baseline consumption at delta=1.
    country_keys = sorted(countries.unique())
    denom = targets.groupby(countries.to_numpy()).sum()
    budget = (2.0 * delta * denom.reindex(country_keys)).to_numpy(dtype=float)

    group_map = xr.DataArray(
        countries.to_numpy(),
        coords={"name": consume_links.index},
        dims="name",
        name="churn_country",
    )
    churn = churn_abs.groupby(group_map).sum()
    budget_xr = xr.DataArray(
        budget, coords={"churn_country": country_keys}, dims="churn_country"
    )

    if diet_cfg["enable_slack"]:
        slack = m.add_variables(
            lower=0,
            coords=[pd.Index(country_keys, name="churn_country")],
            dims=["churn_country"],
            name="diet_regional_budget_slack",
        )
        m.add_constraints(
            churn - slack <= budget_xr,
            name="GlobalConstraint-diet_regional_budget",
        )
        m.objective += slack_marginal_cost * slack.sum()
    else:
        m.add_constraints(
            churn <= budget_xr,
            name="GlobalConstraint-diet_regional_budget",
        )

    n.global_constraints.add(
        [f"diet_regional_budget_{c}" for c in country_keys],
        sense="<=",
        constant=budget,
        type="diet_stability",
    )

    logger.info(
        "Added %d per-country diet churn-budget constraints (delta=%.1f%% of a "
        "full reorganization, slack=%s); total budget %.1f Mt over %d countries",
        len(country_keys),
        delta * 100,
        diet_cfg["enable_slack"],
        float(budget.sum()),
        len(country_keys),
    )


def evaluate_diet_stability_cost(
    n: pypsa.Network,
    matched_baseline: pd.DataFrame,
    dp_cfg: dict,
    slack_marginal_cost: float,
) -> float:
    """Re-evaluate the diet deviation cost from a solved network.

    Used for objective-breakdown bookkeeping. Returns the L1 (or quadratic)
    penalty contribution to the objective in bn USD/yr -- or, in hard mode,
    the churn-budget slack cost -- and 0.0 when disabled.
    """
    if not dp_cfg["enabled"]:
        return 0.0
    diet_cfg = dp_cfg["diet"]
    if not diet_cfg["enabled"]:
        return 0.0
    if matched_baseline is None or matched_baseline.empty:
        return 0.0

    consume_links = n.links.static[n.links.static["carrier"] == "food_consumption"]
    if consume_links.empty:
        return 0.0

    if dp_cfg["penalty_mode"] == "hard":
        # Hard mode carries no penalty term; only the churn-budget slack
        # (when enabled) enters the objective.
        if n.model is None or "diet_regional_budget_slack" not in n.model.variables:
            return 0.0
        sol = n.model.variables["diet_regional_budget_slack"].solution
        return float(slack_marginal_cost) * float(sol.sum())

    targets = (
        matched_baseline.set_index("name")["target_mt"]
        .reindex(consume_links.index)
        .fillna(0.0)
        .astype(float)
        .to_numpy()
    )
    actual = (
        n.links.dynamic["p0"]
        .iloc[-1]
        .reindex(consume_links.index)
        .fillna(0.0)
        .to_numpy()
    )

    min_baseline = float(diet_cfg["min_baseline"])
    deviation_type = dp_cfg["deviation_type"]

    if deviation_type == "relative":
        denominator = np.where(targets > min_baseline, targets, min_baseline)
        dev = (actual - targets) / denominator
    else:
        dev = actual - targets

    penalty_mode = dp_cfg["penalty_mode"]
    if penalty_mode == "l1":
        return float(diet_cfg["l1_cost"]) * float(np.abs(dev).sum())
    if penalty_mode == "quadratic":
        return 0.5 * float(dp_cfg["quadratic_cost"]) * float((dev * dev).sum())
    return 0.0
