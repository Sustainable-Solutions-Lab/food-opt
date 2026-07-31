# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Convex piecewise-linear market response (positive mathematical programming).

The deviation penalty in :mod:`production_stability` anchors production with a
single kink at the baseline: one absolute-value term per link, priced at a
coefficient calibrated so that total deviation lands on a chosen target. Below
the kink the marginal friction is constant, above it too, so the response to a
shock is a dead band followed by a jump to whatever bound stops it.

This module replaces the aggregate part of that friction with a supply curve in
the tradition of positive mathematical programming (Howitt 1995; Heckelei and
Britz 2005; exact calibration following Merel and Bucaram 2010): each production
group gets a convex, increasing marginal-cost curve through its observed
activity, so response is graduated and governed by an exogenous supply
elasticity rather than by a fitted deviation target.

Calibration is the standard two-phase procedure. Phase 1 (``pin_baseline``)
solves the model with each group's activity fixed at its observed value by an
equality constraint; the dual ``lambda_g`` of that constraint is the wedge
between the marginal value of the group's output and its accounting marginal
cost ``c_g`` at the observed point. Phase 2 (``intercepts``) feeds the wedges
back as intercepts on the curves, so each group's marginal cost at the baseline
becomes ``c_g + lambda_g`` -- exactly the marginal value there. The observed
allocation then satisfies the optimality conditions of the *unpinned* model and
is reproduced exactly, with no hard constraints and nothing tuned.

For a group ``g`` with observed activity ``b_g`` and accounting marginal cost
``c_g``, a linear marginal-cost curve reproducing supply elasticity ``eta`` has
slope

    gamma_g = c_g / (eta * b_g),

since equating marginal cost to the marginal value of output gives
``dX/dv = 1 / gamma`` and hence ``eta = v / (gamma * b)`` with ``v = c_g`` at
the calibrated point. The elasticity is stated at the accounting cost (Howitt's
original form) rather than at the calibrated marginal cost ``c_g + lambda_g``
(Merel and Bucaram 2010): under a pinned diet the wedges trace the Ricardian
rent gradient, and rent-poor groups -- in practice a large share, not an edge
case -- have ``c_g + lambda_g <= 0``, where the calibrated-cost form is
undefined. Exactness does not depend on this choice -- only the intercept
enters the optimality condition at the baseline -- and the response to a
policy shock scales as ``eta`` per relative change of the *accounting* cost:
a group with wedge ``lambda_g`` realises elasticity
``eta * (c_g + lambda_g) / c_g`` with respect to its own calibrated marginal
cost, so rent-poor groups respond more stiffly than ``eta`` when measured at
that price. The configured elasticity is therefore a response per relative
accounting-cost change, not an arc elasticity at the curve's equilibrium
price; the demand side is analogous, with the slope stated at the
``slope_basis`` reference price rather than at the fitted wedge level.

The curve enters the objective as its deviation cost, the intercept plus the
integral of the slope from the baseline,

    Psi_g(X) = lambda_g * (X - b_g) + gamma_g / 2 * (X - b_g)^2,

which is convex and, net of the value of output, minimised at ``X = b_g``, so
the observed activity stays optimal at reference prices. ``Psi_g`` enters the
LP as its piecewise-linear interpolant, sampled at ``n_blocks`` breakpoints per
side of the baseline: the deviation from baseline splits into one bounded
segment variable per breakpoint interval whose objective coefficient is the
chord slope over that interval. Convexity makes the segment costs nondecreasing
away from the baseline, so minimisation fills them in merit order; the chords
enter as variable bounds rather than constraint rows. The outermost segment on
each side is unbounded, so beyond the sampled range the curve extrapolates at
the end chords' slopes and imposes no hidden activity bound.

The breakpoint spacing sets the finest response the curve can resolve: between
adjacent breakpoints the marginal cost is a single chord slope. With equal
spacing over a range of half the baseline, four points resolve nothing below
12.5% of baseline -- coarser than the moves most groups make, so the curve
collapses into a flat-rate penalty. ``width_growth`` above 1 narrows the
spacing near the baseline, concentrating resolution where groups sit while the
outer breakpoints still resolve large moves.

Granularity: curves apply to *groups*, whose spatial resolution is
configurable. ``country`` prices how much of a commodity a country produces
(``(crop, country)``, ``(combination, country)`` for multi-cropping,
``(country)`` for grassland, ``(product, country)`` for animals); ``region``
resolves that per optimisation region; ``link`` gives every production link its
own curve. A coarse curve says nothing about which regions, resource classes or
water-supply types activity sits in, so it leaves spatial reallocation to the
per-link deviation penalty; at ``link`` granularity the curves price that
reallocation themselves and replace the deviation penalty outright.

The ``demand`` component applies the same machinery to the consumption side:
each (food, country) consumption link gets a concave marginal-*utility* curve
through its observed intake. The deviation-cost form is identical -- for
consumption ``X`` with baseline ``b``, the objective carries
``lambda * (X - b) + gamma / 2 * (X - b)^2`` with ``gamma > 0`` -- because a
convex deviation cost *is* a diminishing marginal utility: the KKT condition at
the baseline reads ``price + lambda = 0``, so the fitted intercept is minus the
marginal willingness to pay, and the positive slope makes that willingness
decline with quantity.

The demand calibration is *sequential* rather than joint. Pinning production
and consumption in one solve closes every commodity chain, leaving each
chain's price level -- and hence the split of the total wedge between producer
and consumer -- undetermined (the duals park at the pin slack bound). Instead,
the demand pin runs against the already-fitted, elastic supply curves
(``pin_baseline: ["demand"]`` with production intercepts supplied), so the
food-bus prices are unique and the demand wedge is measured against exactly
the supply side deployed solves will carry. The wedge *level* still inherits
the supply fit's value convention, so it cannot set the curve's stiffness; the
slope is instead stated at a reference price fitted by a separate slope-basis
solve -- the same demand pin against *zero-intercept* supply curves, whose
prices sit at the accounting-cost chain (``slope_basis`` column of the
artefact): ``gamma = P_ref / (eta * b)``. Demand curves therefore require the
sequentially fitted artefact; there is no uncalibrated variant.
"""

import logging

import numpy as np
import pandas as pd
import pypsa
import xarray as xr

logger = logging.getLogger(__name__)

# Dimension of the per-side deviation segments of a group's curve.
SEGMENT_DIM = "mr_segment"

# Component -> (carriers, baseline column, commodity-identifying columns). The
# commodity columns say *what* is produced; the granularity setting adds the
# columns saying *where*, so both are needed to form a group key.
COMPONENTS = {
    "crops": (("crop_production",), "baseline_area_mha", ("crop",)),
    "multi_crops": (("crop_production_multi",), "baseline_area_mha", ("combination",)),
    "grassland": (("grassland_production",), "baseline_area_mha", ()),
    "animals": (("animal_production",), "baseline_feed_use_mt_dm", ("product",)),
    "demand": (("food_consumption",), "baseline_consumption_mt", ("food",)),
}

# Components whose activity is consumption rather than production: their curve
# is a marginal utility (slope stated at the fitted willingness to pay, not at
# an accounting cost) and their cost is reported as demand response.
DEMAND_COMPONENTS = frozenset({"demand"})

# Spatial columns appended to a component's commodity key at each granularity.
# ``link`` is the exception: it resolves every production link separately, so it
# replaces the whole key rather than extending it (a link already identifies its
# commodity and its location), and is represented by an empty tuple.
GRANULARITY_COLUMNS = {
    "country": ("country",),
    "region": ("country", "region"),
    "link": None,
}


def _evaluate_curve_cost(n: pypsa.Network, components) -> float:
    m = getattr(n, "model", None)
    if m is None:
        return 0.0
    return sum(
        float(m.variables[name].solution.sum())
        for component in components
        if (name := f"{component}_market_response_cost") in m.variables
    )


def evaluate_supply_response_cost(n: pypsa.Network) -> float:
    """Total deviation cost the production curves contributed to the objective.

    Each component's cost variable enters the objective directly, so the term
    is just the sum of their solutions, and the objective breakdown can report
    it as its own category with an exact identity check.
    """
    return _evaluate_curve_cost(n, COMPONENTS.keys() - DEMAND_COMPONENTS)


def evaluate_demand_response_cost(n: pypsa.Network) -> float:
    """Deviation cost of the demand curves, reported separately from supply."""
    return _evaluate_curve_cost(n, DEMAND_COMPONENTS)


def extract_market_response_intercepts(
    n: pypsa.Network, pin_slack_cost: "float | None" = None
) -> pd.DataFrame:
    """Per-group price wedges from a baseline-pinned solve (PMP phase 1).

    The dual of each group's pinning constraint is the gap between the marginal
    value of the group's output and its marginal cost at the observed activity.
    Returned as a frame with columns ``component``, ``mr_group``, ``intercept``;
    fed back via ``market_response.intercepts``, the wedges make the observed
    allocation the exact optimum of the unpinned model. When ``pin_slack_cost``
    is given, groups whose dual sits at that bound are reported as censored
    even if they used no slack.
    """
    m = n.model
    frames = []
    for component in COMPONENTS:
        name = f"GlobalConstraint-{component}_market_response_balance"
        if name not in m.constraints:
            continue
        dual = m.constraints[name].dual
        slack = (
            m.variables[f"{component}_market_response_pin_slack_up"].solution
            - m.variables[f"{component}_market_response_pin_slack_down"].solution
        )
        frames.append(
            pd.DataFrame(
                {
                    "component": component,
                    "mr_group": dual.coords["mr_group"].values,
                    # The constraint reads activity == baseline, so the dual is
                    # the objective change per unit of *required* activity: the
                    # cost of producing one more unit minus the value gained.
                    # The wedge (value minus cost) is its negation.
                    "intercept": -dual.values,
                    # Net pin slack (activity minus baseline). Nonzero means
                    # the reference point is inconsistent with some other
                    # constraint and the intercept is censored at the slack
                    # price -- worth investigating upstream.
                    "slack": slack.values,
                }
            )
        )
    if not frames:
        raise ValueError(
            "no market-response pinning constraints found; intercepts can only "
            "be extracted from a solve with market_response.pin_baseline: true"
        )
    table = pd.concat(frames, ignore_index=True)
    slacked = table["slack"].abs() > 1e-6
    if slacked.any():
        worst = table.loc[slacked].reindex(
            table.loc[slacked, "slack"].abs().sort_values(ascending=False).index
        )
        logger.warning(
            "%d of %d pinned market-response groups used pin slack (total "
            "|slack| %.3f); their intercepts are censored at the slack price. "
            "Worst: %s",
            int(slacked.sum()),
            len(table),
            float(table["slack"].abs().sum()),
            worst[["component", "mr_group", "slack"]].head(5).to_dict("records"),
        )
    if pin_slack_cost is not None:
        # A dual parked at the slack bound is censored whether or not slack
        # was used: the true wedge lies at or beyond the bound (or the pin is
        # degenerate with another binding constraint).
        at_bound = table["intercept"].abs() >= float(pin_slack_cost) * (1 - 1e-6)
        if at_bound.any():
            counts = table.loc[at_bound, "component"].value_counts().to_dict()
            logger.warning(
                "%d of %d pinned market-response groups have their dual at the "
                "pin slack bound (+/-%.3g); their intercepts are censored "
                "(per component: %s)",
                int(at_bound.sum()),
                len(table),
                float(pin_slack_cost),
                counts,
            )
    return table


def lift_pinned_capacity_bounds(n: pypsa.Network, cfg: dict) -> None:
    """Lift ``p_nom_max`` on links about to be pinned at link granularity.

    Many links sit exactly at their capacity bound (baseline equals
    ``p_nom_max`` after the suitable-area clip), which makes a phase-1 pin
    redundant with the capacity constraint and its dual degenerate: the
    solver then reports the pin slack bound instead of the true wedge.
    Lifting the capacity of the pinned links makes the pin the unique binding
    constraint, so its dual is identified (and includes any scarcity rent,
    which the curve rightly prices on contraction). The pin itself still
    holds every link at baseline, so no other part of the solution changes.

    Must be called BEFORE ``n.optimize.create_model()``: the capacity bounds
    enter the (frozen) linopy model at creation. Only meaningful at ``link``
    granularity, where every pinned group is a single link; at coarser
    granularities the group pin binds the aggregate while individual capacity
    bounds legitimately shape the within-group allocation.
    """
    if not cfg["enabled"] or not cfg["pin_baseline"]:
        return
    if str(cfg["granularity"]) != "link":
        return
    pin_cfg = cfg["pin_baseline"]
    pinned = set(pin_cfg) if isinstance(pin_cfg, list) else set(COMPONENTS)
    carriers = [
        carrier
        for component, (component_carriers, _, _) in COMPONENTS.items()
        if component in pinned
        and component not in DEMAND_COMPONENTS
        and cfg["components"][component]
        for carrier in component_carriers
    ]
    links = n.links.static
    mask = (
        links["carrier"].isin(carriers)
        & links["p_nom_extendable"].astype(bool)
        & np.isfinite(links["p_nom_max"].astype(float))
    )
    if mask.any():
        n.links.static.loc[mask, "p_nom_max"] = np.inf
        logger.info(
            "Lifted capacity bounds of %d pinned production links to "
            "identify pin duals",
            int(mask.sum()),
        )


def add_market_response_curves(n: pypsa.Network, cfg: dict) -> None:
    """Add convex piecewise-linear supply curves to the objective.

    Parameters
    ----------
    n : pypsa.Network
        Network whose ``model`` has already been created.
    cfg : dict
        The resolved ``market_response`` configuration block.
    """
    if not cfg["enabled"]:
        return

    n_blocks = int(cfg["n_blocks"])
    if n_blocks < 1:
        raise ValueError(f"market_response.n_blocks must be >= 1, got {n_blocks}")
    expansion = float(cfg["expansion_range"])
    contraction = float(cfg["contraction_range"])
    if expansion <= 0:
        raise ValueError(
            f"market_response.expansion_range must be > 0, got {expansion}"
        )
    if not 0 < contraction <= 1:
        raise ValueError(
            "market_response.contraction_range must satisfy 0 < r <= 1, got "
            f"{contraction}"
        )
    factor = float(cfg["elasticity_factor"])
    if factor <= 0:
        raise ValueError(f"market_response.elasticity_factor must be > 0, got {factor}")
    width_growth = float(cfg["width_growth"])
    if width_growth < 1:
        raise ValueError(
            f"market_response.width_growth must be >= 1, got {width_growth}; "
            "below 1 the breakpoints would spread toward the baseline, putting "
            "the coarsest resolution where groups actually sit"
        )

    granularity = str(cfg["granularity"])
    if granularity not in GRANULARITY_COLUMNS:
        raise ValueError(
            f"market_response.granularity must be one of "
            f"{sorted(GRANULARITY_COLUMNS)}, got '{granularity}'"
        )
    spatial_cols = GRANULARITY_COLUMNS[granularity]

    # pin_baseline is false (curves everywhere), true (pin every enabled
    # component), or a list of component names to pin while the remaining
    # components carry curves from ``intercepts``. The list form is what the
    # sequential calibration uses: demand wedges are measured against the
    # already-fitted, *elastic* supply side, which keeps the food-bus prices
    # unique where joint pinning would close the chain and degenerate them.
    pin_cfg = cfg["pin_baseline"]
    if isinstance(pin_cfg, list):
        unknown = set(pin_cfg) - set(COMPONENTS)
        if unknown:
            raise ValueError(
                f"market_response.pin_baseline names unknown components "
                f"{sorted(unknown)}; valid components: {sorted(COMPONENTS)}"
            )
        pinned_components = set(pin_cfg)
    else:
        pinned_components = set(COMPONENTS) if pin_cfg else set()
    if pinned_components == set(COMPONENTS) and cfg["intercepts"]:
        raise ValueError(
            "market_response.pin_baseline and market_response.intercepts are "
            "mutually exclusive when every component is pinned: the pinned "
            "solve fits the intercepts against the accounting costs, so it "
            "must not itself carry intercepts"
        )
    intercepts_table = None
    intercepts = None
    if cfg["intercepts"]:
        intercepts_table = pd.read_csv(cfg["intercepts"]).set_index(
            ["component", "mr_group"]
        )
        intercepts = intercepts_table["intercept"]

    for component, (carriers, baseline_col, commodity_cols) in COMPONENTS.items():
        if not cfg["components"][component]:
            continue
        pin = component in pinned_components
        raw_elasticity = cfg["elasticities"][component]
        if component in DEMAND_COMPONENTS:
            # Per-food-group demand elasticity magnitudes.
            elasticity = {
                str(group): float(value) * factor
                for group, value in raw_elasticity.items()
            }
            bad_groups = [g for g, v in elasticity.items() if v <= 0]
            if bad_groups:
                raise ValueError(
                    f"market_response.elasticities.{component} must be > 0 "
                    f"after applying elasticity_factor for every food group; "
                    f"non-positive: {bad_groups}"
                )
        else:
            elasticity = float(raw_elasticity) * factor
            if elasticity <= 0:
                raise ValueError(
                    f"market_response.elasticities.{component} must be > 0 after "
                    f"applying elasticity_factor, got {elasticity}"
                )
        if component in DEMAND_COMPONENTS and not pin and intercepts is None:
            raise ValueError(
                "market_response demand curves require fitted intercepts: the "
                "marginal-utility level comes from the fitted willingness to "
                "pay, so there is no uncalibrated variant. Provide "
                "market_response.intercepts or disable "
                "market_response.components.demand"
            )
        component_intercepts = None
        component_slope_basis = None
        if intercepts is not None and not pin:
            if component not in intercepts.index.get_level_values("component"):
                raise ValueError(
                    f"market_response.intercepts has no entries for component "
                    f"'{component}'; the intercepts must come from a pinned "
                    "solve covering the same components"
                )
            component_intercepts = intercepts.xs(component, level="component")
            if component in DEMAND_COMPONENTS:
                if "slope_basis" not in intercepts_table.columns:
                    raise ValueError(
                        "market_response.intercepts has no slope_basis column; "
                        "demand curves need the reference-price basis from the "
                        "sequential calibration (tools/calibrate "
                        "market_response)"
                    )
                component_slope_basis = intercepts_table.xs(
                    component, level="component"
                )["slope_basis"]
        _add_component_curves(
            n,
            component=component,
            carriers=carriers,
            baseline_col=baseline_col,
            group_cols=() if spatial_cols is None else commodity_cols + spatial_cols,
            elasticity=elasticity,
            n_blocks=n_blocks,
            expansion=expansion,
            contraction=contraction,
            width_growth=width_growth,
            pin=pin,
            pin_slack_cost=float(cfg["pin_slack_cost"]),
            intercepts=component_intercepts,
            slope_basis=component_slope_basis,
        )


def _group_key(links: pd.DataFrame, group_cols: tuple[str, ...]) -> pd.Series:
    """Group label per link: the commodity and spatial columns joined.

    With no columns at all -- the ``link`` granularity, where nothing is
    aggregated -- the link's own name is the key, so each production link gets
    its own curve. Every link then carries the configured elasticity with
    respect to its own margin, and the group aggregate carries it too whenever
    costs within the group are similar.
    """
    if not group_cols:
        return pd.Series(links.index.astype(str), index=links.index)
    null_cols = [c for c in group_cols if links[c].isna().any()]
    if null_cols:
        raise ValueError(
            f"market-response grouping columns {null_cols} contain missing "
            "values; every covered link needs a complete group key"
        )
    return links[list(group_cols)].astype(str).agg("::".join, axis=1)


def _group_table(
    links: pd.DataFrame, baseline_col: str, group_cols: tuple[str, ...]
) -> pd.DataFrame:
    """Aggregate per-link baselines and costs into per-group curve parameters.

    Returns a frame indexed by the group key with the group's baseline activity
    and its baseline-weighted mean marginal cost. Zero-baseline groups are kept
    so callers can constrain their activity explicitly; their cost is undefined
    because there is no observed activity over which to average it.
    """
    work = links.copy()
    work["_baseline"] = work[baseline_col].astype(float)
    work["_group"] = _group_key(work, group_cols)
    if work.empty:
        return pd.DataFrame(columns=["baseline", "cost"])

    negative = work["_baseline"] < 0
    if negative.any():
        examples = work.index[negative][:5].tolist()
        raise ValueError(
            f"{int(negative.sum())} market-response link(s) have a negative "
            f"baseline (e.g. {examples})"
        )

    work["_cost_weighted"] = work["_baseline"] * work["marginal_cost"].astype(float)
    grouped = work.groupby("_group").agg(
        baseline=("_baseline", "sum"), _cost_weighted=("_cost_weighted", "sum")
    )
    positive = grouped["baseline"] > 0
    grouped["cost"] = np.nan
    grouped.loc[positive, "cost"] = (
        grouped.loc[positive, "_cost_weighted"] / grouped.loc[positive, "baseline"]
    )
    return grouped[["baseline", "cost"]].sort_index()


def _add_component_curves(
    n: pypsa.Network,
    *,
    component: str,
    carriers: tuple[str, ...],
    baseline_col: str,
    group_cols: tuple[str, ...],
    elasticity: "float | dict[str, float]",
    n_blocks: int,
    expansion: float,
    contraction: float,
    width_growth: float,
    pin: bool,
    pin_slack_cost: float,
    intercepts: "pd.Series | None",
    slope_basis: "pd.Series | None" = None,
) -> None:
    """Add one component's curves (or its phase-1 pin) to the model."""
    links_df = n.links.static
    links = links_df[links_df["carrier"].isin(carriers)]
    if links.empty:
        logger.info(
            "No %s links with %s; skipping market-response curves",
            "/".join(carriers),
            baseline_col,
        )
        return

    if baseline_col not in links.columns:
        if component in DEMAND_COMPONENTS:
            raise ValueError(
                f"market-response {component} links have no '{baseline_col}' "
                "column; every food-consumption link needs an explicit observed "
                "baseline, including zero"
            )
        logger.info(
            "No %s links with %s; skipping market-response curves",
            "/".join(carriers),
            baseline_col,
        )
        return

    baseline = pd.to_numeric(links[baseline_col], errors="coerce")
    missing = baseline.isna()
    if missing.any():
        examples = links.index[missing][:5].tolist()
        raise ValueError(
            f"{int(missing.sum())} market-response {component} link(s) have no "
            f"explicit baseline in '{baseline_col}' (e.g. {examples}); every "
            "covered link must carry one observed baseline value, using zero "
            "for observed absence"
        )

    groups = _group_table(links, baseline_col, group_cols)
    if groups.empty:
        logger.info("No %s groups; skipping market-response curves", component)
        return

    work = links.copy()
    work["_group"] = _group_key(work, group_cols)

    m = n.model
    link_p = m.variables["Link-p"].sel(snapshot="now")
    group_map = xr.DataArray(
        work["_group"].to_numpy(),
        coords={"name": work.index},
        dims="name",
        name="mr_group",
    )
    activity = link_p.sel(name=work.index).groupby(group_map).sum()
    # groupby orders groups lexicographically, matching the sorted group table.
    group_index = pd.Index(activity.coords["mr_group"].values)
    if not group_index.equals(groups.index):
        groups = groups.reindex(group_index)

    # A zero observed aggregate has no intensive-margin elasticity to identify:
    # gamma = price / (eta * baseline) is undefined. Keep it outside the curve
    # and fix its activity to zero in both calibration and deployment. At coarse
    # granularity this applies only when the whole aggregate is zero; a
    # zero-baseline member of a positive aggregate remains free to participate
    # in that aggregate's reallocation.
    zero_groups = groups.index[groups["baseline"] == 0]
    if len(zero_groups):
        m.add_constraints(
            activity.sel(mr_group=zero_groups.to_numpy()) == 0,
            name=f"GlobalConstraint-{component}_market_response_zero_baseline",
        )
        logger.info(
            "Fixed %d of %d %s market-response groups (%.0f%%) with zero "
            "baseline to zero activity",
            len(zero_groups),
            len(groups),
            component,
            100 * len(zero_groups) / len(groups),
        )

    positive_groups = groups.index[groups["baseline"] > 0]
    if not len(positive_groups):
        logger.info("No %s groups with positive baseline; skipping curves", component)
        return

    groups = groups.loc[positive_groups]
    activity = activity.sel(mr_group=positive_groups.to_numpy())
    group_index = positive_groups
    work = work[work["_group"].isin(positive_groups)]

    if pin:
        # PMP phase 1: hold every group at its observed activity. The dual of
        # this constraint is the group's price wedge, extracted after the solve
        # by extract_market_response_intercepts. The pin is elastic -- slack at
        # a high finite price -- because the reference data is never perfectly
        # consistent with every other constraint (land availability, water,
        # feed balances), and a hard pin would make the whole solve infeasible
        # over the smallest such residual. Where slack is used the dual is
        # censored at the slack price, which flags the group for investigation
        # instead of hiding the inconsistency.
        slack_cost = float(pin_slack_cost)
        if slack_cost <= 0:
            raise ValueError(
                f"market_response.pin_slack_cost must be > 0, got {slack_cost}"
            )
        baseline_da = xr.DataArray(
            groups["baseline"].to_numpy(),
            coords={"mr_group": group_index.to_numpy()},
            dims="mr_group",
        )
        slack_coords = {"mr_group": group_index.to_numpy()}
        slack_up = m.add_variables(
            lower=0,
            coords=slack_coords,
            dims=("mr_group",),
            name=f"{component}_market_response_pin_slack_up",
        )
        slack_down = m.add_variables(
            lower=0,
            coords=slack_coords,
            dims=("mr_group",),
            name=f"{component}_market_response_pin_slack_down",
        )
        m.add_constraints(
            activity - slack_up + slack_down == baseline_da,
            name=f"GlobalConstraint-{component}_market_response_balance",
        )
        m.objective += slack_cost * (slack_up.sum() + slack_down.sum())
        logger.info(
            "Pinned %s activity to baseline: %d groups (baseline=%.1f, "
            "slack cost=%.3g)",
            component,
            len(group_index),
            float(groups["baseline"].sum()),
            slack_cost,
        )
        return

    lam = pd.Series(0.0, index=group_index)
    if intercepts is not None:
        lam = intercepts.reindex(group_index)
        if lam.isna().any():
            missing = lam.index[lam.isna()][:5].tolist()
            raise ValueError(
                f"{int(lam.isna().sum())} {component} group(s) have no entry in "
                f"market_response.intercepts (e.g. {missing}); the intercepts "
                "must come from a pinned solve of the same model at the same "
                "granularity"
            )

    if component in DEMAND_COMPONENTS:
        # Marginal-utility slope, stated at the reference price fitted by the
        # slope-basis solve (demand pinned against zero-intercept supply
        # curves, so food prices sit at the accounting-cost chain). The fitted
        # intercept cannot serve as the basis: its level depends on how the
        # calibration splits value between producer and consumer wedges, which
        # exactness leaves free, while the slope needs a real price scale. A
        # zero basis means supplying the food costs nothing at the margin; the
        # curve degenerates to perfectly elastic demand there, which anchors
        # nothing but also distorts nothing.
        if slope_basis is None:
            raise ValueError(
                f"{component} curves need a slope_basis series; the intercepts "
                "artefact must come from the sequential calibration"
            )
        basis = slope_basis.reindex(group_index)
        if basis.isna().any():
            missing = group_index[basis.isna()][:5].tolist()
            raise ValueError(
                f"{int(basis.isna().sum())} {component} group(s) have no "
                f"slope_basis entry (e.g. {missing})"
            )
        # Per-food-group elasticity: every consumed group must have an entry.
        group_food = (
            work.drop_duplicates("_group")
            .set_index("_group")["food_group"]
            .reindex(group_index)
        )
        eta = group_food.map(elasticity)
        if eta.isna().any():
            missing_groups = sorted(group_food[eta.isna()].dropna().unique())
            raise ValueError(
                f"market_response.elasticities.{component} has no entry for "
                f"food group(s) {missing_groups}"
            )
        slope = basis.abs() / (eta * groups["baseline"])
        flat = slope <= 0
        if flat.any():
            logger.warning(
                "%d %s group(s) have a zero reference price; their demand "
                "curve is flat (perfectly elastic) and does not anchor "
                "consumption. Examples: %s",
                int(flat.sum()),
                component,
                group_index[flat][:5].tolist(),
            )
    else:
        bad = groups["cost"] <= 0
        if bad.any():
            raise ValueError(
                f"{int(bad.sum())} {component} group(s) have a non-positive "
                f"baseline-weighted marginal cost (e.g. "
                f"{group_index[bad][:5].tolist()}); the supply-curve slope "
                "gamma = cost / (elasticity * baseline) is undefined. Marginal "
                "costs come from the cost data and the cost calibration, so a "
                "zero here is an upstream data problem rather than something "
                "to default away."
            )

        # Slope of the marginal-cost curve, in objective units per activity
        # squared. Stated at the accounting cost, not the calibrated marginal
        # cost c + lambda: under a pinned diet the wedges trace the Ricardian
        # rent gradient, and sub-marginal groups at the extensive margin have
        # c + lambda <= 0, where a calibrated-cost slope is undefined. The
        # intercept alone carries exactness -- the observed optimum only needs
        # the marginal cost *level* at the baseline to match the marginal value
        # there -- so the slope choice affects only the stiffness of the
        # response, which for a group with wedge lambda realises elasticity
        # eta * (c + lambda) / c with respect to the market price.
        slope = groups["cost"] / (elasticity * groups["baseline"])

    # Breakpoint offsets from the baseline, as fractions of it: n_blocks points
    # per side, spaced by shares that grow by width_growth away from the
    # baseline. At width_growth == 1 the spacing is uniform over each range;
    # above 1 it concentrates resolution at the baseline where nearly every
    # group sits, while the outer points still resolve large moves.
    shares = width_growth ** np.arange(n_blocks)
    cum = np.cumsum(shares / shares.sum())
    offsets = np.concatenate([-contraction * cum[::-1], [0.0], expansion * cum])

    # The deviation cost lambda * d + slope / 2 * d^2 sampled on the grid; the
    # chords between consecutive samples give each segment's marginal cost.
    # The deviation from baseline is split into one bounded variable per
    # segment; because the curve is convex the segment costs are nondecreasing
    # away from the baseline, so minimisation fills them in merit order and
    # the formulation is exactly the piecewise-linear interpolant of the
    # curve, kinked at the baseline point so the observed activity is strictly
    # optimal. The chords enter as variable bounds and objective coefficients
    # rather than constraint rows -- 2 * n_blocks rows per group would be
    # about half the rows of the whole model. The outermost segment on each
    # side is unbounded, so beyond the sampled range the curve extrapolates at
    # the end chords' slopes, imposing no hidden activity bound.
    baseline_np = groups["baseline"].to_numpy()
    deviation = baseline_np[:, None] * offsets[None, :]
    x_points = baseline_np[:, None] + deviation
    y_points = (
        lam.to_numpy()[:, None] * deviation
        + slope.to_numpy()[:, None] / 2 * deviation**2
    )
    chord_slopes = np.diff(y_points, axis=1) / np.diff(x_points, axis=1)

    seg_coords = {
        "mr_group": group_index.to_numpy(),
        SEGMENT_DIM: np.arange(n_blocks),
    }
    seg_dims = ("mr_group", SEGMENT_DIM)
    # Costs and widths ordered from the baseline outward on each side; a
    # downward segment's cost is the negated chord slope (the curve rises
    # leftward wherever contraction is costly).
    up_slopes = chord_slopes[:, n_blocks:].copy()
    if component in DEMAND_COMPONENTS:
        # A demand group whose marginal utility is still positive at the outer
        # end of the sampled range would otherwise extrapolate that utility
        # forever on the unbounded tail, making it perfectly elastic there
        # with no satiation. Floor the tail's cost slope at zero: expansion
        # beyond the sampled range earns no further utility, while remaining
        # feasible if nutrition constraints force intake up.
        runaway = up_slopes[:, -1] < 0
        if runaway.any():
            up_slopes[runaway, -1] = 0.0
            logger.info(
                "Floored the unbounded expansion tail at zero marginal "
                "utility for %d of %d %s market-response groups",
                int(runaway.sum()),
                len(group_index),
                component,
            )
    up_cost = xr.DataArray(up_slopes, seg_coords, seg_dims)
    down_cost = xr.DataArray(-chord_slopes[:, :n_blocks][:, ::-1], seg_coords, seg_dims)
    widths_up = baseline_np[:, None] * np.diff(np.concatenate([[0.0], expansion * cum]))
    widths_down = baseline_np[:, None] * np.diff(
        np.concatenate([[0.0], contraction * cum])
    )
    widths_up[:, -1] = np.inf
    widths_down[:, -1] = np.inf
    up = m.add_variables(
        lower=0,
        upper=xr.DataArray(widths_up, seg_coords, seg_dims),
        coords=seg_coords,
        dims=seg_dims,
        name=f"{component}_market_response_up",
    )
    down = m.add_variables(
        lower=0,
        upper=xr.DataArray(widths_down, seg_coords, seg_dims),
        coords=seg_coords,
        dims=seg_dims,
        name=f"{component}_market_response_down",
    )
    baseline_da = xr.DataArray(
        baseline_np, coords={"mr_group": group_index.to_numpy()}, dims="mr_group"
    )
    m.add_constraints(
        activity - up.sum(SEGMENT_DIM) + down.sum(SEGMENT_DIM) == baseline_da,
        name=f"GlobalConstraint-{component}_market_response_segments",
    )
    # The per-group deviation cost as an explicit variable: it keeps the
    # objective breakdown an exact identity and presolve folds it away.
    cost = m.add_variables(
        coords={"mr_group": group_index.to_numpy()},
        dims=("mr_group",),
        name=f"{component}_market_response_cost",
    )
    m.add_constraints(
        cost - (up_cost * up).sum(SEGMENT_DIM) - (down_cost * down).sum(SEGMENT_DIM)
        == 0,
        name=f"GlobalConstraint-{component}_market_response_curve",
    )
    m.objective += cost.sum()

    elasticity_label = (
        f"{min(elasticity.values()):.2f}-{max(elasticity.values()):.2f} by group"
        if isinstance(elasticity, dict)
        else f"{elasticity:.3f}"
    )
    logger.info(
        "Added %s market-response curves: %d groups x %d breakpoints/side "
        "(elasticity=%s, baseline=%.1f, slope median=%.4g, intercept "
        "median=%.4g)",
        component,
        len(group_index),
        n_blocks,
        elasticity_label,
        float(groups["baseline"].sum()),
        float(slope.median()),
        float(lam.median()),
    )
