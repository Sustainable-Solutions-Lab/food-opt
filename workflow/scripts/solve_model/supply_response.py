# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Convex piecewise-linear supply response (positive mathematical programming).

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
rent gradient, and sub-marginal groups at the extensive margin have
``c_g + lambda_g <= 0``, where the calibrated-cost form is undefined. Exactness
does not depend on this choice -- only the intercept enters the optimality
condition at the baseline -- but a group with wedge ``lambda_g`` realises
elasticity ``eta * (c_g + lambda_g) / c_g`` with respect to the market price,
so rent-poor groups respond more stiffly than ``eta`` suggests.

The curve enters the objective as its deviation cost, the intercept plus the
integral of the slope from the baseline,

    Psi_g(X) = lambda_g * (X - b_g) + gamma_g / 2 * (X - b_g)^2,

which is convex and, net of the value of output, minimised at ``X = b_g``, so
the observed activity stays optimal at reference prices. ``Psi_g`` enters the
LP through linopy's piecewise machinery: it is sampled at ``n_blocks``
breakpoints per side of the baseline and one free variable per group is bounded
below by the chords between consecutive samples (``linopy.piecewise
.tangent_lines``). Minimisation presses the variable onto the upper envelope of
the chords -- the piecewise-linear interpolant of ``Psi_g`` -- and beyond the
outermost breakpoints the envelope extrapolates at the end chords' slopes, so
the curve imposes no hidden activity bound.

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
"""

import logging
import warnings

from linopy.constants import BREAKPOINT_DIM, EvolvingAPIWarning
from linopy.piecewise import tangent_lines
import numpy as np
import pandas as pd
import pypsa
import xarray as xr

logger = logging.getLogger(__name__)

# Component -> (carriers, baseline column, commodity-identifying columns). The
# commodity columns say *what* is produced; the granularity setting adds the
# columns saying *where*, so both are needed to form a group key.
COMPONENTS = {
    "crops": (("crop_production",), "baseline_area_mha", ("crop",)),
    "multi_crops": (("crop_production_multi",), "baseline_area_mha", ("combination",)),
    "grassland": (("grassland_production",), "baseline_area_mha", ()),
    "animals": (("animal_production",), "baseline_feed_use_mt_dm", ("product",)),
}

# Spatial columns appended to a component's commodity key at each granularity.
# ``link`` is the exception: it resolves every production link separately, so it
# replaces the whole key rather than extending it (a link already identifies its
# commodity and its location), and is represented by an empty tuple.
GRANULARITY_COLUMNS = {
    "country": ("country",),
    "region": ("country", "region"),
    "link": None,
}


def evaluate_supply_response_cost(n: pypsa.Network) -> float:
    """Total deviation cost the curves contributed to the objective.

    Each component's cost variable enters the objective directly, so the term
    is just the sum of their solutions, and the objective breakdown can report
    it as its own category with an exact identity check.
    """
    m = getattr(n, "model", None)
    if m is None:
        return 0.0
    return sum(
        float(m.variables[name].solution.sum())
        for component in COMPONENTS
        if (name := f"{component}_supply_response_cost") in m.variables
    )


def extract_supply_response_intercepts(n: pypsa.Network) -> pd.DataFrame:
    """Per-group price wedges from a baseline-pinned solve (PMP phase 1).

    The dual of each group's pinning constraint is the gap between the marginal
    value of the group's output and its marginal cost at the observed activity.
    Returned as a frame with columns ``component``, ``sr_group``, ``intercept``;
    fed back via ``supply_response.intercepts``, the wedges make the observed
    allocation the exact optimum of the unpinned model.
    """
    m = n.model
    frames = []
    for component in COMPONENTS:
        name = f"GlobalConstraint-{component}_supply_response_balance"
        if name not in m.constraints:
            continue
        dual = m.constraints[name].dual
        slack = (
            m.variables[f"{component}_supply_response_pin_slack_up"].solution
            - m.variables[f"{component}_supply_response_pin_slack_down"].solution
        )
        frames.append(
            pd.DataFrame(
                {
                    "component": component,
                    "sr_group": dual.coords["sr_group"].values,
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
            "no supply-response pinning constraints found; intercepts can only "
            "be extracted from a solve with supply_response.pin_baseline: true"
        )
    table = pd.concat(frames, ignore_index=True)
    slacked = table["slack"].abs() > 1e-6
    if slacked.any():
        worst = table.loc[slacked].reindex(
            table.loc[slacked, "slack"].abs().sort_values(ascending=False).index
        )
        logger.warning(
            "%d of %d pinned supply-response groups used pin slack (total "
            "|slack| %.3f); their intercepts are censored at the slack price. "
            "Worst: %s",
            int(slacked.sum()),
            len(table),
            float(table["slack"].abs().sum()),
            worst[["component", "sr_group", "slack"]].head(5).to_dict("records"),
        )
    return table


def add_supply_response_curves(n: pypsa.Network, cfg: dict) -> None:
    """Add convex piecewise-linear supply curves to the objective.

    Parameters
    ----------
    n : pypsa.Network
        Network whose ``model`` has already been created.
    cfg : dict
        The resolved ``supply_response`` configuration block.
    """
    if not cfg["enabled"]:
        return

    n_blocks = int(cfg["n_blocks"])
    if n_blocks < 1:
        raise ValueError(f"supply_response.n_blocks must be >= 1, got {n_blocks}")
    expansion = float(cfg["expansion_range"])
    contraction = float(cfg["contraction_range"])
    if expansion <= 0:
        raise ValueError(
            f"supply_response.expansion_range must be > 0, got {expansion}"
        )
    if not 0 < contraction <= 1:
        raise ValueError(
            "supply_response.contraction_range must satisfy 0 < r <= 1, got "
            f"{contraction}"
        )
    factor = float(cfg["elasticity_factor"])
    if factor <= 0:
        raise ValueError(f"supply_response.elasticity_factor must be > 0, got {factor}")
    width_growth = float(cfg["width_growth"])
    if width_growth < 1:
        raise ValueError(
            f"supply_response.width_growth must be >= 1, got {width_growth}; "
            "below 1 the breakpoints would spread toward the baseline, putting "
            "the coarsest resolution where groups actually sit"
        )

    granularity = str(cfg["granularity"])
    if granularity not in GRANULARITY_COLUMNS:
        raise ValueError(
            f"supply_response.granularity must be one of "
            f"{sorted(GRANULARITY_COLUMNS)}, got '{granularity}'"
        )
    spatial_cols = GRANULARITY_COLUMNS[granularity]

    pin = bool(cfg["pin_baseline"])
    if pin and cfg["intercepts"]:
        raise ValueError(
            "supply_response.pin_baseline and supply_response.intercepts are "
            "mutually exclusive: the pinned solve fits the intercepts against "
            "the accounting costs, so it must not itself carry intercepts"
        )
    intercepts = None
    if cfg["intercepts"]:
        table = pd.read_csv(cfg["intercepts"])
        intercepts = table.set_index(["component", "sr_group"])["intercept"]

    for component, (carriers, baseline_col, commodity_cols) in COMPONENTS.items():
        if not cfg["components"][component]:
            continue
        elasticity = float(cfg["elasticities"][component]) * factor
        if elasticity <= 0:
            raise ValueError(
                f"supply_response.elasticities.{component} must be > 0 after "
                f"applying elasticity_factor, got {elasticity}"
            )
        component_intercepts = None
        if intercepts is not None:
            if component not in intercepts.index.get_level_values("component"):
                raise ValueError(
                    f"supply_response.intercepts has no entries for component "
                    f"'{component}'; the intercepts must come from a pinned "
                    "solve covering the same components"
                )
            component_intercepts = intercepts.xs(component, level="component")
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
    return links[list(group_cols)].astype(str).agg("::".join, axis=1)


def _group_table(
    links: pd.DataFrame, baseline_col: str, group_cols: tuple[str, ...]
) -> pd.DataFrame:
    """Aggregate per-link baselines and costs into per-group curve parameters.

    Returns a frame indexed by the group key with the group's baseline activity
    and its baseline-weighted mean marginal cost. Groups with no baseline
    activity are dropped: they have no observed point to build a curve through,
    and the growth caps already govern whether an activity absent from the
    baseline may appear at all.
    """
    work = links.copy()
    work["_baseline"] = work[baseline_col].fillna(0.0).astype(float)
    work["_group"] = _group_key(work, group_cols)
    work = work[work["_baseline"] > 0]
    if work.empty:
        return pd.DataFrame(columns=["baseline", "cost"])

    work["_cost_weighted"] = work["_baseline"] * work["marginal_cost"].astype(float)
    grouped = work.groupby("_group").agg(
        baseline=("_baseline", "sum"), _cost_weighted=("_cost_weighted", "sum")
    )
    grouped["cost"] = grouped["_cost_weighted"] / grouped["baseline"]
    return grouped[["baseline", "cost"]].sort_index()


def _add_component_curves(
    n: pypsa.Network,
    *,
    component: str,
    carriers: tuple[str, ...],
    baseline_col: str,
    group_cols: tuple[str, ...],
    elasticity: float,
    n_blocks: int,
    expansion: float,
    contraction: float,
    width_growth: float,
    pin: bool,
    pin_slack_cost: float,
    intercepts: "pd.Series | None",
) -> None:
    """Add one component's curves (or its phase-1 pin) to the model."""
    links_df = n.links.static
    links = links_df[links_df["carrier"].isin(carriers)]
    if links.empty or baseline_col not in links.columns:
        logger.info(
            "No %s links with %s; skipping supply-response curves",
            "/".join(carriers),
            baseline_col,
        )
        return

    groups = _group_table(links, baseline_col, group_cols)
    if groups.empty:
        logger.info("No %s groups with positive baseline; skipping", component)
        return

    # Only links whose group carries a baseline participate; the rest are
    # governed by the growth caps.
    work = links.copy()
    work["_group"] = _group_key(work, group_cols)
    work = work[work["_group"].isin(groups.index)]

    m = n.model
    link_p = m.variables["Link-p"].sel(snapshot="now")
    group_map = xr.DataArray(
        work["_group"].to_numpy(),
        coords={"name": work.index},
        dims="name",
        name="sr_group",
    )
    activity = link_p.sel(name=work.index).groupby(group_map).sum()
    # groupby orders groups lexicographically, matching the sorted group table.
    group_index = pd.Index(activity.coords["sr_group"].values)
    if not group_index.equals(groups.index):
        groups = groups.reindex(group_index)

    if pin:
        # PMP phase 1: hold every group at its observed activity. The dual of
        # this constraint is the group's price wedge, extracted after the solve
        # by extract_supply_response_intercepts. The pin is elastic -- slack at
        # a high finite price -- because the reference data is never perfectly
        # consistent with every other constraint (land availability, water,
        # feed balances), and a hard pin would make the whole solve infeasible
        # over the smallest such residual. Where slack is used the dual is
        # censored at the slack price, which flags the group for investigation
        # instead of hiding the inconsistency.
        slack_cost = float(pin_slack_cost)
        if slack_cost <= 0:
            raise ValueError(
                f"supply_response.pin_slack_cost must be > 0, got {slack_cost}"
            )
        baseline_da = xr.DataArray(
            groups["baseline"].to_numpy(),
            coords={"sr_group": group_index.to_numpy()},
            dims="sr_group",
        )
        slack_coords = {"sr_group": group_index.to_numpy()}
        slack_up = m.add_variables(
            lower=0,
            coords=slack_coords,
            dims=("sr_group",),
            name=f"{component}_supply_response_pin_slack_up",
        )
        slack_down = m.add_variables(
            lower=0,
            coords=slack_coords,
            dims=("sr_group",),
            name=f"{component}_supply_response_pin_slack_down",
        )
        m.add_constraints(
            activity - slack_up + slack_down == baseline_da,
            name=f"GlobalConstraint-{component}_supply_response_balance",
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
                f"supply_response.intercepts (e.g. {missing}); the intercepts "
                "must come from a pinned solve of the same model at the same "
                "granularity"
            )

    bad = groups["cost"] <= 0
    if bad.any():
        raise ValueError(
            f"{int(bad.sum())} {component} group(s) have a non-positive "
            f"baseline-weighted marginal cost (e.g. "
            f"{group_index[bad][:5].tolist()}); the supply-curve slope "
            "gamma = cost / (elasticity * baseline) is undefined. Marginal "
            "costs come from the cost data and the cost calibration, so a zero "
            "here is an upstream data problem rather than something to default "
            "away."
        )

    # Slope of the marginal-cost curve, in objective units per activity
    # squared. Stated at the accounting cost, not the calibrated marginal cost
    # c + lambda: under a pinned diet the wedges trace the Ricardian rent
    # gradient, and sub-marginal groups at the extensive margin have
    # c + lambda <= 0, where a calibrated-cost slope is undefined. The
    # intercept alone carries exactness -- the observed optimum only needs the
    # marginal cost *level* at the baseline to match the marginal value there
    # -- so the slope choice affects only the stiffness of the response, which
    # for a group with wedge lambda realises elasticity
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

    # The deviation cost lambda * d + slope / 2 * d^2 sampled on the grid; one
    # free variable per group bounded below by the chords between consecutive
    # samples. Minimisation presses it onto the chords' upper envelope -- the
    # piecewise-linear interpolant of the curve, kinked at the baseline point
    # so the observed activity is strictly optimal -- and beyond the outermost
    # samples the envelope extrapolates at the end chords' slopes, imposing no
    # hidden activity bound.
    deviation = groups["baseline"].to_numpy()[:, None] * offsets[None, :]
    bp_coords = {
        "sr_group": group_index.to_numpy(),
        BREAKPOINT_DIM: np.arange(len(offsets)),
    }
    bp_dims = ("sr_group", BREAKPOINT_DIM)
    x_points = xr.DataArray(
        groups["baseline"].to_numpy()[:, None] + deviation, bp_coords, bp_dims
    )
    y_points = xr.DataArray(
        lam.to_numpy()[:, None] * deviation
        + slope.to_numpy()[:, None] / 2 * deviation**2,
        bp_coords,
        bp_dims,
    )
    cost = m.add_variables(
        coords={"sr_group": group_index.to_numpy()},
        dims=("sr_group",),
        name=f"{component}_supply_response_cost",
    )
    with warnings.catch_warnings():
        # The piecewise API is marked evolving upstream; our linopy fork is
        # pinned by tag, so any evolution arrives through a deliberate bump.
        warnings.simplefilter("ignore", EvolvingAPIWarning)
        chords = tangent_lines(activity, x_points, y_points)
    m.add_constraints(
        cost >= chords,
        name=f"GlobalConstraint-{component}_supply_response_curve",
    )
    m.objective += cost.sum()

    logger.info(
        "Added %s supply-response curves: %d groups x %d breakpoints/side "
        "(elasticity=%.3f, baseline=%.1f, slope median=%.4g, intercept "
        "median=%.4g)",
        component,
        len(group_index),
        n_blocks,
        elasticity,
        float(groups["baseline"].sum()),
        float(slope.median()),
        float(lam.median()),
    )
