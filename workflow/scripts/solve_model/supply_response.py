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

For a group ``g`` with observed activity ``b_g`` and marginal cost ``c_g`` at
that activity, a linear marginal-cost curve reproducing supply elasticity
``eta`` has slope

    gamma_g = c_g / (eta * b_g),

since equating marginal cost to the marginal value of output gives
``dX/dv = 1 / gamma`` and hence ``eta = v / (gamma * b)`` with ``v = c_g`` at the
calibrated point. The cost calibration already aligns each link's marginal cost
with the marginal value of its output at baseline, so ``c_g`` is read off the
network rather than supplied separately.

The curve enters the objective as its deviation cost, the integral of the slope
from the baseline,

    Psi_g(X) = gamma_g / 2 * (X - b_g)^2,

which is convex and minimised at ``X = b_g``, so the observed activity stays
optimal at reference prices. ``Psi_g`` is approximated by ``n_blocks`` tranches
per direction, each priced at the slope times its midpoint deviation. Tranche
prices increase with distance from the baseline, so a cost-minimising solution
fills the near tranches first and no ordering constraints are needed. The last
tranche in each direction is left unbounded, which extrapolates the curve
linearly at its top price and keeps the formulation feasible whatever the
growth caps allow.

Because activity moves in whole tranches, the tranche width sets the finest
response the curve can express. With equal widths and a range of half the
baseline, four blocks resolve nothing below 12.5% of baseline -- coarser than
the moves most groups make, so the curve collapses into a flat-rate penalty and
loses the graduated response it exists for. ``width_growth`` above 1 narrows the
near tranches so resolution concentrates at the baseline where groups sit, while
the outer tranches still resolve large moves; that is far cheaper than the many
equal tranches the same near-baseline resolution would otherwise need.

Granularity: curves apply to *groups* -- ``(crop, country)``, ``(country)`` for
grassland, ``(product, country)`` for animals -- matching the units the cost
calibration extracts corrections for. A group curve prices changes in a
country's total activity for a commodity; it says nothing about which regions,
resource classes or water-supply types that activity sits in. Reallocation
within a group is a distinct friction (land is heterogeneous, and rotations and
infrastructure are local), so it stays with the per-link deviation penalty,
which the two mechanisms are meant to share: the curve carries the elasticity,
the per-link term carries spatial inertia.

Multi-cropping links are not covered yet. Their duals price a joint cycle
bundle rather than one constituent crop, so they need either their own
``(combination, country)`` curves or a rule for splitting a bundle's activity
across the crop curves it contributes to; until that is settled they keep the
per-link deviation penalty alone.
"""

import logging

import numpy as np
import pandas as pd
import pypsa
import xarray as xr

logger = logging.getLogger(__name__)

# Tranche prices by variable name, kept so the post-solve cost evaluation can
# recover the objective term without rebuilding the curve parameters. Populated
# by add_supply_response_curves and reset on each call, since one process may
# build several models in sequence (the calibration drivers do).
_PRICES: dict[str, "xr.DataArray"] = {}

# Component -> (carriers, baseline column, group-key columns).
COMPONENTS = {
    "crops": (("crop_production",), "baseline_area_mha", ("crop", "country")),
    "grassland": (("grassland_production",), "baseline_area_mha", ("country",)),
    "animals": (
        ("animal_production",),
        "baseline_feed_use_mt_dm",
        ("product", "country"),
    ),
}


def evaluate_supply_response_cost(n: pypsa.Network) -> float:
    """Total deviation cost the curves contributed to the objective.

    Recovers the term from the solved tranche variables and the prices stashed
    when the curves were built, so the objective breakdown can report it as its
    own category and its identity check stays exact.
    """
    m = getattr(n, "model", None)
    if m is None:
        return 0.0
    total = 0.0
    for component in COMPONENTS:
        for direction in ("up", "down"):
            name = f"{component}_supply_response_{direction}"
            if name not in m.variables:
                continue
            prices = _PRICES.get(name)
            if prices is None:
                raise ValueError(
                    f"supply-response variable '{name}' is in the model but its "
                    "tranche prices were not recorded, so its objective "
                    "contribution cannot be evaluated. This means the curves "
                    "were built by a different call than the one being "
                    "evaluated."
                )
            total += float((prices * m.variables[name].solution).sum())
    return total


def add_supply_response_curves(n: pypsa.Network, cfg: dict) -> None:
    """Add convex piecewise-linear supply curves to the objective.

    Parameters
    ----------
    n : pypsa.Network
        Network whose ``model`` has already been created.
    cfg : dict
        The resolved ``supply_response`` configuration block.
    """
    _PRICES.clear()
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
            "below 1 the tranches would widen toward the baseline, putting the "
            "coarsest resolution where groups actually sit"
        )

    for component, (carriers, baseline_col, group_cols) in COMPONENTS.items():
        if not cfg["components"][component]:
            continue
        elasticity = float(cfg["elasticities"][component]) * factor
        if elasticity <= 0:
            raise ValueError(
                f"supply_response.elasticities.{component} must be > 0 after "
                f"applying elasticity_factor, got {elasticity}"
            )
        _add_component_curves(
            n,
            component=component,
            carriers=carriers,
            baseline_col=baseline_col,
            group_cols=group_cols,
            elasticity=elasticity,
            n_blocks=n_blocks,
            expansion=expansion,
            contraction=contraction,
            width_growth=width_growth,
        )


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
    work["_group"] = work[list(group_cols)].astype(str).agg("::".join, axis=1)
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
) -> None:
    """Add the tranche variables, balance rows and objective terms for one component."""
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

    zero_cost = groups["cost"] <= 0
    if zero_cost.any():
        raise ValueError(
            f"{int(zero_cost.sum())} {component} group(s) have a non-positive "
            f"baseline-weighted marginal cost (e.g. "
            f"{groups.index[zero_cost][:5].tolist()}); the supply-curve slope "
            "gamma = cost / (elasticity * baseline) is undefined. Marginal costs "
            "come from the cost data and the cost calibration, so a zero here is "
            "an upstream data problem rather than something to default away."
        )

    # Slope of the marginal-cost curve, in objective units per activity squared.
    slope = groups["cost"] / (elasticity * groups["baseline"])

    # Only links whose group carries a baseline participate; the rest are
    # governed by the growth caps and the per-link deviation penalty.
    work = links.copy()
    work["_group"] = work[list(group_cols)].astype(str).agg("::".join, axis=1)
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
        slope = slope.reindex(group_index)

    coords = {
        "sr_group": group_index.to_numpy(),
        "sr_block": np.arange(1, n_blocks + 1),
    }
    dims = ("sr_group", "sr_block")

    # Relative tranche widths. At width_growth == 1 these are all equal, which
    # spends resolution uniformly over the range; above 1 the near tranches are
    # narrower, concentrating resolution at the baseline where nearly every
    # group sits while the outer tranches still resolve large moves.
    shares = width_growth ** np.arange(n_blocks)
    shares = shares / shares.sum()

    def _tranches(range_fraction: float, direction: str):
        """Build tranche variables and their objective prices for one direction.

        Widths split the group's directional range by ``shares``; each tranche is
        priced at the slope times its midpoint deviation -- the width accumulated
        before it plus half its own -- so the piecewise cost traces the quadratic
        whatever the widths are. Prices rise with distance from the baseline, so
        a cost-minimising solution fills the near tranches first and convexity
        needs no ordering rows. The outermost tranche is unbounded so the curve
        extrapolates at its top price instead of imposing a hidden bound.
        """
        w = range_fraction * groups["baseline"].to_numpy()[:, None] * shares[None, :]
        upper = w.copy()
        upper[:, -1] = np.inf
        midpoints = np.cumsum(w, axis=1) - w / 2
        prices = slope.to_numpy()[:, None] * midpoints
        var = m.add_variables(
            lower=0,
            upper=xr.DataArray(upper, coords=coords, dims=dims),
            coords=coords,
            dims=dims,
            name=f"{component}_supply_response_{direction}",
        )
        price_da = xr.DataArray(prices, coords=coords, dims=dims)
        _PRICES[f"{component}_supply_response_{direction}"] = price_da
        return var, price_da

    up, up_price = _tranches(expansion, "up")
    down, down_price = _tranches(contraction, "down")

    baseline_da = xr.DataArray(
        groups["baseline"].to_numpy(),
        coords={"sr_group": group_index.to_numpy()},
        dims="sr_group",
    )
    m.add_constraints(
        activity - up.sum("sr_block") + down.sum("sr_block") == baseline_da,
        name=f"GlobalConstraint-{component}_supply_response_balance",
    )
    m.objective += (up_price * up).sum() + (down_price * down).sum()

    logger.info(
        "Added %s supply-response curves: %d groups x %d blocks/direction "
        "(elasticity=%.3f, baseline=%.1f, slope median=%.4g, first up-tranche "
        "price median=%.4g)",
        component,
        len(group_index),
        n_blocks,
        elasticity,
        float(groups["baseline"].sum()),
        float(slope.median()),
        float(up_price.isel(sr_block=0).median()),
    )
