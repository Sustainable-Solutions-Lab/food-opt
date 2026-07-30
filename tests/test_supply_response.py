# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Tests for the convex piecewise-linear supply-response curves.

Deliberately formulation-agnostic: each test drives a small LP through a real
HiGHS solve and probes the resulting dispatch, so the suite pins the economic
contract (the curve delivers the configured elasticity) rather than the tranche
bookkeeping that delivers it.

The fixture puts a single crop group in front of a fixed output price. At a
price equal to the links' marginal cost the observed activity is optimal, and
raising the price by a fraction ``s`` should expand activity by ``eta * s`` of
baseline -- that is what an own-price elasticity of ``eta`` means.
"""

import numpy as np
import pypsa
import pytest

from workflow.scripts.solve_model.supply_response import add_supply_response_curves

BASELINE_MHA = (6.0, 4.0)  # two links in one (crop, country) group
COST = 2.0  # bnUSD per Mha, equal on both links


def _cfg(
    *,
    enabled: bool = True,
    # Fine enough that tranche quantisation (activity moves in whole tranches
    # of expansion_range * baseline / n_blocks) stays well inside the
    # tolerances below, so the assertions measure the curve and not the
    # piecewise approximation of it.
    n_blocks: int = 200,
    expansion_range: float = 0.5,
    contraction_range: float = 1.0,
    width_growth: float = 1.0,
    granularity: str = "country",
    crops: float = 0.5,
    elasticity_factor: float = 1.0,
    components: dict | None = None,
) -> dict:
    return {
        "enabled": enabled,
        "n_blocks": n_blocks,
        "expansion_range": expansion_range,
        "contraction_range": contraction_range,
        "width_growth": width_growth,
        "granularity": granularity,
        "elasticities": {"crops": crops, "grassland": 0.3, "animals": 0.4},
        "elasticity_factor": elasticity_factor,
        "components": components
        or {"crops": True, "grassland": False, "animals": False},
    }


def _make_network(
    price: float,
    *,
    baselines: tuple[float, ...] = BASELINE_MHA,
    cost: float = COST,
) -> pypsa.Network:
    """One crop group of ``len(baselines)`` links selling output at ``price``.

    Yield is 1 Mt per Mha so per-Mha and per-Mt quantities coincide and the
    elasticity arithmetic stays readable. Land is abundant, so the supply curve
    is the only thing that bounds expansion.
    """
    n = pypsa.Network()
    n.set_snapshots(["now"])
    n.carriers.add(["land", "crop_wheat", "crop_production", "sale"], unit="Mha")
    n.buses.add(["land:R", "crop:wheat:USA", "sink:USA"], carrier="land")

    n.generators.add("supply:land:R", bus="land:R", p_nom=1e4, marginal_cost=0.0)
    n.stores.add("sink", bus="sink:USA", e_nom=1e6, e_initial=0.0, e_cyclic=False)

    names = [f"produce:wheat:R{i}" for i in range(len(baselines))]
    n.links.add(
        names,
        bus0="land:R",
        bus1="crop:wheat:USA",
        carrier="crop_production",
        efficiency=1.0,
        p_nom=1e4,
        marginal_cost=cost,
        baseline_area_mha=list(baselines),
        crop="wheat",
        country="USA",
    )
    # Negative marginal cost = revenue per Mt sold.
    n.links.add(
        "sell:wheat:USA",
        bus0="crop:wheat:USA",
        bus1="sink:USA",
        carrier="sale",
        efficiency=1.0,
        p_nom=1e4,
        marginal_cost=-price,
    )
    return n


def _solve(n: pypsa.Network, cfg: dict) -> float:
    """Create the model, add curves, solve, and return total group activity (Mha)."""
    n.optimize.create_model(include_objective_constant=False)
    add_supply_response_curves(n, cfg)
    status, condition = n.model.solve(solver_name="highs")
    assert (status, condition) == ("ok", "optimal"), (status, condition)
    sol = n.model.variables["Link-p"].solution.sel(snapshot="now")
    produce = [
        str(v) for v in sol.coords["name"].values if str(v).startswith("produce")
    ]
    return float(sol.sel(name=produce).sum().item())


def test_baseline_is_optimal_at_reference_price():
    """At price == marginal cost the curve's minimum is the observed activity."""
    total = _solve(_make_network(COST), _cfg())
    assert total == pytest.approx(sum(BASELINE_MHA), rel=1e-6)


@pytest.mark.parametrize("eta", [0.25, 0.5, 1.0])
def test_price_shock_delivers_the_configured_elasticity(eta):
    """A price rise of s expands activity by eta * s of baseline."""
    shock = 0.1
    baseline = sum(BASELINE_MHA)
    total = _solve(_make_network(COST * (1 + shock)), _cfg(crops=eta))
    arc_elasticity = ((total - baseline) / baseline) / shock
    assert arc_elasticity == pytest.approx(eta, rel=0.05)


def test_price_drop_contracts_symmetrically():
    """The curve is symmetric about baseline, so a price cut contracts likewise."""
    shock = 0.1
    baseline = sum(BASELINE_MHA)
    total = _solve(_make_network(COST * (1 - shock)), _cfg(crops=0.5))
    arc_elasticity = ((baseline - total) / baseline) / shock
    assert arc_elasticity == pytest.approx(0.5, rel=0.05)


def test_elasticity_factor_scales_the_response():
    """elasticity_factor multiplies the elasticity, so it scales the response."""
    shock = 0.1
    baseline = sum(BASELINE_MHA)
    plain = _solve(_make_network(COST * (1 + shock)), _cfg(crops=0.4))
    scaled = _solve(
        _make_network(COST * (1 + shock)), _cfg(crops=0.4, elasticity_factor=2.0)
    )
    assert (scaled - baseline) == pytest.approx(2.0 * (plain - baseline), rel=0.05)


def test_response_is_monotone_in_the_price():
    """Graduated response: each further price rise expands activity further."""
    totals = [
        _solve(_make_network(COST * (1 + s)), _cfg()) for s in (0.0, 0.05, 0.1, 0.2)
    ]
    assert all(np.diff(totals) > 0)


def test_extreme_price_stays_feasible_beyond_the_expansion_range():
    """The outermost tranche is unbounded, so a large shock extrapolates."""
    baseline = sum(BASELINE_MHA)
    total = _solve(_make_network(COST * 20), _cfg(expansion_range=0.5))
    assert total > baseline * (1 + 0.5)


def test_zero_baseline_group_gets_no_curve():
    """A group with no observed activity has no point to build a curve through."""
    n = _make_network(COST * 2, baselines=(0.0, 0.0))
    n.optimize.create_model(include_objective_constant=False)
    add_supply_response_curves(n, _cfg())
    assert not [k for k in n.model.variables if "supply_response" in k]


def test_disabled_config_is_a_no_op():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    add_supply_response_curves(n, _cfg(enabled=False))
    assert not [k for k in n.model.variables if "supply_response" in k]


def test_non_positive_marginal_cost_raises():
    """The slope is undefined at zero cost; that is an upstream data problem."""
    n = _make_network(COST, cost=0.0)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="non-positive"):
        add_supply_response_curves(n, _cfg())


def test_invalid_block_count_raises():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="n_blocks"):
        add_supply_response_curves(n, _cfg(n_blocks=0))


def test_geometric_widths_reach_the_elasticity_with_far_fewer_blocks():
    """Narrow near tranches resolve a small move that equal widths cannot.

    This is the reason `width_growth` exists. At 4 equal tranches over a range of
    half the baseline, nothing finer than 12.5% of baseline is resolved, so a 5%
    target response is quantised to a whole tranche and the arc elasticity comes
    out badly wrong. Growing the widths puts the resolution near the baseline.
    """
    shock = 0.1
    eta = 0.5
    baseline = sum(BASELINE_MHA)
    uniform = _solve(_make_network(COST * (1 + shock)), _cfg(n_blocks=4, crops=eta))
    geometric = _solve(
        _make_network(COST * (1 + shock)),
        _cfg(n_blocks=4, width_growth=3.0, crops=eta),
    )
    target = eta * shock * baseline
    err_uniform = abs((uniform - baseline) - target)
    err_geometric = abs((geometric - baseline) - target)
    assert err_geometric < err_uniform


def test_width_growth_one_matches_equal_widths():
    """The default reproduces the equal-width formulation exactly."""
    a = _solve(_make_network(COST * 1.1), _cfg(n_blocks=8, width_growth=1.0))
    b = _solve(_make_network(COST * 1.1), _cfg(n_blocks=8))
    assert a == pytest.approx(b, rel=1e-9)


def test_link_granularity_curves_each_link_and_keeps_the_elasticity():
    """Per-link curves still deliver the configured elasticity in aggregate.

    Each link gets slope ``c / (eta * b_link)``, so it carries the elasticity
    with respect to its own margin; with equal costs across the group the
    aggregate carries the same elasticity, which is what lets granularity be
    changed without re-tuning it.
    """
    shock = 0.1
    baseline = sum(BASELINE_MHA)
    total = _solve(_make_network(COST * (1 + shock)), _cfg(granularity="link"))
    arc_elasticity = ((total - baseline) / baseline) / shock
    assert arc_elasticity == pytest.approx(0.5, rel=0.05)


def test_link_granularity_builds_one_group_per_link():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    add_supply_response_curves(n, _cfg(granularity="link"))
    var = n.model.variables["crops_supply_response_up"]
    assert var.sizes["sr_group"] == len(BASELINE_MHA)


def test_unknown_granularity_raises():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="granularity"):
        add_supply_response_curves(n, _cfg(granularity="continent"))


def test_width_growth_below_one_raises():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="width_growth"):
        add_supply_response_curves(n, _cfg(width_growth=0.5))


def test_contraction_range_above_one_raises():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="contraction_range"):
        add_supply_response_curves(n, _cfg(contraction_range=1.5))
