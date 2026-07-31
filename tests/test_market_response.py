# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Tests for the convex piecewise-linear market-response curves.

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
import pandas as pd
import pypsa
import pytest

from workflow.scripts.solve_model.market_response import (
    add_market_response_curves,
    extract_market_response_intercepts,
)

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
    pin_baseline: bool = False,
    pin_slack_cost: float = 10.0,
    intercepts: str | None = None,
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
        "pin_baseline": pin_baseline,
        "pin_slack_cost": pin_slack_cost,
        "intercepts": intercepts,
        "elasticities": {
            "crops": crops,
            "multi_crops": crops,
            "grassland": 0.3,
            "animals": 0.4,
            "demand": 0.4,
        },
        "elasticity_factor": elasticity_factor,
        "components": components
        or {
            "crops": True,
            "multi_crops": False,
            "grassland": False,
            "animals": False,
            "demand": False,
        },
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
    add_market_response_curves(n, cfg)
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
    add_market_response_curves(n, _cfg())
    assert not [k for k in n.model.variables if "market_response" in k]


def test_disabled_config_is_a_no_op():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    add_market_response_curves(n, _cfg(enabled=False))
    assert not [k for k in n.model.variables if "market_response" in k]


def test_non_positive_marginal_cost_raises():
    """The slope is undefined at zero cost; that is an upstream data problem."""
    n = _make_network(COST, cost=0.0)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="non-positive"):
        add_market_response_curves(n, _cfg())


def test_invalid_block_count_raises():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="n_blocks"):
        add_market_response_curves(n, _cfg(n_blocks=0))


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
    add_market_response_curves(n, _cfg(granularity="link"))
    var = n.model.variables["crops_market_response_cost"]
    assert var.sizes["mr_group"] == len(BASELINE_MHA)


def test_unknown_granularity_raises():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="granularity"):
        add_market_response_curves(n, _cfg(granularity="continent"))


def _fit_intercepts(price: float, cfg: dict, tmp_path) -> str:
    """PMP phase 1 on the fixture: pinned solve, wedges written to a CSV."""
    n = _make_network(price)
    n.optimize.create_model(include_objective_constant=False)
    add_market_response_curves(n, {**cfg, "pin_baseline": True, "intercepts": None})
    status, condition = n.model.solve(solver_name="highs")
    assert (status, condition) == ("ok", "optimal"), (status, condition)
    path = tmp_path / "intercepts.csv"
    extract_market_response_intercepts(n).to_csv(path, index=False)
    return str(path)


def test_pinned_solve_holds_activity_at_baseline():
    """Phase 1 fixes the group at its observed activity whatever the price."""
    n = _make_network(COST * 3)
    n.optimize.create_model(include_objective_constant=False)
    add_market_response_curves(n, _cfg(pin_baseline=True))
    status, condition = n.model.solve(solver_name="highs")
    assert (status, condition) == ("ok", "optimal")
    sol = n.model.variables["Link-p"].solution.sel(snapshot="now")
    produce = [
        str(v) for v in sol.coords["name"].values if str(v).startswith("produce")
    ]
    total = float(sol.sel(name=produce).sum().item())
    assert total == pytest.approx(sum(BASELINE_MHA), rel=1e-9)


def test_pinned_duals_measure_the_price_wedge(tmp_path):
    """The extracted intercept is the value-minus-cost gap at the baseline."""
    wedge = 0.4
    path = _fit_intercepts(COST + wedge, _cfg(), tmp_path)
    table = pd.read_csv(path)
    assert set(table["component"]) == {"crops"}
    assert table["intercept"].to_numpy() == pytest.approx(wedge, rel=1e-9)


def test_intercepts_reproduce_the_baseline_exactly(tmp_path):
    """PMP phase 2: with fitted intercepts, the observed activity is the exact
    optimum even though the output price sits well above the accounting cost."""
    price = COST * 1.5
    path = _fit_intercepts(price, _cfg(), tmp_path)
    total = _solve(_make_network(price), _cfg(intercepts=path))
    assert total == pytest.approx(sum(BASELINE_MHA), rel=1e-9)


def test_calibrated_curve_keeps_the_accounting_cost_slope(tmp_path):
    """The slope is stated at the accounting cost, so around the calibrated
    point a price rise of s * cost expands activity by eta * s of baseline."""
    eta, shock = 0.5, 0.1
    calibrated_price = COST * 1.5
    path = _fit_intercepts(calibrated_price, _cfg(crops=eta), tmp_path)
    baseline = sum(BASELINE_MHA)
    total = _solve(
        _make_network(calibrated_price + shock * COST),
        _cfg(crops=eta, intercepts=path),
    )
    arc_elasticity = ((total - baseline) / baseline) / shock
    assert arc_elasticity == pytest.approx(eta, rel=0.05)


def test_elastic_pin_survives_an_infeasible_reference_and_censors_the_dual():
    """A reference the model cannot reach uses pin slack instead of going
    infeasible, and the affected group's intercept is censored at the slack
    price so the inconsistency is visible rather than hidden."""
    slack_cost = 10.0
    n = _make_network(COST)
    # Capacity below the observed activity: the pin cannot be met exactly.
    n.links.static.loc[n.links.static.index.str.startswith("produce"), "p_nom"] = 3.0
    n.optimize.create_model(include_objective_constant=False)
    add_market_response_curves(n, _cfg(pin_baseline=True, pin_slack_cost=slack_cost))
    status, condition = n.model.solve(solver_name="highs")
    assert (status, condition) == ("ok", "optimal")
    table = extract_market_response_intercepts(n)
    assert table["slack"].to_numpy() == pytest.approx(-4.0)  # 6 capped at 3, x2
    # Activity is stuck below the pin, so at the margin the pinned unit is
    # worth less than its cost by at least the slack price: censored at
    # -slack_cost (an unreachable overshoot would censor at +slack_cost).
    assert table["intercept"].to_numpy() == pytest.approx(-slack_cost)


def test_pin_and_intercepts_are_mutually_exclusive(tmp_path):
    path = _fit_intercepts(COST, _cfg(), tmp_path)
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="mutually exclusive"):
        add_market_response_curves(n, _cfg(pin_baseline=True, intercepts=path))


def test_intercepts_from_a_different_granularity_raise(tmp_path):
    """Group keys must match between fit and apply, so a granularity switch
    between phases fails loudly instead of silently zeroing the wedges."""
    path = _fit_intercepts(COST * 1.5, _cfg(granularity="link"), tmp_path)
    n = _make_network(COST * 1.5)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="no entry"):
        add_market_response_curves(n, _cfg(granularity="country", intercepts=path))


def test_width_growth_below_one_raises():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="width_growth"):
        add_market_response_curves(n, _cfg(width_growth=0.5))


def test_contraction_range_above_one_raises():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="contraction_range"):
        add_market_response_curves(n, _cfg(contraction_range=1.5))


# ---------------------------------------------------------------------------
# Demand component: concave marginal-utility curves on consumption links.

DEMAND_BASELINE_MT = 5.0
SUPPLY_PRICE = 1.5  # bnUSD per Mt of food offered to the consumption link


def _demand_cfg(**kwargs) -> dict:
    components = {
        "crops": False,
        "multi_crops": False,
        "grassland": False,
        "animals": False,
        "demand": True,
    }
    return _cfg(components=components, **kwargs)


def _make_demand_network(
    supply_price: float, *, baseline: float = DEMAND_BASELINE_MT
) -> pypsa.Network:
    """One (food, country) consumption link buying food at ``supply_price``.

    Food is available in any quantity at the supply price and nothing else
    values it, so consumption is driven entirely by the demand curve.
    """
    n = pypsa.Network()
    n.set_snapshots(["now"])
    n.carriers.add(["food", "food_consumption", "nutrient"], unit="Mt")
    n.buses.add(["food:apple:USA", "nutrient:cal:USA"], carrier="food")
    n.generators.add(
        "supply:apple:USA",
        bus="food:apple:USA",
        p_nom=1e4,
        marginal_cost=supply_price,
    )
    n.stores.add(
        "intake", bus="nutrient:cal:USA", e_nom=1e6, e_initial=0.0, e_cyclic=False
    )
    n.links.add(
        "consume:apple:USA",
        bus0="food:apple:USA",
        bus1="nutrient:cal:USA",
        carrier="food_consumption",
        efficiency=1.0,
        p_nom=1e4,
        marginal_cost=0.0,
        baseline_consumption_mt=baseline,
        food="apple",
        country="USA",
    )
    return n


def _solve_demand(n: pypsa.Network, cfg: dict) -> float:
    n.optimize.create_model(include_objective_constant=False)
    add_market_response_curves(n, cfg)
    status, condition = n.model.solve(solver_name="highs")
    assert (status, condition) == ("ok", "optimal"), (status, condition)
    sol = n.model.variables["Link-p"].solution.sel(snapshot="now")
    return float(sol.sel(name="consume:apple:USA").item())


def _fit_demand_intercepts(supply_price: float, cfg: dict, tmp_path) -> str:
    n = _make_demand_network(supply_price)
    n.optimize.create_model(include_objective_constant=False)
    add_market_response_curves(n, {**cfg, "pin_baseline": True, "intercepts": None})
    status, condition = n.model.solve(solver_name="highs")
    assert (status, condition) == ("ok", "optimal"), (status, condition)
    path = tmp_path / "demand_intercepts.csv"
    table = extract_market_response_intercepts(n)
    # In this fixture supply is a bare generator at a fixed price, so the
    # slope-basis solve would measure the same wedge; the calibration driver
    # normally stitches this column from a separate zero-intercept solve.
    table["slope_basis"] = table["intercept"].abs()
    table.to_csv(path, index=False)
    return str(path)


def test_demand_pin_holds_consumption_and_measures_willingness_to_pay(tmp_path):
    """Phase 1 pins intake; the intercept is minus the price paid at the
    margin, i.e. minus the willingness to pay the calibration imputes."""
    path = _fit_demand_intercepts(SUPPLY_PRICE, _demand_cfg(), tmp_path)
    table = pd.read_csv(path)
    assert set(table["component"]) == {"demand"}
    assert table["intercept"].to_numpy() == pytest.approx(-SUPPLY_PRICE, rel=1e-9)
    assert table["slack"].to_numpy() == pytest.approx(0.0, abs=1e-9)


def test_demand_intercepts_reproduce_baseline_consumption(tmp_path):
    """With fitted intercepts, observed intake is the exact optimum: the
    marginal utility at the baseline equals the supply price there."""
    path = _fit_demand_intercepts(SUPPLY_PRICE, _demand_cfg(), tmp_path)
    consumed = _solve_demand(
        _make_demand_network(SUPPLY_PRICE), _demand_cfg(intercepts=path)
    )
    assert consumed == pytest.approx(DEMAND_BASELINE_MT, rel=1e-6)


@pytest.mark.parametrize("eta", [0.25, 0.5])
def test_demand_curve_delivers_the_configured_elasticity(eta, tmp_path):
    """A price rise of s contracts intake by eta * s of baseline: the slope is
    stated at the fitted willingness to pay, so the realised own-price demand
    elasticity is eta itself."""
    shock = 0.1
    cfg = _demand_cfg()
    cfg["elasticities"]["demand"] = eta
    path = _fit_demand_intercepts(SUPPLY_PRICE, cfg, tmp_path)
    consumed = _solve_demand(
        _make_demand_network(SUPPLY_PRICE * (1 + shock)), {**cfg, "intercepts": path}
    )
    arc = ((consumed - DEMAND_BASELINE_MT) / DEMAND_BASELINE_MT) / shock
    assert arc == pytest.approx(-eta, rel=0.05)


def test_demand_curves_without_intercepts_raise():
    n = _make_demand_network(SUPPLY_PRICE)
    n.optimize.create_model(include_objective_constant=False)
    with pytest.raises(ValueError, match="demand curves require"):
        add_market_response_curves(n, _demand_cfg())


def test_sequential_pin_list_pins_demand_against_elastic_supply(tmp_path):
    """pin_baseline as a component list: production carries its fitted curves
    while demand alone is pinned, so the demand wedge is measured against the
    elastic (deployed) supply side rather than a closed pinned chain."""
    supply_cfg = _cfg(granularity="link")
    supply_path = _fit_intercepts(COST * 1.5, supply_cfg, tmp_path)

    n = _make_network(COST * 1.5)
    n.optimize.create_model(include_objective_constant=False)
    cfg = _cfg(
        granularity="link",
        intercepts=supply_path,
        components={
            "crops": True,
            "multi_crops": False,
            "grassland": False,
            "animals": False,
            "demand": False,
        },
    )
    cfg["pin_baseline"] = ["demand"]
    # No demand links in this fixture: the point is that the list form applies
    # curves to the unpinned crops component instead of raising the
    # pin-vs-intercepts exclusivity error, and reproduces the baseline.
    add_market_response_curves(n, cfg)
    status, condition = n.model.solve(solver_name="highs")
    assert (status, condition) == ("ok", "optimal")
    sol = n.model.variables["Link-p"].solution.sel(snapshot="now")
    produce = [
        str(v) for v in sol.coords["name"].values if str(v).startswith("produce")
    ]
    assert float(sol.sel(name=produce).sum().item()) == pytest.approx(
        sum(BASELINE_MHA), rel=1e-6
    )


def test_pin_list_with_unknown_component_raises():
    n = _make_network(COST)
    n.optimize.create_model(include_objective_constant=False)
    cfg = _cfg()
    cfg["pin_baseline"] = ["croops"]
    with pytest.raises(ValueError, match="unknown components"):
        add_market_response_curves(n, cfg)
