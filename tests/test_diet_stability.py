# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Unit tests for the diet deviation component.

Covers the post-hoc cost evaluator that the objective-breakdown extraction
relies on, and the hard-mode per-country diet churn budget via a minimal
two-food LP solved with HiGHS.
"""

import pandas as pd
import pypsa
import pytest

from workflow.scripts.solve_model.diet_stability import (
    add_diet_stability_constraints,
    evaluate_diet_stability_cost,
)


def _make_network(
    consumption_mt: list[float], baseline_mt: list[float]
) -> pypsa.Network:
    """Build a minimal pypsa.Network with food_consumption links and dispatch."""
    n = pypsa.Network()
    n.set_snapshots(["now"])
    names = [f"consume:f{i}:CTY" for i in range(len(consumption_mt))]
    n.add("Bus", "food:f0:CTY")
    for nm in names:
        n.links.static.loc[nm, "carrier"] = "food_consumption"
        n.links.static.loc[nm, "bus0"] = "food:f0:CTY"
    p0 = pd.DataFrame(
        {nm: [val] for nm, val in zip(names, consumption_mt)}, index=n.snapshots
    )
    n.links.dynamic["p0"] = p0
    matched = pd.DataFrame({"name": names, "target_mt": baseline_mt})
    return n, matched


def _dp_cfg(
    *,
    enabled: bool = True,
    diet_enabled: bool = True,
    penalty_mode: str = "l1",
    deviation_type: str = "absolute",
    l1_cost: float = 0.0,
    quadratic_cost: float = 0.0,
    max_relative_deviation: float = 1.0,
    enable_slack: bool = False,
    min_baseline: float = 1e-6,
) -> dict:
    return {
        "enabled": enabled,
        "penalty_mode": penalty_mode,
        "deviation_type": deviation_type,
        "quadratic_cost": quadratic_cost,
        "land": {
            "enabled": False,
            "crops": {"enabled": False, "l1_cost": 0.0, "l1_cost_factor": 1.0},
            "grassland": {"enabled": False, "l1_cost": 0.0, "l1_cost_factor": 1.0},
        },
        "feed": {"enabled": False, "l1_cost": 0.0, "l1_cost_factor": 1.0},
        "diet": {
            "enabled": diet_enabled,
            "l1_cost": l1_cost,
            "l1_cost_factor": 1.0,
            "max_relative_deviation": max_relative_deviation,
            "enable_slack": enable_slack,
            "min_baseline": min_baseline,
        },
    }


def test_disabled_returns_zero():
    n, matched = _make_network([1.0, 2.0], [1.0, 2.0])
    assert evaluate_diet_stability_cost(n, matched, _dp_cfg(enabled=False), 0.0) == 0.0
    assert (
        evaluate_diet_stability_cost(n, matched, _dp_cfg(diet_enabled=False), 0.0)
        == 0.0
    )


def test_l1_absolute_no_deviation():
    n, matched = _make_network([1.0, 2.0], [1.0, 2.0])
    cfg = _dp_cfg(penalty_mode="l1", deviation_type="absolute", l1_cost=10.0)
    assert evaluate_diet_stability_cost(n, matched, cfg, 0.0) == 0.0


def test_l1_absolute_symmetric_deviation():
    # Deviations: +0.5, -0.5 => |.|.sum() = 1.0
    n, matched = _make_network([1.5, 1.5], [1.0, 2.0])
    cfg = _dp_cfg(penalty_mode="l1", deviation_type="absolute", l1_cost=10.0)
    assert evaluate_diet_stability_cost(n, matched, cfg, 0.0) == 10.0


def test_l1_relative_uses_min_baseline_floor():
    # Baseline 0 with min_baseline floor 1e-3, consumption 0.01 -> rel dev = 10
    n, matched = _make_network([0.01], [0.0])
    cfg = _dp_cfg(
        penalty_mode="l1",
        deviation_type="relative",
        l1_cost=1.0,
        min_baseline=1e-3,
    )
    assert evaluate_diet_stability_cost(n, matched, cfg, 0.0) == 10.0


def test_quadratic_absolute():
    # Deviations: +1, -1 => (.)^2 .sum = 2 => 0.5 * 5 * 2 = 5
    n, matched = _make_network([2.0, 1.0], [1.0, 2.0])
    cfg = _dp_cfg(
        penalty_mode="quadratic",
        deviation_type="absolute",
        quadratic_cost=5.0,
    )
    assert evaluate_diet_stability_cost(n, matched, cfg, 5.0) == 5.0


def test_hard_without_model_returns_zero():
    n, matched = _make_network([1.0, 2.0], [1.0, 2.0])
    cfg = _dp_cfg(penalty_mode="hard")
    assert evaluate_diet_stability_cost(n, matched, cfg, 500.0) == 0.0


def _make_solvable_network() -> tuple[pypsa.Network, pd.DataFrame]:
    """Two foods, one country: baseline (f0, f1) = (2, 1) Mt, total demand 3 Mt.

    f1 is much cheaper to supply, so unconstrained the model sources all 3 Mt
    via f1; the diet churn budget limits how far consumption may shift.
    """
    n = pypsa.Network()
    for food, cost in [("f0", 10.0), ("f1", 1.0)]:
        n.add("Bus", f"food:{food}:CTY")
        n.add(
            "Generator",
            f"supply:{food}",
            bus=f"food:{food}:CTY",
            p_nom=10.0,
            marginal_cost=cost,
        )
        n.add(
            "Link",
            f"consume:{food}:CTY",
            bus0=f"food:{food}:CTY",
            bus1="sink",
            p_nom=10.0,
            carrier="food_consumption",
        )
    n.add("Bus", "sink")
    n.add("Load", "demand", bus="sink", p_set=3.0)
    n.links.static["country"] = "CTY"
    matched = pd.DataFrame(
        {"name": ["consume:f0:CTY", "consume:f1:CTY"], "target_mt": [2.0, 1.0]}
    )
    return n, matched


def test_hard_churn_budget_binds():
    # Budget = 2 * 0.25 * 3 = 1.5 Mt churn => shift 0.75 Mt from f0 to f1.
    n, matched = _make_solvable_network()
    n.optimize.create_model()
    cfg = _dp_cfg(penalty_mode="hard", max_relative_deviation=0.25)
    add_diet_stability_constraints(n, matched, cfg, 0.0)
    status, _ = n.optimize.solve_model(solver_name="highs")
    assert status == "ok"
    p0 = n.links.dynamic["p0"].iloc[-1]
    assert p0["consume:f0:CTY"] == pytest.approx(1.25)
    assert p0["consume:f1:CTY"] == pytest.approx(1.75)
    assert evaluate_diet_stability_cost(n, matched, cfg, 500.0) == 0.0


def test_hard_full_flexibility_is_unconstrained():
    # delta=1 allows a complete reorganization: everything moves to cheap f1.
    n, matched = _make_solvable_network()
    n.optimize.create_model()
    cfg = _dp_cfg(penalty_mode="hard", max_relative_deviation=1.0)
    add_diet_stability_constraints(n, matched, cfg, 0.0)
    status, _ = n.optimize.solve_model(solver_name="highs")
    assert status == "ok"
    p0 = n.links.dynamic["p0"].iloc[-1]
    assert p0["consume:f0:CTY"] == pytest.approx(0.0)
    assert p0["consume:f1:CTY"] == pytest.approx(3.0)


def test_hard_slack_priced_into_evaluation():
    # delta=0 freezes the diet; cheap slack (1 bn/Mt vs a 9 bn/Mt supply-cost
    # saving) absorbs the full 4 Mt churn of a complete shift to f1.
    n, matched = _make_solvable_network()
    n.optimize.create_model()
    cfg = _dp_cfg(penalty_mode="hard", max_relative_deviation=0.0, enable_slack=True)
    add_diet_stability_constraints(n, matched, cfg, 1.0)
    status, _ = n.optimize.solve_model(solver_name="highs")
    assert status == "ok"
    p0 = n.links.dynamic["p0"].iloc[-1]
    assert p0["consume:f1:CTY"] == pytest.approx(3.0)
    assert evaluate_diet_stability_cost(n, matched, cfg, 1.0) == pytest.approx(4.0)
