# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Tests for the health Stage 2 piecewise formulation."""

import linopy
import numpy as np
import pandas as pd
import pytest

from workflow.scripts.solve_model.health import _add_stage2_piecewise


def _build_stage2_model(log_total: float) -> tuple[linopy.Model, object]:
    m = linopy.Model()
    store_names = pd.Index(["store-0", "store-1", "store-2"], name="name")
    store = m.add_variables(lower=-np.inf, coords=[store_names], name="store")
    log_0 = m.add_variables(lower=-10, upper=10, name="log-0")
    log_1 = m.add_variables(lower=-10, upper=10, name="log-1")
    log_2 = m.add_variables(lower=-10, upper=10, name="log-2")
    m.add_constraints(log_0 == log_total, name="fix-log-0")
    m.add_constraints(log_1 == log_total, name="fix-log-1")
    m.add_constraints(log_2 == log_total, name="fix-log-2")

    pairs = [(0, "cause"), (1, "cause"), (2, "cause")]
    health_stores = pd.DataFrame(
        {"name": store_names},
        index=pd.MultiIndex.from_tuples(pairs, names=["health_cluster", "cause"]),
    )
    cluster_cause_data = {
        (0, "cause"): {
            "yll_total": 1e6,
            "rr_baseline": 1.0,
            "rr_ref": 1.0,
        },
        (1, "cause"): {
            "yll_total": 2e6,
            "rr_baseline": 1.0,
            "rr_ref": 1.0,
        },
        (2, "cause"): {
            "yll_total": 0.0,
            "rr_baseline": 1.0,
            "rr_ref": 1.0,
        },
    }
    count = _add_stage2_piecewise(
        m=m,
        log_rr_totals={
            (0, "cause"): log_0,
            (1, "cause"): log_1,
            (2, "cause"): log_2,
        },
        cluster_cause_pairs=pairs,
        cluster_cause_data=cluster_cause_data,
        health_stores=health_stores,
        store_level_var=store,
        log_pts=np.array([0.0, 1.0, 2.0]),
        rr_pts=np.array([1.0, 3.0, 8.0]),
    )
    assert count == 3
    m.objective = store.sum()
    return m, store


@pytest.mark.parametrize(
    "log_total, expected",
    [
        (0.5, [1.0, 2.0, 0.0]),
        (1.0, [2.0, 4.0, 0.0]),
    ],
)
def test_stage2_matches_chord_at_interior_and_knot(log_total, expected):
    m, store = _build_stage2_model(log_total)

    status, condition = m.solve(solver_name="highs")

    assert (status, condition) == ("ok", "optimal")
    assert store.solution.values.tolist() == pytest.approx(expected)
    assert set(m.variables) == {"store", "log-0", "log-1", "log-2"}


def test_stage2_enforces_breakpoint_domain():
    m, _ = _build_stage2_model(-0.5)

    _, condition = m.solve(solver_name="highs")

    assert condition == "infeasible"
