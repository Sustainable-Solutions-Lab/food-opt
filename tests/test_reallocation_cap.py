# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Unit tests for endpoint-normalized production reallocation caps."""

import pytest

from workflow.scripts.solve_model.reallocation_cap import _component_weights
from workflow.scripts.solve_namespace import reallocation_reference_inputs


def _deviation_config() -> dict:
    return {
        "penalty_mode": "l1",
        "deviation_type": "absolute",
        "land": {
            "crops": {"l1_cost": 2.0},
            "grassland": {"l1_cost": 0.5},
        },
        "feed": {"l1_cost": 0.1},
    }


def test_component_weights_use_resolved_l1_coefficients():
    assert _component_weights(_deviation_config()) == {
        "cropland": 2.0,
        "grassland": 0.5,
        "feed": 0.1,
    }


def test_component_weights_reject_incomparable_geometry():
    config = _deviation_config()
    config["deviation_type"] = "relative"
    with pytest.raises(ValueError, match="deviation_type='absolute'"):
        _component_weights(config)


def test_reallocation_reference_path_uses_both_endpoints():
    disabled = {
        "enabled": False,
        "fraction": 1.0,
        "baseline_scenario": None,
        "endpoint_scenario": None,
    }
    enabled = {
        "enabled": True,
        "fraction": 0.2,
        "baseline_scenario": "base",
        "endpoint_scenario": "end",
    }
    scenario_defs = {"base": {}, "end": {}, "cap": {"reallocation_cap": enabled}}
    inputs = reallocation_reference_inputs(
        enabled,
        "cap",
        {"reallocation_cap": disabled},
        scenario_defs,
    )
    assert inputs == {
        "reallocation_reference": (
            "<results>/{name}/reallocation_reference/baseline-base/endpoint-end.parquet"
        )
    }


def test_reallocation_reference_rejects_undefined_endpoint():
    enabled = {
        "enabled": True,
        "fraction": 0.2,
        "baseline_scenario": "base",
        "endpoint_scenario": "missing",
    }
    with pytest.raises(ValueError, match="undefined reallocation-cap reference"):
        reallocation_reference_inputs(
            enabled,
            "cap",
            {"reallocation_cap": {"enabled": False}},
            {"base": {}, "cap": {}},
        )
