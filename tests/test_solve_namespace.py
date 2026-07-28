# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Tests for solve_namespace helpers shared between Snakemake and the cluster
manifest exporter."""

import pytest
import yaml

from workflow.scripts.solve_namespace import (
    get_effective_config,
    validate_scenario_config_schemas,
    validate_scenario_overrides,
)


class TestValidateScenarioOverrides:
    def test_accepts_solve_time_only_overrides(self):
        defs = {
            "low_ghg": {"emissions": {"ghg_price": 50.0}},
            "no_health": {"health": {"enabled": False}},
            "stability": {"deviation_penalty": {"diet": {"enabled": True}}},
        }
        validate_scenario_overrides(defs)

    def test_rejects_structural_topology_override(self):
        """Scenario overriding 'countries' must fail: same built network is
        reused across scenarios, so changing topology silently mismatches it."""
        defs = {"bad": {"countries": ["USA", "FRA"]}}
        with pytest.raises(ValueError, match="structural key 'countries'"):
            validate_scenario_overrides(defs)

    def test_rejects_structural_residue_override(self):
        defs = {"bad": {"residues": {"max_feed_fraction": 0.5}}}
        with pytest.raises(ValueError, match="structural key"):
            validate_scenario_overrides(defs)

    def test_collects_multiple_errors(self):
        defs = {
            "bad1": {"countries": ["USA"]},
            "bad2": {"residues": {"max_feed_fraction": 0.5}},
        }
        with pytest.raises(ValueError) as info:
            validate_scenario_overrides(defs)
        msg = str(info.value)
        assert "bad1" in msg and "bad2" in msg


class TestValidateScenarioConfigSchemas:
    @pytest.fixture(scope="class")
    def base_config(self):
        with open("config/default.yaml") as f:
            return yaml.safe_load(f)

    def test_accepts_current_structure(self, base_config):
        defs = {
            "ok": {
                "deviation_penalty": {"land": {"crops": {"l1_cost_factor": 0.31623}}},
                "land": {"reforestation_cap": {"max_fraction": 0.5}},
            },
        }
        validate_scenario_config_schemas(base_config, defs, ".")

    def test_rejects_pre_split_deviation_penalty_structure(self, base_config):
        """The old flat land.l1_cost_factor layout must fail loudly instead of
        merging in silently and being ignored at solve time."""
        defs = {"stale": {"deviation_penalty": {"land": {"l1_cost_factor": 0.3}}}}
        with pytest.raises(ValueError, match="stale"):
            validate_scenario_config_schemas(base_config, defs, ".")

    def test_rejects_unknown_sensitivity_key(self, base_config):
        defs = {"stale": {"sensitivity": {"max_reforestation_fraction": 0.5}}}
        with pytest.raises(ValueError, match="stale"):
            validate_scenario_config_schemas(base_config, defs, ".")

    @pytest.mark.parametrize(
        "override",
        [
            {"solving": {"inline_analysis": True}},
            {"remote_solve": {"enabled": True}},
            {"plotting": {"comparison_scenarios": "all"}},
            {"food_groups": {"max_per_capita": {"grain": 100.0}}},
            {"biomass": {"marginal_values_usd_per_tonne": 100.0}},
        ],
    )
    def test_rejects_parse_or_build_time_override(self, override):
        with pytest.raises(ValueError, match="structural key"):
            validate_scenario_overrides({"ignored": override})

    def test_accepts_solver_resource_override(self):
        validate_scenario_overrides(
            {
                "parallel": {
                    "solving": {
                        "threads": 4,
                        "runtime": 20,
                        "mem_mb": 8000,
                    }
                }
            }
        )

    def test_validates_one_representative_per_structure(self, base_config, monkeypatch):
        """Thousands of same-template samples must cost one validation."""
        import workflow.validation.config_schema as cs

        calls = []
        monkeypatch.setattr(
            cs, "validate_config_schema", lambda cfg, root: calls.append(1)
        )
        defs = {f"gsa_{i}": {"emissions": {"ghg_price": float(i)}} for i in range(50)}
        validate_scenario_config_schemas(base_config, defs, ".")
        assert len(calls) == 1


def test_get_effective_config_applies_nested_overrides_without_mutating_base() -> None:
    base = {"solving": {"solver": "highs", "threads": 1}}
    scenarios = {"parallel": {"solving": {"threads": 4}}}

    effective = get_effective_config(base, "parallel", scenarios)

    assert effective == {"solving": {"solver": "highs", "threads": 4}}
    assert base["solving"]["threads"] == 1


def test_get_effective_config_rejects_unknown_scenario() -> None:
    with pytest.raises(KeyError, match="missing"):
        get_effective_config({}, "missing", {})
