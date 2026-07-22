# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Tests for solve_namespace helpers shared between Snakemake and the cluster
manifest exporter."""

import pytest
import yaml

from workflow.scripts.solve_namespace import (
    affordability_reference_inputs,
    barrier_reference_inputs,
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


class TestAffordabilityReferenceInputs:
    @staticmethod
    def _base_config(cost_cap):
        return {"affordability": {"cost_cap": cost_cap}}

    def test_absolute_cap_needs_no_references(self):
        cap = {"enabled": True, "max_cost_bnusd": 3300.0}
        assert (
            affordability_reference_inputs(cap, "s", self._base_config(cap), {}) == {}
        )

    def test_disabled_cap_needs_no_references(self):
        spec = {"reference_scenario": "a"}
        cap = {"enabled": False, "max_cost_bnusd": spec}
        assert (
            affordability_reference_inputs(cap, "s", self._base_config(cap), {}) == {}
        )

    def test_reference_cap_returns_reference_path(self):
        spec = {"reference_scenario": "a"}
        cap = {"enabled": True, "max_cost_bnusd": spec}
        base = self._base_config({"enabled": False, "max_cost_bnusd": 1.0})
        defs = {"a": {}}
        inputs = affordability_reference_inputs(cap, "s", base, defs)
        assert inputs == {
            "cost_reference": "<results>/{name}/analysis/scen-a/production_cost.parquet"
        }

    def test_undefined_reference_scenario_fails(self):
        spec = {"reference_scenario": "a"}
        cap = {"enabled": True, "max_cost_bnusd": spec}
        base = self._base_config({"enabled": False, "max_cost_bnusd": 1.0})
        with pytest.raises(ValueError, match="not a defined scenario"):
            affordability_reference_inputs(cap, "s", base, {})

    def test_reference_based_reference_scenario_fails(self):
        """Reference chains must not recurse: a reference scenario whose own
        cap is reference-based is rejected at input-building time."""
        spec = {"reference_scenario": "a"}
        cap = {"enabled": True, "max_cost_bnusd": spec}
        base = self._base_config({"enabled": False, "max_cost_bnusd": 1.0})
        defs = {
            "a": {
                "affordability": {
                    "cost_cap": {"enabled": True, "max_cost_bnusd": dict(spec)}
                }
            },
        }
        with pytest.raises(ValueError, match="itself uses a reference-based"):
            affordability_reference_inputs(cap, "s", base, defs)


class TestBarrierReferenceInputs:
    @pytest.fixture
    def base_config(self):
        with open("config/default.yaml") as f:
            return yaml.safe_load(f)

    def test_shared_reference_returns_one_artifact(self, base_config):
        defs = {"frozen": {}}
        scenario = yaml.safe_load(yaml.safe_dump(base_config))
        scenario["food_energy"]["floor"].update(
            {"enabled": True, "reference_scenario": "frozen"}
        )
        scenario["biodiversity"]["cap"].update(
            {
                "enabled": True,
                "max_conversion_mha": {"reference_scenario": "frozen"},
            }
        )
        assert barrier_reference_inputs(scenario, "sample", base_config, defs) == {
            "barrier_reference": (
                "<results>/{name}/analysis/scen-frozen/barrier_reference.parquet"
            )
        }

    def test_mixed_references_fail(self, base_config):
        defs = {"frozen": {}, "other": {}}
        scenario = yaml.safe_load(yaml.safe_dump(base_config))
        scenario["food_energy"]["floor"].update(
            {"enabled": True, "reference_scenario": "frozen"}
        )
        scenario["protein"]["floor"].update(
            {"enabled": True, "reference_scenario": "other"}
        )
        with pytest.raises(ValueError, match="must share one reference"):
            barrier_reference_inputs(scenario, "sample", base_config, defs)
