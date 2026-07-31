# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later


import pytest
import yaml

from workflow.validation.calibration import validate_calibration


@pytest.fixture
def default_config() -> dict:
    with open("config/default.yaml") as f:
        return yaml.safe_load(f)


def test_schema_requires_fixed_default_sections(default_config) -> None:
    with open("config/schemas/config.schema.yaml") as f:
        schema = yaml.safe_load(f)

    assert set(default_config) <= set(schema["required"])
    for section in ("numerics", "diet", "health", "solving"):
        section_schema = schema["properties"][section]
        assert set(default_config[section]) <= set(section_schema["required"])


class TestStabilityYamlExistenceGate:
    """The legacy stability yaml is required only when a solve would open it."""

    @staticmethod
    def _config(dp_enabled: bool, l1_cost) -> dict:
        disabled = {"enabled": False, "generate": False, "scenario": "calibration"}
        return {
            "scenarios": {"default": {}},
            "grazing": {"grassland_forage_calibration": dict(disabled)},
            "exogenous_feed_calibration": dict(disabled),
            "food_loss_waste_calibration": dict(disabled),
            "food_demand_calibration": dict(disabled),
            "cost_calibration": dict(disabled),
            "deviation_penalty": {
                "enabled": dp_enabled,
                "penalty_mode": "l1",
                "deviation_type": "absolute",
                "land": {
                    "crops": {"l1_cost": l1_cost},
                    "grassland": {"l1_cost": l1_cost},
                },
                "feed": {"l1_cost": l1_cost},
                "diet": {"l1_cost": 0.0},
                "calibration": {
                    "enabled": True,
                    "generate": False,
                    "calibrated_yaml": "missing/deviation_penalty.yaml",
                },
            },
        }

    def test_missing_yaml_accepted_when_penalty_is_off(self, tmp_path):
        validate_calibration(self._config(False, "calibrated"), tmp_path)

    def test_missing_yaml_accepted_without_calibrated_sentinel(self, tmp_path):
        validate_calibration(self._config(True, 0.1), tmp_path)

    def test_missing_yaml_rejected_when_sentinel_resolves(self, tmp_path):
        with pytest.raises(ValueError, match="calibrated_yaml"):
            validate_calibration(self._config(True, "calibrated"), tmp_path)
