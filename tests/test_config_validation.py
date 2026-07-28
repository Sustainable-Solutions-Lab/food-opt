# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

from copy import deepcopy

import pytest
import yaml

from workflow.validation.consumer_values import validate_consumer_values


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


def test_disabled_piecewise_utility_needs_no_baseline(default_config) -> None:
    validate_consumer_values(default_config)


def test_piecewise_utility_requires_configured_baseline(default_config) -> None:
    config = deepcopy(default_config)
    config["food_utility_piecewise"]["enabled"] = True

    with pytest.raises(ValueError, match="baseline scenario 'baseline' is not defined"):
        validate_consumer_values(config)


def test_piecewise_utility_requires_enforced_baseline(default_config) -> None:
    config = deepcopy(default_config)
    config["scenarios"]["baseline"] = {}
    config["food_utility_piecewise"]["enabled"] = True

    with pytest.raises(ValueError, match="enforce_baseline_diet=true"):
        validate_consumer_values(config)


def test_scenario_piecewise_utility_accepts_enforced_baseline(default_config) -> None:
    config = deepcopy(default_config)
    config["scenarios"] = {
        "baseline": {"validation": {"enforce_baseline_diet": True}},
        "utility": {"food_utility_piecewise": {"enabled": True}},
    }

    validate_consumer_values(config)
