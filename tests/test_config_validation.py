# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later


import pytest
import yaml


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
