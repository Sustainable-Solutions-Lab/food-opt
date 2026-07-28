# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Validation checks for consumer values configuration."""

from workflow.scenario_generators import expand_scenario_defs
from workflow.scripts.solve_namespace import get_effective_config


def validate_consumer_values(config: dict, _project_root=None) -> None:
    """Validate piecewise food utility against its configured baseline."""
    if config["deviation_penalty"]["calibration"]["generate"]:
        # The Broyden iteration synthesizes a baseline and a main scenario per
        # iteration in-process, so the configured baseline names one of those
        # rather than a scenario the workflow declares.
        return

    scenario_defs = expand_scenario_defs(config["scenarios"])
    effective_configs = [("base config", config)]
    effective_configs.extend(
        (f"scenario '{name}'", get_effective_config(config, name, scenario_defs))
        for name in scenario_defs
    )

    for label, effective in effective_configs:
        if not effective["food_utility_piecewise"]["enabled"]:
            continue
        baseline = effective["consumer_values"]["baseline_scenario"]
        if baseline not in scenario_defs:
            raise ValueError(
                f"{label} enables food_utility_piecewise but configured consumer "
                f"values baseline scenario '{baseline}' is not defined"
            )
        baseline_config = get_effective_config(config, baseline, scenario_defs)
        if not baseline_config["validation"]["enforce_baseline_diet"]:
            raise ValueError(
                f"consumer values baseline scenario '{baseline}' must set "
                "validation.enforce_baseline_diet=true"
            )
