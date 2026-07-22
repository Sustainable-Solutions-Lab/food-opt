# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Tests for affordability cost-cap resolution."""

from types import SimpleNamespace

import pandas as pd
import pytest

from workflow.scripts.solve_model.affordability import resolve_cost_cap


def _write_cost(path, total: float) -> None:
    pd.DataFrame(
        {
            "component": ["Link", "Generator"],
            "carrier": ["crop_production", "fertilizer"],
            "cost_bnusd": [total - 40.0, 40.0],
        }
    ).to_parquet(path)


class TestResolveCostCap:
    def test_absolute_cap_passes_through(self):
        cap_cfg = {"enabled": True, "max_cost_bnusd": 3300.0}
        assert resolve_cost_cap(cap_cfg, SimpleNamespace()) == 3300.0

    def test_reference_cap_uses_realized_cost(self, tmp_path):
        reference_path = tmp_path / "reference.parquet"
        _write_cost(reference_path, 4200.0)
        cap_cfg = {
            "enabled": True,
            "max_cost_bnusd": {"reference_scenario": "frozen"},
        }
        inputs = SimpleNamespace(cost_reference=str(reference_path))
        assert resolve_cost_cap(cap_cfg, inputs) == pytest.approx(4200.0)

    def test_empty_reference_fails(self, tmp_path):
        reference_path = tmp_path / "reference.parquet"
        pd.DataFrame(columns=["component", "carrier", "cost_bnusd"]).to_parquet(
            reference_path
        )
        cap_cfg = {
            "enabled": True,
            "max_cost_bnusd": {"reference_scenario": "frozen"},
        }
        inputs = SimpleNamespace(cost_reference=str(reference_path))
        with pytest.raises(ValueError, match="empty"):
            resolve_cost_cap(cap_cfg, inputs)
