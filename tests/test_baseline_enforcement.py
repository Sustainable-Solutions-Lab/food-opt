# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Unit tests for food-level baseline enforcement."""

import pandas as pd
import pypsa
import pytest

from workflow.scripts.solve_model.core import (
    _build_ratios_from_baseline,
    _prepare_baseline_diet_for_food_constraints,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def baseline_df():
    """Minimal baseline diet DataFrame."""
    return pd.DataFrame(
        {
            "food": [
                "wheat",
                "rice",
                "maize",
                "beef",
                "poultry",
                "lentils",
            ],
            "country": ["USA", "USA", "USA", "USA", "USA", "USA"],
            "food_group": [
                "grain",
                "grain",
                "grain",
                "red_meat",
                "poultry_meat",
                "legumes",
            ],
            "consumption_g_per_day": [150.0, 50.0, 100.0, 75.0, 60.0, 20.0],
        }
    )


@pytest.fixture
def baseline_df_multi_country():
    """Baseline diet with two countries."""
    return pd.DataFrame(
        {
            "food": ["wheat", "rice", "wheat", "rice"],
            "country": ["USA", "USA", "IND", "IND"],
            "food_group": ["grain", "grain", "grain", "grain"],
            "consumption_g_per_day": [150.0, 50.0, 30.0, 170.0],
        }
    )


@pytest.fixture
def food_network():
    """Minimal PyPSA network with food consumption links."""
    n = pypsa.Network()

    n.carriers.add("food_consumption", unit="Mt")

    # Buses for foods and food groups
    n.buses.add(
        [
            "food:wheat:USA",
            "food:rice:USA",
            "food:beef:USA",
            "food:wheat:IND",
            "food:rice:IND",
            "group:grain:USA",
            "group:red_meat:USA",
            "group:grain:IND",
        ],
        carrier=[
            "food_wheat",
            "food_rice",
            "food_beef",
            "food_wheat",
            "food_rice",
            "group_grain",
            "group_red_meat",
            "group_grain",
        ],
    )

    # Food consumption links
    n.links.add(
        [
            "consume:wheat:USA",
            "consume:rice:USA",
            "consume:beef:USA",
            "consume:wheat:IND",
            "consume:rice:IND",
        ],
        bus0=[
            "food:wheat:USA",
            "food:rice:USA",
            "food:beef:USA",
            "food:wheat:IND",
            "food:rice:IND",
        ],
        bus1=[
            "group:grain:USA",
            "group:grain:USA",
            "group:red_meat:USA",
            "group:grain:IND",
            "group:grain:IND",
        ],
        carrier="food_consumption",
        marginal_cost=[0.01, 0.01, 0.01, 0.01, 0.01],
        food=["wheat", "rice", "beef", "wheat", "rice"],
        country=["USA", "USA", "USA", "IND", "IND"],
        food_group=["grain", "grain", "red_meat", "grain", "grain"],
        flw_multiplier=[1.0, 1.0, 1.0, 1.0, 1.0],
    )

    return n


# ---------------------------------------------------------------------------
# Test _build_ratios_from_baseline
# ---------------------------------------------------------------------------


class TestBuildRatiosFromBaseline:
    def test_ratios_sum_to_one_per_group(self, baseline_df):
        result = _build_ratios_from_baseline(baseline_df)

        grain = result[result["food_group"] == "grain"]
        assert grain["ratio"].sum() == pytest.approx(1.0)

    def test_single_food_group_gets_ratio_one(self, baseline_df):
        result = _build_ratios_from_baseline(baseline_df)

        beef = result[result["food"] == "beef"]
        assert beef["ratio"].values[0] == pytest.approx(1.0)

        lentils = result[result["food"] == "lentils"]
        assert lentils["ratio"].values[0] == pytest.approx(1.0)

    def test_correct_within_group_proportions(self, baseline_df):
        result = _build_ratios_from_baseline(baseline_df)

        grain = result[result["food_group"] == "grain"].set_index("food")
        # wheat=150, rice=50, maize=100 → total=300
        assert grain.at["wheat", "ratio"] == pytest.approx(150 / 300)
        assert grain.at["rice", "ratio"] == pytest.approx(50 / 300)
        assert grain.at["maize", "ratio"] == pytest.approx(100 / 300)

    def test_country_codes_uppercased(self, baseline_df_multi_country):
        result = _build_ratios_from_baseline(baseline_df_multi_country)

        assert (result["country"] == result["country"].str.upper()).all()

    def test_per_country_ratios_independent(self, baseline_df_multi_country):
        result = _build_ratios_from_baseline(baseline_df_multi_country)

        usa = result[result["country"] == "USA"].set_index("food")
        ind = result[result["country"] == "IND"].set_index("food")

        # USA: wheat=150, rice=50
        assert usa.at["wheat", "ratio"] == pytest.approx(0.75)
        assert usa.at["rice", "ratio"] == pytest.approx(0.25)

        # IND: wheat=30, rice=170
        assert ind.at["wheat", "ratio"] == pytest.approx(30 / 200)
        assert ind.at["rice", "ratio"] == pytest.approx(170 / 200)

    def test_zero_consumption_gives_zero_ratio(self):
        df = pd.DataFrame(
            {
                "food": ["wheat", "rice"],
                "country": ["USA", "USA"],
                "food_group": ["grain", "grain"],
                "consumption_g_per_day": [0.0, 0.0],
            }
        )
        result = _build_ratios_from_baseline(df)

        assert (result["ratio"] == 0.0).all()

    def test_output_columns(self, baseline_df):
        result = _build_ratios_from_baseline(baseline_df)

        assert list(result.columns) == ["country", "food_group", "food", "ratio"]


# ---------------------------------------------------------------------------
# Test _prepare_baseline_diet_for_food_constraints
# ---------------------------------------------------------------------------


class TestPrepareBaselineDietForFoodConstraints:
    def test_true_zero_preserved(self):
        """A genuine zero baseline must not be lifted to the minimum floor."""
        baseline = pd.DataFrame(
            {
                "food": ["wheat", "rice"],
                "country": ["USA", "USA"],
                "consumption_g_per_day_intake": [120.0, 0.0],
            }
        )
        consume_links = pd.DataFrame(
            {"food": ["wheat", "rice"], "country": ["USA", "USA"]}
        )

        result = _prepare_baseline_diet_for_food_constraints(baseline, consume_links)

        rice = result.loc[result["food"] == "rice", "consumption_g_per_day"].iloc[0]
        wheat = result.loc[result["food"] == "wheat", "consumption_g_per_day"].iloc[0]
        assert rice == pytest.approx(0.0)
        assert wheat == pytest.approx(120.0)

    def test_tiny_positive_lifted_to_floor(self):
        """Tiny-but-positive values are still lifted to the configured floor."""
        baseline = pd.DataFrame(
            {
                "food": ["wheat"],
                "country": ["USA"],
                "consumption_g_per_day_intake": [0.03],
            }
        )
        consume_links = pd.DataFrame({"food": ["wheat"], "country": ["USA"]})

        result = _prepare_baseline_diet_for_food_constraints(
            baseline, consume_links, min_consumption_g_per_day=0.1
        )
        assert result["consumption_g_per_day"].iloc[0] == pytest.approx(0.1)
