# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Tests for GDD-IA source selection."""

import pandas as pd
import pytest

from workflow.scripts.diet.gdd_ia import closest_gdd_ia_release_year
from workflow.scripts.prepare_gdd_ia_dietary_intake import (
    UNIT_BY_GROUP,
    _materialize_sparse_zeros,
)


@pytest.mark.parametrize(
    ("reference_year", "expected"),
    [
        (2015, 2015),
        (2017, 2015),
        (2018, 2020),
    ],
)
def test_closest_gdd_ia_release_year(reference_year, expected):
    assert closest_gdd_ia_release_year(reference_year) == expected


def test_sparse_country_group_cells_become_explicit_zeros():
    intake = pd.DataFrame(
        {
            "country": ["AAA", "BBB", "CCC"],
            "group": ["grain", "grain", "legumes"],
            "kcal": [100.0, 80.0, 12.0],
            "value": [30.0, 25.0, 5.0],
        }
    )

    result = _materialize_sparse_zeros(
        intake,
        required_countries={"AAA", "BBB"},
        included_groups={"grain", "legumes"},
    ).set_index(["country", "group"])

    assert result.loc[("AAA", "grain"), "value"] == 30.0
    assert result.loc[("AAA", "legumes"), ["kcal", "value"]].tolist() == [0.0, 0.0]
    assert result.loc[("BBB", "legumes"), ["kcal", "value"]].tolist() == [0.0, 0.0]
    assert result.loc[("CCC", "legumes"), "value"] == 5.0


def test_all_gdd_groups_have_output_units():
    assert UNIT_BY_GROUP["animal_fat"] == "g/day (fresh wt)"


def test_sparse_zero_materialization_can_leave_documented_fallback_group():
    intake = pd.DataFrame(
        {"country": ["AAA"], "group": ["grain"], "kcal": [100.0], "value": [30.0]}
    )

    result = _materialize_sparse_zeros(
        intake,
        required_countries={"AAA"},
        included_groups={"grain", "animal_fat"},
        excluded_groups={"animal_fat"},
    )

    assert set(result["group"]) == {"grain"}
