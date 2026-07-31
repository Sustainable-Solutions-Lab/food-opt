# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Tests for merging baseline dietary-intake sources."""

import pandas as pd

from workflow.scripts.merge_dietary_sources import _faostat_animal_fat_fallback


def test_animal_fat_fallback_preserves_existing_and_fills_missing():
    existing = pd.DataFrame(
        {"country": ["AAA"], "item": ["animal_fat"], "value": [0.0]}
    )
    supply = pd.DataFrame(
        {
            "country": ["AAA", "BBB"],
            "item": ["animal_fat", "animal_fat"],
            "value": [10.0, 20.0],
        }
    )

    result = _faostat_animal_fat_fallback(supply, existing)

    assert result[["country", "value"]].to_dict("records") == [
        {"country": "BBB", "value": 17.0}
    ]
