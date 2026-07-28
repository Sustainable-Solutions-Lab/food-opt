# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

from workflow.scripts.ihme import country_name_to_iso3


def test_country_name_to_iso3_prefers_ihme_override() -> None:
    assert country_name_to_iso3("Niger") == "NER"


def test_country_name_to_iso3_uses_pycountry() -> None:
    assert country_name_to_iso3("Canada") == "CAN"


def test_country_name_to_iso3_returns_none_for_unknown_name() -> None:
    assert country_name_to_iso3("Not a country") is None
