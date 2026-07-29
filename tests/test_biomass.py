# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Tests for biomass-routing helpers in `workflow.scripts.build_model.biomass`.

Pins the skip-on-missing-supply behaviour of `add_biofuel_links` so
regional configs without all global producers (e.g. Europe-only without
sugarcane) do not silently introduce infeasible fixed-demand links.
Also exercises the food-bus DM deflation introduced for biomass routing.
"""

import pandas as pd
import pypsa

from workflow.scripts.build_model import biomass


def _build_network() -> pypsa.Network:
    n = pypsa.Network()
    n.set_snapshots(["now"])
    n.add(
        "Bus",
        ["food:wheat:USA", "crop:wheat:USA", "biomass:USA", "land:USA"],
        carrier=["food", "crop", "biomass", "land"],
    )
    # Only wheat has a crop_production link; sugarcane and maize do not.
    n.add(
        "Link",
        "produce:wheat",
        bus0="land:USA",
        bus1="crop:wheat:USA",
        carrier="crop_production",
        crop="wheat",
    )
    return n


def test_biofuel_skips_crops_without_production(caplog):
    n = _build_network()
    df = pd.DataFrame(
        {
            "source_item": ["wheat", "sugarcane"],
            "crop": ["wheat", "sugarcane"],
            "country": ["USA", "USA"],
            "bus_type": ["food", "food"],
            "demand_mt": [10.0, 5.0],
            "dm_mt": [8.7, 1.5],
        }
    )
    # bus0 for sugarcane does not exist either, but the no-supply check
    # fires first so we test the *no-supply* skip path independently
    # of the bus-not-found path.
    n.add("Bus", "food:sugarcane:USA", carrier="food")

    with caplog.at_level("WARNING"):
        biomass.add_biofuel_links(n, df)

    names = n.links.static.index.tolist()
    assert any("biofuel:wheat:" in s for s in names)
    assert not any("biofuel:sugarcane:" in s for s in names)
    assert any("no production anywhere" in rec.message for rec in caplog.records)


def test_biofuel_draws_demand_and_delivers_dry_matter():
    """The link draws demand_mt from the source bus and delivers dm_mt."""
    n = _build_network()
    df = pd.DataFrame(
        {
            "source_item": ["wheat"],
            "crop": ["wheat"],
            "country": ["USA"],
            "bus_type": ["food"],
            "demand_mt": [10.0],
            "dm_mt": [8.7],
        }
    )
    biomass.add_biofuel_links(n, df)

    link = n.links.static.loc["biofuel:wheat:USA"]
    assert abs(link["p_nom"] - 10.0) < 1e-9
    assert abs(link["p_nom"] * link["efficiency"] - 8.7) < 1e-9


def test_biofuel_commercial_basis_item_is_not_over_drawn():
    """An item already in commercial units draws exactly its reported demand.

    This is the palm-oil case: the demand is reported in oil, which is what the
    food bus carries and is dry matter throughout, so the link must draw the
    reported quantity rather than inflating it by the source crop's moisture
    (oil palm's fresh fruit bunches are 60% water; the oil is not).
    """
    n = _build_network()
    df = pd.DataFrame(
        {
            "source_item": ["wheat"],
            "crop": ["wheat"],
            "country": ["USA"],
            "bus_type": ["food"],
            "demand_mt": [38.13],
            "dm_mt": [38.13],
        }
    )
    biomass.add_biofuel_links(n, df)

    link = n.links.static.loc["biofuel:wheat:USA"]
    assert abs(link["p_nom"] - 38.13) < 1e-9
    assert link["efficiency"] == 1.0
