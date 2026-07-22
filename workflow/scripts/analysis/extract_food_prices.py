# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Extract food prices from bus marginal prices in solved networks.

Food bus shadow prices (marginal prices on the nodal balance constraint)
represent the marginal system cost of delivering one additional unit of
food at a given location. They naturally incorporate all upstream costs
(production, trade, processing) plus any externality pricing (GHG,
health). Units are bnUSD/Mt = USD/kg.

Outputs:
- food_prices.parquet: Per-food, per-country prices (USD/kg) with
  consumption-weighted averages and per-capita daily diet cost.
"""

import logging

import numpy as np
import pandas as pd
import pypsa

from workflow.scripts.constants import DAYS_PER_YEAR
from workflow.scripts.population import get_country_population

logger = logging.getLogger(__name__)


def extract_food_prices(n: pypsa.Network) -> pd.DataFrame:
    """Extract food prices from bus marginal prices.

    Parameters
    ----------
    n : pypsa.Network
        Solved network with marginal prices computed.

    A food bus whose marginal price is set by an active positive food-slack
    generator (unmet fixed-diet demand priced at ``slack_marginal_cost``) does
    not carry a meaningful market price -- the whole trade-linked flow of that
    food equalises to the slack penalty (e.g. apple at ~500 USD/kg when its
    baseline demand exceeds achievable production). Such rows are marked
    ``is_slack_pinned`` and their price/cost columns are set to NaN so they are
    excluded from diet-cost aggregates rather than corrupting them.

    Returns
    -------
    pd.DataFrame
        Columns: food, food_group, country, price_usd_per_kg,
        consumption_mt, cost_bnusd, cost_usd_per_person_per_day,
        is_slack_pinned.
    """
    columns = [
        "food",
        "food_group",
        "country",
        "price_usd_per_kg",
        "consumption_mt",
        "cost_bnusd",
        "cost_usd_per_person_per_day",
        "is_slack_pinned",
    ]

    links = n.links.static
    consume_links = links[links["carrier"] == "food_consumption"]

    if consume_links.empty:
        return pd.DataFrame(columns=columns)

    # Get marginal prices on food buses (bnUSD/Mt = USD/kg)
    if "marginal_price" not in n.buses.dynamic:
        logger.warning("No marginal prices found in network; returning empty prices")
        return pd.DataFrame(columns=columns)

    marginal_price = n.buses.dynamic.marginal_price.iloc[0]

    # Get food consumption flows (p0 on consume links)
    p0 = n.links.dynamic.p0
    snapshot = n.snapshots[-1]
    consumption = p0.loc[snapshot].reindex(consume_links.index).fillna(0.0)

    population = pd.Series(get_country_population(n), dtype=float)

    df = consume_links[["food", "food_group", "country", "bus0"]].copy()
    df = df.astype({"food": str, "food_group": str, "country": str})
    df["consumption_mt"] = consumption.to_numpy()
    df = df[df["consumption_mt"] >= 1e-12]
    if df.empty:
        return pd.DataFrame(columns=columns)

    df["price_usd_per_kg"] = df["bus0"].map(marginal_price).fillna(0.0)
    df["cost_bnusd"] = df["price_usd_per_kg"] * df["consumption_mt"]
    pop = df["country"].map(population)
    df["cost_usd_per_person_per_day"] = np.where(
        pop > 0,
        (df["cost_bnusd"] * 1e9) / (pop.to_numpy() * DAYS_PER_YEAR),
        0.0,
    )

    # Flag foods whose price is set by an active positive slack generator
    # (unmet fixed-diet demand priced at slack_marginal_cost). The slack
    # generator dispatches on only a few source buses, but trade equalises its
    # penalty price across every bus of that food, so flag a (food, country)
    # row when its food slacks *somewhere* and its price sits at the slack
    # ceiling. Null those rows' price/cost so aggregates skip the artifact.
    gens = n.generators.static
    slack_gens = gens[gens["carrier"] == "slack_positive_food"]
    df["is_slack_pinned"] = False
    if not slack_gens.empty and "p" in n.generators.dynamic:
        gp = n.generators.dynamic.p.loc[snapshot].reindex(slack_gens.index).fillna(0.0)
        active_buses = slack_gens.loc[gp.to_numpy() > 1e-6, "bus"]
        bus_to_food = n.buses.static["food"]
        active_foods = set(bus_to_food.reindex(active_buses).dropna())
        slack_cost = float(slack_gens["marginal_cost"].max())
        df["is_slack_pinned"] = df["food"].isin(active_foods) & (
            df["price_usd_per_kg"] >= 0.5 * slack_cost
        )
    pinned = df["is_slack_pinned"].to_numpy()
    for col in ("price_usd_per_kg", "cost_bnusd", "cost_usd_per_person_per_day"):
        df.loc[pinned, col] = np.nan
    if pinned.any():
        logger.info(
            "Food prices: %d (food, country) rows slack-pinned (nulled): %s",
            int(pinned.sum()),
            ", ".join(sorted(df.loc[pinned, "food"].unique())[:8]),
        )
    df = df.drop(columns="bus0")

    # If duplicates exist across (food, food_group, country), the price
    # should be identical (it comes from the dual on the same food bus)
    # while flow-quantity columns sum. agg("first") on price would
    # silently drop divergent prices; assert uniqueness instead.
    group_cols = ["food", "food_group", "country"]
    if df.duplicated(subset=group_cols).any():
        prices_per_group = df.groupby(group_cols)["price_usd_per_kg"].nunique()
        if (prices_per_group > 1).any():
            offenders = prices_per_group[prices_per_group > 1].head().to_dict()
            raise ValueError(
                "Duplicate (food, food_group, country) rows with divergent "
                f"prices in food-price extraction: {offenders}"
            )
    df = df.groupby(group_cols, as_index=False).agg(
        {
            "price_usd_per_kg": "first",
            "consumption_mt": "sum",
            "cost_bnusd": "sum",
            "cost_usd_per_person_per_day": "sum",
            "is_slack_pinned": "max",
        }
    )

    return df.sort_values(["country", "food"]).reset_index(drop=True)
