# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Extract the realized production cost of the food system.

Serialises the affordability guardrail's capped quantity -- the sum of
``marginal_cost * p`` over the production-cost links and generators (see
``workflow/scripts/solve_model/affordability.py`` for the carrier
whitelist and what it deliberately excludes) -- as one row per
(component, carrier) in bnUSD. The total over all rows is, by
construction, the left-hand side of ``affordability.cost_cap``;
reference-based cap specifications resolve against this output, so it
can never drift from the constraint's own cost definition.
"""

import logging

import pandas as pd
import pypsa

from workflow.scripts.solve_model.affordability import (
    production_cost_generators,
    production_cost_links,
)

logger = logging.getLogger(__name__)


def extract_production_cost(n: pypsa.Network) -> pd.DataFrame:
    """Realized production cost by (component, carrier), in bnUSD."""
    parts = []
    link_mc = production_cost_links(n)
    if not link_mc.empty:
        parts.append(
            pd.DataFrame(
                {
                    "component": "Link",
                    "carrier": n.links.static.loc[link_mc.index, "carrier"],
                    "cost_bnusd": link_mc
                    * n.links.dynamic.p0.loc["now", link_mc.index],
                }
            )
        )
    gen_mc = production_cost_generators(n)
    if not gen_mc.empty:
        parts.append(
            pd.DataFrame(
                {
                    "component": "Generator",
                    "carrier": n.generators.static.loc[gen_mc.index, "carrier"],
                    "cost_bnusd": gen_mc
                    * n.generators.dynamic.p.loc["now", gen_mc.index],
                }
            )
        )
    if not parts:
        logger.info("No production-cost components found; empty cost output")
        return pd.DataFrame(columns=["component", "carrier", "cost_bnusd"])
    out = (
        pd.concat(parts)
        .groupby(["component", "carrier"], as_index=False)["cost_bnusd"]
        .sum()
    )
    logger.info("Realized production cost: %.1f bnUSD", out["cost_bnusd"].sum())
    return out
