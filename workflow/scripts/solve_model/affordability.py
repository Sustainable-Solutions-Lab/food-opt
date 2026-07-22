# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Affordability guardrail: a cap on total food-system production cost.

Water-scarcity relief reorganises production toward water-rich but not
necessarily cost-optimal regions, which can raise what it costs to produce
the (fixed) diet. This module caps the total marginal production cost from
above at solve time, so relief cannot be bought with an unbounded increase
in the food system's cost.

The capped quantity is the sum of ``marginal_cost * p`` over the
production-cost carriers (in bnUSD): crop, grassland and animal
production, food processing (incl. marketing), feed conversion, trade,
non-renewable groundwater pumping, synthetic fertilizer supply, and the
exogenous feed backstops (priced at grazing cost, so they carry a real
feed-production cost). This is the recurring resource cost of producing
the diet, and deliberately nothing else:

- Externality prices (water scarcity, groundwater depletion, GHG,
  health) are marginal costs on *stores* and never enter.
- Slack (baseline-diet, feed, land, water) sits on generators outside
  the whitelist; churn (deviation) penalties and production-stability
  slacks are linopy-level objective terms, likewise excluded.
- The one-off land-conversion investment cost on ``land_conversion`` /
  ``new_to_pasture`` links is a cost of *change*, not of production --
  it belongs to the reallocation-feasibility barrier and is excluded
  here by the carrier whitelist. So are the degeneracy-breaking
  regularizers (water merit-order epsilon on renewable supply rows,
  the consume-link tie-break cost) and biomass export revenue.

The dual is the marginal value of relaxing the cost cap -- how much
objective (relief) one more bnUSD of allowed cost would buy.

This is a global cap, so it bounds aggregate cost, not its distribution; a
per-country diet-cost cap would be the distributionally faithful "access"
variant.
"""

import logging

import pandas as pd
import pypsa

from workflow.scripts.solve_model.guardrails import (
    add_scalar_dispatch_cap,
    assign_scalar_cap_dual,
)

logger = logging.getLogger(__name__)

_LABEL = "production_cost_cap"

# Link carriers whose marginal cost is a recurring production cost. Land
# conversion (a one-off cost of change), regularizers (water merit-order
# epsilon, consume-link tie-break) and any future abstract link cost stay
# out of the capped sum.
_PRODUCTION_COST_LINK_CARRIERS = [
    "animal_production",
    "crop_production",
    "crop_production_multi",
    "feed_conversion",
    "food_processing",
    "grassland_production",
    "trade_crop",
    "trade_feed",
    "trade_food",
]

# Generator carriers whose marginal cost is a recurring production cost:
# synthetic fertilizer and the exogenous feed backstops (build-time
# exogenous_feed plus the solve-time calibration forage/roughage/protein
# supplies, all priced at grazing cost). Slack generators and the biomass
# export sinks (revenue) stay out.
_PRODUCTION_COST_GEN_CARRIERS = [
    "exogenous_feed",
    "exogenous_forage_cal",
    "exogenous_protein_cal",
    "exogenous_roughage_cal",
    "fertilizer",
]


def _cost_bearing(static: pd.DataFrame, keep: pd.Series, kind: str) -> pd.Series:
    """Nonzero marginal costs among ``keep`` rows, logging dropped carriers."""
    mc = pd.to_numeric(static["marginal_cost"], errors="coerce").fillna(0.0)
    dropped = sorted(static.loc[(mc != 0.0) & ~keep, "carrier"].unique())
    if dropped:
        logger.info(
            "Cost-bearing %s carriers outside the production-cost cap: %s",
            kind,
            ", ".join(dropped),
        )
    return mc[keep & (mc != 0.0)]


def production_cost_links(n: pypsa.Network) -> pd.Series:
    """Per-link marginal cost (bnUSD per unit dispatch) for production-cost links.

    Selects the production-cost carriers plus non-renewable groundwater
    pumping (the renewable ``water_supply`` rows carry only the merit-order
    regularizer). Logs any other cost-bearing link carriers it drops.
    """
    links = n.links.static
    keep = links["carrier"].isin(_PRODUCTION_COST_LINK_CARRIERS) | (
        (links["carrier"] == "water_supply")
        & (links["source"] == "groundwater_nonrenewable")
    )
    return _cost_bearing(links, keep, "link")


def production_cost_generators(n: pypsa.Network) -> pd.Series:
    """Per-generator marginal cost (bnUSD/Mt) for production-cost generators.

    Selects fertilizer supply and the exogenous feed backstops. Logs any
    other cost-bearing generator carriers it drops (slack, biomass sinks).
    """
    gens = n.generators.static
    keep = gens["carrier"].isin(_PRODUCTION_COST_GEN_CARRIERS)
    return _cost_bearing(gens, keep, "generator")


def _reference_cost(path: str) -> float:
    """Total realized production cost (bnUSD) from a production_cost.parquet."""
    df = pd.read_parquet(path)
    if df.empty:
        raise ValueError(f"Reference production cost output is empty: {path}")
    return float(df["cost_bnusd"].sum())


def resolve_cost_cap(cap_cfg: dict, inputs) -> float:
    """Resolve the configured cap to bnUSD.

    ``max_cost_bnusd`` is either a number (absolute cap) or a mapping
    ``{reference_scenario}``: the cap is the reference scenario's realized
    production cost, read from the ``cost_reference`` input that the workflow
    wires up (see ``affordability_reference_inputs`` in
    workflow/scripts/solve_namespace.py).
    """
    spec = cap_cfg["max_cost_bnusd"]
    if not isinstance(spec, dict):
        return float(spec)
    cap = _reference_cost(inputs.cost_reference)
    logger.info(
        "Resolved production cost cap from %s: %.1f bnUSD",
        spec["reference_scenario"],
        cap,
    )
    return cap


def add_production_cost_cap(n: pypsa.Network, cap_cfg: dict, inputs) -> None:
    """Cap total marginal production cost at solve time.

    ``cap_cfg`` is the ``affordability.cost_cap`` block with keys ``enabled``
    (bool) and ``max_cost_bnusd`` (absolute bnUSD or the reference-based
    mapping resolved by ``resolve_cost_cap``, which reads the reference cost
    file from ``inputs``). Bounds the sum of ``marginal_cost * p`` over
    production-cost links and generators at the resolved cap. Like the other
    caps, baseline dispatch is not automatically feasible: a cap below the
    baseline production cost forces reallocation even at zero water price.
    """
    if not cap_cfg["enabled"]:
        return
    cap = resolve_cost_cap(cap_cfg, inputs)
    if cap < 0.0:
        raise ValueError(
            f"affordability.cost_cap resolved to a negative cap ({cap} bnUSD)"
        )
    link_mc = production_cost_links(n)
    gen_mc = production_cost_generators(n)
    if link_mc.empty and gen_mc.empty:
        logger.info("No cost-bearing components found; skipping affordability cap")
        return
    terms = {}
    if not link_mc.empty:
        terms["Link"] = (link_mc.index, link_mc.to_numpy())
    if not gen_mc.empty:
        terms["Generator"] = (gen_mc.index, gen_mc.to_numpy())
    add_scalar_dispatch_cap(n, terms, cap, _LABEL)
    logger.info(
        "Added production cost cap at %.1f bnUSD over %d links + %d generators",
        cap,
        len(link_mc),
        len(gen_mc),
    )


def assign_production_cost_cap_duals(n: pypsa.Network) -> None:
    """Persist the production-cost-cap shadow price onto the registry (``mu``)."""
    assign_scalar_cap_dual(n, _LABEL)
