# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Emissions guardrail: a cap on total net food-system GHG emissions.

Water-scarcity relief reorganises production across space, and relocating
output to water-rich regions can require converting natural land (a land-use-
change CO2 pulse) and shifting the livestock mix (CH4 / N2O). This module
bounds the total net GHG accumulated on ``store:emission:ghg`` (MtCO2eq, net
of spared-land sequestration credits) from above at solve time, so relief
cannot be bought with an unbounded increase in emissions.

The capped quantity is the model's own aggregate emissions store: every CO2,
CH4 and N2O source (and the negative spared-land credit) flows into it, so the
bound is on *net* CO2-equivalent -- the quantity a "keep emissions from rising"
policy actually targets. It is distinct from the biodiversity area cap (which
bounds converted hectares regardless of their carbon) and from GHG pricing
(which charges emissions in the objective rather than bounding them).

The bound is a one-sided linopy constraint on the store's energy level::

    Store-e[now, store:emission:ghg]  <=  cap

added after ``create_model``. A one-sided upper bound is required rather than
the store-capacity (``e_nom_max``) trick used for the water-scarcity cap: the
emissions store is bidirectional (``e_min_pu = -1``), so an ``e_nom_max`` bound
would symmetrically cap sequestration (the negative side) too. The dual (bnUSD
of objective per MtCO2eq of headroom -- the endogenous carbon shadow price) is
persisted to ``n.global_constraints`` for analysis.
"""

import logging

import pandas as pd
import pypsa

from workflow.scripts.solve_model.guardrails import assign_scalar_cap_dual

logger = logging.getLogger(__name__)

_STORE = "store:emission:ghg"
_LABEL = "emissions_cap"


def _reference_net_emissions(path: str) -> float:
    """Total net GHG (MtCO2eq) from a reference scenario's net_emissions.parquet."""
    df = pd.read_parquet(path)
    if df.empty:
        raise ValueError(f"Reference net-emissions output is empty: {path}")
    return float(df["mtco2eq"].sum())


def resolve_emissions_cap(cap_cfg: dict, inputs) -> float:
    """Resolve the configured cap to MtCO2eq.

    ``max_mtco2eq`` is either a number (absolute cap) or a mapping
    ``{relaxation, reference_scenario}``: the cap then anchors on the reference
    scenario's realized net emissions ``E_ref`` (read from its
    ``net_emissions.parquet`` analysis output via the ``emissions_reference``
    input) and scales as ``cap = E_ref * (1 + relaxation)``. Relaxation 0 pins
    the cap at the reference emissions (non-increasing: the reference state is
    exactly feasible) and relaxation 1 allows twice the reference. Because the
    reference (e.g. frozen 2020) state emits ``E_ref`` at zero conversion, it
    satisfies the cap at every relaxation, keeping every level feasible -- the
    property a joint barrier sweep relies on.
    """
    spec = cap_cfg["max_mtco2eq"]
    if not isinstance(spec, dict):
        return float(spec)
    relaxation = float(spec["relaxation"])
    if relaxation < 0.0:
        raise ValueError(
            f"emissions.cap.max_mtco2eq.relaxation must be >= 0, got {relaxation}"
        )
    reference = _reference_net_emissions(inputs.emissions_reference)
    cap = reference * (1.0 + relaxation)
    logger.info(
        "Resolved emissions cap: reference %.1f (%s), relaxation %.2f -> cap %.1f MtCO2eq",
        reference,
        spec["reference_scenario"],
        relaxation,
        cap,
    )
    return cap


def add_emissions_cap(n: pypsa.Network, cap_cfg: dict, inputs) -> None:
    """Cap total net GHG emissions at solve time (one-sided epsilon-constraint).

    ``cap_cfg`` is the ``emissions.cap`` block with keys ``enabled`` (bool) and
    ``max_mtco2eq`` (absolute MtCO2eq or the reference-relative mapping resolved
    by :func:`resolve_emissions_cap`, which reads the reference net-emissions
    file from ``inputs``). Bounds the ``store:emission:ghg`` energy level from
    above, so its accumulated net CO2-equivalent cannot exceed the cap while the
    spared-land sequestration (negative) side stays unconstrained.

    Must be called AFTER ``n.optimize.create_model()`` so the ``Store-e``
    variable exists. Registers a one-row global constraint for dual bookkeeping.
    """
    if not cap_cfg["enabled"]:
        return
    cap = resolve_emissions_cap(cap_cfg, inputs)
    m = n.model
    store_e = m.variables["Store-e"].sel(snapshot="now").sel(name=_STORE)
    m.add_constraints(store_e <= cap, name=f"GlobalConstraint-{_LABEL}")
    n.global_constraints.add(_LABEL, sense="<=", constant=cap, type=_LABEL)
    logger.info("Added net-emissions cap at %.1f MtCO2eq on %s", cap, _STORE)


def assign_emissions_cap_dual(n: pypsa.Network) -> None:
    """Persist the emissions cap's shadow price onto the registry (``mu``)."""
    assign_scalar_cap_dual(n, _LABEL)
