# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Biodiversity guardrail: a cap on natural-land conversion to agriculture.

Water-scarcity relief reorganises production across space, and relocating
output to water-rich regions requires converting natural land to cropland
and pasture there. This module bounds the total converted area at solve
time, so the water-relief optimisation cannot buy scarcity reductions with
unlimited habitat loss.

The cap sums dispatch over the conversion links -- ``land_conversion``
(natural land -> cropland) and ``new_to_pasture`` (natural land ->
pasture), whose ``Link-p`` is the converted area in Mha -- and bounds it
from above. It is distinct from GHG pricing, which charges the land-use-
change CO2 on the same links' efficiency: a pure *area* cap also protects
biodiverse-but-low-carbon systems (savanna, steppe, drylands) that a
carbon price undervalues.

The dual (bnUSD of objective per Mha of conversion headroom) is the
marginal system cost of protecting one more Mha of natural land, persisted
to ``n.global_constraints`` for analysis.
"""

import logging

import numpy as np
import pypsa

from workflow.scripts.solve_model.guardrails import (
    add_scalar_dispatch_cap,
    assign_scalar_cap_dual,
    load_barrier_reference,
)

logger = logging.getLogger(__name__)

CONVERSION_CARRIERS = ("land_conversion", "new_to_pasture")

_LABEL = "biodiversity_conversion_cap"


def resolve_conversion_cap(cap_cfg: dict, reference_path: str | None = None) -> float:
    """Resolve the configured cap to Mha.

    ``max_conversion_mha`` is the absolute cap. If ``relaxation`` (a 0-1
    dial) is set, ``max_conversion_mha`` instead anchors the loose end and
    the cap scales linearly, ``cap = relaxation * anchor``: relaxation 0
    forbids all new conversion (pinning the axis at the 2020 land state,
    which needs none) and relaxation 1 is the loose anchor (off when the
    anchor exceeds the achievable conversion). The zero-conversion point is
    always feasible, so a relaxation dial keeps every level feasible -- the
    property a joint barrier sweep relies on.
    """
    spec = cap_cfg["max_conversion_mha"]
    if isinstance(spec, dict):
        if reference_path is None:
            raise ValueError("biodiversity cap requires a barrier_reference input")
        reference = load_barrier_reference(reference_path, _LABEL)
        anchor = float(reference.loc["global"]) * (1.0 + 1e-9) + 1e-9
    else:
        anchor = float(spec)
    relaxation = cap_cfg["relaxation"]
    if relaxation is None:
        return anchor
    r = float(relaxation)
    if not 0.0 <= r <= 1.0:
        raise ValueError(f"biodiversity.cap.relaxation must be in [0, 1], got {r}")
    return r * anchor


def add_conversion_cap_constraints(
    n: pypsa.Network, cap_cfg: dict, reference_path: str | None = None
) -> None:
    """Cap total natural-land conversion to agriculture at solve time.

    ``cap_cfg`` is the ``biodiversity.cap`` config block with keys
    ``enabled`` (bool), ``max_conversion_mha`` (float), and ``relaxation``
    (a 0-1 dial or null; see :func:`resolve_conversion_cap`). The bound is
    global::

        sum over links L with carrier in CONVERSION_CARRIERS of p[L]
            <= cap

    In the absolute form, baseline dispatch is not automatically feasible: a
    cap below the baseline conversion level forces reallocation even at zero
    water price, which is intentional (the guardrail bites). The relaxation
    form instead anchors the tight end at zero conversion, which the 2020
    land state always satisfies.
    """
    if not cap_cfg["enabled"]:
        return
    cap = resolve_conversion_cap(cap_cfg, reference_path)
    if cap < 0.0:
        raise ValueError(f"biodiversity.cap resolved to a negative cap ({cap} Mha)")

    links_df = n.links.static
    conv_links = links_df[links_df["carrier"].isin(CONVERSION_CARRIERS)]
    if conv_links.empty:
        logger.info("No conversion links found; skipping biodiversity cap")
        return

    add_scalar_dispatch_cap(
        n, {"Link": (conv_links.index, np.ones(len(conv_links)))}, cap, _LABEL
    )
    logger.info(
        "Added biodiversity conversion cap at %.1f Mha over %d conversion links",
        cap,
        len(conv_links),
    )


def assign_conversion_cap_duals(n: pypsa.Network) -> None:
    """Persist the conversion-cap shadow price onto the registry (``mu``)."""
    assign_scalar_cap_dual(n, _LABEL)
