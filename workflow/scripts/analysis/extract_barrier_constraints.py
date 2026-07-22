# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Extract active guardrail constraints and their shadow prices.

Serialises the solved network's ``global_constraints`` registry -- the
uniform "barrier" object of the barriers study. Every solve-time guardrail
(production-value and food-energy floors, the biodiversity conversion cap,
the per-crop production-concentration caps, the reforestation cap) records
its bound in ``constant`` and, after dual assignment, its shadow price in
``mu``. One row per constraint, so relief-vs-tightness and cost-of-
protection can be read for all barriers from a single output.
"""

import logging

import pandas as pd
import pypsa

logger = logging.getLogger(__name__)

# Columns of interest from the registry; kept if present. ``max_share`` is the
# per-crop concentration cap (the constraint's ``constant`` is 0 for the
# rearranged max-share form, so the resolved share lives in its own column).
_COLS = ["type", "sense", "constant", "mu", "crop", "country", "carrier", "max_share"]


def extract_barrier_constraints(n: pypsa.Network) -> pd.DataFrame:
    """One row per solve-time guardrail constraint with its bound and dual.

    Columns: name, plus whichever of type/sense/constant/mu/crop/country/
    carrier the registry carries. Empty when no guardrails were active.
    """
    gc = n.global_constraints.static
    if gc.empty:
        logger.info("No global constraints registered; empty barrier output")
        return pd.DataFrame()

    cols = [c for c in _COLS if c in gc.columns]
    out = gc[cols].copy()
    out.insert(0, "name", gc.index.astype(str))
    if "mu" not in out.columns:
        out["mu"] = pd.NA
    return out.reset_index(drop=True)
