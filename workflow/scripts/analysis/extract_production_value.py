# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Extract per-country agricultural production value from a solved network.

Uses the ``value_bnusd_per_unit`` coefficients stamped on production links
at solve time (see solve_model/production_value.py), so no price inputs
are needed here. Emits one row per (country, carrier, item) with realized
and baseline gross production value, plus the per-country floor level and
shadow price when the production-value floor constraint was active.

The ``item`` column is the crop for single-crop links, the combination
label for multi-cropping links (their per-cycle value is folded into one
coefficient at solve time), and the product for animal links.
"""

import logging

import pandas as pd
import pypsa

logger = logging.getLogger(__name__)


def extract_production_value(n: pypsa.Network) -> pd.DataFrame:
    """Per-(country, carrier, item) production value with floor metadata.

    Columns: country, carrier, item, realized_value_bnusd,
    baseline_value_bnusd, floor_bnusd, floor_mu. The floor columns are
    per-country values replicated across that country's rows; NaN when no
    floor constraint was active for the country.
    """
    links_df = n.links.static
    if "value_bnusd_per_unit" not in links_df.columns:
        logger.info(
            "Network carries no production value coefficients; "
            "writing empty production_value output"
        )
        return pd.DataFrame()

    coeff = pd.to_numeric(links_df["value_bnusd_per_unit"], errors="coerce")
    valued = links_df[coeff.notna()].copy()
    if valued.empty:
        return pd.DataFrame()
    valued["value_bnusd_per_unit"] = coeff[valued.index]

    p0 = n.links.dynamic.p0.loc["now"].reindex(valued.index).fillna(0.0)
    is_animal = valued["carrier"] == "animal_production"

    def _numeric(col: str) -> pd.Series:
        if col not in valued.columns:
            return pd.Series(0.0, index=valued.index)
        return pd.to_numeric(valued[col], errors="coerce").fillna(0.0)

    baseline_qty = _numeric("baseline_area_mha").where(
        ~is_animal, _numeric("baseline_feed_use_mt_dm")
    )
    valued["realized_value_bnusd"] = valued["value_bnusd_per_unit"] * p0
    valued["baseline_value_bnusd"] = valued["value_bnusd_per_unit"] * baseline_qty
    valued["item"] = valued["crop"].where(~is_animal, valued["product"])

    out = (
        valued.groupby(["country", "carrier", "item"], observed=True)[
            ["realized_value_bnusd", "baseline_value_bnusd"]
        ]
        .sum()
        .reset_index()
    )

    # Floor level and shadow price from the global-constraints registry
    # (populated by add_production_value_floor_constraints and
    # assign_production_value_floor_duals; absent when the floor is off).
    gc = n.global_constraints.static
    floor_rows = gc[gc["type"] == "production_value_floor"]
    if not floor_rows.empty:
        floors = floor_rows.set_index(floor_rows["country"].astype(str))
        out["floor_bnusd"] = out["country"].map(floors["constant"])
        mu = floors["mu"] if "mu" in floors.columns else pd.Series(dtype=float)
        out["floor_mu"] = out["country"].map(mu)
    else:
        out["floor_bnusd"] = pd.NA
        out["floor_mu"] = pd.NA

    return out.sort_values(["country", "carrier", "item"]).reset_index(drop=True)
