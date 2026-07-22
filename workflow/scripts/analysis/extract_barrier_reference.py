# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Extract realized guardrail left-hand sides for a reference scenario."""

import pandas as pd
import pypsa

from workflow.scripts.solve_model.biodiversity import CONVERSION_CARRIERS
from workflow.scripts.solve_model.production_concentration import _crop_output_terms
from workflow.scripts.solve_model.production_value import output_floor_baselines

_FLOOR_COLUMNS = {
    "food_energy_floor": "energy_tkcal_per_unit",
    "protein_floor": "protein_mt_per_unit",
    "production_value_floor": "value_bnusd_per_unit",
}


def extract_barrier_reference(n: pypsa.Network) -> pd.DataFrame:
    """Return realized frozen-reference bounds for all non-scalar guardrails.

    Rows are keyed by ``barrier`` and ``key``. ``value`` is the realized
    guardrail left-hand side and ``baseline`` records the corresponding
    build-time quantity for auditing. The solve applies a tiny directional
    tolerance when consuming these values.
    """
    dispatch = n.links_t.p.loc["now"]
    rows: list[dict] = []

    conversion = n.links.static["carrier"].isin(CONVERSION_CARRIERS)
    rows.append(
        {
            "barrier": "biodiversity_conversion_cap",
            "key": "global",
            "value": float(dispatch[conversion].sum()),
            "baseline": float("nan"),
        }
    )

    terms = _crop_output_terms(n)
    if not terms.empty:
        terms["realized_mt"] = terms["eff"] * terms["link"].map(dispatch)
        terms["baseline_mt"] = terms["eff"] * terms["baseline_mha"]
        for crop, sub in terms.groupby("crop"):
            realized = sub.groupby("country")["realized_mt"].sum()
            baseline = sub.groupby("country")["baseline_mt"].sum()
            if realized.sum() <= 0.0 or baseline.sum() <= 0.0:
                continue
            rows.append(
                {
                    "barrier": "production_concentration",
                    "key": str(crop),
                    "value": float(realized.max() / realized.sum()),
                    "baseline": float(baseline.max() / baseline.sum()),
                }
            )

    links = n.links.static
    for label, coeff_column in _FLOOR_COLUMNS.items():
        coeff = links[coeff_column]
        covered = coeff.notna()
        countries = links.loc[covered, "country"].astype(str)
        realized = (coeff[covered] * dispatch[covered]).groupby(countries).sum()
        baseline = output_floor_baselines(n, coeff_column)
        baseline = baseline.groupby(
            links.loc[baseline.index, "country"].astype(str)
        ).sum()
        for country in baseline.index.union(realized.index):
            rows.append(
                {
                    "barrier": label,
                    "key": str(country),
                    "value": float(realized.get(country, 0.0)),
                    "baseline": float(baseline.get(country, 0.0)),
                }
            )

    return pd.DataFrame(rows).sort_values(["barrier", "key"]).reset_index(drop=True)
