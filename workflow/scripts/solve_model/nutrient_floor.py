# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Generic per-country nutrient-production accounting and floor machinery.

A nutrient self-sufficiency floor bounds each country's primal production of
a nutrient (calories, protein, ...) from below at a fraction of its baseline.
Nutrient densities come from nutrition.csv (per-food, in nutrient-unit per
100 g); a crop's density is the maximum over its food-processing pathways of
``pathway efficiency * food nutrient density`` -- the nutrient obtainable
from a tonne of crop in its best food use. Crops with no food pathway carry
zero and drop out.

Primal convention (shared with the value floor in ``production_value.py``):
nutrient content counts where a crop is grown regardless of destination, and
animal products count at their own density, so feed nutrients are
double-counted. The floor binds on where the nutrient is grown, not how it is
routed. The food-energy (calorie) and protein floors are thin wrappers around
this module with different nutrient rows.
"""

import logging

import pandas as pd
import pypsa

from workflow.scripts.solve_model.production_value import (
    add_output_floor_constraints,
    assign_output_floor_duals,
    fold_link_output_coefficients,
    output_floor_baselines,
)

logger = logging.getLogger(__name__)

# nutrition.csv values are per 100 g; a per-100g density in X/100g equals
# X-mass per Mt of food * 1e-2 once X and food masses share their Mt/t scale
# (Tkcal per Mt for kcal, Mt per Mt for grams-of-protein-in-Mt-basis).
_PER_100G_TO_PER_MT = 1.0 / 100.0


def nutrient_food_densities(
    nutrition_path: str, nutrient: str, expected_unit: str
) -> pd.Series:
    """Per-food nutrient density from nutrition.csv, scaled to a per-Mt basis."""
    nutrition = pd.read_csv(nutrition_path)
    rows = nutrition[nutrition["nutrient"] == nutrient]
    if rows.empty:
        raise ValueError(f"No '{nutrient}' rows in {nutrition_path}")
    bad_units = set(rows["unit"].unique()) - {expected_unit}
    if bad_units:
        raise ValueError(
            f"Unexpected units for nutrient '{nutrient}': {bad_units} "
            f"(expected {expected_unit})"
        )
    return rows.set_index("food")["value"].astype(float) * _PER_100G_TO_PER_MT


def crop_nutrient_densities(n: pypsa.Network, e_food: pd.Series) -> pd.Series:
    """Per-crop nutrient density (per Mt of crop-bus DM) via best food use.

    For each crop, the maximum over its food-processing output ports of
    ``efficiency * food nutrient density`` (pathway efficiencies are
    country-neutral, so the max over links is the max over pathways). Crops
    without any food pathway get no entry.
    """
    links_df = n.links.static
    proc = links_df[links_df["carrier"] == "food_processing"]
    if proc.empty:
        return pd.Series(dtype=float)
    bus_to_crop = n.buses.static["crop"]
    bus_to_food = n.buses.static["food"]
    crop_of_link = proc["bus0"].map(bus_to_crop)

    best = pd.Series(dtype=float)
    port_cols = [("bus1", "efficiency")] + [
        (c, f"efficiency{c[3:]}")
        for c in proc.columns
        if c.startswith("bus") and c not in ("bus0", "bus1") and c[3:].isdigit()
    ]
    for bus_col, eff_col in port_cols:
        if bus_col not in proc.columns or eff_col not in proc.columns:
            continue
        port_food = proc[bus_col].map(bus_to_food)
        valid = port_food.notna() & (port_food != "") & crop_of_link.notna()
        if not valid.any():
            continue
        density = (
            pd.to_numeric(proc.loc[valid, eff_col], errors="coerce").fillna(0.0)
            * port_food[valid].map(e_food).fillna(0.0).to_numpy()
        )
        port_best = density.groupby(crop_of_link[valid].to_numpy()).max()
        best = best.combine(port_best, max, fill_value=0.0)
    return best[best > 0]


def stamp_nutrient_coefficients(
    n: pypsa.Network,
    nutrition_path: str,
    *,
    nutrient: str,
    expected_unit: str,
    link_col: str,
    bus_col: str,
) -> None:
    """Stamp per-link and per-bus nutrient-production coefficients.

    ``link_col`` receives the per-Mha (crop) / per-Mt-feed (animal) nutrient
    output; ``bus_col`` the per-Mt bus density (so analysis can value trade and
    consumption flows from the network alone). Stamped unconditionally, so the
    solved network supports the accounting whether or not the floor is active.
    """
    e_food = nutrient_food_densities(nutrition_path, nutrient, expected_unit)
    e_crop = crop_nutrient_densities(n, e_food)
    coeff = fold_link_output_coefficients(n, e_crop, e_food, on_missing_crop="zero")
    n.links.static[link_col] = coeff
    buses = n.buses.static
    density = buses["food"].map(e_food)
    if "crop" in buses.columns:
        density = density.fillna(buses["crop"].map(e_crop))
    n.buses.static[bus_col] = density
    logger.info(
        "Stamped %s coefficients on %d links (%d crops with a food pathway, "
        "global baseline %.0f)",
        nutrient,
        int(coeff.notna().sum()),
        len(e_crop),
        float(output_floor_baselines(n, link_col).sum()),
    )


def add_nutrient_floor(
    n: pypsa.Network,
    floor_cfg: dict,
    *,
    link_col: str,
    label: str,
    min_baseline_key: str,
    unit: str,
    reference_path: str | None = None,
) -> None:
    """Per-country lower bound on primal production of a stamped nutrient."""
    add_output_floor_constraints(
        n,
        floor_cfg,
        coeff_column=link_col,
        label=label,
        min_baseline_key=min_baseline_key,
        unit=unit,
        reference_path=reference_path,
    )


def assign_nutrient_floor_dual(n: pypsa.Network, label: str) -> None:
    """Persist a nutrient floor's shadow prices for analysis."""
    assign_output_floor_duals(n, label)
