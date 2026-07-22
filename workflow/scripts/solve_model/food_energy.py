# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Per-country food-energy (calorie) production accounting and floor.

The calorie instance of the generic nutrient floor (see
``solve_model/nutrient_floor.py``): each country's primal production of
human-edible food energy is bounded from below at a fraction of its
baseline -- the production-capacity notion of calorie self-sufficiency (the
SSR numerator). With fixed diets this caps how much more import-dependent a
country can become without attributing pooled hub-trade flows.

Energy densities come from the ``cal`` rows of nutrition.csv (kcal/100g),
scaled to Tkcal (1e12 kcal) per Mt of bus output.
"""

import pypsa

from workflow.scripts.solve_model.nutrient_floor import (
    add_nutrient_floor,
    assign_nutrient_floor_dual,
    stamp_nutrient_coefficients,
)

_LINK_COL = "energy_tkcal_per_unit"
_BUS_COL = "energy_tkcal_per_mt"
_LABEL = "food_energy_floor"


def stamp_food_energy_coefficients(n: pypsa.Network, nutrition_path: str) -> None:
    """Stamp Tkcal-per-unit link and per-bus food-energy coefficients."""
    stamp_nutrient_coefficients(
        n,
        nutrition_path,
        nutrient="cal",
        expected_unit="kcal/100g",
        link_col=_LINK_COL,
        bus_col=_BUS_COL,
    )


def add_food_energy_floor_constraints(
    n: pypsa.Network, floor_cfg: dict, reference_path: str | None = None
) -> None:
    """Per-country lower bounds on primal food-energy production."""
    add_nutrient_floor(
        n,
        floor_cfg,
        link_col=_LINK_COL,
        label=_LABEL,
        min_baseline_key="min_baseline_tkcal",
        unit="Tkcal",
        reference_path=reference_path,
    )


def assign_food_energy_floor_duals(n: pypsa.Network) -> None:
    """Persist food-energy floor shadow prices for analysis."""
    assign_nutrient_floor_dual(n, _LABEL)
