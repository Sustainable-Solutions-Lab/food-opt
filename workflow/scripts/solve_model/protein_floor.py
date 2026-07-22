# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Per-country protein production accounting and floor.

The protein instance of the generic nutrient floor (see
``solve_model/nutrient_floor.py``): a diet-quality twin of the calorie
floor. Calorie self-sufficiency is cheap to protect because staples are
water-flexible; protein self-sufficiency is the sterner test, since the
protein-dense production (pulses, dairy, animal products) is more
water-intensive. Densities come from the ``protein`` rows of nutrition.csv
(g/100g), scaled to Mt of protein per Mt of bus output.
"""

import pypsa

from workflow.scripts.solve_model.nutrient_floor import (
    add_nutrient_floor,
    assign_nutrient_floor_dual,
    stamp_nutrient_coefficients,
)

_LINK_COL = "protein_mt_per_unit"
_BUS_COL = "protein_mt_per_mt"
_LABEL = "protein_floor"


def stamp_protein_coefficients(n: pypsa.Network, nutrition_path: str) -> None:
    """Stamp Mt-protein-per-unit link and per-bus protein coefficients."""
    stamp_nutrient_coefficients(
        n,
        nutrition_path,
        nutrient="protein",
        expected_unit="g/100g",
        link_col=_LINK_COL,
        bus_col=_BUS_COL,
    )


def add_protein_floor_constraints(
    n: pypsa.Network, floor_cfg: dict, reference_path: str | None = None
) -> None:
    """Per-country lower bounds on primal protein production."""
    add_nutrient_floor(
        n,
        floor_cfg,
        link_col=_LINK_COL,
        label=_LABEL,
        min_baseline_key="min_baseline_mt",
        unit="Mt protein",
        reference_path=reference_path,
    )


def assign_protein_floor_duals(n: pypsa.Network) -> None:
    """Persist protein floor shadow prices for analysis."""
    assign_nutrient_floor_dual(n, _LABEL)
