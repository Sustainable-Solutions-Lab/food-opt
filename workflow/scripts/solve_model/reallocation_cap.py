# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Endpoint-normalized cap on L1-weighted production reallocation."""

import numpy as np
import pandas as pd
import pypsa
import xarray as xr

from workflow.scripts.solve_model.production_stability import add_churn_abs_split

COMPONENT_CARRIERS = {
    "cropland": ("crop_production", "crop_production_multi"),
    "grassland": ("grassland_production",),
    "feed": ("animal_production",),
}


def _component_weights(dp_cfg: dict) -> dict[str, float]:
    if dp_cfg["penalty_mode"] != "l1" or dp_cfg["deviation_type"] != "absolute":
        raise ValueError(
            "reallocation_cap requires deviation_penalty penalty_mode='l1' "
            "and deviation_type='absolute'"
        )
    values = {
        "cropland": dp_cfg["land"]["crops"]["l1_cost"],
        "grassland": dp_cfg["land"]["grassland"]["l1_cost"],
        "feed": dp_cfg["feed"]["l1_cost"],
    }
    if any(value is None or not np.isscalar(value) for value in values.values()):
        raise ValueError("reallocation_cap requires explicit resolved L1 costs")
    return {component: float(value) for component, value in values.items()}


def add_reallocation_cap(
    n: pypsa.Network,
    cap_cfg: dict,
    dp_cfg: dict,
    reference_path: str | None,
) -> None:
    """Limit reallocation to a fraction of the realized endpoint distance.

    The distance is the weighted sum of absolute changes in crop and grassland
    area (Mha) and animal feed use (Mt DM). Component weights are the resolved
    L1 coefficients already used by the production-stability objective.
    """
    if not cap_cfg["enabled"]:
        return
    if reference_path is None:
        raise ValueError("Enabled reallocation_cap requires a reference input")

    fraction = float(cap_cfg["fraction"])
    if not 0.0 <= fraction <= 1.0:
        raise ValueError(f"reallocation_cap fraction must be in [0, 1], got {fraction}")
    reference = pd.read_parquet(reference_path).set_index(["component", "link"])
    if reference.index.has_duplicates:
        raise ValueError(
            "Reallocation reference contains duplicate component/link rows"
        )

    weights = _component_weights(dp_cfg)
    link_p = n.model.variables["Link-p"].sel(snapshot="now")
    links = n.links.static
    weighted_distance = 0
    endpoint_distance = 0.0
    for component, carriers in COMPONENT_CARRIERS.items():
        names = links.index[links["carrier"].isin(carriers)]
        expected = pd.MultiIndex.from_product([[component], names])
        missing = expected.difference(reference.index)
        if not missing.empty:
            raise ValueError(
                f"Reallocation reference lacks {len(missing)} {component} links"
            )
        component_reference = reference.reindex(expected)
        baseline = xr.DataArray(
            component_reference["baseline_dispatch"].to_numpy(),
            coords={"name": names},
            dims="name",
        )
        endpoint = component_reference["endpoint_dispatch"].to_numpy()
        churn = add_churn_abs_split(
            n.model,
            names,
            link_p.sel(name=names) - baseline,
            f"reallocation_cap_{component}",
        )
        weight = weights[component]
        weighted_distance += weight * churn.sum()
        endpoint_distance += weight * float(
            np.abs(endpoint - baseline.to_numpy()).sum()
        )

    if endpoint_distance <= 0.0:
        raise ValueError("Reallocation reference endpoint has zero weighted distance")
    budget = fraction * endpoint_distance
    n.model.add_constraints(
        weighted_distance <= budget,
        name="GlobalConstraint-reallocation_cap",
    )
    n.global_constraints.add(
        "reallocation_cap",
        sense="<=",
        constant=budget,
        type="reallocation_cap",
    )
