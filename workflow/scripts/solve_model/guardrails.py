# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Shared machinery for scalar "guardrail" caps on a weighted dispatch aggregate.

A guardrail bounds one linear function of dispatch from above::

    sum over components c of weight[c] * p[c]  <=  cap

and records its dual (the marginal system cost of one more unit of headroom)
on ``n.global_constraints`` so it survives into analysis. The biodiversity
conversion cap (weight = 1 over conversion links) and the affordability
production-cost cap (weight = per-component marginal cost over links and
generators) are both instances. The per-country output *floors* (production
value, food energy, protein) use the separate lower-bound machinery in
``production_value.py``.
"""

import logging

import numpy as np
import pandas as pd
import pypsa
import xarray as xr

logger = logging.getLogger(__name__)


def load_barrier_reference(path: str, barrier: str) -> pd.Series:
    """Load one barrier's realized reference values keyed by scope."""
    data = pd.read_parquet(path)
    required = {"barrier", "key", "value"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Barrier reference is missing columns: {sorted(missing)}")
    rows = data[data["barrier"] == barrier]
    if rows.empty:
        raise ValueError(f"Barrier reference has no rows for '{barrier}'")
    if rows["key"].duplicated().any():
        raise ValueError(f"Barrier reference has duplicate keys for '{barrier}'")
    return rows.set_index("key")["value"].astype(float)


def add_scalar_dispatch_cap(
    n: pypsa.Network,
    terms: dict[str, tuple[pd.Index, np.ndarray]],
    cap: float,
    label: str,
) -> None:
    """Add a single upper bound on a weighted sum of dispatch.

    ``terms`` maps a component name (``"Link"`` or ``"Generator"``) to a
    ``(index, weights)`` pair selecting that component's assets and their
    coefficients; the constraint bounds the sum of ``weights * p`` over all
    terms at ``cap``. Registers a one-row global constraint named ``label``
    for dual bookkeeping.
    """
    m = n.model
    expr = None
    for component, (index, weights) in terms.items():
        p = m.variables[f"{component}-p"].sel(snapshot="now")
        w = xr.DataArray(
            np.asarray(weights, dtype=float),
            coords={"name": np.asarray(index)},
            dims="name",
        )
        term = (p.sel(name=np.asarray(index)) * w).sum()
        expr = term if expr is None else expr + term
    m.add_constraints(expr <= cap, name=f"GlobalConstraint-{label}")
    n.global_constraints.add(label, sense="<=", constant=cap, type=label)


def assign_scalar_cap_dual(n: pypsa.Network, label: str) -> None:
    """Persist a scalar cap's shadow price onto the registry (``mu``)."""
    cname = f"GlobalConstraint-{label}"
    if n.model is None or cname not in n.model.constraints:
        return
    dual = float(np.asarray(n.model.constraints[cname].dual).reshape(-1)[0])
    n.global_constraints.static.loc[label, "mu"] = dual
    logger.info("Assigned %s dual (mu = %.4f)", label, dual)
