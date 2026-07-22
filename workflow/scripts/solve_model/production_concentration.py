# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Stability guardrail: a per-crop cap on production concentration.

Water-scarcity relief can concentrate the world's production of a staple
into fewer countries; a shock or export ban then hits harder. This module
caps, for each broadly-grown crop, the share of global production any
single country may hold -- the linear "max-share" form of a concentration
limit (the literal Herfindahl is convex-quadratic and would leave the LP
regime):

    x[c, j]  <=  s * sum_k x[c, k]      for every producing country j

where ``x[c, j]`` is country j's output of crop c in Mt, summed over its
``crop_production`` / ``crop_production_multi`` links and their crop-output
ports (``efficiency_port * p[link]``). Rearranged to
``x[c,j] - s * X[c] <= 0`` the constraint is fully linear even with the
global total ``X[c]`` endogenous. The dual per (crop, country) is the
marginal system cost of forcing that crop's production to spread -- the
cost of resilience -- persisted to ``n.global_constraints`` for analysis.

Only crops that are grown in enough countries and above a global-tonnage
floor are capped: a crop producible in very few places cannot satisfy a
share cap below its structural minimum, so capping it would be infeasible
rather than informative.

Trade is hub-based (k-means), so bilateral import-source concentration is
not representable; production concentration sidesteps the hub layer and is
the cleaner GLADE-native fragility metric.
"""

import logging

import numpy as np
import pandas as pd
import pypsa
import xarray as xr

from workflow.scripts.solve_model.guardrails import load_barrier_reference
from workflow.scripts.solve_model.production_value import VALUED_CROP_CARRIERS

logger = logging.getLogger(__name__)

_LABEL = "production_concentration"
_DIM = "concentration_country"


def _crop_output_terms(n: pypsa.Network) -> pd.DataFrame:
    """Per-(link, crop, country) physical crop-output coefficient (Mt/Mha).

    One row per crop-output port of every crop-production link: the port's
    ``efficiency`` (Mt of crop per Mha of land dispatch) tagged with the
    crop it outputs (resolved via the buses table, so multi-cropping links
    contribute one row per crop they grow) and the producing country. Rows
    for the same (link, crop) are summed defensively.
    """
    links_df = n.links.static
    crop_links = links_df[links_df["carrier"].isin(VALUED_CROP_CARRIERS)]
    if crop_links.empty:
        return pd.DataFrame(columns=["link", "crop", "country", "eff", "baseline_mha"])
    bus_to_crop = n.buses.static["crop"]
    port_cols = [("bus1", "efficiency")] + [
        (c, f"efficiency{c[3:]}")
        for c in crop_links.columns
        if c.startswith("bus") and c not in ("bus0", "bus1") and c[3:].isdigit()
    ]
    baseline_mha = pd.to_numeric(
        crop_links.get("baseline_area_mha", 0.0), errors="coerce"
    ).fillna(0.0)
    rows = []
    for bus_col, eff_col in port_cols:
        if bus_col not in crop_links.columns or eff_col not in crop_links.columns:
            continue
        port_crop = crop_links[bus_col].map(bus_to_crop)
        has_crop = port_crop.notna() & (port_crop != "")
        if not has_crop.any():
            continue
        sub = crop_links[has_crop]
        eff = pd.to_numeric(sub[eff_col], errors="coerce").fillna(0.0)
        keep = eff > 0
        if not keep.any():
            continue
        rows.append(
            pd.DataFrame(
                {
                    "link": sub.index[keep],
                    "crop": port_crop[has_crop][keep].astype(str).to_numpy(),
                    "country": sub["country"][keep].astype(str).to_numpy(),
                    "eff": eff[keep].to_numpy(),
                    "baseline_mha": baseline_mha[sub.index[keep]].to_numpy(),
                }
            )
        )
    if not rows:
        return pd.DataFrame(columns=["link", "crop", "country", "eff", "baseline_mha"])
    terms = pd.concat(rows, ignore_index=True)
    return terms.groupby(["link", "crop", "country"], as_index=False).agg(
        eff=("eff", "sum"), baseline_mha=("baseline_mha", "first")
    )


def _crop_cap_share(sub: pd.DataFrame, cap_cfg: dict) -> float:
    """Resolve the max-share cap for one crop.

    In the absolute form (``relaxation`` null) every crop shares the scalar
    ``max_share``. With ``relaxation`` set (a 0-1 dial) the cap self-anchors
    on this crop's baseline concentration ``b`` (the largest producer's
    share of baseline output): ``cap = b + relaxation * (1 - b)``, so
    relaxation 0 pins the crop at its baseline concentration and relaxation
    1 leaves it uncapped (share 1). The baseline dispatch satisfies the
    former with equality, keeping every level feasible.
    """
    relaxation = cap_cfg["relaxation"]
    if relaxation is None:
        return float(cap_cfg["max_share"])
    r = float(relaxation)
    if not 0.0 <= r <= 1.0:
        raise ValueError(
            f"production_concentration.cap.relaxation must be in [0, 1], got {r}"
        )
    global_mt = float(sub["baseline_mt"].sum())
    baseline_share = (
        float(sub.groupby("country")["baseline_mt"].sum().max()) / global_mt
        if global_mt > 0.0
        else 1.0
    )
    return baseline_share + r * (1.0 - baseline_share)


def add_concentration_cap_constraints(
    n: pypsa.Network, cap_cfg: dict, reference_path: str | None = None
) -> None:
    """Add per-crop max-share caps on production concentration.

    ``cap_cfg`` is the ``production_concentration.cap`` block with keys
    ``enabled`` (bool), ``max_share`` (float in (0, 1]), ``relaxation`` (a
    0-1 dial or null; see :func:`_crop_cap_share`), ``min_global_mt``
    (exempt trace crops), and ``min_producing_countries`` (skip crops too
    concentrated by geography to satisfy the cap).
    """
    if not cap_cfg["enabled"]:
        return
    if cap_cfg["relaxation"] is None:
        s_abs = float(cap_cfg["max_share"])
        if not 0.0 < s_abs <= 1.0:
            raise ValueError(
                f"production_concentration.cap.max_share must be in (0, 1], got {s_abs}"
            )
        if s_abs >= 1.0:
            logger.info(
                "Concentration max_share %.3f >= 1.0; non-binding, skipping", s_abs
            )
            return
    min_global_mt = float(cap_cfg["min_global_mt"])
    min_countries = int(cap_cfg["min_producing_countries"])
    reference_scenario = cap_cfg.get("reference_scenario")
    reference = None
    if reference_scenario is not None:
        if reference_path is None:
            raise ValueError("concentration cap requires a barrier_reference input")
        reference = load_barrier_reference(reference_path, _LABEL)

    terms = _crop_output_terms(n)
    if terms.empty:
        logger.info("No crop-output terms found; skipping concentration cap")
        return
    terms["baseline_mt"] = terms["eff"] * terms["baseline_mha"]

    m = n.model
    link_p = m.variables["Link-p"].sel(snapshot="now")

    reg_names, reg_crops, reg_countries, cap_shares = [], [], [], []
    n_capped = 0
    for crop, sub in terms.groupby("crop"):
        global_mt = float(sub["baseline_mt"].sum())
        countries = sub["country"].unique()
        if global_mt < min_global_mt or len(countries) < min_countries:
            continue
        if reference is None:
            s = _crop_cap_share(sub, cap_cfg)
        else:
            if crop not in reference.index:
                raise ValueError(f"Concentration reference is missing crop '{crop}'")
            s = min(1.0, float(reference.loc[crop]) * (1.0 + 1e-9) + 1e-9)
        if s >= 1.0:
            # Self-anchored cap at or above 1 (relaxation 1, or a crop with a
            # degenerate single producer): non-binding, so skip it.
            continue
        weights = xr.DataArray(
            sub["eff"].to_numpy(dtype=float),
            coords={"name": sub["link"].to_numpy()},
            dims="name",
        )
        country_map = xr.DataArray(
            sub["country"].to_numpy(),
            coords={"name": sub["link"].to_numpy()},
            dims="name",
            name=_DIM,
        )
        expr = (
            (link_p.sel(name=sub["link"].to_numpy()) * weights)
            .groupby(country_map)
            .sum()
        )
        # A per-crop aggregate variable X[c] carries the global-total term so
        # the per-country cap rows stay sparse (each references only that
        # country's links plus X[c]); writing s * expr.sum() inline instead
        # puts every one of the crop's links into every country row, a dense
        # block that dominates the barrier factorization.
        total_var = m.add_variables(lower=0, name=f"concentration_total_{crop}")
        m.add_constraints(
            total_var - expr.sum() == 0.0,
            name=f"GlobalConstraint-{_LABEL}_total_{crop}",
        )
        m.add_constraints(
            expr - s * total_var <= 0.0, name=f"GlobalConstraint-{_LABEL}_{crop}"
        )
        con_countries = [str(c) for c in expr.coords[_DIM].values]
        reg_names.extend(f"{_LABEL}_{crop}_{c}" for c in con_countries)
        reg_crops.extend([crop] * len(con_countries))
        reg_countries.extend(con_countries)
        cap_shares.extend([s] * len(con_countries))
        n_capped += 1

    if not reg_names:
        logger.info(
            "Concentration cap: no crops met min_global_mt/min_producing_countries; "
            "skipping"
        )
        return

    n.global_constraints.add(
        reg_names,
        sense="<=",
        constant=0.0,
        type=_LABEL,
        crop=reg_crops,
        country=reg_countries,
        max_share=cap_shares,
    )
    mode = (
        f"reference {reference_scenario}"
        if reference_scenario is not None
        else f"relaxation {float(cap_cfg['relaxation']):.2f}"
        if cap_cfg["relaxation"] is not None
        else f"max_share {float(cap_cfg['max_share']):.3f}"
    )
    logger.info(
        "Added production concentration caps (%s; per-crop cap %.3f-%.3f) on "
        "%d crops, %d (crop, country) rows",
        mode,
        min(cap_shares),
        max(cap_shares),
        n_capped,
        len(reg_names),
    )


def assign_concentration_cap_duals(n: pypsa.Network) -> None:
    """Persist per-(crop, country) concentration shadow prices to ``mu``."""
    if n.model is None:
        return
    gc = n.global_constraints.static
    rows = gc[gc["type"] == _LABEL]
    if rows.empty:
        return
    max_mu = 0.0
    for crop, grp in rows.groupby("crop"):
        cname = f"GlobalConstraint-{_LABEL}_{crop}"
        if cname not in n.model.constraints:
            continue
        dual = n.model.constraints[cname].dual
        country_dual = {
            str(c): float(v) for c, v in zip(dual.coords[_DIM].values, np.asarray(dual))
        }
        for name, country in zip(grp.index, grp["country"].astype(str)):
            if country in country_dual:
                n.global_constraints.static.loc[name, "mu"] = country_dual[country]
                max_mu = max(max_mu, abs(country_dual[country]))
    logger.info(
        "Assigned concentration cap duals for %d rows (max |mu| = %.4f)",
        len(rows),
        max_mu,
    )
