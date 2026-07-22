# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Per-country agricultural production value accounting and floor constraint.

Production value is the gross production value of primary agriculture at
fixed producer prices (the FAOSTAT GPV convention): every marketed crop
output (single and multi-cropping links) and every animal product output
is valued at its per-(item, country) producer price; grassland output and
the fodder crops in ``production_value.exclude_crops`` are unmarketed
intermediates and are not valued. This intentionally double-counts feed
crops (a feed maize tonne is valued, and so is the resulting meat) --
gross output measures sector size, not value added.

The floor constraint bounds each country's production value from below at
a fraction of its baseline value, where the baseline uses the same per-link
coefficients evaluated at the baseline quantities (``baseline_area_mha``
for crops, ``baseline_feed_use_mt_dm`` for animals). Because floor and
realized value share prices and model conventions (post-seed, post-loss
output), level errors in prices largely cancel: the constraint binds on
the price-weighted quantity mix, not on absolute price levels.

This module also hosts the generic per-country output-floor machinery
(coefficient folding over link output ports, baseline evaluation,
constraint and dual assignment) shared with the food-energy floor in
``solve_model/food_energy.py``.

Interpretation caveats (see docs/production_value.rst):

- Observed producer prices include support-price distortions (e.g. Indian
  MSP), so the floor partially protects politically constructed revenue;
  the ``uniform`` price basis isolates the quantity-mix effect.
- The objective satisfies the floor with the cheapest available value,
  where "cheapest" includes cost-calibration corrections; check the value
  composition ex post before interpreting.
- Livestock fed on (imported) feed generates gross value with little land
  or water, so a binding floor can be met by expanding animal production;
  this is economically real but should be read jointly with feed trade.
"""

import logging

import numpy as np
import pandas as pd
import pypsa
import xarray as xr

from workflow.scripts import constants
from workflow.scripts.solve_model.guardrails import load_barrier_reference

logger = logging.getLogger(__name__)

VALUED_CROP_CARRIERS = ("crop_production", "crop_production_multi")

# USD/t on a bus in Mt: USD/t / (t/Mt) -> bnUSD/Mt
_USD_PER_T_TO_BNUSD_PER_MT = constants.USD_TO_BNUSD / constants.TONNE_TO_MEGATONNE


def _load_bus_prices(prices_path: str, price_basis: str) -> pd.Series:
    """Load producer prices on the bus basis, keyed by (item, country).

    ``price_basis`` selects between observed per-country producer prices
    and a uniform world price per item.
    """
    if price_basis not in ("observed", "uniform"):
        raise ValueError(
            f"production_value.floor.price_basis must be 'observed' or "
            f"'uniform', got '{price_basis}'"
        )
    column = (
        "price_usd_per_t_bus"
        if price_basis == "observed"
        else "price_usd_per_t_bus_uniform"
    )
    df = pd.read_csv(prices_path)
    return df.set_index(["item", "country"])[column].astype(float)


def _lookup(coeff: pd.Series, items: pd.Series, countries: pd.Series) -> pd.Series:
    """Look up per-Mt output coefficients keyed by item or (item, country)."""
    if coeff.index.nlevels == 2:
        keys = pd.MultiIndex.from_arrays([items.astype(str), countries.astype(str)])
    else:
        keys = items.astype(str)
    return pd.Series(coeff.reindex(keys).to_numpy(), index=items.index)


def fold_link_output_coefficients(
    n: pypsa.Network,
    crop_coeff: pd.Series,
    product_coeff: pd.Series,
    *,
    exclude_crops: set[str] | None = None,
    on_missing_crop: str = "raise",
) -> pd.Series:
    """Per-link output coefficient over marketable output ports.

    Sums ``efficiency * coeff`` over the crop-output ports of
    ``crop_production`` / ``crop_production_multi`` links (ports resolved
    to crops via the buses table, so multi-cropping links fold their
    per-cycle outputs into a single scalar per Mha of dispatch) and takes
    ``efficiency * coeff`` on the product output of ``animal_production``
    links (per Mt of feed DM). Coefficient series are per Mt of bus output,
    keyed by item or (item, country).

    ``on_missing_crop`` controls crops absent from ``crop_coeff``:
    ``"raise"`` (producer prices must cover every crop) or ``"zero"``
    (e.g. no human-edible energy). Missing animal products always raise.
    Returns a Series over all links, NaN for non-covered carriers.
    """
    excluded = set(exclude_crops or [])
    links_df = n.links.static
    coeff = pd.Series(0.0, index=links_df.index)

    crop_mask = links_df["carrier"].isin(VALUED_CROP_CARRIERS)
    crop_links = links_df[crop_mask]
    bus_to_crop = n.buses.static["crop"]
    port_cols = [("bus1", "efficiency")] + [
        (c, f"efficiency{c[3:]}")
        for c in links_df.columns
        if c.startswith("bus") and c != "bus0" and c != "bus1" and c[3:].isdigit()
    ]
    missing_crops: set[str] = set()
    for bus_col, eff_col in port_cols:
        if bus_col not in crop_links.columns or eff_col not in crop_links.columns:
            continue
        port_crop = crop_links[bus_col].map(bus_to_crop)
        has_crop = port_crop.notna() & (port_crop != "") & ~port_crop.isin(excluded)
        if not has_crop.any():
            continue
        sub = crop_links[has_crop]
        port_coeff = _lookup(crop_coeff, port_crop[has_crop], sub["country"])
        unmatched = port_coeff.isna()
        if unmatched.any():
            missing_crops.update(port_crop[has_crop][unmatched].unique())
            port_coeff = port_coeff.fillna(0.0)
        coeff.loc[sub.index] += (
            pd.to_numeric(sub[eff_col], errors="coerce").fillna(0.0) * port_coeff
        )
    if missing_crops and on_missing_crop == "raise":
        raise ValueError(
            f"Output coefficients missing for crops: {sorted(missing_crops)}"
        )

    animal_mask = links_df["carrier"] == "animal_production"
    animal_links = links_df[animal_mask]
    if not animal_links.empty:
        product_c = _lookup(
            product_coeff, animal_links["product"], animal_links["country"]
        )
        if product_c.isna().any():
            missing = sorted(animal_links["product"][product_c.isna()].unique())
            raise ValueError(f"Output coefficients missing for products: {missing}")
        coeff.loc[animal_links.index] = (
            pd.to_numeric(animal_links["efficiency"], errors="coerce").fillna(0.0)
            * product_c
        )

    return coeff.where(crop_mask | animal_mask)


def output_floor_baselines(n: pypsa.Network, coeff_column: str) -> pd.Series:
    """Baseline output (value / energy) per link from a stamped coefficient.

    Evaluates the coefficients at the baseline quantities:
    ``baseline_area_mha`` for crop links, ``baseline_feed_use_mt_dm`` for
    animal links. Links without a baseline contribute zero.
    """
    links_df = n.links.static
    coeff = links_df[coeff_column]
    valued = links_df[coeff.notna()]
    is_animal = (valued["carrier"] == "animal_production").to_numpy()

    def _qty(col: str) -> np.ndarray:
        if col not in valued.columns:
            return np.zeros(len(valued))
        return pd.to_numeric(valued[col], errors="coerce").fillna(0.0).to_numpy()

    qty = np.where(
        is_animal, _qty("baseline_feed_use_mt_dm"), _qty("baseline_area_mha")
    )
    return pd.Series(coeff[valued.index].to_numpy() * qty, index=valued.index)


def resolve_floor_fraction(floor_cfg: dict) -> float:
    """Effective floor fraction from the absolute or relaxation form.

    ``fraction`` is the direct fraction-of-baseline form. If ``relaxation``
    (a 0-1 dial) is set it overrides ``fraction`` and self-anchors the floor
    at baseline: ``fraction = 1 - relaxation``, so relaxation 0 pins output
    at baseline (fraction 1, always feasible) and relaxation 1 removes the
    floor (fraction 0). This mirrors the biodiversity/concentration caps'
    relaxation dial so every guardrail runs baseline (0) -> unconstrained (1).
    """
    relaxation = floor_cfg["relaxation"]
    if relaxation is None:
        return float(floor_cfg["fraction"])
    r = float(relaxation)
    if not 0.0 <= r <= 1.0:
        raise ValueError(f"floor relaxation must be in [0, 1], got {r}")
    return 1.0 - r


def add_output_floor_constraints(
    n: pypsa.Network,
    floor_cfg: dict,
    *,
    coeff_column: str,
    label: str,
    min_baseline_key: str,
    unit: str,
    reference_path: str | None = None,
) -> None:
    """Add per-country lower bounds on a stamped link-output aggregate.

    For each country c with baseline aggregate above the configured
    threshold::

        sum over links L in c of coeff[L] * p[L] >= fraction * baseline[c]

    where ``fraction`` comes from :func:`resolve_floor_fraction` (absolute
    or relaxation form). At baseline dispatch the two sides are equal by
    construction, so any fraction <= 1 is feasible whenever the baseline
    itself is feasible; joint infeasibility with tight water/emission caps
    at high fractions is informative (the country cannot retain that share
    within the caps), not a bug.
    """
    if not floor_cfg["enabled"]:
        return
    reference_scenario = floor_cfg.get("reference_scenario")
    if reference_scenario is None:
        fraction = resolve_floor_fraction(floor_cfg)
        if fraction < 0.0:
            raise ValueError(f"{label} fraction must be >= 0, got {fraction}")
        if fraction == 0.0:
            logger.info("%s fraction 0; non-binding, skipping", label)
            return
    elif reference_path is None:
        raise ValueError(f"{label} requires a barrier_reference input")
    min_baseline = float(floor_cfg[min_baseline_key])

    links_df = n.links.static
    if coeff_column not in links_df.columns:
        raise RuntimeError(
            f"{coeff_column} not stamped; the coefficient stamping must run "
            f"before add_output_floor_constraints({label})"
        )

    baseline = output_floor_baselines(n, coeff_column)
    countries = links_df.loc[baseline.index, "country"].astype(str)
    baseline_per_country = baseline.groupby(countries.to_numpy()).sum()
    eligible = baseline_per_country[baseline_per_country >= min_baseline].sort_index()
    skipped = baseline_per_country.index.difference(eligible.index)
    if len(skipped):
        logger.info(
            "%s: exempting %d countries below %.3f %s baseline (%s)",
            label,
            len(skipped),
            min_baseline,
            unit,
            ", ".join(skipped[:8]),
        )
    if eligible.empty:
        logger.info("%s: no countries above the baseline threshold; skipping", label)
        return

    coeff = links_df[coeff_column]
    sel = links_df[
        coeff.notna() & (coeff > 0) & links_df["country"].isin(eligible.index)
    ]
    if reference_scenario is None:
        floors = fraction * eligible
        mode = f"fraction {fraction:.3f}"
    else:
        reference = load_barrier_reference(reference_path, label)
        missing = eligible.index.difference(reference.index)
        if len(missing):
            raise ValueError(f"{label} reference is missing countries: {list(missing)}")
        floors = reference.reindex(eligible.index).clip(lower=0.0) * (1.0 - 1e-9)
        mode = f"reference {reference_scenario}"
    dim = f"{label}_country"

    m = n.model
    link_p = m.variables["Link-p"].sel(snapshot="now")
    weights = xr.DataArray(
        coeff[sel.index].to_numpy(dtype=float),
        coords={"name": sel.index},
        dims="name",
    )
    country_map = xr.DataArray(
        sel["country"].astype(str).to_numpy(),
        coords={"name": sel.index},
        dims="name",
        name=dim,
    )
    expr = (link_p.sel(name=sel.index) * weights).groupby(country_map).sum()
    expr_countries = pd.Index(expr.coords[dim].values)
    lower_bounds = xr.DataArray(
        floors.reindex(expr_countries).to_numpy(dtype=float),
        coords={dim: expr_countries.to_numpy()},
        dims=dim,
    )
    if lower_bounds.isnull().any():
        missing = expr_countries[np.isnan(lower_bounds.to_numpy())]
        raise RuntimeError(
            f"Countries with covered links but no baseline floor: {list(missing)}"
        )

    m.add_constraints(expr >= lower_bounds, name=f"GlobalConstraint-{label}")

    n.global_constraints.add(
        [f"{label}_{c}" for c in expr_countries],
        sense=">=",
        constant=lower_bounds.to_numpy(),
        type=label,
        country=list(expr_countries),
    )

    logger.info(
        "Added %s (%s) over %d countries; baseline total %.1f %s, "
        "floored total %.1f %s",
        label,
        mode,
        len(eligible),
        float(baseline_per_country.sum()),
        unit,
        float(floors.sum()),
        unit,
    )


def assign_output_floor_duals(n: pypsa.Network, label: str) -> None:
    """Persist floor shadow prices onto the global-constraints registry.

    Runs after ``assign_duals``; writes the per-country dual (bnUSD of
    objective per floored unit -- the marginal system cost of protecting
    one unit in that country) into the ``mu`` column so it survives the
    netCDF round trip into analysis.
    """
    cname = f"GlobalConstraint-{label}"
    if n.model is None or cname not in n.model.constraints:
        return
    dual = n.model.constraints[cname].dual
    countries = [str(c) for c in dual.coords[f"{label}_country"].values]
    names = [f"{label}_{c}" for c in countries]
    n.global_constraints.static.loc[names, "mu"] = dual.to_numpy()
    logger.info(
        "Assigned %s duals for %d countries (max |mu| = %.4f)",
        label,
        len(names),
        float(np.abs(dual.to_numpy()).max()) if len(names) else 0.0,
    )


def stamp_production_value_coefficients(
    n: pypsa.Network,
    prices_path: str,
    price_basis: str,
    exclude_crops: list[str] | None = None,
) -> None:
    """Stamp ``value_bnusd_per_unit`` on crop and animal production links.

    The coefficient converts link dispatch into gross production value:
    bnUSD per Mha of land dispatch for crop links and bnUSD per Mt of feed
    DM for animal links. Stamped unconditionally so the solved network
    supports value accounting in analysis even when no floor constraint is
    active.

    Crops in ``exclude_crops`` (unmarketed fodder intermediates, priced
    only by proxy) contribute no value: valuing them lets endogenous
    fodder expansion masquerade as production value and satisfy the
    floor. Animal co-product ports (rendered fats, bus5+) carry no
    producer price and are likewise not valued.
    """
    prices = _load_bus_prices(prices_path, price_basis) * _USD_PER_T_TO_BNUSD_PER_MT
    coeff = fold_link_output_coefficients(
        n,
        prices,
        prices,
        exclude_crops=set(exclude_crops or []),
        on_missing_crop="raise",
    )
    n.links.static["value_bnusd_per_unit"] = coeff
    logger.info(
        "Stamped production value coefficients (%s prices) on %d links",
        price_basis,
        int(coeff.notna().sum()),
    )


def production_value_baselines(n: pypsa.Network) -> pd.Series:
    """Baseline production value per link (bnUSD)."""
    return output_floor_baselines(n, "value_bnusd_per_unit").rename(
        "baseline_value_bnusd"
    )


def add_production_value_floor_constraints(
    n: pypsa.Network, floor_cfg: dict, reference_path: str | None = None
) -> None:
    """Per-country lower bounds on gross agricultural production value."""
    add_output_floor_constraints(
        n,
        floor_cfg,
        coeff_column="value_bnusd_per_unit",
        label="production_value_floor",
        min_baseline_key="min_baseline_bnusd",
        unit="bnUSD",
        reference_path=reference_path,
    )


def assign_production_value_floor_duals(n: pypsa.Network) -> None:
    """Persist production value floor shadow prices for analysis."""
    assign_output_floor_duals(n, "production_value_floor")
