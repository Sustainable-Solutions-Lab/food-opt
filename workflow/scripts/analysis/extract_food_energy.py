# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Extract per-country food-energy production, trade and floor metrics.

Uses the ``energy_tkcal_per_unit`` link coefficients and the per-bus
``energy_tkcal_per_mt`` densities stamped at solve time (see
solve_model/food_energy.py), so the solved network file is self-contained.
One row per country with:

- ``production_tkcal`` / ``baseline_production_tkcal``: primal human-edible
  energy production (realized / at baseline quantities).
- ``floor_tkcal`` / ``floor_mu``: food-energy floor level and shadow price
  when the constraint was active.
- ``import_tkcal`` / ``export_tkcal`` / ``net_import_tkcal``: gross and net
  hub-trade flows of crops and foods at human-edible energy (crops at
  their best-food-use density; the two-tier convention -- feed trade is
  reported separately in mass terms as ``feed_import_mt_dm`` /
  ``feed_export_mt_dm`` since feed calories are not human-edible in kind).
- ``consumption_tkcal``: final food energy consumed, the natural
  denominator for import dependence.

All energies in Tkcal = 1e12 kcal (~ the annual intake of ~1.2 million
people at 2300 kcal/day).
"""

import logging

import pandas as pd
import pypsa

logger = logging.getLogger(__name__)

TRADE_CARRIERS = ("trade_crop", "trade_food")


def extract_food_energy(n: pypsa.Network) -> pd.DataFrame:
    """Per-country food-energy metrics; empty if coefficients are absent."""
    links_df = n.links.static
    buses_df = n.buses.static
    if (
        "energy_tkcal_per_unit" not in links_df.columns
        or "energy_tkcal_per_mt" not in buses_df.columns
    ):
        logger.info("Network carries no food energy coefficients; writing empty output")
        return pd.DataFrame()

    p0 = n.links.dynamic.p0.loc["now"]
    coeff = pd.to_numeric(links_df["energy_tkcal_per_unit"], errors="coerce")
    covered = links_df[coeff.notna()]
    is_animal = covered["carrier"] == "animal_production"

    def _numeric(col: str) -> pd.Series:
        if col not in covered.columns:
            return pd.Series(0.0, index=covered.index)
        return pd.to_numeric(covered[col], errors="coerce").fillna(0.0)

    baseline_qty = _numeric("baseline_area_mha").where(
        ~is_animal, _numeric("baseline_feed_use_mt_dm")
    )
    prod = coeff[covered.index] * p0.reindex(covered.index).fillna(0.0)
    base = coeff[covered.index] * baseline_qty
    countries = covered["country"].astype(str)
    out = pd.DataFrame(
        {
            "production_tkcal": prod.groupby(countries.to_numpy()).sum(),
            "baseline_production_tkcal": base.groupby(countries.to_numpy()).sum(),
        }
    )
    out.index.name = "country"

    bus_country = buses_df["country"].astype(object)
    bus_country = bus_country.where(
        bus_country.notna() & (bus_country != ""), other=pd.NA
    )
    bus_density = pd.to_numeric(buses_df["energy_tkcal_per_mt"], errors="coerce")

    def _trade_flows(carriers, density_of_bus) -> tuple[pd.Series, pd.Series]:
        trade = links_df[links_df["carrier"].isin(carriers)]
        if trade.empty:
            zero = pd.Series(0.0, index=out.index)
            return zero, zero
        flows = p0.reindex(trade.index).fillna(0.0)
        c0 = trade["bus0"].map(bus_country)
        c1 = trade["bus1"].map(bus_country)
        # The country-side bus carries the density; hub buses do not.
        density = (
            trade["bus0"]
            .map(density_of_bus)
            .fillna(trade["bus1"].map(density_of_bus))
            .fillna(0.0)
        )
        weighted = density * flows
        exports = weighted[c0.notna()].groupby(c0[c0.notna()].to_numpy()).sum()
        imports = weighted[c1.notna()].groupby(c1[c1.notna()].to_numpy()).sum()
        return (
            imports.reindex(out.index).fillna(0.0),
            exports.reindex(out.index).fillna(0.0),
        )

    imports, exports = _trade_flows(TRADE_CARRIERS, bus_density)
    out["import_tkcal"] = imports
    out["export_tkcal"] = exports
    out["net_import_tkcal"] = imports - exports

    unit_density = pd.Series(1.0, index=buses_df.index)
    feed_imp, feed_exp = _trade_flows(("trade_feed",), unit_density)
    out["feed_import_mt_dm"] = feed_imp
    out["feed_export_mt_dm"] = feed_exp

    consume = links_df[links_df["carrier"] == "food_consumption"]
    if not consume.empty:
        cons_energy = consume["bus0"].map(bus_density).fillna(0.0) * p0.reindex(
            consume.index
        ).fillna(0.0)
        out["consumption_tkcal"] = (
            cons_energy.groupby(consume["country"].astype(str).to_numpy())
            .sum()
            .reindex(out.index)
            .fillna(0.0)
        )
    else:
        out["consumption_tkcal"] = 0.0

    gc = n.global_constraints.static
    floor_rows = gc[gc["type"] == "food_energy_floor"]
    if not floor_rows.empty:
        floors = floor_rows.set_index(floor_rows["country"].astype(str))
        out["floor_tkcal"] = floors["constant"].reindex(out.index)
        mu = (
            pd.to_numeric(floors["mu"], errors="coerce")
            if "mu" in floors.columns
            else pd.Series(dtype=float)
        )
        out["floor_mu"] = mu.reindex(out.index)
    else:
        out["floor_tkcal"] = pd.NA
        out["floor_mu"] = pd.NA

    return out.reset_index().sort_values("country").reset_index(drop=True)
