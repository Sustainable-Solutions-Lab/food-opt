# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Prepare per-(item, country) producer prices for crops and animal products.

Emits the price table consumed by the production-value accounting at solve
time (see workflow/scripts/solve_model/production_value.py). Prices come
from the FAOSTAT PP domain (element 5532, Producer Price USD/tonne),
CPI-deflated to the configured base year and averaged over the configured
period -- the same source, deflation and averaging used for crop costs in
prepare_faostat_crop_costs.py, so value and cost are in consistent real USD.

Two price bases are emitted per row:

- ``price_usd_per_t_bus``: observed producer price converted to the units
  of the model bus the item is produced onto. Crop buses carry dry matter,
  so crop prices are divided by (1 - moisture_fraction); the whole-commodity
  moisture inversion is used (NOT the edible-portion factor) because
  producer prices are quoted per tonne of marketed commodity. Animal
  product buses carry fresh weight, so animal prices pass through
  unchanged. Note the meat basis mismatch: FAOSTAT prices carcass-basis
  commodities ("with the bone") while the model's meat buses carry retail
  weight, so meat value is conservatively underestimated by roughly the
  retail-to-carcass factor. This cancels between baseline and realized
  value in the production-value floor, but shifts weight from livestock
  to crops in cross-commodity value composition.
- ``price_usd_per_t_bus_uniform``: a single world price per item (median
  of observed, non-fallback country prices on the bus basis), broadcast to
  all countries. Used for the uniform price basis of the production-value
  floor, which removes cross-country price distortions (e.g. support
  prices) from the value accounting.

Gap-filling follows the crop-cost conventions: per-item median across
configured countries, then a global (all-country) median for items absent
from the configured-country subset; filled rows are flagged
``is_fallback``. Per-item upper winsorization at the configured quantile
bounds greenhouse / data-quality outliers, flagged ``is_capped``.

Inputs
------
- PP.parquet : FAOSTAT Prices bulk
- faostat_crop_item_map.csv : model crop -> FAOSTAT item mapping
- faostat_animal_price_item_map.csv : model animal product -> FAOSTAT item mapping
- faostat_cost_proxies.yaml : proxy mappings for unmapped crops
- crop_moisture_content.csv : moisture fractions for the DM bus conversion
- M49-codes.csv, cpi_annual.csv

Output
------
- producer_prices.csv : columns (kind, item, country, price_usd_per_t,
  price_usd_per_t_bus, price_usd_per_t_bus_uniform, n_years, is_fallback,
  is_capped) where kind is "crop" or "animal" and item is the model crop /
  animal product name.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from workflow.scripts.faostat_bulk import (
    add_iso3_column,
    filter_bulk,
    get_item_map,
    load_bulk,
    load_m49_to_iso3,
)
from workflow.scripts.logging_config import setup_script_logging
from workflow.scripts.prepare_faostat_crop_costs import _deflate_to_base_year

logger = logging.getLogger(__name__)


def _compute_prices(
    pp_bulk: pd.DataFrame,
    *,
    items: list[str],
    item_to_codes: dict[str, list[int]],
    cpi_df: pd.DataFrame,
    base_year: int,
    years: list[int],
    price_element: int,
    iso3_codes: list[str] | None,
) -> pd.DataFrame:
    """Average deflated producer price per (item, country).

    Pass ``iso3_codes=None`` to compute across all countries (used to
    derive the global median fallback). Returns columns ``item``,
    ``country``, ``price_usd_per_t``, ``n_years``.
    """
    all_codes = sorted({code for codes in item_to_codes.values() for code in codes})
    pp_df = filter_bulk(
        pp_bulk,
        element_codes=[price_element],
        item_codes=all_codes,
        years=years,
        iso3_codes=iso3_codes,
    )
    pp_df = pp_df.dropna(subset=["Value"])
    pp_df["country"] = pp_df["iso3"].str.upper()
    pp_df["year"] = pd.to_numeric(pp_df["Year"], errors="coerce").astype(int)
    pp_df["value"] = pp_df["Value"].astype(float)

    results = []
    for item in items:
        codes = item_to_codes.get(item)
        if not codes:
            continue
        item_pp = pp_df[pp_df["Item Code"].isin(codes)]
        if item_pp.empty:
            continue
        # Average across mapped FAOSTAT items first (e.g. sheep + goat meat),
        # then deflate and average across years.
        agg = item_pp.groupby(["country", "year"])["value"].mean().reset_index()
        agg = _deflate_to_base_year(agg, cpi_df, base_year)
        agg = agg.dropna(subset=["value"])
        agg = agg[agg["value"] > 0]
        if agg.empty:
            continue
        avg = (
            agg.groupby("country")
            .agg(price_usd_per_t=("value", "mean"), n_years=("year", "nunique"))
            .reset_index()
        )
        avg["item"] = item
        results.append(avg)

    if not results:
        return pd.DataFrame(columns=["item", "country", "price_usd_per_t", "n_years"])
    return pd.concat(results, ignore_index=True)


def _map_items_to_codes(
    mapping: pd.DataFrame,
    item_col: str,
    pp_item_map: dict[str, int],
    proxies: dict | None = None,
) -> dict[str, list[int]]:
    """Resolve FAOSTAT item labels to PP item codes per model item.

    Model items without a direct mapping fall back to their proxy item's
    codes when ``proxies`` is given.
    """
    mapping = mapping.copy()
    mapping["item_code"] = mapping["faostat_item"].map(pp_item_map)
    missing = mapping[mapping["item_code"].isna()]["faostat_item"].unique()
    if len(missing) > 0:
        logger.warning(
            "FAOSTAT items missing from PP data: %s",
            ", ".join(str(m) for m in missing),
        )
    mapping = mapping.dropna(subset=["item_code"])
    mapping["item_code"] = mapping["item_code"].astype(int)
    item_to_codes = mapping.groupby(item_col)["item_code"].apply(list).to_dict()
    for item, proxy in (proxies or {}).items():
        if item not in item_to_codes and proxy in item_to_codes:
            item_to_codes[item] = item_to_codes[proxy]
    return item_to_codes


def main() -> None:
    pp_path = snakemake.input.pp_parquet
    crop_mapping_path = snakemake.input.crop_mapping
    animal_mapping_path = snakemake.input.animal_mapping
    m49_path = snakemake.input.m49_codes
    cpi_path = snakemake.input.cpi
    proxies_path = snakemake.input.proxies
    moisture_path = snakemake.input.moisture
    output_path = Path(snakemake.output[0])

    countries = [str(c).upper() for c in snakemake.params.countries]
    crops = list(snakemake.params.crops)
    animal_products = list(snakemake.params.animal_products)
    base_year = int(snakemake.params.currency_base_year)
    avg_start = int(snakemake.params.averaging_period["start_year"])
    avg_end = int(snakemake.params.averaging_period["end_year"])
    price_element = int(snakemake.params.price_element_code)
    cap_quantile_raw = snakemake.params.outlier_cap_quantile
    cap_quantile = float(cap_quantile_raw) if cap_quantile_raw is not None else None

    years = list(range(avg_start, avg_end + 1))

    crop_mapping = pd.read_csv(crop_mapping_path)
    crop_mapping = crop_mapping.loc[crop_mapping["faostat_item"].notna()]
    animal_mapping = pd.read_csv(animal_mapping_path)
    m49_to_iso3 = load_m49_to_iso3(m49_path)
    cpi_df = pd.read_csv(cpi_path)
    with open(proxies_path) as f:
        proxies = yaml.safe_load(f)

    logger.info("Loading FAOSTAT PP bulk data")
    pp_bulk = load_bulk(pp_path)
    pp_bulk = add_iso3_column(pp_bulk, m49_to_iso3)
    pp_item_map = get_item_map(pp_bulk)

    item_to_codes = _map_items_to_codes(
        crop_mapping.rename(columns={"crop": "item"}), "item", pp_item_map, proxies
    )
    item_to_codes.update(
        _map_items_to_codes(
            animal_mapping.rename(columns={"product": "item"}), "item", pp_item_map
        )
    )
    kind_by_item = dict.fromkeys(crops, "crop") | dict.fromkeys(
        animal_products, "animal"
    )
    all_items = crops + animal_products

    prices_df = _compute_prices(
        pp_bulk,
        items=all_items,
        item_to_codes=item_to_codes,
        cpi_df=cpi_df,
        base_year=base_year,
        years=years,
        price_element=price_element,
        iso3_codes=countries,
    )
    if prices_df.empty:
        raise RuntimeError("No valid FAOSTAT producer price data produced")
    prices_df["is_fallback"] = False
    logger.info(
        "Configured-country FAOSTAT producer prices: %d (item, country) rows",
        len(prices_df),
    )

    # Upper winsorization per item, mirroring the crop-cost pipeline
    # rationale: greenhouse / protected-cultivation producers report
    # producer prices far above field-production levels.
    prices_df["is_capped"] = False
    if cap_quantile is not None:
        caps = prices_df.groupby("item")["price_usd_per_t"].quantile(cap_quantile)
        for item, cap in caps.items():
            if not np.isfinite(cap):
                continue
            mask = (prices_df["item"] == item) & (prices_df["price_usd_per_t"] > cap)
            if not mask.any():
                continue
            prices_df.loc[mask, "price_usd_per_t"] = cap
            prices_df.loc[mask, "is_capped"] = True
        logger.info(
            "Outlier cap (q=%.2f): clipped %d of %d rows",
            cap_quantile,
            int(prices_df["is_capped"].sum()),
            len(prices_df),
        )

    # Per-item median fallback across configured countries; global median
    # for items with no configured-country data at all.
    item_medians = prices_df.groupby("item")["price_usd_per_t"].median()
    items_needing_global = [
        i for i in all_items if not np.isfinite(item_medians.get(i, np.nan))
    ]
    if items_needing_global:
        logger.info(
            "Computing global FAOSTAT median price for: %s",
            ", ".join(sorted(items_needing_global)),
        )
        global_df = _compute_prices(
            pp_bulk,
            items=items_needing_global,
            item_to_codes=item_to_codes,
            cpi_df=cpi_df,
            base_year=base_year,
            years=years,
            price_element=price_element,
            iso3_codes=None,
        )
        if not global_df.empty:
            item_medians = item_medians.combine_first(
                global_df.groupby("item")["price_usd_per_t"].median()
            )

    fallback_rows = []
    for item in all_items:
        existing = set(prices_df.loc[prices_df["item"] == item, "country"])
        missing_countries = sorted(set(countries) - existing)
        if not missing_countries:
            continue
        median_price = item_medians.get(item)
        if median_price is None or not np.isfinite(median_price):
            raise ValueError(
                f"No producer price available for item '{item}' (neither "
                "configured-country subset nor global FAOSTAT has data); "
                f"{len(missing_countries)} countries would inherit no price. "
                "Add a proxy mapping or check the PP parquet."
            )
        for country in missing_countries:
            fallback_rows.append(
                {
                    "item": item,
                    "country": country,
                    "price_usd_per_t": median_price,
                    "n_years": 0,
                    "is_fallback": True,
                    "is_capped": False,
                }
            )
    if fallback_rows:
        prices_df = pd.concat(
            [prices_df, pd.DataFrame(fallback_rows)], ignore_index=True
        )
    logger.info(
        "Producer prices: %d (item, country) pairs, %d fallbacks",
        len(prices_df),
        len(fallback_rows),
    )

    prices_df["kind"] = prices_df["item"].map(kind_by_item)

    # Convert to bus basis: crop buses carry dry matter, so divide the
    # fresh-commodity price by (1 - moisture); animal buses carry fresh
    # weight already.
    moisture = pd.read_csv(moisture_path, comment="#").set_index("crop")[
        "moisture_fraction"
    ]
    missing_moisture = sorted(set(crops) - set(moisture.index))
    if missing_moisture:
        raise ValueError(
            f"Missing crop_moisture_content.csv rows for: {missing_moisture}"
        )
    dm_factor = (
        prices_df["item"]
        .map(1.0 / (1.0 - moisture))
        .where(prices_df["kind"] == "crop", 1.0)
    )
    prices_df["price_usd_per_t_bus"] = prices_df["price_usd_per_t"] * dm_factor

    # Uniform world price per item: median of observed (non-fallback)
    # bus-basis prices; items known only through fallback keep it.
    observed = prices_df[~prices_df["is_fallback"]]
    uniform = observed.groupby("item")["price_usd_per_t_bus"].median()
    uniform = uniform.combine_first(
        prices_df.groupby("item")["price_usd_per_t_bus"].median()
    )
    prices_df["price_usd_per_t_bus_uniform"] = prices_df["item"].map(uniform)

    output_cols = [
        "kind",
        "item",
        "country",
        "price_usd_per_t",
        "price_usd_per_t_bus",
        "price_usd_per_t_bus_uniform",
        "n_years",
        "is_fallback",
        "is_capped",
    ]
    prices_df = (
        prices_df[output_cols]
        .sort_values(["kind", "item", "country"])
        .reset_index(drop=True)
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    prices_df.to_csv(output_path, index=False)
    logger.info("Wrote %d rows to %s", len(prices_df), output_path)


if __name__ == "__main__":
    logger = setup_script_logging(log_file=snakemake.log[0] if snakemake.log else None)
    main()
