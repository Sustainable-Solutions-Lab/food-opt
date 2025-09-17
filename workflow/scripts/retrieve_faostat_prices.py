# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Retrieve FAOSTAT producer prices (USD/tonne) for configured crops and
compute a 10-year global average (2015–2024) using the `faostat` Python package.

Inputs
- snakemake.params.crops: list of crop names from config
- snakemake.input.mapping: CSV mapping of model crop -> FAOSTAT item label

Output
- snakemake.output.prices: CSV with columns:
    crop,faostat_item,n_obs,price_usd_per_tonne

Notes
- Strictly uses the `faostat` package to obtain ready-to-use DataFrames.
- If any crop lacks data, the script writes NaN for its price and logs a warning.
"""

import logging
from typing import Iterable

import pandas as pd
import faostat

logger = logging.getLogger(__name__)


YEARS = list(range(2015, 2025))  # inclusive [2015, 2024]
DATASET = "PP"  # Producer Prices
ELEMENT_LABEL = "Producer Price (USD/tonne)"


def _unique_preserve(seq: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for s in seq:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def _fetch_with_faostat(items: list[str]) -> pd.DataFrame:
    """Retrieve data using the `faostat` package only.

    Returns a DataFrame with columns: Item, Element, Year, Area, Value, Unit
    """
    # Convert labels to parameter codes
    elem_map = faostat.get_par(DATASET, "Element")  # label -> code
    item_map = faostat.get_par(DATASET, "Item")  # label -> code
    year_codes = [str(y) for y in YEARS]

    if ELEMENT_LABEL not in elem_map:
        raise RuntimeError(f"Element '{ELEMENT_LABEL}' not found in FAOSTAT")
    element_code = elem_map[ELEMENT_LABEL]

    missing = [lab for lab in items if lab not in item_map]
    if missing:
        logger.warning(
            "Items not found in FAOSTAT and will be skipped: %s", ", ".join(missing)
        )
    item_codes = [item_map[lab] for lab in items if lab in item_map]
    if not item_codes:
        return pd.DataFrame(
            columns=["Item", "Element", "Year", "Area", "Value", "Unit"]
        )

    pars = {
        "element": element_code,
        "item": item_codes,
        "year": year_codes,
        # omit area to include all areas
    }
    df = faostat.get_data_df(DATASET, pars=pars, strval=True)

    # Keep the essential columns and return
    cols = ["Item", "Element", "Year", "Area", "Value", "Unit"]
    # Some datasets may have different capitalization; be lenient
    for c in cols:
        if c not in df.columns:
            raise KeyError(
                f"Expected column '{c}' not in FAOSTAT response columns {list(df.columns)}"
            )
    return df[cols].copy()


# No REST fallback by design: we rely on the faostat client per requirements.


def fetch_prices(items: list[str]) -> pd.DataFrame:
    """Fetch FAOSTAT USD/tonne producer prices for given item labels.

    Returns tidy DataFrame with columns: Item, Year, Area, Value (USD/tonne)
    """
    # Retrieve using faostat package only
    df = _fetch_with_faostat(items)
    logger.info("Retrieved data via faostat package (%d rows)", len(df))

    # Filter strictly to the element and years of interest
    df = df[df["Element"] == ELEMENT_LABEL]
    df = df[df["Year"].astype(int).isin(YEARS)]
    # Coerce numeric
    df["Value"] = pd.to_numeric(df["Value"], errors="coerce")
    df = df.dropna(subset=["Value"]).copy()
    return df


def main():
    crops: list[str] = list(snakemake.params.crops)  # type: ignore[name-defined]
    mapping_path: str = str(snakemake.input.mapping)  # type: ignore[name-defined]
    out_path: str = str(snakemake.output.prices)  # type: ignore[name-defined]

    map_df = pd.read_csv(mapping_path)
    if not {"crop", "faostat_item"}.issubset(map_df.columns):
        raise ValueError("Mapping file must have columns: crop, faostat_item")

    # Keep only crops of interest, preserve input order
    map_df = map_df.set_index("crop").reindex(crops)

    missing = map_df[map_df["faostat_item"].isna()].index.tolist()
    if missing:
        logger.warning("Missing FAOSTAT item mapping for crops: %s", ", ".join(missing))

    items = _unique_preserve([i for i in map_df["faostat_item"].dropna().tolist()])
    if not items:
        # Nothing to fetch; produce empty with NaNs
        out = pd.DataFrame(
            {
                "crop": crops,
                "faostat_item": [pd.NA] * len(crops),
                "n_obs": [0] * len(crops),
                "price_usd_per_tonne": [pd.NA] * len(crops),
            }
        )
        out.to_csv(out_path, index=False)
        return

    df = fetch_prices(items)

    # Compute unweighted global average across countries and years
    grp = df.groupby("Item", as_index=False).agg(
        n_obs=("Value", "count"), price_usd_per_tonne=("Value", "mean")
    )

    # Map back to crops
    item_to_price = grp.set_index("Item").to_dict(orient="index")
    rows = []
    for crop in crops:
        item = map_df.at[crop, "faostat_item"] if crop in map_df.index else pd.NA
        if pd.isna(item):
            rows.append(
                {
                    "crop": crop,
                    "faostat_item": pd.NA,
                    "n_obs": 0,
                    "price_usd_per_tonne": pd.NA,
                }
            )
            continue
        rec = item_to_price.get(item)
        if rec is None:
            logger.warning(
                "No FAOSTAT price found for item '%s' (crop '%s')", item, crop
            )
            rows.append(
                {
                    "crop": crop,
                    "faostat_item": item,
                    "n_obs": 0,
                    "price_usd_per_tonne": pd.NA,
                }
            )
        else:
            rows.append(
                {
                    "crop": crop,
                    "faostat_item": item,
                    "n_obs": int(rec["n_obs"]),
                    "price_usd_per_tonne": float(rec["price_usd_per_tonne"]),
                }
            )

    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_path, index=False)


if __name__ == "__main__":
    main()
