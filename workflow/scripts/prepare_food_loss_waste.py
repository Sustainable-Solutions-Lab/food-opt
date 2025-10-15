#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Retrieve and process SDG 12.3.1 food loss and waste data from UNSD API.

Queries the UN Statistics Division API for:
- SDG 12.3.1(a): Food loss percentage by country and product type
- SDG 12.3.1(b): Food waste per capita by country and sector

Maps UN SDG food categories to internal model food groups, and converts
food waste from kg/person/year to fractions relative to food supply.

Input:
    - M49 codes (for regional mapping)
    - Countries list from config
    - Food groups list from config
    - Reference year from config (for FAOSTAT food supply data)

Output:
    - CSV with columns: country, food_group, loss_fraction, waste_fraction
"""

import logging
import sys

import faostat
import pandas as pd
import pycountry
import requests

logger = logging.getLogger(__name__)

FALLBACK_FOOD_SUPPLY: dict[str, list[str]] = {
    # Map territories or small countries to a proxy with similar dietary patterns.
    "ASM": ["USA"],  # American Samoa -> United States
    "BRN": ["MYS", "SGP"],  # Brunei -> Malaysia / Singapore
    "BTN": ["NPL", "IND"],  # Bhutan -> Nepal / India
    "ERI": ["ETH"],  # Eritrea -> Ethiopia
    "GNQ": ["GAB", "CMR"],  # Equatorial Guinea -> Gabon / Cameroon
    "GUF": ["GUY", "SUR"],  # French Guiana -> Guyana / Suriname
    "PRI": ["USA"],  # Puerto Rico -> United States
    "PSE": ["ISR", "JOR"],  # Palestine -> Israel / Jordan
    "SSD": ["SDN", "ETH"],  # South Sudan -> Sudan / Ethiopia
}

UN_TO_MODEL_FOOD_GROUPS: dict[str, list[str] | None] = {
    "CRL_PUL": ["grain", "whole_grains", "legumes"],
    "FRT_VGT": ["fruits", "vegetables"],
    "RT_TBR": ["starchy_vegetable", "oil", "nuts_seeds"],
    "ANMPROD": ["red_meat", "poultry", "dairy", "eggs"],
    "FSH_FSHPROD": ["fish"],
    "ALP": None,
}

MODEL_GROUP_TO_PRODUCT: dict[str, str] = {
    group: product
    for product, groups in UN_TO_MODEL_FOOD_GROUPS.items()
    if groups
    for group in groups
}


def load_m49_regions(m49_file: str) -> dict[str, dict]:
    """Load UN M49 region mappings.

    Args:
        m49_file: Path to M49 CSV file

    Returns:
        Dict mapping ISO3 code to region info (region_code, subregion_code, etc.)
    """
    df = pd.read_csv(m49_file, sep=";", encoding="utf-8-sig", comment="#")

    # Build mapping from ISO3 to region info
    mapping = {}
    for _, row in df.iterrows():
        iso3 = row["ISO-alpha3 Code"]
        if pd.notna(iso3):
            # Convert numeric codes to integers first to remove .0 suffix
            region_code = None
            if pd.notna(row["Region Code"]):
                region_code = str(int(float(row["Region Code"])))

            subregion_code = None
            if pd.notna(row["Sub-region Code"]):
                subregion_code = str(int(float(row["Sub-region Code"])))

            mapping[iso3] = {
                "m49_code": str(int(float(row["M49 Code"]))),
                "region_code": region_code,
                "region_name": row["Region Name"]
                if pd.notna(row["Region Name"])
                else None,
                "subregion_code": subregion_code,
                "subregion_name": row["Sub-region Name"]
                if pd.notna(row["Sub-region Name"])
                else None,
            }

    return mapping


def iso3_to_m49(iso3: str) -> str | None:
    """Convert ISO3 country code to M49 numeric code.

    Args:
        iso3: ISO 3166-1 alpha-3 country code (e.g., "USA")

    Returns:
        M49 numeric code as string, or None if not found
    """
    try:
        country = pycountry.countries.get(alpha_3=iso3)
        return country.numeric if country else None
    except (KeyError, AttributeError):
        return None


def query_unsd_series(series_code: str) -> pd.DataFrame:
    """Query UNSD SDG API for a given series code with pagination support.

    Args:
        series_code: UNSD series code (e.g., "AG_FLS_PCT")

    Returns:
        DataFrame with all observations for the series
    """
    url = "https://unstats.un.org/sdgapi/v1/sdg/Series/Data"
    all_data = []
    page = 1

    logger.info("Querying UNSD API: %s", series_code)

    while True:
        params = {"seriesCode": series_code, "page": page}

        try:
            response = requests.get(url, params=params, timeout=60)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            logger.error("Failed to query UNSD API for %s: %s", series_code, e)
            sys.exit(1)

        data = response.json()

        if "data" not in data or not data["data"]:
            break

        all_data.extend(data["data"])

        # Check if we have all pages
        total_pages = data.get("totalPages", 1)
        if page >= total_pages:
            break

        page += 1

    if not all_data:
        logger.error("No data returned for series %s", series_code)
        sys.exit(1)

    df = pd.DataFrame(all_data)
    logger.info("Retrieved %d observations for %s", len(df), series_code)

    return df


def map_un_food_category_to_model_groups(un_category_code: str) -> list[str]:
    """Map UN SDG food categories to internal model food groups.

    UN Categories (Type of product dimension):
    - CRL_PUL: Cereals and pulses
    - FRT_VGT: Fruits and vegetables
    - RT_TBR: Roots, tubers and oil-bearing crops
    - ANMPROD: Animal products
    - FSH_FSHPROD: Fish and fish products
    - ALP: Total or no breakdown

    Args:
        un_category_code: UN SDG product type code

    Returns:
        List of internal food group names
    """
    groups = UN_TO_MODEL_FOOD_GROUPS.get(un_category_code)
    if groups is None:
        return []
    return list(groups)


def process_food_loss_data(
    df: pd.DataFrame,
    countries: list[str],
    food_groups: list[str],
    m49_regions: dict[str, dict],
) -> pd.DataFrame:
    """Process food loss percentage data using regional aggregates.

    Food loss data is only available at regional/sub-regional level, not country level.
    We map each country to its UN M49 sub-region and use the regional average.

    Args:
        df: Raw UNSD data for AG_FLS_PCT series
        countries: List of ISO3 country codes
        food_groups: List of internal food group names
        m49_regions: Mapping of ISO3 to region info

    Returns:
        DataFrame with columns: country, food_group, loss_fraction, year
    """
    # Extract product type from dimensions
    df["product_type"] = df["dimensions"].apply(
        lambda x: x.get("Type of product") if isinstance(x, dict) else None
    )

    # Parse year and value
    df["year"] = pd.to_numeric(df["timePeriodStart"], errors="coerce")
    df["value_numeric"] = pd.to_numeric(df["value"], errors="coerce")

    # Build lookup of available regional data: {region_code: {product_type: {year: value}}}
    regional_data = {}
    for _, row in df.iterrows():
        region_code = str(row["geoAreaCode"])
        product_type = row["product_type"]
        year = row["year"]
        value = row["value_numeric"]

        if pd.isna(value) or pd.isna(year):
            continue

        if region_code not in regional_data:
            regional_data[region_code] = {}
        if product_type not in regional_data[region_code]:
            regional_data[region_code][product_type] = {}

        regional_data[region_code][product_type][year] = value

    logger.info("Found food loss data for %d regions", len(regional_data))

    product_type_summary = {
        region: {ptype for ptype in types if ptype}
        for region, types in regional_data.items()
    }
    available_types = {
        ptype for types in product_type_summary.values() for ptype in types if ptype
    }
    if available_types:
        logger.info(
            "Food loss product types reported: %s",
            ", ".join(sorted(available_types)),
        )
    else:
        logger.warning("No product-specific food loss breakdown found in UNSD response")

    regions_with_breakdown = sum(
        1
        for types in product_type_summary.values()
        if any(ptype and ptype != "ALP" for ptype in types)
    )
    total_regions = len(product_type_summary)
    logger.info(
        "Regions with product-specific breakdown: %d/%d",
        regions_with_breakdown,
        total_regions,
    )

    # Derive world-level product correction factors for disaggregation
    global_shares: dict[str, float] = {}
    world_alp_value: float | None = None
    world_data = df[df["geoAreaCode"].astype(str) == "1"]
    if not world_data.empty:
        latest_world_year = world_data["year"].max()
        world_latest = world_data[world_data["year"] == latest_world_year]
        product_values: dict[str, float] = {}
        for _, row in world_latest.iterrows():
            ptype = row["product_type"]
            value = row["value_numeric"]
            if ptype == "ALP" and pd.notna(value):
                world_alp_value = value
            if pd.isna(value) or ptype is None or ptype == "ALP":
                continue
            product_values[ptype] = value

        if world_alp_value and world_alp_value > 0:
            global_shares = {
                ptype: v / world_alp_value
                for ptype, v in product_values.items()
                if v > 0 and ptype
            }
            logger.info(
                "Using world food loss shares (%d) for fallback disaggregation: %s",
                int(latest_world_year),
                ", ".join(
                    f"{ptype}={share:.2f}" for ptype, share in global_shares.items()
                ),
            )
    if not global_shares:
        logger.warning(
            "World food loss shares unavailable; ALP totals will remain un-disaggregated"
        )

    global_group_corrections: dict[str, float] = {}
    if global_shares:
        for group, product in MODEL_GROUP_TO_PRODUCT.items():
            ratio = global_shares.get(product)
            if ratio is not None:
                global_group_corrections[group] = ratio
    for group in food_groups:
        global_group_corrections.setdefault(group, 1.0)

    if global_group_corrections:
        logger.info(
            "Applying global loss correction factors per food group: %s",
            ", ".join(
                f"{group}={factor:.2f}"
                for group, factor in sorted(global_group_corrections.items())
            ),
        )

    results = []

    for country_code in countries:
        # Get region info for this country
        region_info = m49_regions.get(country_code)
        if not region_info:
            logger.warning("No M49 region info for %s", country_code)
            continue

        # Try sub-region first, then region
        region_code = None
        if (
            region_info["subregion_code"]
            and region_info["subregion_code"] in regional_data
        ):
            region_code = region_info["subregion_code"]
        elif region_info["region_code"] and region_info["region_code"] in regional_data:
            region_code = region_info["region_code"]

        if not region_code:
            # No regional data available for this country
            continue
        region_entries = regional_data[region_code]

        alp_year_data = region_entries.get("ALP")
        if not alp_year_data:
            continue

        latest_year = max(alp_year_data.keys())
        loss_pct = alp_year_data[latest_year]
        if pd.isna(loss_pct):
            continue

        base_loss_fraction = loss_pct / 100.0

        for food_group in food_groups:
            correction = global_group_corrections.get(food_group, 1.0)
            loss_fraction = base_loss_fraction * correction
            results.append(
                {
                    "country": country_code,
                    "food_group": food_group,
                    "loss_fraction": loss_fraction,
                    "year": int(latest_year),
                }
            )

    return pd.DataFrame(results)


def fetch_faostat_food_supply(
    countries: list[str], reference_year: int
) -> pd.DataFrame:
    """Fetch total food supply per capita from FAOSTAT Food Balance Sheets.

    Args:
        countries: List of ISO3 country codes
        reference_year: Year for which to retrieve data

    Returns:
        DataFrame with columns: country (ISO3), food_supply_g_day
    """
    logger.info("Fetching FAOSTAT food supply data for %d", reference_year)

    # Use FBS (Food Balance Sheets) dataset
    # Element: "Food supply quantity (kg/capita/yr)"
    DATASET = "FBS"
    ELEMENT_LABEL = "Food supply quantity (kg/capita/yr)"

    # Get parameter mappings
    elem_map = faostat.get_par(DATASET, "Element")

    if ELEMENT_LABEL not in elem_map:
        raise RuntimeError(f"Element '{ELEMENT_LABEL}' not found in FAOSTAT FBS")
    element_code = elem_map[ELEMENT_LABEL]

    # Query FAOSTAT for the requested year.
    # The API does not expose the aggregate "Grand Total" alongside this element,
    # so we pull the per-item values and aggregate locally.
    logger.info("Querying FAOSTAT for total food supply")
    target_countries = [code.upper() for code in countries]
    pars = {
        "element": element_code,
        "year": str(reference_year),
    }

    try:
        df = faostat.get_data_df(
            DATASET, pars=pars, coding={"area": "ISO3"}, strval=True
        )
    except Exception as exc:  # pragma: no cover - network/runtime errors
        logger.error("Failed to fetch FAOSTAT data for %d: %s", reference_year, exc)
        return pd.DataFrame(columns=["country", "food_supply_g_day"])

    if df.empty:
        logger.warning("FAOSTAT returned no food supply data for requested selection")
        return pd.DataFrame(columns=["country", "food_supply_g_day"])

    # Process results
    df["Value"] = pd.to_numeric(df["Value"], errors="coerce")
    df = df.dropna(subset=["Value"])
    if df.empty:
        logger.warning(
            "FAOSTAT food supply data contains no numeric values after parsing"
        )
        return pd.DataFrame(columns=["country", "food_supply_g_day"])

    # Identify the ISO3 column (returned via coding system)
    iso_col = None
    for candidate in df.columns:
        normalized = candidate.strip().lower().replace(" ", "_")
        if normalized in {"area_code", "area_code_(iso3)", "area_code_iso3"}:
            iso_col = candidate
            break
    if iso_col is None and "Area Code" in df.columns:
        iso_col = "Area Code"

    if iso_col is None:
        logger.warning("Could not identify ISO3 area column in FAOSTAT response")
        return pd.DataFrame(columns=["country", "food_supply_g_day"])

    df["country"] = df[iso_col].astype(str).str.upper()
    df = df[df["country"].isin(target_countries)]
    df = df[df["Unit"].str.lower() == "kg/cap"]

    if df.empty:
        logger.warning(
            "FAOSTAT returned no records for the requested ISO3 country list"
        )
        return pd.DataFrame(columns=["country", "food_supply_g_day"])

    # Convert from kg/year to g/day: kg/yr * 1000 / 365. Sum across items to obtain total supply.
    total_supply = (
        df.groupby("country", as_index=False)["Value"]
        .sum()
        .rename(columns={"Value": "food_supply_kg_year"})
    )
    total_supply["food_supply_g_day"] = (
        total_supply["food_supply_kg_year"] * 1000.0
    ) / 365.0
    total_supply = total_supply.drop(columns=["food_supply_kg_year"])

    missing_countries = sorted(set(target_countries) - set(total_supply["country"]))
    fallback_rows: list[dict[str, float]] = []
    applied_fallbacks: list[tuple[str, str]] = []

    if missing_countries:
        supply_lookup = dict(
            zip(total_supply["country"], total_supply["food_supply_g_day"])
        )
        for iso_code in missing_countries:
            for candidate in FALLBACK_FOOD_SUPPLY.get(iso_code, []):
                value = supply_lookup.get(candidate.upper())
                if value is not None:
                    fallback_rows.append(
                        {"country": iso_code, "food_supply_g_day": value}
                    )
                    supply_lookup[iso_code] = value
                    applied_fallbacks.append((iso_code, candidate.upper()))
                    break

    if fallback_rows:
        total_supply = pd.concat(
            [total_supply, pd.DataFrame(fallback_rows)],
            ignore_index=True,
        )
        logger.info(
            "Applied food supply fallbacks: %s",
            "; ".join(f"{iso}->{proxy}" for iso, proxy in applied_fallbacks),
        )

    remaining_missing = sorted(set(target_countries) - set(total_supply["country"]))
    if remaining_missing:
        logger.warning(
            "Missing FAOSTAT food supply data for: %s",
            ", ".join(remaining_missing),
        )

    logger.info("Retrieved food supply for %d countries", len(total_supply))
    logger.info(
        "Mean food supply: %.1f g/day",
        total_supply["food_supply_g_day"].mean(),
    )

    return total_supply


def process_food_waste_data(
    df: pd.DataFrame,
    food_supply: pd.DataFrame,
    countries: list[str],
    food_groups: list[str],
    m49_regions: dict[str, dict],
    reference_year: int,
) -> pd.DataFrame:
    """Process food waste per capita data and convert to fractions.

    Args:
        df: Raw UNSD data for AG_FOOD_WST_PC series
        food_supply: FAOSTAT food supply per capita (g/day)
        countries: List of ISO3 country codes
        food_groups: List of internal food group names
        m49_regions: Mapping of ISO3 to region metadata (for fallbacks)
        reference_year: Reference model year for fallback rows

    Returns:
        DataFrame with columns: country, food_group, waste_fraction, year
    """
    # Note: Data is classified as "Global" or "Estimated" reporting type but contains country-level estimates
    # We'll match by geoAreaCode instead of filtering by Reporting Type

    # Extract sector from dimensions
    df["sector"] = df["dimensions"].apply(
        lambda x: x.get("Food Waste Sector") if isinstance(x, dict) else None
    )

    # Filter to "ALL" sector totals
    df_total = df[df["sector"] == "ALL"].copy()

    df_total["year"] = pd.to_numeric(df_total["timePeriodStart"], errors="coerce")
    df_total["waste_kg_year"] = pd.to_numeric(df_total["value"], errors="coerce")

    results = []

    for country_code in countries:
        # Convert ISO3 to M49 numeric code
        m49_code = iso3_to_m49(country_code)
        if not m49_code:
            continue

        # Get waste data for this country using M49 code
        country_waste = df_total[df_total["geoAreaCode"] == m49_code]

        if country_waste.empty:
            continue

        # Get latest year
        latest = country_waste.sort_values("year").iloc[-1]

        if pd.isna(latest["waste_kg_year"]):
            continue

        waste_kg_year = latest["waste_kg_year"]
        year = int(latest["year"])

        # Get food supply for this country
        country_supply = food_supply[food_supply["country"] == country_code]

        if country_supply.empty:
            logger.warning(
                "No food supply data for %s, skipping waste calculation",
                country_code,
            )
            continue

        food_supply_g_day = country_supply["food_supply_g_day"].iloc[0]

        # Convert waste kg/year to g/day, then to fraction
        # waste_fraction = waste / food_supply
        waste_g_day = (waste_kg_year / 365.0) * 1000.0
        waste_fraction = waste_g_day / food_supply_g_day if food_supply_g_day > 0 else 0

        # Apply to all food groups (waste is not broken down by food type)
        for food_group in food_groups:
            results.append(
                {
                    "country": country_code,
                    "food_group": food_group,
                    "waste_fraction": waste_fraction,
                    "year": year,
                }
            )

    waste_df = pd.DataFrame(results)

    if waste_df.empty:
        logger.warning(
            "No country-level food waste data available; fallbacks not applied"
        )
        return waste_df

    # Attach regional metadata for fallback calculations
    m49_meta = (
        pd.DataFrame.from_dict(m49_regions, orient="index")[
            ["subregion_code", "region_code"]
        ]
        if m49_regions
        else pd.DataFrame()
    )
    if not m49_meta.empty:
        waste_df = waste_df.merge(
            m49_meta,
            left_on="country",
            right_index=True,
            how="left",
        )
    else:
        waste_df["subregion_code"] = None
        waste_df["region_code"] = None

    subregion_avg = (
        waste_df.dropna(subset=["subregion_code"])
        .groupby(["subregion_code", "food_group"])["waste_fraction"]
        .mean()
    )
    region_avg = (
        waste_df.dropna(subset=["region_code"])
        .groupby(["region_code", "food_group"])["waste_fraction"]
        .mean()
    )
    global_avg = waste_df.groupby("food_group")["waste_fraction"].mean()

    year_candidates = pd.to_numeric(waste_df["year"], errors="coerce").dropna()
    fallback_year = (
        int(year_candidates.median()) if not year_candidates.empty else reference_year
    )

    countries_with_data = set(waste_df["country"].unique())
    missing_countries = [iso for iso in countries if iso not in countries_with_data]

    fallback_rows: list[dict] = []
    fallback_sources: list[str] = []

    for iso_code in missing_countries:
        region_info = m49_regions.get(iso_code, {})
        fallback_series = None
        source_label = None

        subregion_code = region_info.get("subregion_code")
        if subregion_code:
            try:
                fallback_series = subregion_avg.xs(subregion_code)
                source_label = f"subregion {subregion_code}"
            except KeyError:
                fallback_series = None

        if fallback_series is None:
            region_code = region_info.get("region_code")
            if region_code:
                try:
                    fallback_series = region_avg.xs(region_code)
                    source_label = f"region {region_code}"
                except KeyError:
                    fallback_series = None

        if fallback_series is None and not global_avg.empty:
            fallback_series = global_avg
            source_label = "global average"

        if fallback_series is None:
            logger.warning("Unable to determine fallback food waste for %s", iso_code)
            continue

        fallback_values = fallback_series.reindex(food_groups)
        fallback_values = fallback_values.fillna(global_avg.reindex(food_groups))

        if fallback_values.isna().all():
            logger.warning("Fallback food waste values remain NaN for %s", iso_code)
            continue

        for food_group, waste_fraction in fallback_values.items():
            if pd.isna(waste_fraction):
                continue
            fallback_rows.append(
                {
                    "country": iso_code,
                    "food_group": food_group,
                    "waste_fraction": float(waste_fraction),
                    "year": fallback_year,
                    "subregion_code": region_info.get("subregion_code"),
                    "region_code": region_info.get("region_code"),
                }
            )
        fallback_sources.append(f"{iso_code}->{source_label}")

    if fallback_rows:
        waste_df = pd.concat([waste_df, pd.DataFrame(fallback_rows)], ignore_index=True)
        logger.info(
            "Applied food waste fallbacks for %d countries: %s",
            len(fallback_sources),
            ", ".join(fallback_sources),
        )

    return waste_df.drop(columns=["subregion_code", "region_code"], errors="ignore")


def main():
    m49_file = snakemake.input["m49"]
    output_file = snakemake.output["food_loss_waste"]
    countries = snakemake.params["countries"]
    food_groups = snakemake.params["food_groups"]
    reference_year = snakemake.params["health_reference_year"]

    logger.info("Processing food loss and waste data")
    logger.info("Countries: %d", len(countries))
    logger.info("Food groups: %s", food_groups)
    logger.info("Reference year: %d", reference_year)

    # Load M49 region mappings
    m49_regions = load_m49_regions(m49_file)
    logger.info("Loaded M49 regions for %d countries", len(m49_regions))

    # Fetch FAOSTAT food supply data
    food_supply = fetch_faostat_food_supply(countries, reference_year)

    # Query UNSD API
    loss_data = query_unsd_series("AG_FLS_PCT")
    waste_data = query_unsd_series("AG_FOOD_WST_PC")

    # Process food loss
    loss_df = process_food_loss_data(loss_data, countries, food_groups, m49_regions)
    logger.info("Processed %d food loss observations", len(loss_df))

    # Process food waste
    waste_df = process_food_waste_data(
        waste_data,
        food_supply,
        countries,
        food_groups,
        m49_regions,
        reference_year,
    )
    logger.info("Processed %d food waste observations", len(waste_df))

    # Merge loss and waste data
    if not loss_df.empty and not waste_df.empty:
        result = pd.merge(
            loss_df[["country", "food_group", "loss_fraction"]],
            waste_df[["country", "food_group", "waste_fraction"]],
            on=["country", "food_group"],
            how="outer",
        )
    elif not loss_df.empty:
        result = loss_df[["country", "food_group", "loss_fraction"]].copy()
        result["waste_fraction"] = None
    elif not waste_df.empty:
        result = waste_df[["country", "food_group", "waste_fraction"]].copy()
        result["loss_fraction"] = None
    else:
        logger.error("No food loss or waste data retrieved")
        sys.exit(1)

    # Fill missing values with 0 (no data = assume no loss/waste)
    result = result.infer_objects(copy=False)
    result["loss_fraction"] = result["loss_fraction"].fillna(0.0)
    result["waste_fraction"] = result["waste_fraction"].fillna(0.0)

    # Sort for readability
    result = result.sort_values(["country", "food_group"]).reset_index(drop=True)

    logger.info("Final output: %d rows", len(result))
    logger.info("Countries with data: %d", result["country"].nunique())
    logger.info("Mean loss fraction: %.3f", result["loss_fraction"].mean())
    logger.info("Mean waste fraction: %.3f", result["waste_fraction"].mean())

    # Write output
    result.to_csv(output_file, index=False)
    logger.info("Wrote output to %s", output_file)


if __name__ == "__main__":
    main()
