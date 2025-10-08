#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Process Global Dietary Database (GDD) country-level dietary intake data.

Reads multiple v*_cnty.csv files from GDD, extracts national-level baseline
intake estimates (aggregated across age/sex/urban/education strata), maps
GDD food variables to model food groups, and outputs a consolidated
baseline diet file.

Input:
    - GDD directory with Country-level estimates/*.csv files
    - Reference year from config
    - Food groups from config

Output:
    - CSV with columns: scenario,unit,item,country,year,value
"""

import sys
from pathlib import Path

import pandas as pd


def main():
    gdd_dir = Path(snakemake.input["gdd_dir"])
    reference_year = snakemake.params["reference_year"]
    food_groups = snakemake.params["food_groups"]
    output_file = snakemake.output["diet"]

    # Map GDD variable codes (vXX) to model food groups
    # Based on GDD codebook and canonical food_groups.csv
    gdd_to_model_items = {
        "v01": "fruits",  # Fruits (all types)
        "v02": "vegetables",  # Non-starchy vegetables
        "v03": "starchy_vegetable",  # Potatoes
        "v04": "starchy_vegetable",  # Other starchy vegetables (yam, cassava, etc.)
        "v05": "legumes",  # Beans and legumes
        "v06": "nuts_seeds",  # Nuts and seeds
        "v07": "grain",  # Refined grains (white flour, white rice)
        "v08": "whole_grains",  # Whole grains
        "v09": "prc_meat",  # Total processed meats
        "v10": "red_meat",  # Unprocessed red meats (cattle, pig)
        "v11": "fish",  # Total seafoods (fish + shellfish)
        "v12": "eggs",  # Eggs
        "v57": "dairy",  # Total Milk (includes milk equivalents from all dairy)
        "v15": None,  # Sugar-sweetened beverages (not tracked as food group)
        "v16": "fruits",  # Fruit juices (count as fruits)
        "v17": None,  # Coffee (not tracked as food group)
        "v18": None,  # Tea (not tracked as food group)
        # Note: We use v57 "Total Milk" for dairy, which aligns with the GBD dairy
        # risk factor definition and includes milk equivalents from all dairy products.
        # Individual components (v13 cheese, v14 yogurt) are not used separately
        # to avoid double-counting.
    }

    # Filter to only food groups that are in the config
    # Multiple GDD variables may map to the same food group
    requested_food_groups = set(food_groups)
    food_group_vars = {}
    for varcode, item in gdd_to_model_items.items():
        if item is not None and item in requested_food_groups:
            if item not in food_group_vars:
                food_group_vars[item] = []
            food_group_vars[item].append(varcode)

    print(f"[prepare_gdd_dietary_intake] Processing GDD data for year {reference_year}")
    print(f"[prepare_gdd_dietary_intake] Food groups: {sorted(food_group_vars.keys())}")
    for item, varcodes in sorted(food_group_vars.items()):
        print(f"  {item}: {varcodes}")

    country_estimates_dir = gdd_dir / "Country-level estimates"
    if not country_estimates_dir.exists():
        print(
            f"ERROR: Country-level estimates directory not found: {country_estimates_dir}",
            file=sys.stderr,
        )
        sys.exit(1)

    all_data = []

    # Process each food group (which may aggregate multiple GDD variables)
    for model_item, varcodes in food_group_vars.items():
        item_data = []

        for varcode in varcodes:
            csv_file = country_estimates_dir / f"{varcode}_cnty.csv"
            if not csv_file.exists():
                print(f"WARNING: File not found: {csv_file}", file=sys.stderr)
                continue

            print(
                f"[prepare_gdd_dietary_intake] Reading {csv_file.name} ({model_item})..."
            )
            df = pd.read_csv(csv_file)

            # Filter to reference year
            df_year = df[df["year"] == reference_year].copy()

            if df_year.empty:
                print(
                    f"WARNING: No data for year {reference_year} in {csv_file.name}. Trying nearest year...",
                    file=sys.stderr,
                )
                # Find nearest year
                available_years = sorted(df["year"].unique())
                nearest_year = min(
                    available_years, key=lambda y: abs(y - reference_year)
                )
                print(f"  Using nearest year: {nearest_year}", file=sys.stderr)
                df_year = df[df["year"] == nearest_year].copy()
                if df_year.empty:
                    print(f"ERROR: Still no data for {csv_file.name}", file=sys.stderr)
                    continue

            # Aggregate to country level (weighted mean across strata)
            # GDD stratifies by age, sex, urban, education
            # We want population-weighted national averages
            # The 'median' column is the mean intake (50th percentile of modeled simulations)

            # For national-level estimates, we want rows where strata are "999" (all ages/all)
            # Check if such aggregated rows exist
            if "age" in df_year.columns:
                natl = df_year[df_year["age"] == 999].copy()
                if natl.empty:
                    # If no pre-aggregated national data, compute weighted average across strata
                    # This is complex without population weights, so we'll take simple mean
                    print(
                        f"  No pre-aggregated national data (age==999) for {csv_file.name}. "
                        "Computing unweighted mean across strata.",
                        file=sys.stderr,
                    )
                    natl = (
                        df_year.groupby("iso3")["median"]
                        .mean()
                        .reset_index()
                        .rename(columns={"median": "value"})
                    )
                else:
                    # Use pre-aggregated national data
                    # Some files have multiple national rows (e.g., by sex/urban), so still average
                    natl = (
                        natl.groupby("iso3")["median"]
                        .mean()
                        .reset_index()
                        .rename(columns={"median": "value"})
                    )
            else:
                # Simpler structure: already aggregated
                natl = (
                    df_year.groupby("iso3")["median"]
                    .mean()
                    .reset_index()
                    .rename(columns={"median": "value"})
                )

            natl["varcode"] = varcode
            item_data.append(natl)

        # Aggregate multiple GDD variables into the same food group
        # For example, dairy = cheese + yogurt, fruits = fruits + fruit juices, etc.
        if item_data:
            if len(item_data) == 1:
                # Single variable maps to this food group
                combined = item_data[0].copy()
                combined["item"] = model_item
                combined = combined[["iso3", "value", "item"]]
            else:
                # Multiple variables map to this food group - sum them
                print(
                    f"[prepare_gdd_dietary_intake] Aggregating {len(item_data)} GDD variables for {model_item}"
                )
                combined = pd.concat(item_data, ignore_index=True)
                combined = combined.groupby("iso3")["value"].sum().reset_index()
                combined["item"] = model_item

            combined["year"] = reference_year
            all_data.append(combined)

    if not all_data:
        print("ERROR: No data collected from GDD files", file=sys.stderr)
        sys.exit(1)

    # Concatenate all food groups
    result = pd.concat(all_data, ignore_index=True)

    # Add scenario and unit columns
    result["scenario"] = "BMK"  # Baseline scenario
    result["unit"] = "g/d_w"  # grams per day per person (weighted)
    result = result.rename(columns={"iso3": "country"})

    # Reorder columns to standard format
    result = result[["scenario", "unit", "item", "country", "year", "value"]]

    # Sort by country and item for readability
    result = result.sort_values(["country", "item"]).reset_index(drop=True)

    print(
        f"[prepare_gdd_dietary_intake] Processed {len(result)} country-item pairs "
        f"for {result['country'].nunique()} countries"
    )
    print(
        f"[prepare_gdd_dietary_intake] Food groups: {sorted(result['item'].unique())}"
    )

    # Fill in missing countries using proxy data from similar countries
    # This is for territories/dependencies that don't have separate GDD data
    COUNTRY_PROXIES = {
        "ASM": "WSM",  # American Samoa -> Samoa
        "GUF": "FRA",  # French Guiana -> France
        "PRI": "USA",  # Puerto Rico -> USA
        "SOM": "ETH",  # Somalia -> Ethiopia (similar region, data available)
    }

    required_countries = set(snakemake.params["countries"])
    requested_food_groups = set(snakemake.params["food_groups"])
    output_countries = set(result["country"].unique())
    output_food_groups = set(result["item"].unique())

    # Only validate food groups that we actually tried to extract from GDD
    # Some food groups (like oil, poultry) are not tracked in dietary surveys
    expected_food_groups = set(food_group_vars.keys())

    missing_countries = required_countries - output_countries
    if missing_countries:
        filled = []
        still_missing = []
        for missing in sorted(missing_countries):
            if missing in COUNTRY_PROXIES:
                proxy = COUNTRY_PROXIES[missing]
                if proxy in output_countries:
                    # Duplicate proxy country's data for the missing country
                    proxy_data = result[result["country"] == proxy].copy()
                    proxy_data["country"] = missing
                    result = pd.concat([result, proxy_data], ignore_index=True)
                    filled.append(f"{missing} (using {proxy} data)")
                else:
                    still_missing.append(missing)
            else:
                still_missing.append(missing)

        if filled:
            print(
                f"[prepare_gdd_dietary_intake] Filled {len(filled)} missing countries using proxies:"
            )
            for entry in filled:
                print(f"  - {entry}")

        # Update missing list after filling
        output_countries = set(result["country"].unique())
        missing_countries = required_countries - output_countries

    # Validate that we have all required countries and food groups
    if missing_countries:
        raise ValueError(
            f"[prepare_gdd_dietary_intake] ERROR: GDD dietary data is missing {len(missing_countries)} required countries: "
            f"{sorted(list(missing_countries))[:20]}{'...' if len(missing_countries) > 20 else ''}. "
            f"Please ensure the GDD download includes all countries listed in config."
        )

    missing_food_groups = expected_food_groups - output_food_groups
    if missing_food_groups:
        raise ValueError(
            f"[prepare_gdd_dietary_intake] ERROR: GDD dietary data is missing {len(missing_food_groups)} expected food groups: "
            f"{sorted(missing_food_groups)}. Available: {sorted(output_food_groups)}. "
            f"Please ensure the GDD download includes all necessary data."
        )

    # Report food groups not available in GDD
    unavailable_in_gdd = requested_food_groups - expected_food_groups
    if unavailable_in_gdd:
        print(
            f"[prepare_gdd_dietary_intake] Note: {len(unavailable_in_gdd)} food groups not available in GDD: "
            f"{sorted(unavailable_in_gdd)}"
        )

    print(
        f"[prepare_gdd_dietary_intake] ✓ Validation passed: all required countries and "
        f"{len(expected_food_groups)} GDD-tracked food groups present"
    )

    # Write output
    result.to_csv(output_file, index=False)
    print(f"[prepare_gdd_dietary_intake] Wrote output to {output_file}")


if __name__ == "__main__":
    main()
