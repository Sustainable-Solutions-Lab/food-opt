#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Process Global Dietary Database (GDD) country-level dietary intake data.

Reads multiple v*_cnty.csv files from GDD, extracts national-level baseline
intake estimates (aggregated across age/sex/urban/education strata), maps
GDD food variables to model dietary risk factors, and outputs a consolidated
baseline diet file matching the DIA format.

Input:
    - GDD directory with Country-level estimates/*.csv files
    - Reference year from config

Output:
    - CSV with columns: scenario,unit,item,country,year,value
    - Format compatible with prepare_health_costs.py expectations
"""

import sys
from pathlib import Path

import pandas as pd


def main():
    gdd_dir = Path(snakemake.input["gdd_dir"])
    reference_year = snakemake.params["reference_year"]
    output_file = snakemake.output["diet"]

    # Map GDD variable codes (vXX) to model dietary items
    # Based on codebook and model's item_to_risk mapping in prepare_health_costs.py
    gdd_to_model_items = {
        "v01": "fruits",  # Fruits → fruits (aggregated from temp/trop/starch)
        "v02": "vegetables",  # Non-starchy vegetables → vegetables
        "v03": None,  # Potatoes → not used (starchy, not a risk factor)
        "v04": None,  # Other starchy vegetables → not used
        "v05": "legumes",  # Beans and legumes → legumes
        "v06": "nuts_seeds",  # Nuts and seeds → nuts_seeds
        "v07": None,  # Refined grains → not used
        "v08": "whole_grains",  # Whole grains → whole_grains
        "v09": "prc_meat",  # Total processed meats → prc_meat
        "v10": "red_meat",  # Unprocessed red meats → red_meat
        "v11": "fish",  # Total seafoods → fish (aggregated from types + shellfish)
        "v12": None,  # Eggs → not a GBD dietary risk factor
        "v13": None,  # Cheese → not used
        "v14": None,  # Yoghurt → not used
        "v15": None,  # Sugar-sweetened beverages → not used
        "v16": None,  # Fruit juices → not used
        "v17": None,  # Coffee → not used
        "v18": None,  # Tea → not used
    }

    # Filter to only risk-factor variables
    risk_factor_vars = {
        varcode: item
        for varcode, item in gdd_to_model_items.items()
        if item is not None
    }

    print(f"[prepare_gdd_dietary_intake] Processing GDD data for year {reference_year}")
    print(
        f"[prepare_gdd_dietary_intake] Risk factor variables: {list(risk_factor_vars.keys())}"
    )

    country_estimates_dir = gdd_dir / "Country-level estimates"
    if not country_estimates_dir.exists():
        print(
            f"ERROR: Country-level estimates directory not found: {country_estimates_dir}",
            file=sys.stderr,
        )
        sys.exit(1)

    all_data = []

    for varcode, model_item in risk_factor_vars.items():
        csv_file = country_estimates_dir / f"{varcode}_cnty.csv"
        if not csv_file.exists():
            print(f"WARNING: File not found: {csv_file}", file=sys.stderr)
            continue

        print(f"[prepare_gdd_dietary_intake] Reading {csv_file.name} ({model_item})...")
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
            nearest_year = min(available_years, key=lambda y: abs(y - reference_year))
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

        natl["item"] = model_item
        natl["year"] = reference_year
        all_data.append(natl)

    if not all_data:
        print("ERROR: No data collected from GDD files", file=sys.stderr)
        sys.exit(1)

    # Concatenate all risk factors
    result = pd.concat(all_data, ignore_index=True)

    # Add scenario and unit columns to match DIA format
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
        f"[prepare_gdd_dietary_intake] Risk factors: {sorted(result['item'].unique())}"
    )

    # Write output
    result.to_csv(output_file, index=False)
    print(f"[prepare_gdd_dietary_intake] Wrote output to {output_file}")


if __name__ == "__main__":
    main()
