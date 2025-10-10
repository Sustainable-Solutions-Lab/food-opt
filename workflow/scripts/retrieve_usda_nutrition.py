# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Retrieve nutritional data from USDA FoodData Central API.

This script fetches nutrition data for foods listed in data/usda_food_mapping.csv
and outputs in the format expected by the model (data/nutrition.csv).

Requires network access to call the USDA API.
"""

import time
from pathlib import Path

import pandas as pd
import requests


BASE_URL = "https://api.nal.usda.gov/fdc/v1"


def get_food_nutrients(fdc_id: int, api_key: str) -> dict[str, tuple[float, str]]:
    """
    Retrieve nutrient data for a specific FDC ID.

    Parameters
    ----------
    fdc_id : int
        USDA FoodData Central food ID
    api_key : str
        USDA API key

    Returns
    -------
    dict[str, tuple[float, str]]
        Mapping from nutrient name to (amount, unit) tuple
    """
    url = f"{BASE_URL}/food/{fdc_id}"
    params = {"api_key": api_key}

    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()
    data = response.json()

    nutrients = {}
    for nutrient_entry in data.get("foodNutrients", []):
        name = nutrient_entry["nutrient"]["name"]
        amount = nutrient_entry.get("amount", 0.0)
        unit = nutrient_entry["nutrient"]["unitName"]
        nutrients[name] = (amount, unit)

    return nutrients


def main():
    # Get inputs from Snakemake or defaults
    try:
        mapping_path = snakemake.input.mapping  # noqa: F821
        output_path = snakemake.output[0]  # noqa: F821
        api_key = snakemake.config["data"]["usda"]["api_key"]  # noqa: F821
        nutrient_name_map = snakemake.config["data"]["usda"]["nutrients"]  # noqa: F821
    except NameError:
        # Fallback for manual testing
        mapping_path = "data/usda_food_mapping.csv"
        output_path = "data/nutrition.csv"
        api_key = "NZESBYEsg9Utlh3OmIsHu7pgEI2AD3jN76SxgCq6"
        nutrient_name_map = {
            "protein": "Protein",
            "carb": "Carbohydrate, by difference",
            "fat": "Total lipid (fat)",
            "kcal": "Energy",
        }

    # Read mapping file
    mapping_df = pd.read_csv(mapping_path, comment="#")

    # Prepare output data
    output_rows = []

    # Map USDA units to our internal units
    unit_map = {
        "G": "g/100g",
        "KCAL": "kcal/100g",
        "KJ": "kcal/100g",  # Will convert from kJ to kcal
    }

    print("Retrieving nutrition data from USDA FoodData Central...")
    print(f"Processing {len(mapping_df)} foods...")

    for idx, row in mapping_df.iterrows():
        food_name = row["food"]
        fdc_id = int(row["fdc_id"])

        print(f"  [{idx + 1}/{len(mapping_df)}] {food_name} (FDC {fdc_id})")

        try:
            nutrients = get_food_nutrients(fdc_id, api_key)

            # Extract requested nutrients
            for internal_name, usda_nutrient_name in nutrient_name_map.items():
                if usda_nutrient_name in nutrients:
                    amount, usda_unit = nutrients[usda_nutrient_name]

                    # Convert units if necessary
                    if usda_unit.upper() == "KJ":
                        # Convert kilojoules to kilocalories
                        amount = amount / 4.184
                        internal_unit = "kcal/100g"
                    else:
                        # Map USDA units to internal units
                        internal_unit = unit_map.get(
                            usda_unit.upper(), usda_unit.lower()
                        )

                    output_rows.append(
                        {
                            "food": food_name,
                            "nutrient": internal_name,
                            "unit": internal_unit,
                            "value": round(amount, 2),
                        }
                    )
                else:
                    print(
                        f"    Warning: {usda_nutrient_name} not found for {food_name}"
                    )

            # Be polite to the API
            time.sleep(0.5)

        except Exception as e:
            print(f"    Error retrieving data for {food_name}: {e}")
            continue

    # Write output
    output_df = pd.DataFrame(output_rows)
    output_df = output_df.sort_values(["food", "nutrient"])

    # Ensure output directory exists
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # Write with comment header
    with open(output_path, "w", newline="") as f:
        f.write("# Nutritional data retrieved from USDA FoodData Central\n")
        f.write("# Source: https://fdc.nal.usda.gov/\n")
        output_df.to_csv(f, index=False)

    print(f"\nWrote {len(output_rows)} nutrient entries to {output_path}")


if __name__ == "__main__":
    main()
