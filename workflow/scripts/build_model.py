# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import pypsa


def add_carriers_and_buses(n: pypsa.Network, config: dict) -> None:
    """Add all carriers and their corresponding buses to the network."""
    # Primary resources
    n.add("Carrier", name="land", unit="ha")  # Changed to ha to match crops data
    n.add("Carrier", "water", unit="m^3")
    n.add("Carrier", "fertilizer", unit="kg")
    n.add("Carrier", "co2", unit="kg")
    n.add("Carrier", "ch4", unit="kg")

    # Crops
    for crop in config["crops"]:
        n.add("Carrier", f"crop_{crop}", unit="t")

    # Foods (we'll add these dynamically based on foods.csv)

    # Nutritional values (macronutrients)
    n.add("Carrier", "carb", unit="t")
    n.add("Carrier", "protein", unit="t")
    n.add("Carrier", "fat", unit="t")

    # Add buses for each carrier
    for carrier in n.carriers.index:
        n.add("Bus", carrier, carrier=carrier)


def add_primary_resources(n: pypsa.Network, config: dict) -> None:
    """Add stores for primary resources with their limits."""
    # Add stores for primary resources
    for carrier in ["land", "water", "fertilizer"]:
        n.add("Store", carrier, bus=carrier, carrier=carrier)
        n.add("Generator", carrier, bus=carrier, carrier=carrier, p_nom_extendable=True)

    # Add stores for emissions with costs to create objective
    n.add("Store", "co2", bus="co2", carrier="co2", e_nom_extendable=True)
    n.add("Store", "ch4", bus="ch4", carrier="ch4", e_nom_extendable=True)

    # Set resource limits from config
    for resource in ["land", "water", "fertilizer"]:
        if resource in config["primary"]:
            n.stores.at[resource, "e_nom_max"] = float(
                config["primary"][resource]["limit"]
            )


def add_food_group_buses_and_stores(
    n: pypsa.Network, food_groups: pd.DataFrame
) -> None:
    """Add carriers, buses, and stores for food groups defined in the CSV."""
    # Get unique food groups from CSV (excluding empty groups)
    unique_groups = food_groups[
        food_groups["group"].notna() & (food_groups["group"] != "")
    ]["group"].unique()

    # Add carriers, buses, and stores for food groups
    for group in unique_groups:
        # Add carrier if it doesn't exist
        if group not in n.carriers.index:
            n.add("Carrier", group, unit="t")
        # Add bus if it doesn't exist
        if group not in n.buses.index:
            n.add("Bus", group, carrier=group)
        # Add store if it doesn't exist
        if group not in n.stores.index:
            n.add("Store", group, bus=group, carrier=group, e_nom_extendable=True)


def add_crop_production_links(
    n: pypsa.Network, crops: pd.DataFrame, config: dict
) -> None:
    """Add links for crop production (land/water/fertilizer -> crop + emissions)."""
    for crop in config["crops"]:
        if crop not in crops.index.get_level_values(0):
            continue

        crop_data = crops.loc[crop]

        # Get production coefficients
        yield_value = crop_data.loc["yield", "value"]  # t/ha
        water_use = crop_data.loc["water", "value"]  # m^3/t
        fert_use = crop_data.loc["fertilizer", "value"]  # kg/t

        # Get emission coefficients (if they exist)
        co2_emission = 0
        ch4_emission = 0
        if "co2" in crop_data.index:
            co2_emission = crop_data.loc["co2", "value"]  # kg/t
        if "ch4" in crop_data.index:
            ch4_emission = crop_data.loc["ch4", "value"]  # kg/t

        # Add crop production link
        # bus0 is the primary input (land), bus1 is the main output (crop)
        # Other buses are additional inputs (negative efficiency) or outputs (positive)
        link_name = f"produce_{crop}"
        crop_bus = f"crop_{crop}"

        # Start with basic parameters - add small cost to make optimization work
        link_params = {
            "bus0": "land",  # Primary input
            "bus1": crop_bus,  # Main output
            "efficiency": yield_value,  # Efficiency from land to crop (t crop per ha land)
            "bus2": "water",
            "efficiency2": (
                -water_use / yield_value
            ),  # Additional input per unit output (m^3 per t crop)
            "bus3": "fertilizer",
            "efficiency3": (
                -fert_use / yield_value
            ),  # Additional input per unit output (kg per t crop)
            "marginal_cost": 0.01,  # Small cost to enable optimization
        }

        # Add emission outputs if they exist
        bus_idx = 4
        if co2_emission > 0:
            link_params[f"bus{bus_idx}"] = "co2"
            link_params[f"efficiency{bus_idx}"] = (
                co2_emission / yield_value  # Output per unit crop (kg per t crop)
            )
            bus_idx += 1

        if ch4_emission > 0:
            link_params[f"bus{bus_idx}"] = "ch4"
            link_params[f"efficiency{bus_idx}"] = (
                ch4_emission / yield_value  # Output per unit crop (kg per t crop)
            )

        n.add("Link", link_name, p_nom_extendable=True, **link_params)


def add_food_conversion_links(n: pypsa.Network, foods: pd.DataFrame) -> None:
    """Add links for converting crops to foods."""
    # First add food carriers and buses
    unique_foods = foods["food"].unique()
    for food in unique_foods:
        food_bus = f"food_{food}"
        if food_bus not in n.carriers.index:
            n.add("Carrier", food_bus, unit="t")
            n.add("Bus", food_bus, carrier=food_bus)

    # Add conversion links
    for _, row in foods.iterrows():
        crop = row["crop"]
        food = row["food"]
        factor = row["factor"]

        crop_bus = f"crop_{crop}"
        food_bus = f"food_{food}"

        # Only add link if the crop is in our network
        if crop_bus in n.buses.index:
            link_name = f"convert_{crop}_to_{food.replace(' ', '_').replace('(', '').replace(')', '')}"
            n.add(
                "Link",
                link_name,
                bus0=crop_bus,  # Input: crop (primary input)
                bus1=food_bus,  # Output: food
                efficiency=factor,  # efficiency from crop to food
                marginal_cost=0.01,
                p_nom_extendable=True,
            )  # Small cost to enable optimization


def add_macronutrient_stores(n: pypsa.Network) -> None:
    """Add stores for macronutrients."""
    for nutrient in ["carb", "protein", "fat"]:
        if nutrient in n.buses.index:
            n.add(
                "Store", nutrient, bus=nutrient, carrier=nutrient, e_nom_extendable=True
            )


def add_food_nutrition_links(
    n: pypsa.Network,
    foods: pd.DataFrame,
    food_groups: pd.DataFrame,
    nutrition: pd.DataFrame,
) -> None:
    """Add multilinks for converting foods to food groups and macronutrients."""

    # Create nutrition lookup for faster access
    nutrition_dict = {}
    for _, row in nutrition.iterrows():
        food = row["food"]
        nutrient = row["nutrient"]
        value = row["value"]  # g/100g
        if food not in nutrition_dict:
            nutrition_dict[food] = {}
        nutrition_dict[food][nutrient] = value / 100.0  # Convert to fraction

    # Add multilinks for each food from foods.csv
    unique_foods = foods["food"].unique()
    for food in unique_foods:
        food_bus = f"food_{food}"
        # Only process foods that exist in our network
        if food_bus not in n.buses.index:
            continue

        # Find food group for this food
        food_group = None
        food_group_row = food_groups[food_groups["food"] == food]
        if not food_group_row.empty:
            group_val = food_group_row.iloc[0]["group"]
            if pd.notna(group_val) and group_val != "":
                food_group = group_val

        # Get nutrition data for this food - must exist for all foods
        if food not in nutrition_dict:
            raise ValueError(f"Nutritional values not defined for food: {food}")
        nutrition_data = nutrition_dict[food]

        # Check that all required macronutrients are present
        required_nutrients = ["carb", "protein", "fat"]
        missing_nutrients = [n for n in required_nutrients if n not in nutrition_data]
        if missing_nutrients:
            raise ValueError(
                f"Missing nutritional values for {food}: {missing_nutrients}"
            )

        link_name = (
            f"consume_{food.replace(' ', '_').replace('(', '').replace(')', '')}"
        )

        # Start with basic parameters
        link_params = {
            "bus0": food_bus,  # Input: food (primary input)
            "marginal_cost": 0.01,  # Small cost to enable optimization
        }

        bus_idx = 1

        # Add macronutrient outputs (all foods must have these)
        for nutrient in required_nutrients:
            link_params[f"bus{bus_idx}"] = nutrient
            eff_param_name = "efficiency" if bus_idx == 1 else f"efficiency{bus_idx}"
            link_params[eff_param_name] = nutrition_data[nutrient]
            bus_idx += 1

        # Add food group output if applicable
        if food_group:
            link_params[f"bus{bus_idx}"] = food_group
            link_params[f"efficiency{bus_idx}"] = 1.0

        n.add("Link", link_name, p_nom_extendable=True, **link_params)


def build_network(
    config: dict,
    crops: pd.DataFrame,
    foods: pd.DataFrame,
    food_groups: pd.DataFrame,
    nutrition: pd.DataFrame,
) -> pypsa.Network:
    """Build the complete food systems optimization network."""
    n = pypsa.Network()
    n.set_snapshots(["now"])

    # Build network step by step
    add_carriers_and_buses(n, config)
    add_primary_resources(n, config)
    add_crop_production_links(n, crops, config)
    add_food_conversion_links(n, foods)
    add_food_group_buses_and_stores(n, food_groups)
    add_macronutrient_stores(n)
    add_food_nutrition_links(n, foods, food_groups, nutrition)

    return n


if __name__ == "__main__":
    # Read crop data
    crops = pd.read_csv(snakemake.input.crops, index_col=["crop", "param"])

    # Read food conversion data
    foods = pd.read_csv(snakemake.input.foods)

    # Read food groups data
    food_groups = pd.read_csv(snakemake.input.food_groups)

    # Read nutrition data
    nutrition = pd.read_csv(snakemake.input.nutrition)

    print("Crops data:")
    print(crops.head(10))
    print("\nFoods data:")
    print(foods.head())
    print("\nFood groups data:")
    print(food_groups.head())
    print("\nNutrition data:")
    print(nutrition.head())

    # Build the network
    n = build_network(snakemake.config, crops, foods, food_groups, nutrition)

    print("\nNetwork summary:")
    print(f"Carriers: {len(n.carriers)}")
    print(f"Buses: {len(n.buses)}")
    print(f"Stores: {len(n.stores)}")
    print(f"Links: {len(n.links)}")

    n.export_to_netcdf(snakemake.output.network)
