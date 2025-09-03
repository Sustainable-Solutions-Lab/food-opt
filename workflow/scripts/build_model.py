# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import geopandas as gpd
import pypsa
import logging

logger = logging.getLogger(__name__)


def add_carriers_and_buses(
    n: pypsa.Network,
    crop_list: list,
    food_list: list,
    food_group_list: list,
    regions: list,
    yields_data: dict,
) -> None:
    """Add all carriers and their corresponding buses to the network."""
    # Regional land carriers and buses
    n.add("Carrier", "land", unit="ha")
    n.add("Bus", [f"land_{r}" for r in regions], carrier="land")

    # Crops
    for crop in crop_list:
        crop_bus = f"crop_{crop}"
        n.add("Carrier", crop_bus, unit="t")
        n.add("Bus", crop_bus, carrier=crop_bus)

    # Foods
    for food in food_list:
        food_bus = f"food_{food}"
        n.add("Carrier", food_bus, unit="t")
        n.add("Bus", food_bus, carrier=food_bus)

    # Food groups
    for group in food_group_list:
        group_bus = f"group_{group}"
        n.add("Carrier", group_bus, unit="t")
        n.add("Bus", group_bus, carrier=group_bus)

    # Primary resources and nutrients
    primary_and_nutrients = [
        ("water", "m^3"),
        ("fertilizer", "kg"),
        ("co2", "kg"),
        ("ch4", "kg"),
        ("carb", "t"),
        ("protein", "t"),
        ("fat", "t"),
    ]

    # Add buses for non-regional carriers
    for carrier, unit in primary_and_nutrients:
        n.add("Carrier", carrier, unit=unit)
        n.add("Bus", carrier, carrier=carrier)


def add_primary_resources(
    n: pypsa.Network, config: dict, region_crop_areas: pd.Series
) -> None:
    """Add stores for primary resources with their limits."""
    # Add stores for global resources
    for carrier in ["water", "fertilizer"]:
        n.add("Store", carrier, bus=carrier, carrier=carrier)
        n.add("Generator", carrier, bus=carrier, carrier=carrier, p_nom_extendable=True)

    land_carrier = "land_" + region_crop_areas.index
    n.add(
        "Generator",
        land_carrier,
        bus=land_carrier,
        carrier="land",
        p_nom_extendable=True,
        p_nom_max=region_crop_areas.values,
    )

    # Add stores for emissions with costs to create objective
    n.add("Store", "co2", bus="co2", carrier="co2", e_nom_extendable=True)
    n.add("Store", "ch4", bus="ch4", carrier="ch4", e_nom_extendable=True)

    # Set resource limits from config
    for resource in ["water", "fertilizer"]:
        if resource in config["primary"]:
            n.stores.at[resource, "e_nom_max"] = float(
                config["primary"][resource]["limit"]
            )


def add_regional_crop_production_links(
    n: pypsa.Network, crop_list: list, crops: pd.DataFrame, yields_data: dict
) -> None:
    """Add links for crop production per region and resource class."""
    for crop in crop_list:
        if crop not in crops.index.get_level_values(0):
            logger.warning("Crop '%s' not found in crops data, skipping", crop)
            continue

        crop_data = crops.loc[crop]
        crop_bus = f"crop_{crop}"

        # Get global production coefficients
        water_use = crop_data.loc["water", "value"]  # m^3/t
        fert_use = crop_data.loc["fertilizer", "value"]  # kg/t

        # Get emission coefficients (if they exist)
        co2_emission = 0
        ch4_emission = 0
        if "co2" in crop_data.index:
            co2_emission = crop_data.loc["co2", "value"]  # kg/t
        if "ch4" in crop_data.index:
            ch4_emission = crop_data.loc["ch4", "value"]  # kg/t

        # Get regional yields data for this crop
        crop_yields = yields_data[f"{crop}_yield"]

        # Add a "name" column to crop_yields following this: f"produce_{crop}_{region}_class{resource_class}"
        crop_yields["name"] = crop_yields.index.map(
            lambda x: f"produce_{crop}_{x[0]}_class{x[1]}"
        )

        # Make index levels columns
        df = crop_yields.reset_index()

        # Set index to "name"
        df.set_index("name", inplace=True)
        df.index.name = None

        # Filter out rows with zero suitable area
        df = df[(df["suitable_area"] > 0) & (df["yield"] > 0)]

        # Add links
        link_params = {
            "name": df.index,
            "carrier": "crop_production",
            "bus0": df["region"].apply(lambda x: f"land_{x}"),
            "bus1": crop_bus,
            "efficiency": df["yield"],
            "bus2": "water",
            "efficiency2": -water_use / df["yield"],
            "bus3": "fertilizer",
            "efficiency3": -fert_use / df["yield"],
            "marginal_cost": 0.01,
            "p_nom_max": df["suitable_area"],
            "p_nom_extendable": True,
        }

        # Add emission outputs if they exist
        if co2_emission > 0:
            link_params["bus4"] = "co2"
            link_params["efficiency4"] = co2_emission / df["yield"]

        if ch4_emission > 0:
            bus_idx = 5 if co2_emission > 0 else 4
            link_params[f"bus{bus_idx}"] = "ch4"
            link_params[f"efficiency{bus_idx}"] = ch4_emission / df["yield"]

        n.add("Link", **link_params)


def add_food_conversion_links(
    n: pypsa.Network, food_list: list, foods: pd.DataFrame
) -> None:
    """Add links for converting crops to foods."""
    for _, row in foods.iterrows():
        if row["food"] not in food_list:
            continue
        crop = row["crop"]
        food = row["food"]
        factor = row["factor"]
        link_name = f"convert_{crop}_to_{food.replace(' ', '_').replace('(', '').replace(')', '')}"
        n.add(
            "Link",
            link_name,
            bus0=f"crop_{crop}",  # Input: crop
            bus1=f"food_{food}",  # Output: food
            efficiency=factor,  # efficiency from crop to food
            marginal_cost=0.01,
            p_nom_extendable=True,
        )  # Small cost to enable optimization


def add_food_group_buses_and_loads(
    n: pypsa.Network, food_group_list: list, food_groups: pd.DataFrame, config: dict
) -> None:
    """Add carriers, buses, and loads for food groups defined in the CSV."""
    # Add loads for food groups with requirements
    if "food_groups" in config:
        logger.info("Adding food group loads based on nutrition requirements...")
        for group in food_group_list:
            if group in config["food_groups"]:
                group_config = config["food_groups"][group]
                if "min_per_person_per_day" in group_config:
                    # Calculate total annual requirement
                    population = config.get("population", 1000000)
                    days_per_year = 365
                    min_per_person_per_day = float(
                        group_config["min_per_person_per_day"]
                    )  # g/person/day
                    total_annual_requirement = (
                        min_per_person_per_day
                        * population
                        * days_per_year
                        / 1000000  # Convert g to tonnes
                    )

                    n.add(
                        "Load",
                        group,
                        bus=f"group_{group}",
                        carrier=group,
                        p_set=total_annual_requirement,
                    )

                    # Add an extensible store at the food group bus to store excess food
                    n.add(
                        "Store",
                        f"store_{group}",
                        bus=f"group_{group}",
                        carrier=f"group_{group}",
                        e_nom_extendable=True,
                    )


def add_macronutrient_loads(n: pypsa.Network, config: dict) -> None:
    """Add loads for macronutrients based on minimum requirements."""
    if "macronutrients" in config:
        logger.info("Adding macronutrient loads based on nutrition requirements...")
        for nutrient in ["carb", "protein", "fat"]:
            if nutrient in n.buses.index and nutrient in config["macronutrients"]:
                nutrient_config = config["macronutrients"][nutrient]
                if "min_per_person_per_day" in nutrient_config:
                    # Calculate total annual requirement
                    population = config.get("population", 1000000)
                    days_per_year = 365
                    min_per_person_per_day = float(
                        nutrient_config["min_per_person_per_day"]
                    )  # g/person/day
                    total_annual_requirement = (
                        min_per_person_per_day
                        * population
                        * days_per_year
                        / 1000000  # Convert g to tonnes
                    )

                    n.add(
                        "Load",
                        nutrient,
                        bus=nutrient,
                        carrier=nutrient,
                        p_set=total_annual_requirement,
                    )

                    # Add an extensible store at the macronutrient bus to store excess nutrients
                    n.add(
                        "Store",
                        f"store_{nutrient}",
                        bus=nutrient,
                        carrier=nutrient,
                        e_nom_extendable=True,
                    )


def add_food_nutrition_links(
    n: pypsa.Network,
    food_list: list,
    foods: pd.DataFrame,
    food_groups: pd.DataFrame,
    nutrition: pd.DataFrame,
) -> None:
    """Add multilinks for converting foods to food groups and macronutrients."""
    for food in food_list:
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
                food_group = f"group_{group_val}"

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
        for nutrient in nutrition.index.get_level_values("nutrient").unique():
            link_params[f"bus{bus_idx}"] = nutrient
            eff_param_name = "efficiency" if bus_idx == 1 else f"efficiency{bus_idx}"
            link_params[eff_param_name] = nutrition.loc[(food, nutrient), "value"]
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
    yields_data: dict,
    regions: list,
    region_crop_areas: pd.Series,
) -> pypsa.Network:
    """Build the complete food systems optimization network."""
    n = pypsa.Network()
    n.set_snapshots(["now"])

    # Read in crops from config, and propagate find which foods and food groups we have
    crop_list = config["crops"]
    food_list = foods.loc[foods["crop"].isin(crop_list), "food"].unique()
    food_group_list = food_groups.loc[
        food_groups["food"].isin(food_list), "group"
    ].unique()

    # Build network step by step
    add_carriers_and_buses(
        n, crop_list, food_list, food_group_list, regions, yields_data
    )
    add_primary_resources(n, config, region_crop_areas)
    add_regional_crop_production_links(
        n,
        crop_list,
        crops,
        yields_data,
    )
    add_food_conversion_links(n, food_list, foods)
    add_food_group_buses_and_loads(n, food_group_list, food_groups, config)
    add_macronutrient_loads(n, config)
    add_food_nutrition_links(n, food_list, foods, food_groups, nutrition)

    return n


if __name__ == "__main__":
    # Read crop data
    crops = pd.read_csv(snakemake.input.crops, index_col=["crop", "param"])

    # Read food conversion data
    foods = pd.read_csv(snakemake.input.foods)

    # Read food groups data
    food_groups = pd.read_csv(snakemake.input.food_groups)

    # Read nutrition data
    nutrition = pd.read_csv(snakemake.input.nutrition, index_col=["food", "nutrient"])

    # Read yields data for each crop
    yields_data = {}
    for crop in snakemake.config["crops"]:
        yields_key = f"{crop}_yield"
        yields_df = pd.read_csv(
            snakemake.input[yields_key], index_col=["region", "resource_class"]
        )
        yields_data[yields_key] = yields_df
        logger.info("Loaded yields for %s: %d regions/classes", crop, len(yields_df))

    # Read regions
    regions_df = gpd.read_file(snakemake.input.regions)
    regions = regions_df["region"].tolist()

    # Load region crop areas
    region_crop_areas_df = pd.read_csv(snakemake.input.region_crop_areas)
    region_crop_areas = region_crop_areas_df.set_index("region")["cropland_area_ha"]

    logger.debug("Crops data:\n%s", crops.head(10))
    logger.debug("Foods data:\n%s", foods.head())
    logger.debug("Food groups data:\n%s", food_groups.head())
    logger.debug("Nutrition data:\n%s", nutrition.head())

    # Build the network
    n = build_network(
        snakemake.config,
        crops,
        foods,
        food_groups,
        nutrition,
        yields_data,
        regions,
        region_crop_areas,
    )

    logger.info("Network summary:")
    logger.info("Carriers: %d", len(n.carriers))
    logger.info("Buses: %d", len(n.buses))
    logger.info("Stores: %d", len(n.stores))
    logger.info("Links: %d", len(n.links))

    n.export_to_netcdf(snakemake.output.network)
