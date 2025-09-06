# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import geopandas as gpd
import numpy as np
from sklearn.cluster import KMeans
import pypsa
import logging

logger = logging.getLogger(__name__)


def add_carriers_and_buses(
    n: pypsa.Network,
    crop_list: list,
    food_list: list,
    food_group_list: list,
    regions: list,
    countries: list,
) -> None:
    """Add all carriers and their corresponding buses to the network.

    - Regional land buses remain per-region.
    - Crops, foods, food groups, and macronutrients are created per-country.
    - Primary resources (water, fertilizer) and emissions (co2, ch4) stay global.
    """
    # Regional land carriers and buses
    n.add("Carrier", "land", unit="ha")
    n.add("Bus", [f"land_{r}" for r in regions], carrier="land")

    # Crops per country
    crop_buses = [
        f"crop_{crop}_{country}" for country in countries for crop in crop_list
    ]
    crop_carriers = [f"crop_{crop}" for country in countries for crop in crop_list]
    if crop_buses:
        n.add("Carrier", sorted(set(f"crop_{crop}" for crop in crop_list)), unit="t")
        n.add("Bus", crop_buses, carrier=crop_carriers)

    # Foods per country
    food_buses = [
        f"food_{food}_{country}" for country in countries for food in food_list
    ]
    food_carriers = [f"food_{food}" for country in countries for food in food_list]
    if food_buses:
        n.add("Carrier", sorted(set(f"food_{food}" for food in food_list)), unit="t")
        n.add("Bus", food_buses, carrier=food_carriers)

    # Food groups per country
    group_buses = [
        f"group_{group}_{country}" for country in countries for group in food_group_list
    ]
    group_carriers = [
        f"group_{group}" for country in countries for group in food_group_list
    ]
    if group_buses:
        n.add(
            "Carrier",
            sorted(set(f"group_{group}" for group in food_group_list)),
            unit="t",
        )
        n.add("Bus", group_buses, carrier=group_carriers)

    # Macronutrients per country
    nutrient_buses = [
        f"{nut}_{country}"
        for country in countries
        for nut in ["carb", "protein", "fat"]
    ]
    nutrient_carriers = [
        nut for country in countries for nut in ["carb", "protein", "fat"]
    ]
    if nutrient_buses:
        n.add("Carrier", ["carb", "protein", "fat"], unit="t")
        n.add("Bus", nutrient_buses, carrier=nutrient_carriers)

    # Primary resources and emissions remain global
    for carrier, unit in [
        ("water", "m^3"),
        ("fertilizer", "kg"),
        ("co2", "kg"),
        ("ch4", "kg"),
    ]:
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

    land_carrier = "land_" + region_crop_areas.index.astype(str)
    n.add(
        "Generator",
        land_carrier,
        bus=land_carrier,
        carrier="land",
        p_nom_extendable=True,
        p_nom_max=(
            config["primary"]["land"]["regional_limit"] * region_crop_areas.values
        ),
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
    n: pypsa.Network,
    crop_list: list,
    crops: pd.DataFrame,
    yields_data: dict,
    region_to_country: pd.Series,
    allowed_countries: set,
) -> None:
    """Add links for crop production per region and resource class."""
    for crop in crop_list:
        if crop not in crops.index.get_level_values(0):
            logger.warning("Crop '%s' not found in crops data, skipping", crop)
            continue

        crop_data = crops.loc[crop]

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

        # Map regions to countries and filter to allowed countries
        df["country"] = df["region"].map(region_to_country)
        df = df[df["country"].isin(allowed_countries)]

        # Add links
        link_params = {
            "name": df.index,
            "carrier": "crop_production",
            "bus0": df["region"].apply(lambda x: f"land_{x}"),
            "bus1": df["country"].apply(lambda c: f"crop_{crop}_{c}"),
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
    n: pypsa.Network, food_list: list, foods: pd.DataFrame, countries: list
) -> None:
    """Add links for converting crops to foods."""
    for _, row in foods.iterrows():
        if row["food"] not in food_list:
            continue
        crop = row["crop"]
        food = row["food"]
        factor = float(row["factor"]) if pd.notna(row["factor"]) else 1.0
        names = [
            f"convert_{crop}_to_{food.replace(' ', '_').replace('(', '').replace(')', '')}_{c}"
            for c in countries
        ]
        bus0 = [f"crop_{crop}_{c}" for c in countries]
        bus1 = [f"food_{food}_{c}" for c in countries]
        n.add(
            "Link",
            names,
            bus0=bus0,
            bus1=bus1,
            efficiency=[factor] * len(countries),
            marginal_cost=[0.01] * len(countries),
            p_nom_extendable=[True] * len(countries),
        )


def add_food_group_buses_and_loads(
    n: pypsa.Network,
    food_group_list: list,
    food_groups: pd.DataFrame,
    config: dict,
    countries: list,
    population: pd.Series,
) -> None:
    """Add carriers, buses, and loads for food groups defined in the CSV."""
    # Add loads for food groups with requirements
    if "food_groups" in config:
        logger.info("Adding food group loads based on nutrition requirements...")
        for group in food_group_list:
            if group in config["food_groups"]:
                group_config = config["food_groups"][group]
                if "min_per_person_per_day" in group_config:
                    days_per_year = 365
                    min_per_person_per_day = float(
                        group_config["min_per_person_per_day"]
                    )  # g/person/day
                    names = [f"{group}_{c}" for c in countries]
                    buses = [f"group_{group}_{c}" for c in countries]
                    carriers = [f"group_{group}"] * len(countries)
                    p_set = [
                        min_per_person_per_day
                        * float(population[c])
                        * days_per_year
                        / 1_000_000.0
                        for c in countries
                    ]
                    n.add("Load", names, bus=buses, carrier=carriers, p_set=p_set)

                    store_names = [f"store_{group}_{c}" for c in countries]
                    n.add(
                        "Store",
                        store_names,
                        bus=buses,
                        carrier=carriers,
                        e_nom_extendable=[True] * len(countries),
                    )


def add_macronutrient_loads(
    n: pypsa.Network, config: dict, countries: list, population: pd.Series
) -> None:
    """Add per-country loads for macronutrients based on minimum requirements."""
    if "macronutrients" in config:
        logger.info("Adding macronutrient loads per country based on requirements...")
        for nutrient in ["carb", "protein", "fat"]:
            if nutrient in config["macronutrients"]:
                nutrient_config = config["macronutrients"][nutrient]
                if "min_per_person_per_day" in nutrient_config:
                    days_per_year = 365
                    min_per_person_per_day = float(
                        nutrient_config["min_per_person_per_day"]
                    )  # g/person/day
                    names = [f"{nutrient}_{c}" for c in countries]
                    buses = names
                    carriers = [nutrient] * len(countries)
                    p_set = [
                        min_per_person_per_day
                        * float(population[c])
                        * days_per_year
                        / 1_000_000.0
                        for c in countries
                    ]
                    n.add("Load", names, bus=buses, carrier=carriers, p_set=p_set)
                    store_names = [f"store_{nutrient}_{c}" for c in countries]
                    n.add(
                        "Store",
                        store_names,
                        bus=buses,
                        carrier=carriers,
                        e_nom_extendable=[True] * len(countries),
                    )


def add_food_nutrition_links(
    n: pypsa.Network,
    food_list: list,
    foods: pd.DataFrame,
    food_groups: pd.DataFrame,
    nutrition: pd.DataFrame,
    countries: list,
) -> None:
    """Add multilinks per country for converting foods to groups and macronutrients."""
    # Pre-index food_groups for lookup
    food_to_group = food_groups.set_index("food")["group"].to_dict()

    nutrients = list(nutrition.index.get_level_values("nutrient").unique())
    for food in food_list:
        group_val = food_to_group.get(food, None)
        names = [
            f"consume_{food.replace(' ', '_').replace('(', '').replace(')', '')}_{c}"
            for c in countries
        ]
        bus0 = [f"food_{food}_{c}" for c in countries]

        # macronutrient outputs
        out_bus_lists = []
        eff_lists = []
        for i, nutrient in enumerate(nutrients, start=1):
            out_bus_lists.append([f"{nutrient}_{c}" for c in countries])
            eff_val = (
                float(nutrition.loc[(food, nutrient), "value"])
                if (food, nutrient) in nutrition.index
                else 0.0
            )
            eff_lists.append([eff_val] * len(countries))

        params = {"bus0": bus0, "marginal_cost": [0.01] * len(countries)}
        for i, (buses, effs) in enumerate(zip(out_bus_lists, eff_lists), start=1):
            params[f"bus{i}"] = buses
            params["efficiency" if i == 1 else f"efficiency{i}"] = effs

        # optional food group output as last leg
        if group_val is not None and pd.notna(group_val):
            idx = len(nutrients) + 1
            params[f"bus{idx}"] = [f"group_{group_val}_{c}" for c in countries]
            params[f"efficiency{idx}"] = [1.0] * len(countries)

        n.add("Link", names, p_nom_extendable=[True] * len(countries), **params)


def add_crop_trade_hubs_and_links(
    n: pypsa.Network,
    config: dict,
    regions_gdf: gpd.GeoDataFrame,
    countries: list,
    crop_list: list,
) -> None:
    """Add crop trading hubs and connect crop buses via hubs.

    - Cluster all model regions' centroids (in Equal Earth projection) into K hubs.
    - Add a per-crop bus at each hub centroid (to match crop carriers).
    - Connect each country's crop bus to its nearest hub (bidirectional, extendable links).
      Nearest hub is computed from the centroid of the dissolved country's regions.
    - Fully connect the hub graph (bidirectional, extendable links).

    Marginal cost is distance-dependent and controlled by config key
    config["trade"]["crop_trade_marginal_cost_per_km"] with default 1e-4 per km.
    Number of hubs from config["trade"]["crop_hubs"] with default 20.
    """
    trade_cfg = config.get("trade", {})
    n_hubs = int(trade_cfg.get("crop_hubs", 20))
    cost_per_km = float(trade_cfg.get("crop_trade_marginal_cost_per_km", 1e-4))

    if len(regions_gdf) == 0 or len(crop_list) == 0 or len(countries) == 0:
        logger.info("Skipping trade hubs: no regions/crops/countries available")
        return

    # Ensure CRS and project to an equal-area/earth projection for distance in meters
    gdf = regions_gdf.copy()
    gdf_ee = gdf.to_crs(6933)

    # Region centroids (in meters) for clustering
    cent = gdf_ee.geometry.centroid
    X = np.column_stack([cent.x.values, cent.y.values])
    k = min(max(1, n_hubs), len(X))
    if k < n_hubs:
        logger.info("Reducing hub count from %d to %d (regions=%d)", n_hubs, k, len(X))
        n_hubs = k

    # K-means clustering to get hub centers
    km = KMeans(n_clusters=n_hubs, n_init=10, random_state=0)
    km.fit_predict(X)
    centers = km.cluster_centers_  # shape (n_hubs, 2)

    # Add per-crop hub buses at these centers
    hub_ids = list(range(n_hubs))
    for crop in crop_list:
        hub_bus_names = [f"hub_{h}_{crop}" for h in hub_ids]
        hub_carriers = [f"crop_{crop}"] * n_hubs
        n.add("Bus", hub_bus_names, carrier=hub_carriers)

    # Compute per-country centroid from dissolved regions (projected)
    gdf_countries = gdf_ee[gdf_ee["country"].isin(countries)].dissolve(
        by="country", as_index=True
    )
    ccent = gdf_countries.geometry.centroid
    C = np.column_stack([ccent.x.values, ccent.y.values])  # meters
    # Map country -> nearest hub index
    # Compute distances to all hubs (vectorized)
    dch = ((C[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2) ** 0.5  # meters
    nearest_hub_idx = dch.argmin(axis=1)
    nearest_hub_dist_km = dch[np.arange(len(C)), nearest_hub_idx] / 1000.0

    country_index = gdf_countries.index.to_list()
    country_to_hub = {c: int(h) for c, h in zip(country_index, nearest_hub_idx)}
    country_to_dist_km = {
        c: float(d) for c, d in zip(country_index, nearest_hub_dist_km)
    }

    # Connect each country's crop bus to its nearest hub with a single
    # bidirectional link (p_nom_min = -inf)
    valid_countries = [c for c in countries if c in country_to_hub]
    for crop in crop_list:
        if not valid_countries:
            continue
        names_from_c = [
            f"trade_{crop}_{c}_hub{country_to_hub[c]}" for c in valid_countries
        ]
        names_from_hub = [
            f"trade_{crop}_hub{country_to_hub[c]}_{c}" for c in valid_countries
        ]
        bus0 = [f"crop_{crop}_{c}" for c in valid_countries]
        bus1 = [f"hub_{country_to_hub[c]}_{crop}" for c in valid_countries]
        costs = [country_to_dist_km[c] * cost_per_km for c in valid_countries]
        n.add(
            "Link",
            names_from_c,
            bus0=bus0,
            bus1=bus1,
            marginal_cost=costs,
            p_nom_extendable=True,
        )
        n.add(
            "Link",
            names_from_hub,
            bus0=bus1,
            bus1=bus0,
            marginal_cost=costs,
            p_nom_extendable=True,
        )

    # Fully connect hubs (complete graph), per crop, add two directed links per pair
    if n_hubs >= 2:
        H = centers
        # Pairwise distances (km) for all ordered hub pairs (i != j)
        D = np.sqrt(((H[:, None, :] - H[None, :, :]) ** 2).sum(axis=2)) / 1000.0
        ii, jj = np.where(~np.eye(n_hubs, dtype=bool))
        dists_km = D[ii, jj].tolist()

        for crop in crop_list:
            if len(ii) == 0:
                continue
            names = [f"trade_{crop}_hub{i}_to_hub{j}" for i, j in zip(ii, jj)]
            bus0 = [f"hub_{i}_{crop}" for i in ii]
            bus1 = [f"hub_{j}_{crop}" for j in jj]
            costs = [d * cost_per_km for d in dists_km]
            n.add(
                "Link",
                names,
                bus0=bus0,
                bus1=bus1,
                marginal_cost=costs,
                p_nom_extendable=True,
            )


def build_network(
    config: dict,
    crops: pd.DataFrame,
    foods: pd.DataFrame,
    food_groups: pd.DataFrame,
    nutrition: pd.DataFrame,
    yields_data: dict,
    regions: list,
    regions_gdf: gpd.GeoDataFrame,
    region_to_country: pd.Series,
    countries: list,
    region_crop_areas: pd.Series,
    population: pd.Series,
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
    add_carriers_and_buses(n, crop_list, food_list, food_group_list, regions, countries)
    add_primary_resources(n, config, region_crop_areas)
    add_regional_crop_production_links(
        n,
        crop_list,
        crops,
        yields_data,
        region_to_country,
        set(countries),
    )
    add_food_conversion_links(n, food_list, foods, countries)
    add_food_group_buses_and_loads(
        n, food_group_list, food_groups, config, countries, population
    )
    add_macronutrient_loads(n, config, countries, population)
    add_food_nutrition_links(n, food_list, foods, food_groups, nutrition, countries)

    # Add crop trading hubs and links (hierarchical trade network)
    add_crop_trade_hubs_and_links(n, config, regions_gdf, countries, list(crop_list))

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

    # Load population per country for planning horizon
    population_df = pd.read_csv(snakemake.input.population)
    # Expect columns: iso3, country, year, population
    # Select only configured countries and validate coverage
    cfg_countries = list(snakemake.config.get("countries", []))
    pop_map = population_df.set_index("iso3")["population"].reindex(cfg_countries)
    missing = pop_map[pop_map.isna()].index.tolist()
    if missing:
        raise ValueError("Missing population for countries: " + ", ".join(missing))
    # population series indexed by country code (ISO3)
    population = pop_map.astype(float)

    region_to_country = regions_df.set_index("region")["country"]
    # Warn if any configured countries are missing from regions
    present_countries = set(region_to_country.unique())
    missing_in_regions = [c for c in cfg_countries if c not in present_countries]
    if missing_in_regions:
        logger.warning(
            "Configured countries missing from regions and may be disconnected: %s",
            ", ".join(sorted(missing_in_regions)),
        )
    # Keep only regions whose country is in configured countries
    region_to_country = region_to_country[region_to_country.isin(cfg_countries)]

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
        regions_df,
        region_to_country,
        cfg_countries,
        region_crop_areas,
        population,
    )

    logger.info("Network summary:")
    logger.info("Carriers: %d", len(n.carriers))
    logger.info("Buses: %d", len(n.buses))
    logger.info("Stores: %d", len(n.stores))
    logger.info("Links: %d", len(n.links))

    n.export_to_netcdf(snakemake.output.network)
