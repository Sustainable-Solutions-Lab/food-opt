# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import functools
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
    countries: list,
) -> None:
    """Add all carriers and their corresponding buses to the network.

    - Regional land buses remain per-region.
    - Crops, foods, food groups, and macronutrients are created per-country.
    - Primary resources (water, fertilizer) and emissions (co2, ch4) stay global.
    """
    # Land carrier (class-level buses are added later)
    n.add("Carrier", "land", unit="ha")

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

    # Feed carriers per country (ruminant-wide and monogastric-only pools)
    feed_types = ["ruminant", "monogastric"]
    feed_buses = [f"feed_{ft}_{country}" for country in countries for ft in feed_types]
    feed_carriers = [f"feed_{ft}" for country in countries for ft in feed_types]
    if feed_buses:
        n.add("Carrier", sorted(set(feed_carriers)), unit="t")
        n.add("Bus", feed_buses, carrier=feed_carriers)

    # Primary resources and emissions remain global
    for carrier, unit in [
        ("water", "m^3"),
        ("fertilizer", "kg"),
        ("co2", "kg"),
        ("ch4", "kg"),
    ]:
        n.add("Carrier", carrier, unit=unit)
        n.add("Bus", carrier, carrier=carrier)


def add_primary_resources(n: pypsa.Network, primary_config: dict) -> None:
    """Add stores for primary resources with their limits."""
    # Add stores for global resources
    for carrier in ["water", "fertilizer"]:
        n.add("Store", carrier, bus=carrier, carrier=carrier)
        n.add("Generator", carrier, bus=carrier, carrier=carrier, p_nom_extendable=True)

    # Region-level land limit removed; replaced by class-level land pools

    # Add stores for emissions with costs to create objective
    n.add("Store", "co2", bus="co2", carrier="co2", e_nom_extendable=True)
    n.add("Store", "ch4", bus="ch4", carrier="ch4", e_nom_extendable=True)

    # Set resource limits from config
    for resource in ["water", "fertilizer"]:
        n.stores.at[resource, "e_nom_max"] = float(primary_config[resource]["limit"])


def add_regional_crop_production_links(
    n: pypsa.Network,
    crop_list: list,
    crops: pd.DataFrame,
    yields_data: dict,
    region_to_country: pd.Series,
    allowed_countries: set,
    crop_prices_usd_per_t: pd.Series,
) -> None:
    """Add links for crop production per region/resource class and water supply.

    Expects yields_data to contain per-crop DataFrames for both water supplies
    with keys f"{crop}_yield_i" and f"{crop}_yield_r".
    Output links produce into the same crop bus per country; link names encode
    supply type (i/r) and resource class.
    """
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

        # Process both water supplies
        for ws in ["i", "r"]:
            key = f"{crop}_yield_{ws}"
            if key not in yields_data:
                logger.warning("Missing yields for %s (%s), skipping", crop, ws)
                continue
            crop_yields = yields_data[key].copy()

            # Add a unique name per link including water supply and class
            crop_yields["name"] = crop_yields.index.map(
                lambda x: f"produce_{crop}_{'irrigated' if ws == 'i' else 'rainfed'}_{x[0]}_class{x[1]}"
            )

            # Make index levels columns
            df = crop_yields.reset_index()

            # Set index to "name"
            df.set_index("name", inplace=True)
            df.index.name = None

            # Filter out rows with zero suitable area or zero yield
            df = df[(df["suitable_area"] > 0) & (df["yield"] > 0)]

            # Map regions to countries and filter to allowed countries
            df["country"] = df["region"].map(region_to_country)
            df = df[df["country"].isin(allowed_countries)]

            if df.empty:
                continue

            # Price for this crop (USD/tonne); if missing, warn and use 0
            price = float(crop_prices_usd_per_t.get(crop, float("nan")))
            if not np.isfinite(price):
                logger.warning(
                    "No FAOSTAT price for crop '%s'; defaulting marginal_cost to 0",
                    crop,
                )
                price = 0.0

            # Add links
            # Connect to class-level land bus per region/resource class and water supply
            link_params = {
                "name": df.index,
                # Use the crop's own carrier so no extra carrier is needed
                "carrier": f"crop_{crop}",
                "bus0": df.apply(
                    lambda r: f"land_{r['region']}_class{int(r['resource_class'])}_{'i' if ws == 'i' else 'r'}",
                    axis=1,
                ),
                "bus1": df["country"].apply(lambda c: f"crop_{crop}_{c}"),
                "efficiency": df["yield"],
                "bus2": "water",
                "efficiency2": -water_use / df["yield"],
                "bus3": "fertilizer",
                "efficiency3": -fert_use / df["yield"],
                # Link marginal_cost is per unit of bus0 flow (ha). To apply a
                # cost per tonne on bus1, multiply by efficiency (t/ha).
                "marginal_cost": price * df["yield"],
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


def add_grassland_feed_links(
    n: pypsa.Network,
    grassland: pd.DataFrame,
    land_rainfed: pd.DataFrame,
    region_to_country: pd.Series,
    allowed_countries: set,
) -> None:
    """Add links supplying ruminant feed directly from rainfed land."""

    df = grassland.copy()
    df = df[np.isfinite(df["yield"]) & (df["yield"] > 0)]
    if df.empty:
        logger.info("Grassland yields contain no positive entries; skipping")
        return

    df = df.reset_index()
    df["resource_class"] = df["resource_class"].astype(int)
    df = df.set_index(["region", "resource_class"])

    merged = df.join(
        land_rainfed[["area_ha"]].rename(columns={"area_ha": "land_area"}),
        how="inner",
    )
    if merged.empty:
        logger.info(
            "No overlap between grassland yields and rainfed land areas; skipping"
        )
        return

    candidate_area = merged["suitable_area"].fillna(merged["land_area"])
    available_area = np.minimum(
        candidate_area.to_numpy(), merged["land_area"].to_numpy()
    )
    merged["available_area"] = available_area
    merged = merged[merged["available_area"] > 0]
    if merged.empty:
        logger.info("Grassland entries have zero available area; skipping")
        return

    merged = merged.reset_index()
    merged["country"] = merged["region"].map(region_to_country)
    merged = merged[merged["country"].isin(allowed_countries)]
    merged = merged.dropna(subset=["country"])
    if merged.empty:
        logger.info("No grassland regions map to configured countries; skipping")
        return

    merged["name"] = merged.apply(
        lambda r: f"graze_{r['region']}_class{int(r['resource_class'])}", axis=1
    )
    merged["bus0"] = merged.apply(
        lambda r: f"land_{r['region']}_class{int(r['resource_class'])}_r", axis=1
    )
    merged["bus1"] = merged["country"].apply(lambda c: f"feed_ruminant_{c}")

    n.add(
        "Link",
        merged["name"].tolist(),
        carrier=["feed_ruminant"] * len(merged),
        bus0=merged["bus0"].tolist(),
        bus1=merged["bus1"].tolist(),
        efficiency=merged["yield"].to_numpy(),
        marginal_cost=[0.0] * len(merged),
        p_nom_max=merged["available_area"].to_numpy(),
        p_nom_extendable=[True] * len(merged),
    )


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


def add_crop_to_feed_links(
    n: pypsa.Network,
    crop_list: list,
    feed_conversion: pd.DataFrame,
    countries: list,
) -> None:
    """Add links that turn crops into feed for ruminants and monogastrics."""

    if feed_conversion.empty:
        logger.info("No feed conversion data provided; skipping crop→feed links")
        return

    df = feed_conversion.copy()
    if "feed_type" not in df.columns or "efficiency" not in df.columns:
        raise ValueError(
            "feed_conversion must contain feed_type and efficiency columns"
        )

    df["feed_type"] = df["feed_type"].str.strip().str.lower()
    valid_types = {"ruminant", "monogastric"}
    invalid_types = sorted(set(df["feed_type"]) - valid_types)
    if invalid_types:
        logger.warning(
            "Ignoring unsupported feed types in feed_conversion: %s",
            ", ".join(invalid_types),
        )
        df = df[df["feed_type"].isin(valid_types)]

    for crop in crop_list:
        crop_rows = df[df["crop"] == crop]
        if crop_rows.empty:
            logger.warning("No feed conversion data for crop '%s', skipping", crop)
            continue

        for _, row in crop_rows.iterrows():
            efficiency = float(row["efficiency"])
            if not np.isfinite(efficiency) or efficiency <= 0:
                logger.warning(
                    "Invalid feed efficiency %.3f for crop '%s' (%s), skipping",
                    efficiency,
                    crop,
                    row["feed_type"],
                )
                continue

            feed_suffix = (
                "ruminant" if row["feed_type"] == "ruminant" else "monogastric"
            )
            names = [
                f"convert_{crop}_to_{feed_suffix}_feed_{country}"
                for country in countries
            ]
            bus0 = [f"crop_{crop}_{country}" for country in countries]
            bus1 = [f"feed_{feed_suffix}_{country}" for country in countries]
            n.add(
                "Link",
                names,
                bus0=bus0,
                bus1=bus1,
                efficiency=[efficiency] * len(countries),
                marginal_cost=[0.01] * len(countries),
                p_nom_extendable=[True] * len(countries),
            )


def add_feed_to_animal_product_links(
    n: pypsa.Network,
    animal_products: list,
    feed_requirements: pd.DataFrame,
    countries: list,
) -> None:
    """Add links that convert feed pools into animal products."""

    if not animal_products:
        logger.info("No animal products configured; skipping feed→animal links")
        return

    if feed_requirements.empty:
        logger.warning("Animal products configured but feed requirement data is empty")
        return

    df = feed_requirements.copy()
    if "feed_type" not in df.columns or "efficiency" not in df.columns:
        raise ValueError(
            "feed_to_animal_products must contain feed_type and efficiency columns"
        )

    df["feed_type"] = df["feed_type"].str.strip().str.lower()
    df["product"] = df["product"].str.strip()
    df = df[df["product"].isin(animal_products)]
    if df.empty:
        logger.warning(
            "No feed requirement data for configured animal products: %s",
            ", ".join(sorted(animal_products)),
        )
        return

    valid_types = {"ruminant", "monogastric"}
    invalid_types = sorted(set(df["feed_type"]) - valid_types)
    if invalid_types:
        logger.warning(
            "Ignoring unsupported feed types in feed_to_animal_products: %s",
            ", ".join(invalid_types),
        )
        df = df[df["feed_type"].isin(valid_types)]

    for _, row in df.iterrows():
        efficiency = float(row["efficiency"])
        if not np.isfinite(efficiency) or efficiency <= 0:
            logger.warning(
                "Invalid feed efficiency %.3f for product '%s' (%s), skipping",
                efficiency,
                row["product"],
                row["feed_type"],
            )
            continue

        feed_suffix = "ruminant" if row["feed_type"] == "ruminant" else "monogastric"
        product = row["product"]
        clean_name = product.replace(" ", "_").replace("(", "").replace(")", "")
        names = [
            f"produce_{clean_name}_from_{feed_suffix}_{country}"
            for country in countries
        ]
        bus0 = [f"feed_{feed_suffix}_{country}" for country in countries]
        bus1 = [f"food_{product}_{country}" for country in countries]
        n.add(
            "Link",
            names,
            bus0=bus0,
            bus1=bus1,
            efficiency=[efficiency] * len(countries),
            marginal_cost=[0.02] * len(countries),
            p_nom_extendable=[True] * len(countries),
        )


def add_food_group_buses_and_loads(
    n: pypsa.Network,
    food_group_list: list,
    food_groups: pd.DataFrame,
    food_groups_config: dict,
    countries: list,
    population: pd.Series,
) -> None:
    """Add carriers, buses, and loads for food groups defined in the CSV."""
    # Add loads for food groups with requirements
    if food_groups_config:
        logger.info("Adding food group loads based on nutrition requirements...")
        for group in food_group_list:
            if group in food_groups_config:
                group_config = food_groups_config[group]
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
    n: pypsa.Network,
    macronutrients_config: dict,
    countries: list,
    population: pd.Series,
) -> None:
    """Add per-country loads for macronutrients based on minimum requirements."""
    if macronutrients_config:
        logger.info("Adding macronutrient loads per country based on requirements...")
        for nutrient in ["carb", "protein", "fat"]:
            if nutrient in macronutrients_config:
                nutrient_config = macronutrients_config[nutrient]
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
    trade_config: dict,
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

    Marginal cost is distance-dependent and controlled by
    config["trade"]["crop_trade_marginal_cost_per_km"].
    Number of hubs from config["trade"]["crop_hubs"].
    """
    n_hubs = int(trade_config["crop_hubs"])
    cost_per_km = float(trade_config["crop_trade_marginal_cost_per_km"])

    if len(regions_gdf) == 0 or len(countries) == 0:
        logger.info("Skipping trade hubs: no regions/countries available")
        return

    if len(crop_list) == 0:
        logger.info("Skipping trade hubs: no crops configured")
        return

    non_tradable = {str(crop) for crop in trade_config.get("non_tradable_crops", [])}
    tradable_crops = [crop for crop in crop_list if crop not in non_tradable]
    skipped = sorted(non_tradable.intersection(set(crop_list)))
    if skipped:
        logger.info(
            "Skipping trade network for non-tradable crops: %s",
            ", ".join(skipped),
        )

    if not tradable_crops:
        logger.info("Skipping trade hubs: no tradable crops available")
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
    for crop in tradable_crops:
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
    for crop in tradable_crops:
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

        for crop in tradable_crops:
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


if __name__ == "__main__":
    read_csv = functools.partial(pd.read_csv, comment="#")

    # Read crop data
    crops = read_csv(snakemake.input.crops, index_col=["crop", "param"])

    # Read food conversion data
    foods = read_csv(snakemake.input.foods)

    # Read food groups data
    food_groups = read_csv(snakemake.input.food_groups)

    # Read nutrition data
    nutrition = read_csv(snakemake.input.nutrition, index_col=["food", "nutrient"])

    # Read feed conversion data (crop -> feed pools)
    feed_conversion = read_csv(snakemake.input.feed_conversion)

    # Read feed requirements for animal products (feed pools -> foods)
    feed_to_products = read_csv(snakemake.input.feed_to_products)

    # Read yields data for each configured crop and water supply (skip missing)
    yields_data = {}
    for crop in snakemake.params.crops:
        for ws in ["i", "r"]:
            yields_key = f"{crop}_yield_{ws}"
            try:
                path = snakemake.input[yields_key]
            except AttributeError:
                logger.warning(
                    "Missing yields input for %s (%s), skipping",
                    crop,
                    "irrigated" if ws == "i" else "rainfed",
                )
                continue
            yields_df = read_csv(path, index_col=["region", "resource_class"])
            yields_data[yields_key] = yields_df
            logger.info(
                "Loaded yields for %s (%s): %d rows",
                crop,
                "irrigated" if ws == "i" else "rainfed",
                len(yields_df),
            )

    # Read regions
    regions_df = gpd.read_file(snakemake.input.regions)

    # Load class-level land areas
    land_class_df = read_csv(snakemake.input.land_area_by_class)
    # Expect columns: region, water_supply, resource_class, area_ha
    land_class_df = land_class_df.set_index(
        ["region", "water_supply", "resource_class"]
    ).sort_index()

    land_rainfed_df = land_class_df.xs("r", level="water_supply").copy()
    grassland_df = pd.DataFrame()
    if snakemake.params.grazing["enabled"]:
        grassland_df = read_csv(
            snakemake.input.grassland_yields, index_col=["region", "resource_class"]
        ).sort_index()

    # Load population per country for planning horizon
    population_df = read_csv(snakemake.input.population)
    # Expect columns: iso3, country, year, population
    # Select only configured countries and validate coverage
    cfg_countries = list(snakemake.params.countries)
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

    # Read FAOSTAT prices (USD/tonne) and build crop->price mapping
    prices_df = read_csv(snakemake.input.prices)
    # Expect columns: crop, faostat_item, n_obs, price_usd_per_tonne
    crop_prices = prices_df.set_index("crop")["price_usd_per_tonne"].astype(float)

    # Build the network (inlined)
    n = pypsa.Network()
    n.set_snapshots(["now"])

    crop_list = snakemake.params.crops
    animal_products_cfg = snakemake.params.animal_products
    animal_product_list = list(animal_products_cfg["include"])

    base_food_list = foods.loc[foods["crop"].isin(crop_list), "food"].unique().tolist()
    food_list = sorted(set(base_food_list).union(animal_product_list))
    food_group_list = food_groups.loc[
        food_groups["food"].isin(food_list), "group"
    ].unique()

    add_carriers_and_buses(n, crop_list, food_list, food_group_list, cfg_countries)
    add_primary_resources(n, snakemake.params.primary)

    # Add class-level land buses and generators (shared pools), replacing region-level caps
    # Apply same regional_limit factor per class pool
    reg_limit = float(snakemake.params.primary["land"]["regional_limit"])
    # Build all unique class buses
    bus_names = [f"land_{r}_class{int(k)}_{ws}" for (r, ws, k) in land_class_df.index]
    n.add("Bus", bus_names, carrier=["land"] * len(bus_names))
    n.add(
        "Generator",
        bus_names,
        bus=bus_names,
        carrier=["land"] * len(bus_names),
        p_nom_extendable=[True] * len(bus_names),
        p_nom_max=(reg_limit * land_class_df["area_ha"]).values,
    )
    add_regional_crop_production_links(
        n,
        crop_list,
        crops,
        yields_data,
        region_to_country,
        set(cfg_countries),
        crop_prices,
    )
    if snakemake.params.grazing.get("enabled", False):
        add_grassland_feed_links(
            n,
            grassland_df,
            land_rainfed_df,
            region_to_country,
            set(cfg_countries),
        )
    add_crop_to_feed_links(n, crop_list, feed_conversion, cfg_countries)
    add_food_conversion_links(n, food_list, foods, cfg_countries)
    add_feed_to_animal_product_links(
        n, animal_product_list, feed_to_products, cfg_countries
    )
    add_food_group_buses_and_loads(
        n,
        food_group_list,
        food_groups,
        snakemake.params.food_groups,
        cfg_countries,
        population,
    )
    add_macronutrient_loads(
        n, snakemake.params.macronutrients, cfg_countries, population
    )
    add_food_nutrition_links(n, food_list, foods, food_groups, nutrition, cfg_countries)

    # Add crop trading hubs and links (hierarchical trade network)
    add_crop_trade_hubs_and_links(
        n, snakemake.params.trade, regions_df, cfg_countries, list(crop_list)
    )

    logger.info("Network summary:")
    logger.info("Carriers: %d", len(n.carriers))
    logger.info("Buses: %d", len(n.buses))
    logger.info("Stores: %d", len(n.stores))
    logger.info("Links: %d", len(n.links))

    n.export_to_netcdf(snakemake.output.network)
