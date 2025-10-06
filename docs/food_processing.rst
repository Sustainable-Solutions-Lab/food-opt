.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Food Processing & Trade
========================

Food Processing
---------------

Overview
~~~~~~~~

The food processing module converts raw agricultural products (crops and animal products) into final food products consumed by the population. This captures:

* **Mass losses**: Processing inefficiencies (e.g., wheat → flour loses bran/germ)
* **Multiple inputs**: Some foods require several crops (e.g., bread = flour + oil + salt)
* **Nutritional transformations**: Processing changes nutritional content

Processing is represented in the model as PyPSA links with crop buses as inputs and food buses as outputs.

Data Files
~~~~~~~~~~

**data/foods.csv** (currently mock data)
  Defines food products with columns:

  * ``food``: Food name (e.g., "bread", "pasta", "beef")
  * ``food_group``: Category (e.g., "whole grain", "grain", "animal protein")
  * Processing loss factors
  * Additional metadata

**data/food_groups.csv** (currently mock data)
  Maps foods to food groups for dietary constraint aggregation.

**data/nutrition.csv** (currently mock data)
  Nutritional composition of foods (macronutrients, micronutrients). Columns:

  * ``food``: Food name
  * Nutrient columns (e.g., ``protein_g_per_100g``, ``fat_g_per_100g``, ``kcal_per_100g``)

**Important**: These CSV files contain placeholder values and must be replaced with sourced data (e.g., USDA FoodData Central, FAO food composition databases) before analysis.

Processing Chains
~~~~~~~~~~~~~~~~~

Example: Wheat → Bread

1. **Milling**: Wheat grain → Flour (80% yield)
2. **Baking**: Flour + Water + Oil → Bread (90% yield)
3. **Net**: 1 kg wheat → ~0.72 kg bread

In the model:

* Crop bus ``crop_wheat_{country}`` → Processing link → Food bus ``food_bread_{country}``
* Efficiency: 0.72 (accounting for combined losses)
* Optional additional inputs (oil, salt) via multi-bus links

Example: Soybean → Tofu + Oil

1. **Crushing**: Soybeans → Soy meal (80%) + Soy oil (20%)
2. **Curdling**: Soy meal → Tofu (60% yield)
3. **Net**: 1 kg soybeans → 0.48 kg tofu + 0.20 kg oil

Multiple output processing links allow capturing co-products.

Trade
-----

Overview
~~~~~~~~

The trade module enables inter-regional flows of crops and food products, subject to transport costs. This is essential for:

* **Specialization**: Regions produce what they're best suited for
* **Resource constraints**: Water-scarce regions import water-intensive crops
* **Nutritional diversity**: Small regions import foods they can't produce locally

Hub Network Structure
~~~~~~~~~~~~~~~~~~~~~

To avoid a combinatorial explosion of region-to-region links, the model uses a **hub-based topology**:

1. **Regional buses**: Each region has local crop/food buses
2. **Hub buses**: A small number of hub nodes (configured count)
3. **Hub connections**: Regions connect to nearest hubs; hubs connect to each other

This reduces links from O(n²) to O(n × h + h²), where n = regions and h = hubs.

Configuration
~~~~~~~~~~~~~

.. code-block:: yaml

   trade:
     crop_hubs: 20                          # Number of crop trade hubs
     crop_default_trade_cost_per_km: 1e-2   # Default transport cost (USD/t/km)

     crop_trade_cost_categories:
       bulk_dry_goods:
         cost_per_km: 6e-3                  # Cheaper for bulk (wheat, maize)
         crops:
           - wheat
           - maize
           - soybean
       perishable_high_value:
         cost_per_km: 2.2e-2                # Expensive for fresh produce
         crops:
           - tomato
           - banana

     non_tradable_crops:                    # Local-only commodities
       - alfalfa                            # Fodder crops
       - biomass-sorghum

     animal_product_hubs: 20
     animal_product_default_trade_cost_per_km: 2.1e-2
     animal_product_trade_cost_categories:
       chilled_meat:
         cost_per_km: 2.8e-2
         products:
           - cattle meat
           - pig meat

     non_tradable_animal_products: []

Trade Cost Categories
~~~~~~~~~~~~~~~~~~~~~

Transport costs differentiate by commodity handling requirements:

* **Bulk dry goods** (6e-3 USD/t/km): Cereals, legumes in containers/bulk carriers
* **Bulky fresh** (1.4e-2 USD/t/km): Potatoes, cassava, sugar beets
* **Perishable high-value** (2.2e-2 USD/t/km): Fruits, vegetables, sugarcane requiring refrigeration
* **Chilled meat** (2.8e-2 USD/t/km): Temperature-controlled meat transport

These costs are based on typical freight rates accounting for:

* Handling difficulty (bulk vs. palletized vs. refrigerated)
* Shelf life (perishability)
* Value density (affects insurance, urgency)

Hub Location
~~~~~~~~~~~~

Hub positions are determined by k-means clustering on region centroids:

1. Compute population-weighted centroid for each region
2. Run k-means with k = configured hub count
3. Assign each region to nearest hub
4. Create hub-hub distance matrix for hub-to-hub transport

This ensures hubs are spatially distributed to minimize total transport distance.

Trade Links
~~~~~~~~~~~

Three types of trade links:

1. **Region → Hub**: Local transport from production region to export hub

   * Cost: distance × cost_per_km
   * Efficiency: 1.0 (no loss)

2. **Hub → Hub**: Long-distance transport between hubs

   * Cost: great-circle distance × cost_per_km
   * Efficiency: 1.0

3. **Hub → Region**: Import from hub to consumption region

   * Cost: distance × cost_per_km
   * Efficiency: 1.0

Total transport cost for a region A → region B trade:

.. math::

   \text{Cost} = (\text{dist}(A, \text{hub}_A) + \text{dist}(\text{hub}_A, \text{hub}_B) + \text{dist}(\text{hub}_B, B)) \times c_{\text{km}}

Non-Tradable Commodities
~~~~~~~~~~~~~~~~~~~~~~~~

Certain products are designated non-tradable:

* **Fodder crops** (alfalfa, biomass sorghum): Too bulky/low-value to transport
* **Perishables** (optional): Can restrict local consumption of fragile goods

Non-tradable crops must be consumed (as food or feed) within their production region.

Model Implementation
--------------------

Trade links are created in ``workflow/scripts/build_model.py``:

.. code-block:: python

   # Pseudocode
   for crop in tradable_crops:
       for region in regions:
           hub = nearest_hub(region)
           n.add("Link",
                 f"trade_{crop}_{region}_to_{hub}",
                 bus0=f"crop_{crop}_{region}",
                 bus1=f"crop_{crop}_hub{hub}",
                 p_nom=inf,  # No capacity limit
                 marginal_cost=distance * cost_per_km)

       for hub_i, hub_j in hub_pairs:
           n.add("Link",
                 f"trade_{crop}_hub{hub_i}_to_hub{hub_j}",
                 bus0=f"crop_{crop}_hub{hub_i}",
                 bus1=f"crop_{crop}_hub{hub_j}",
                 p_nom=inf,
                 marginal_cost=hub_distance * cost_per_km)

Similar structure for animal products.

Trade Flow Analysis
-------------------

After solving, trade flows can be analyzed:

* **Net exports/imports by region**: Regions with water/land advantages export; constrained regions import
* **Hub utilization**: Which hubs serve as major trade centers
* **Commodity-specific patterns**: Water-intensive crops flow from water-rich regions

This reveals the spatial structure of optimal food systems under resource constraints.

Scenario Exploration
--------------------

Trade parameters enable exploring policy questions:

**Globalized vs. Localized Food Systems**

* High trade costs → more local production, less specialization
* Low trade costs → comparative advantage drives specialization

**Food Security vs. Efficiency**

* Non-tradable constraints → regions must be self-sufficient (resilient but inefficient)
* Free trade → efficient but vulnerable to disruptions

**Infrastructure Development**

* Increasing crop_hubs or reducing costs → better market integration
* Capturing effects of transport infrastructure investment

Visualization
-------------

Trade flow results can be visualized (future enhancement):

* **Flow maps**: Arrows showing commodity flows between regions
* **Net trade balances**: Choropleth maps of exports (green) vs. imports (red)
* **Hub networks**: Graph visualization of hub connectivity

