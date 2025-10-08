.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Results & Visualization
========================

Overview
--------

After solving, the model produces results in three formats:

1. **PyPSA network** (``results/{name}/solved/model.nc``): Complete optimization results in NetCDF format
2. **Visualizations** (``results/{name}/plots/*.pdf``): Publication-quality plots and maps
3. **Data tables** (``results/{name}/plots/*.csv``): Tabular exports for custom analysis

This page describes how to access and interpret results.

Results Directory Structure
----------------------------

::

    results/{name}/
    ├── build/
    │   └── model.nc           # Built model before solving
    ├── solved/
    │   └── model.nc           # Solved model with optimal values
    └── plots/
        ├── *.pdf              # Visualizations
        └── *.csv              # Data exports

PyPSA Network Results
---------------------

The solved network (``results/{name}/solved/model.nc``) is an xarray Dataset containing:

Accessing in Python
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pypsa

   n = pypsa.Network("results/toy/solved/model.nc")

   # Access component data
   links_df = n.links  # All links (production, processing, trade)
   buses_df = n.buses  # All buses (crops, foods, nutrients, land)
   stores_df = n.stores  # Resource availability (land, water)

   # Optimal flows
   link_flows = n.links_t.p0  # Power/flow on each link (time series if multi-period)

   # Shadow prices (marginal costs)
   bus_prices = n.buses_t.marginal_price  # Marginal value of each commodity

Key Data Structures
~~~~~~~~~~~~~~~~~~~

**n.links**
  Production, processing, and trade links with optimal flows. Columns:

  * ``bus0``, ``bus1``, ``bus2``, ...: Connected buses
  * ``p_nom_opt``: Optimal capacity (usually same as ``p_nom``)
  * ``p0``: Optimal flow (from ``bus0`` to ``bus1``)
  * ``efficiency``, ``efficiency2``, ...: Conversion factors
  * ``marginal_cost``: Cost per unit flow (USD/t, USD/Mcal, etc.)

**n.buses**
  Commodity buses (crops, foods, nutrients) with prices. Columns:

  * ``carrier``: Commodity type (e.g., ``crop_wheat``, ``nutrient_protein``)
  * ``marginal_price``: Shadow price (USD/unit) — economic value of one more unit

**n.stores**
  Resource stores (land, water, fertilizer) with usage. Columns:

  * ``e_nom``: Total capacity (Mha for land, km³ for water, kg for fertilizer)
  * ``e_initial``: Available amount
  * ``e``: Amount used (after solving)

**n.global_constraints**
  System-wide limits (total fertilizer, emissions caps, nutritional requirements).

Example Analysis
~~~~~~~~~~~~~~~~

**Total crop production by crop**:

.. code-block:: python

   crop_production = (
       n.links[n.links.index.str.contains("production_crop")]
       .groupby(lambda x: x.split("_")[2])  # Extract crop name
       ["p0"]
       .sum()
   )
   print(crop_production)

**Regional net trade** (exports - imports):

.. code-block:: python

   trade_links = n.links[n.links.index.str.contains("trade")]
   # Sum outflows (exports) and inflows (imports) per region
   # (Requires parsing link names to extract regions)

**Water shadow prices** (economic value of water):

.. code-block:: python

   water_buses = n.buses[n.buses.carrier == "water"]
   water_prices = water_buses["marginal_price"]
   print(water_prices)

Visualization Outputs
---------------------

Production and Resource Use
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**crop_production.pdf / .csv**
  Bar chart of total crop production (tonnes) by crop, split by rainfed/irrigated.

**food_production.csv**
  Total food production (tonnes) by food product.

**resource_usage.pdf**
  Multi-panel plot showing:

  * Total land use (Mha) by region/class
  * Water use (km³) by region
  * Fertilizer use (Mt) globally

**objective_breakdown.pdf / .csv**
  Stacked bar chart decomposing total objective value into:

  * Production costs
  * Trade costs
  * Environmental costs (emissions × carbon price)
  * Health costs (YLL × value_per_yll)

Spatial Maps
~~~~~~~~~~~~

**regions_map.pdf**
  Choropleth map of optimization regions (colored by region ID or some metric).

**resource_classes_map.pdf**
  Map showing spatial distribution of resource classes (color-coded by class number).

**crop_production_map.pdf**
  Map of total crop production (tonnes) per region.

**crop_land_use_map.pdf**
  Map of total cropland area (hectares) per region.

**cropland_fraction_map.pdf**
  Map showing cropland as fraction of total land area (%).

**irrigated_cropland_fraction_map.pdf**
  Map showing irrigated cropland as fraction of total cropland (%).

**water_value_map.pdf**
  Map of water shadow prices (USD/m³) — regions with high values face water scarcity.

Dietary and Health Outcomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**food_consumption.pdf**
  Per-capita consumption (g/person/day) by food group, compared to requirements.

**health_risk_map.pdf**
  Map of dietary risk-attributable DALYs per capita by region.

**health_baseline_map.pdf**
  Map of baseline (pre-optimization) health burden for comparison.

**health_risk_by_region.csv / health_baseline_by_region.csv**
  Tabular exports of health outcomes by region.

Crop Use Breakdown
~~~~~~~~~~~~~~~~~~

**crop_use_breakdown.pdf / .csv**
  Stacked bar chat showing how crops are allocated:

  * Direct food consumption
  * Animal feed
