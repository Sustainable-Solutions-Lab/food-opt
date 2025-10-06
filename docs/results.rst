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
  * Health costs (DALYs × VSLY)

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
  Sankey-style or stacked bar showing how crops are allocated:

  * Direct food consumption
  * Animal feed
  * Processing losses / waste
  * Exports

This reveals inefficiencies or opportunities for better utilization.

Customizing Visualizations
---------------------------

Visualization scripts (``workflow/scripts/plot_*.py``) can be modified:

**Change color schemes**: Edit ``config.yaml`` ``plotting.colors`` section

**Add annotations**: Modify scripts to highlight specific regions/crops

**Export formats**: Change ``.savefig()`` calls to use ``.png``, ``.svg``, etc.

**Custom plots**: Create new scripts following existing patterns::

    # workflow/scripts/plot_my_analysis.py
    import pypsa
    import matplotlib.pyplot as plt

    n = pypsa.Network(snakemake.input.network)

    # Custom analysis code

    plt.savefig(snakemake.output.plot)

Then add a rule in ``workflow/rules/plotting.smk``:

.. code-block:: python

   rule plot_my_analysis:
       input:
           network="results/{name}/solved/model.nc"
       output:
           plot="results/{name}/plots/my_analysis.pdf"
       script:
           "../scripts/plot_my_analysis.py"

Interpreting Key Results
-------------------------

Production Patterns
~~~~~~~~~~~~~~~~~~~

**Questions to explore**:

* Which regions dominate production of key crops?
* Are high-value crops (fruits, vegetables) concentrated in specific regions?
* How much production is rainfed vs. irrigated?

**Interpretation**:

* Concentration suggests comparative advantage (climate, water, soil)
* Diversification indicates resilience but possibly higher costs
* High irrigation share → water dependency

Trade Flows
~~~~~~~~~~~

**Questions to explore**:

* Which regions are net exporters vs. importers?
* Which crops are traded most intensively?
* Do trade patterns match intuition (water-rich exporting to water-scarce)?

**Interpretation**:

* Net exports indicate surplus capacity
* High trade volumes for staples (wheat, maize) suggest specialization
* Virtual water trade alleviates local water scarcity

Resource Constraints
~~~~~~~~~~~~~~~~~~~~

**Questions to explore**:

* Are water constraints binding (high shadow prices)?
* Is the fertilizer limit reached?
* How much of available land is used?

**Interpretation**:

* High water shadow prices → water scarcity drives decisions
* Fertilizer limit binding → may want to increase cap or add organic options
* Low land use → land not the limiting factor (other constraints dominate)

Environmental Outcomes
~~~~~~~~~~~~~~~~~~~~~~

**Questions to explore**:

* What is the total CO₂-equivalent emissions?
* How much comes from livestock vs. crops vs. land use change?
* How do emissions compare to dietary quality (health outcomes)?

**Interpretation**:

* High livestock emissions → opportunity to reduce ruminant meat
* Land use change emissions → pressure to use existing cropland more efficiently
* Trade-off curves (emissions vs. health) reveal policy options

Dietary Quality
~~~~~~~~~~~~~~~

**Questions to explore**:

* Do optimized diets meet food group requirements?
* How much animal protein vs. plant protein?
* Are diets similar to current patterns or radically different?

**Interpretation**:

* Meeting all requirements → feasibility of healthy diets within constraints
* High plant protein → lower environmental impact
* Divergence from current diets → adjustment challenges

Scenario Comparison
-------------------

To compare scenarios (e.g., baseline vs. high carbon price):

1. **Run multiple scenarios**::

       # Scenario 1
       # Edit config: name: "baseline", ghg_price: 0
       tools/smk -j4 all

       # Scenario 2
       # Edit config: name: "high_carbon", ghg_price: 500
       tools/smk -j4 all

2. **Load results**:

   .. code-block:: python

      n_baseline = pypsa.Network("results/baseline/solved/model.nc")
      n_carbon = pypsa.Network("results/high_carbon/solved/model.nc")

3. **Compare metrics**:

   .. code-block:: python

      # Total emissions
      emissions_baseline = (n_baseline.links["p0"] * n_baseline.links["efficiency4"]).sum()
      emissions_carbon = (n_carbon.links["p0"] * n_carbon.links["efficiency4"]).sum()

      print(f"Emission reduction: {(emissions_baseline - emissions_carbon) / emissions_baseline:.1%}")

4. **Visualize differences**:

   .. code-block:: python

      import pandas as pd
      import matplotlib.pyplot as plt

      comparison = pd.DataFrame({
          "baseline": extract_metrics(n_baseline),
          "high_carbon": extract_metrics(n_carbon),
      })

      comparison.T.plot(kind="bar")
      plt.ylabel("Emissions (tCO2-eq)")
      plt.title("Emissions by Scenario")
      plt.show()

Advanced Analysis
-----------------

Sensitivity Analysis
~~~~~~~~~~~~~~~~~~~~

Vary parameters systematically:

* Water availability: ±20%
* Fertilizer limit: 150 Mt, 200 Mt, 250 Mt
* Carbon price: 0, 100, 200, 500 USD/tCO₂-eq

Plot outcome metrics (emissions, costs, land use) vs. parameter values to identify tipping points.

Regional Decomposition
~~~~~~~~~~~~~~~~~~~~~~

Analyze within-country variation:

* Identify high-productivity regions for investment
* Locate water-stressed regions for infrastructure planning
* Map cropland expansion potential

Multi-Criteria Analysis
~~~~~~~~~~~~~~~~~~~~~~~~

Generate Pareto frontiers:

* Emissions vs. Health DALYs
* Emissions vs. Total Costs
* Land Use vs. Water Use

By solving with different objective weight combinations.

Exporting for External Tools
-----------------------------

**GIS analysis**: Export region shapes with attributes::

    import geopandas as gpd

    regions = gpd.read_file("processing/toy/regions.geojson")
    # Add production/emission data as columns
    regions.to_file("results/toy/regions_with_data.geojson")

**Spreadsheet analysis**: Export summary tables::

    crop_production.to_csv("results/toy/exports/crop_production.csv")

**Visualization tools** (Tableau, Power BI): Use exported CSVs

