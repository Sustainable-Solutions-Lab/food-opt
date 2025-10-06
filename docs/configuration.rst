.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Configuration System
====================

Overview
--------

The food-opt model is configuration-driven: all scenario parameters, crop selections, constraints, and solver options are defined in YAML configuration files under ``config/``. This allows exploring different scenarios without modifying code.

Main Configuration File
-----------------------

The primary configuration is ``config/config.yaml``, structured into thematic sections.

Scenario Identification
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

   name: "toy"  # Scenario name (determines output directory)

Results are saved under ``results/{name}/``, allowing multiple scenarios to coexist.

Planning Horizon
~~~~~~~~~~~~~~~~

.. code-block:: yaml

   planning_horizon: 2030  # Target year for population and climate projections

Matches UN WPP population year and GAEZ climate period.

Download Options
~~~~~~~~~~~~~~~~

.. code-block:: yaml

   downloads:
     show_progress: true  # Display progress bars for dataset downloads

Crop Selection
~~~~~~~~~~~~~~

.. code-block:: yaml

   crops:
     - wheat
     - maize
     - soybean
     # ... ~70 total crops

See :doc:`crop_production` for full list. Add/remove crops to explore specialized vs. diversified production systems.

Country Coverage
~~~~~~~~~~~~~~~~

.. code-block:: yaml

   countries:
     - AFG
     - AGO
     - ALB
     # ... ~190 countries (ISO 3166-1 alpha-3 codes)

Include countries to model; exclude to reduce problem size. Microstate and countries without level-1 GADM data are commented out.

Spatial Aggregation
-------------------

Controls regional resolution and land classification.

.. code-block:: yaml

   aggregation:
     regions:
       type: "cluster"                # Clustering method
       target_count: 400              # Number of optimization regions
       allow_cross_border: false      # Permit cross-border regions
       method: "kmeans"               # k-means or other clustering algorithms
     simplify_tolerance_km: 5         # Geometry simplification tolerance
     simplify_min_area_km: 25         # Remove enclaves smaller than this
     resource_class_quantiles:        # Yield heterogeneity breakpoints
       - 0.25
       - 0.5
       - 0.75
     land_limit_dataset: "suitability"  # "suitability" or "irrigated"

**Trade-offs**:
  * More regions → higher spatial resolution, longer solve time
  * Fewer resource classes → faster solving, less yield heterogeneity
  * ``land_limit_dataset: "irrigated"`` → uniform land base, simpler but less realistic

Primary Resource Constraints
----------------------------

Limits on land, water, and fertilizer availability.

.. code-block:: yaml

   primary:
     land:
       regional_limit: 0.7            # Fraction of potential cropland available (70%)
     fertilizer:
       limit: 2e11                    # kg NPK (200 Mt globally)

Tightening these constraints forces more efficient resource use or extensification.

GAEZ Data Parameters
--------------------

Configures which GAEZ v5 climate scenario and input level to use.

.. code-block:: yaml

   data:
     gaez:
       climate_model: "GFDL-ESM4"         # GCM: GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL, or ENSEMBLE
       period: "FP2140"                   # FP2140 (2021-2040), FP4160 (2041-2060), etc.
       scenario: "SSP126"                 # SSP126 (low), SSP370 (med), SSP585 (high), HIST (historical)
       input_level: "H"                   # H (high inputs), L (low inputs)
       yield_var: "RES05-YCX"            # Attainable yield variable
       water_requirement_var: "RES05-WDC" # Irrigation water requirement
       suitability_var: "RES05-SX1"      # Suitability fraction

**Scenarios**:
  * SSP126: Strong mitigation (1.5-2°C warming)
  * SSP370: Moderate emissions (~3°C)
  * SSP585: High emissions (~4-5°C)

**Input Levels**:
  * H: Modern agriculture (fertilizer, irrigation, pest control)
  * L: Subsistence farming (minimal external inputs)

Irrigation
----------

.. code-block:: yaml

   irrigation:
     irrigated_crops: "all"  # "all" or list of specific crops

Restrict irrigation to water-scarce scenarios or explore rainfed-only production.

Nutritional Requirements
------------------------

Macronutrients
~~~~~~~~~~~~~~

.. code-block:: yaml

   macronutrients:
     carb:
       min: 250         # g/person/day
     protein:
       min: 50          # g/person/day
     fat:
       min: 50          # g/person/day
     kcal:
       equal: 2400      # kcal/person/day (EAT-Lancet 2025)

Use ``min``, ``max``, or ``equal`` constraints.

Food Groups
~~~~~~~~~~~

.. code-block:: yaml

   food_groups:
     whole grain:
       min_per_person_per_day: 50    # g/person/day
     fruit:
       min_per_person_per_day: 50
     vegetable:
       min_per_person_per_day: 50
     animal protein:
       min_per_person_per_day: 30

Increase to promote healthier diets; decrease to relax constraints for faster solving.

Animal Products
---------------

.. code-block:: yaml

   animal_products:
     include:
       - cattle meat
       - pig meat
       - chicken meat
       - dairy
       - eggs

   grazing:
     enabled: true  # Allow grazing-based livestock production

Disable grazing to force intensive feed-based systems.

Trade Configuration
-------------------

.. code-block:: yaml

   trade:
     crop_hubs: 20                              # Number of crop trade hubs
     crop_default_trade_cost_per_km: 1e-2       # USD/t/km

     crop_trade_cost_categories:
       bulk_dry_goods:
         cost_per_km: 6e-3
         crops: [wheat, maize, soybean, ...]
       perishable_high_value:
         cost_per_km: 2.2e-2
         crops: [tomato, banana, ...]

     non_tradable_crops:
       - alfalfa         # Fodder stays local
       - biomass-sorghum

     animal_product_hubs: 20
     animal_product_default_trade_cost_per_km: 2.1e-2

Increase trade costs to explore localized food systems; decrease for globalized trade.

Emissions Pricing
-----------------

.. code-block:: yaml

   emissions:
     ghg_price: 200  # USD/tCO₂-eq

**Values**:
  * 0: No carbon price (baseline)
  * 50-100: Current market prices
  * 200-300: Social cost of carbon
  * 500+: Stringent climate policy

Health Configuration
--------------------

.. code-block:: yaml

   health:
     region_clusters: 30                    # Health cluster count
     reference_year: 2018                   # Baseline health data year
     intake_grid_step: 10                   # g/day resolution for dose-response
     log_rr_points: 10                      # Linearization points for log(RR)
     value_of_statistical_life: 3_500_000   # USD or "regional"
     risk_factors:
       - fruits
       - vegetables
       - nuts_seeds
       - legumes
       - fish
       - red_meat
       - prc_meat
       - whole_grains

Reduce ``region_clusters`` or ``log_rr_points`` to speed up solving.

Solver Configuration
--------------------

.. code-block:: yaml

   solving:
     solver: highs  # or "gurobi"

     options_highs:
       solver: "ipm"          # Interior-point method
       mip_rel_gap: 0.001     # 0.1% optimality gap

     options_gurobi:
       LogToConsole: 0
       OutputFlag: 1
       Method: 2              # Barrier method
       MIPGap: 0.001

**Solver choice**:
  * **HiGHS**: Open-source, fast, good for most problems
  * **Gurobi**: Commercial, often faster for very large problems, requires license

Plotting Configuration
----------------------

.. code-block:: yaml

   plotting:
     colors:
       crops:
         wheat: "#C58E2D"
         maize: "#F1C232"
         soybean: "#7B4F2A"
         # ... color for each crop
     fallback_cmaps:
       crops: "Set3"  # Matplotlib colormap for unconfigured crops

Customize visualization colors for publication-quality plots.

Configuration Workflow
----------------------

Typical workflow for defining a new scenario:

1. **Copy base config**::

       cp config/config.yaml config/my_scenario.yaml

2. **Edit parameters**: Modify crops, constraints, solver options

3. **Update scenario name**:

   .. code-block:: yaml

      name: "my_scenario"

4. **Run workflow**::

       tools/smk -j4 all

5. **Results appear in**: ``results/my_scenario/``

Multiple scenarios can be run in parallel (different terminal sessions) or sequentially by changing the config file.

Configuration Validation
------------------------

The model performs basic validation:

* Missing required keys → error
* Invalid crop names (not in GAEZ mapping) → error during data retrieval
* Inconsistent constraints (e.g., ``kcal.min > kcal.max``) → solver infeasibility

More sophisticated validation (e.g., checking that food groups sum correctly) is future work.

Configuration Best Practices
----------------------------

**Start small**: Use toy config (400 regions, relaxed constraints) for testing

**Scale up gradually**: Increase regions/crops/constraints incrementally

**Document changes**: Comment your config file with scenario rationale

**Version control**: Track config files in Git to reproduce results

**Compare scenarios**: Use consistent naming (``baseline``, ``high_carbon_price``, etc.)

