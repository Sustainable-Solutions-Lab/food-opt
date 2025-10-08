.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Current Diets
=============

Overview
--------

The model uses empirical dietary intake data from the Global Dietary Database (GDD) to represent current consumption patterns. This baseline data serves multiple purposes:

* **Health impact assessment**: Calculating disease burden attributable to current dietary patterns
* **Baseline reference**: Comparing optimized diets against current consumption
* **Model constraints**: (Future) Constraining optimization to remain near current diets

Data Source
-----------

**Global Dietary Database (GDD)**
  * **Provider**: Tufts University Friedman School of Nutrition Science and Policy
  * **Coverage**: 185 countries, individual-level dietary surveys (1990-2018)
  * **Variables**: 54 dietary factors including foods, beverages, and nutrients
  * **Download**: Requires free registration at https://globaldietarydatabase.org/data-download
  * **Citation**: Global Dietary Database. Dietary intake data by country, 2018. https://www.globaldietarydatabase.org/

The GDD compiles and harmonizes national dietary surveys from around the world using standardized protocols. Data are stratified by age, sex, urban/rural residence, and education level, then aggregated to national-level estimates using population weights.

GDD to Food Group Mapping
--------------------------

The model maps GDD dietary variables to the food groups defined in ``config/food_groups``. This mapping is implemented in ``workflow/scripts/prepare_gdd_dietary_intake.py``.

Food Groups with GDD Data
~~~~~~~~~~~~~~~~~~~~~~~~~~

The following food groups are populated from GDD variables:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Food Group
     - GDD Code
     - Description
   * - ``fruits``
     - v01, v16
     - Total fruits + fruit juices (aggregated)
   * - ``vegetables``
     - v02
     - Non-starchy vegetables
   * - ``starchy_vegetable``
     - v03, v04
     - Potatoes + other starchy vegetables (aggregated)
   * - ``legumes``
     - v05
     - Beans and legumes
   * - ``nuts_seeds``
     - v06
     - Nuts and seeds
   * - ``grain``
     - v07
     - Refined grains (white flour, white rice)
   * - ``whole_grains``
     - v08
     - Whole grains
   * - ``red_meat``
     - v10
     - Unprocessed red meats (cattle, pig)
   * - ``prc_meat``
     - v09
     - Total processed meats
   * - ``fish``
     - v11
     - Total seafoods (fish + shellfish)
   * - ``eggs``
     - v12
     - Eggs
   * - ``dairy``
     - v57
     - Total Milk (includes milk equivalents from all dairy products)

**Notes:**

* Multiple GDD variables can map to a single food group (e.g., fruits = v01 + v16)
* When aggregating, values are summed (e.g., starchy_vegetable = potatoes + other starchy veg)
* The ``dairy`` food group uses v57 "Total Milk", which represents milk equivalents from all dairy consumption including liquid milk, cheese, yogurt, and other dairy products
* This aligns with the GBD dairy risk factor definition used in health impact modeling

Food Groups Without GDD Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some food groups in the model do not have direct GDD mappings:

* ``oil``: Not tracked as a dietary intake in GDD (it's an ingredient/processed product)
* ``poultry``: Not tracked separately in GDD (tracked as part of general meat categories)

These food groups rely on model production and trade without baseline dietary constraints.

Data Processing
---------------

The GDD data processing pipeline (``workflow/scripts/prepare_gdd_dietary_intake.py``) performs the following steps:

1. **Load GDD files**: Read country-level CSV files (``v*_cnty.csv``) for each dietary variable
2. **Filter to reference year**: Extract data for ``config.health.reference_year`` (default: 2018)
3. **Aggregate strata**: Compute national averages across age/sex/education/urban-rural groups
4. **Map to food groups**: Apply the GDD-to-food-group mapping defined in the script
5. **Aggregate variables**: Sum multiple GDD variables that map to the same food group
6. **Handle missing countries**: Apply proxies for territories without separate GDD data
7. **Validate completeness**: Ensure all required countries and food groups are present
8. **Output**: Write ``processing/{name}/gdd_dietary_intake.csv``

Output Format
~~~~~~~~~~~~~

The processed dietary intake file has the following structure:

.. code-block:: csv

   scenario,unit,item,country,year,value
   BMK,g/d_w,dairy,USA,2018,187.1
   BMK,g/d_w,fruits,USA,2018,145.2
   BMK,g/d_w,vegetables,USA,2018,185.3
   ...

Where:

* ``scenario``: "BMK" (baseline)
* ``unit``: "g/d_w" (grams per day, population-weighted)
* ``item``: Food group name
* ``country``: ISO 3166-1 alpha-3 country code
* ``year``: Reference year
* ``value``: Mean daily intake in grams per person

Country Coverage
----------------

The GDD dataset covers 185 countries. For a small number of territories without separate dietary surveys, the model uses proxy data from similar countries:

* **American Samoa (ASM)**: Uses Samoa (WSM) data
* **French Guiana (GUF)**: Uses France (FRA) data
* **Puerto Rico (PRI)**: Uses USA data
* **Somalia (SOM)**: Uses Ethiopia (ETH) data

These proxies are defined in the ``COUNTRY_PROXIES`` dictionary in ``prepare_gdd_dietary_intake.py``.

Integration with Health Module
-------------------------------

Current dietary intake data is essential for calculating baseline health burden:

1. **Baseline risk assessment**: GDD provides current intake levels for each dietary risk factor
2. **Relative risk calculation**: Current intake is compared to optimal intake using dose-response curves
3. **Attributable burden**: Disease burden attributable to suboptimal current diet is quantified
4. **Health gains**: Optimization can reduce burden by shifting toward healthier dietary patterns

See :doc:`health` for details on how dietary intake translates to health outcomes.

Example: Dairy Consumption
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GDD "Total Milk" variable (v57) represents total dairy consumption in milk equivalents:

* **USA**: 187.1 g/day per capita
* **France**: 327.8 g/day per capita (high cheese/yogurt consumption)
* **India**: 82.3 g/day per capita
* **China**: 271.7 g/day per capita

This "Total Milk" metric includes liquid milk, cheese, yogurt, and other dairy products converted to milk equivalents, providing a comprehensive measure of dairy consumption that aligns with the GBD dairy risk factor.

Workflow Integration
--------------------

**Snakemake rule**: ``prepare_gdd_dietary_intake``

**Input**:
  * ``data/manually_downloaded/GDD-dietary-intake/Country-level estimates/*.csv``

**Configuration parameters**:
  * ``config.countries``: List of countries to process
  * ``config.food_groups``: Food group definitions (keys used to filter GDD data)
  * ``config.health.reference_year``: Year for dietary intake data

**Output**:
  * ``processing/{name}/gdd_dietary_intake.csv``

**Script**: ``workflow/scripts/prepare_gdd_dietary_intake.py``

To regenerate dietary intake data:

.. code-block:: bash

   tools/smk --configfile config/default.yaml -- processing/default/gdd_dietary_intake.csv

Validation
----------

The processing script validates:

1. **Country coverage**: All countries in ``config.countries`` must have data (or use proxies)
2. **Food group coverage**: All food groups with GDD mappings must have complete data
3. **Data completeness**: Each country must have values for all mapped food groups

Missing data triggers an error with details about which countries or food groups are incomplete.

Future Extensions
-----------------

Planned enhancements for current diet integration:

**Dietary transition constraints**
  * Limit how far optimized diets can deviate from current patterns
  * Model feasibility of large-scale dietary shifts
  * Account for cultural food preferences and acceptance

**Temporal dynamics**
  * Track dietary trends over time using GDD historical data (1990-2018)
  * Project future dietary patterns under different scenarios
  * Model gradual dietary transitions rather than instantaneous shifts

**Subnational detail**
  * Use GDD stratification (urban/rural, education) for within-country heterogeneity
  * Model dietary inequality and access disparities
  * Target interventions to specific population groups

**Food waste**
  * Distinguish between intake and production (accounting for waste)
  * Use FAO food balance sheets to calibrate waste factors
  * Optimize supply chain efficiency alongside dietary patterns

References
----------

.. [Miller2022] Miller V, Reedy J, Cudhea F, et al. Global, regional, and national consumption of animal-source foods between 1990 and 2018: findings from the Global Dietary Database. *The Lancet Planetary Health*, 2022;6(3):e243-e256. doi:10.1016/S2542-5196(21)00352-1

.. [GDD2024] Global Dietary Database. Dietary intake data by country, 2018. Tufts University Friedman School of Nutrition Science and Policy. https://www.globaldietarydatabase.org/ (accessed 2025)
