.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Current Diets
=============

Overview
--------

The model uses empirical dietary intake data from the Global Dietary Database (GDD) [GDD2024]_ [Miller2021]_ to represent current consumption patterns. This baseline data serves multiple purposes:

* **Health impact assessment**: Calculating disease burden attributable to current dietary patterns
* **Baseline reference**: Comparing optimized diets against current consumption
* **Model constraints**: Optionally constrain the optimization to remain near current diets

Data Source
-----------

**Global Dietary Database (GDD)**
  * **Provider**: Tufts University Friedman School of Nutrition Science and Policy
  * **Coverage**: 185 countries, individual-level dietary surveys (1990-2018)
  * **Variables**: 54 dietary factors including foods, beverages, and nutrients
  * **Download**: Requires free registration at https://globaldietarydatabase.org/data-download
  * **Citation**: [GDD2024]_

The GDD compiles and harmonizes national dietary surveys from around the world using standardized protocols. Data are stratified by age, sex, urban/rural residence, and education level, then aggregated to national-level estimates using population weights.

Weight Conventions
~~~~~~~~~~~~~~~~~~

GDD reports all dietary intake values in **grams per day using "as consumed" weights** [Miller2021]_. This means:

* **Fresh vegetables and fruits**: Reported in fresh weight (e.g., a raw apple, fresh tomato)
* **Grains**: Reported in cooked weight (e.g., cooked rice, prepared bread)
* **Dairy**: Reported as **total milk equivalents**, which includes milk, yogurt, cheese and other dairy products converted to their milk equivalent weight
* **Meats**: Reported in cooked/prepared weight

The model preserves these conventions in the processed output files. Units in the output CSV distinguish between general fresh weight (``g/day (fresh wt)``) and dairy milk equivalents (``g/day (milk equiv)``).

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
     - v01
     - Total fruits (whole fruits only, excluding juices)
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

* Multiple GDD variables can map to a single food group (e.g., starchy_vegetable = v03 potatoes + v04 other starchy veg)
* When aggregating, values are summed within each food group
* The ``dairy`` food group uses v57 "Total Milk", which represents milk equivalents from all dairy consumption including liquid milk, cheese, yogurt, and other dairy products
* The ``fruits`` food group uses only v01 (whole fruits), excluding v16 (fruit juices), to align with the GBD fruit risk factor definition used in health impact modeling

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
3. **Map age groups**: Convert GDD age midpoints to GBD-compatible age buckets (0-1, 1-2, 2-5, 6-10, 11-74, 75+ years)
4. **Aggregate strata**: Compute national averages by age group across sex/education/urban-rural strata
5. **Map to food groups**: Apply the GDD-to-food-group mapping defined in the script
6. **Aggregate variables**: Sum multiple GDD variables that map to the same food group (preserving age stratification)
7. **Handle missing countries**: Apply proxies for territories without separate GDD data
8. **Validate completeness**: Ensure all required countries and food groups are present
9. **Output**: Write ``processing/{name}/gdd_dietary_intake.csv`` with age-stratified data

Output Format
~~~~~~~~~~~~~

The processed dietary intake file has the following structure:

.. code-block:: none

   unit,item,country,age,year,value
   g/day (milk equiv),dairy,USA,0-1 years,2018,252.3
   g/day (milk equiv),dairy,USA,1-2 years,2018,258.3
   g/day (milk equiv),dairy,USA,11-74 years,2018,174.6
   g/day (milk equiv),dairy,USA,All ages,2018,187.1
   g/day (fresh wt),fruits,USA,11-74 years,2018,145.2
   ...

Where:

* ``unit``: Weight convention specific to the food group

  * ``g/day (fresh wt)``: Fresh/cooked "as consumed" weight for most foods
  * ``g/day (milk equiv)``: Total milk equivalents for dairy

* ``item``: Food group name
* ``country``: ISO 3166-1 alpha-3 country code
* ``age``: Age group using GBD-compatible naming

  * ``0-1 years``: Infants under 1 year
  * ``1-2 years``: Toddlers 1-2 years
  * ``2-5 years``: Early childhood 2-5 years
  * ``6-10 years``: Middle childhood 6-10 years
  * ``11-74 years``: Adults 11-74 years
  * ``75+ years``: Elderly 75+ years
  * ``All ages``: Population-weighted average across all age groups

* ``year``: Reference year
* ``value``: Mean daily intake in grams per person for the specified age group

Country Coverage
----------------

The GDD dataset covers 185 countries. For a small number of territories without separate dietary surveys, the model uses proxy data from similar countries:

* **American Samoa (ASM)**: Uses Samoa (WSM) data
* **French Guiana (GUF)**: Uses France (FRA) data
* **Puerto Rico (PRI)**: Uses USA data
* **Somalia (SOM)**: Uses Ethiopia (ETH) data

These proxies are defined in the ``COUNTRY_PROXIES`` dictionary in ``prepare_gdd_dietary_intake.py``.

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

Baseline diet enforcement in the optimization can be toggled via
``config.diet.enforce_gdd_baseline``. When enabled, the builder reads
``processing/{name}/gdd_dietary_intake.csv`` (``All ages`` by default) and adds
per-country equality loads for matching food groups, forcing the solution to
replicate observed intake. ``baseline_age`` and ``baseline_reference_year``
override which cohort/year slice the model locks to.

References
----------

.. [GDD2024] Global Dietary Database. Dietary intake data by country, 2018. Tufts University Friedman School of Nutrition Science and Policy. https://www.globaldietarydatabase.org/ (accessed 2025)

.. [Miller2021] Miller V, Singh GM, Onopa J, et al. Global Dietary Database 2017: Data Availability and Gaps on 54 Major Foods, Beverages and Nutrients among 5.6 Million Children and Adults from 1220 Surveys Worldwide. *BMJ Global Health*, 2021;6(2):e003585. https://doi.org/10.1136/bmjgh-2020-003585
