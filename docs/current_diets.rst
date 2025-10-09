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
* **Model constraints**: (Future) Constraining optimization to remain near current diets

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

Age Stratification
------------------

The processed data preserves age stratification from the GDD source, providing dietary intake estimates for seven age groups. This stratification serves multiple purposes:

**Variation across life stages**
  Dietary patterns differ substantially across age groups. For example, dairy consumption is typically highest in early childhood (250-265 g/day for ages 0-10) and lower in adulthood (175 g/day for ages 11-74), reflecting both nutritional needs and cultural feeding practices.

**Energy adjustment**
  The GDD applies age-specific energy adjustment to normalize intakes (700 kcal/day for infants to 2000 kcal/day for adults). This ensures that reported intake values reflect dietary patterns after accounting for differences in total energy consumption across ages.

**Health burden calculation**
  Age-stratified data enables more accurate baseline health burden estimates, as disease risks and mortality rates vary substantially by age. The health module can weight dietary risks appropriately across the age distribution.

**Future extensions**
  Age-stratified baseline data supports planned model features such as age-specific dietary constraints, life-course health dynamics, and demographic transition scenarios.

The ``All ages`` rows provide population-weighted averages useful for simple comparisons and validation against aggregate statistics.

Integration with Health Module
-------------------------------

Current dietary intake data is essential for calculating baseline health burden:

1. **Baseline risk assessment**: GDD provides current intake levels for each dietary risk factor
2. **Relative risk calculation**: Current intake is compared to optimal intake using dose-response curves
3. **Attributable burden**: Disease burden attributable to suboptimal current diet is quantified
4. **Health gains**: Optimization can reduce burden by shifting toward healthier dietary patterns

See :doc:`health` for details on how dietary intake translates to health outcomes.

**Current implementation note**: The health module currently uses the ``All ages`` population-weighted aggregate from the GDD data. Full age-specific matching of dietary intake with age-stratified mortality and morbidity data is planned for future development. The age-stratified dietary data is preserved in the processed output to support this enhancement.

Example: Dairy Consumption
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GDD "Total Milk" variable (v57) represents total dairy consumption in milk equivalents. Age-stratified data shows substantial variation across life stages:

**USA (2018)**
  * 0-10 years: 250-265 g/day (high consumption in childhood)
  * 11-74 years: 175 g/day (adult average)
  * 75+ years: 206 g/day (moderate elderly consumption)
  * All ages: 187 g/day (population average)

**France (2018)**
  * All ages: 328 g/day (high cheese/yogurt consumption culture)

**India (2018)**
  * All ages: 82 g/day (lower but culturally significant)

**China (2018)**
  * All ages: 272 g/day (increasing with economic development)

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

Baseline diet enforcement in the optimization can be toggled via
``config.diet.enforce_gdd_baseline``. When enabled, the builder reads
``processing/{name}/gdd_dietary_intake.csv`` (``All ages`` by default) and adds
per-country equality loads for matching food groups, forcing the solution to
replicate observed intake. ``baseline_age`` and ``baseline_reference_year``
override which cohort/year slice the model locks to.

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

**Age-specific health modeling**
  * Match age-stratified dietary intake with age-specific mortality and morbidity rates
  * Compute age-weighted health burdens rather than using population aggregates
  * Enable life-course health impact analysis and age-targeted interventions
  * Currently: health module uses ``All ages`` aggregate; age-stratified data available for future use

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

.. [GDD2024] Global Dietary Database. Dietary intake data by country, 2018. Tufts University Friedman School of Nutrition Science and Policy. https://www.globaldietarydatabase.org/ (accessed 2025)

.. [Miller2021] Miller V, Singh GM, Onopa J, et al. Global Dietary Database 2017: Data Availability and Gaps on 54 Major Foods, Beverages and Nutrients among 5.6 Million Children and Adults from 1220 Surveys Worldwide. *BMJ Global Health*, 2021;6(2):e003585. https://doi.org/10.1136/bmjgh-2020-003585
