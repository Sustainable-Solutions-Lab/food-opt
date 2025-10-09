.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Data Sources
============

Overview
--------

The model integrates multiple global datasets covering agricultural production, climate, population, health, and water resources. This page documents the key datasets, their licenses, and how to obtain them.

For comprehensive documentation of all datasets, see ``data/DATASETS.md`` in the repository.

.. _manual-download-checklist:

Manual Download Checklist
-------------------------

Several licensed datasets cannot be fetched automatically. Keep the contents of ``data/manually_downloaded`` in sync with the guidance below and the README in that directory.

1. Create an account (free) with IHME and download ``IHME-GBD_2021-dealth-rates.csv`` as described in :ref:`ihme-gbd-mortality`.
2. Download the IHME 2019 relative risk workbook ``IHME_GBD_2019_RELATIVE_RISKS_Y2020M10D15.XLSX`` (:ref:`ihme-relative-risks`).
3. Register at the Global Dietary Database portal and download the dataset, placed locally as the directory ``GDD-dietary-intake`` (:ref:`gdd-dietary-intake`).


Agricultural Production Data
----------------------------

GAEZ (Global Agro-Ecological Zones) v5
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Provider**: FAO/IIASA

**Description**: Global crop suitability and attainable yield estimates under various climate and management scenarios.

**Access**: https://data.apps.fao.org/gaez/; bulk downloads through a Google Cloud Storage interface.

**License**: Creative Commons Attribution 4.0 International (CC BY 4.0) + FAO database terms

**Citation**: FAO/IIASA (2025). Global Agro-Ecological Zones v5 (GAEZ v5).

**Workflow retrieval**: Automatic via Snakemake rules in ``workflow/rules/retrieve.smk``

CROPGRIDS v1.08
~~~~~~~~~~~~~~~

**Provider**: Tang et al., FAO

**Description**: Global harvested and physical crop area maps for 173 crops around 2020 at 0.05° resolution.

**Resolution**: 0.05° × 0.05° (~5.6 km)

**Access**: https://figshare.com/articles/dataset/CROPGRIDS/22491997

**License**: Creative Commons Attribution 4.0 International (CC BY 4.0)

**Citation**: Tang, H., Nguyen, C., Conchedda, G., Casse, L., Tubiello, F. N., & Maggi, F. (2023). CROPGRIDS. *Scientific Data*, 10(1), 1-16.

**Usage**: Yield gap analysis (comparing attainable vs. actual yields)

FAOSTAT Producer Prices
~~~~~~~~~~~~~~~~~~~~~~~~

**Provider**: FAO Statistics Division

**Description**: Crop producer prices by country (2015-2024) in USD/tonne.

**Access**: https://www.fao.org/faostat/en/ (PP domain)

**License**: CC BY 4.0 + FAO database terms

**Retrieval**: Via ``faostat`` Python package (``workflow/scripts/retrieve_faostat_prices.py``)

**Usage**: Calibrating production costs in the objective function

Grassland Yield Data
~~~~~~~~~~~~~~~~~~~~

**Provider**: ISIMIP (Inter-Sectoral Impact Model Intercomparison Project)

**Description**: Historical managed grassland yields from LPJmL model (above-ground dry matter production).

**Resolution**: 0.5° × 0.5°

**Access**: ISIMIP data portal

**Usage**: Grazing-based livestock production potential

Spatial and Administrative Data
--------------------------------

GADM (Global Administrative Areas) v4.1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Provider**: GADM project

**Description**: Global administrative boundary polygons (ADM_0 to ADM_5 levels).

**Format**: GeoPackage with multiple layers

**Access**: https://gadm.org/

**License**: Free for academic/non-commercial use with attribution; redistribution not allowed; commercial use requires permission

**Citation**: GADM (2024). Global Administrative Areas, version 4.1. https://gadm.org/

**Usage**: Building optimization regions via clustering of ADM_1 (states/provinces)

Population Data
---------------

UN World Population Prospects (WPP) 2024
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Provider**: UN DESA Population Division

**Description**: Official UN population estimates and projections by country, age, and sex.

**Variant**: Medium variant projection

**Access**: https://population.un.org/wpp/

**License**: Creative Commons Attribution 3.0 IGO (CC BY 3.0 IGO)

**Files used**:
  * ``WPP2024_TotalPopulationBySex.csv.gz``
  * ``WPP2024_Life_Table_Abridged_Medium_2024-2100.csv.gz``

**Usage**:
  * Scaling per-capita dietary requirements to total demand
  * Age-structured population for health burden calculations
  * Global life expectancy schedule for health loss valuation

Health and Epidemiology Data
-----------------------------

.. _ihme-gbd-mortality:

IHME GBD 2021 — Mortality Rates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Provider**: Institute for Health Metrics and Evaluation (IHME)

**Description**: Cause-specific mortality rates by country, age, and sex from the Global Burden of Disease Study 2021. Used to calculate baseline disease burden attributable to dietary risk factors.

**Query parameters**:
  * Measure: Deaths (Rate per 100,000 population)
  * Causes: Ischemic heart disease, Stroke, Diabetes mellitus, Colon and rectum cancer, Chronic respiratory diseases, All causes
  * Age groups: <1 year, 12-23 months, 2-4 years, 5-9 years, ..., 95+ years (individual age bins)
  * Sex: Both
  * Year: 2021

**License**: Free for non-commercial use with attribution (IHME Free-of-Charge Non-commercial User Agreement)

**Citation**: Global Burden of Disease Collaborative Network. Global Burden of Disease Study 2021 (GBD 2021) Results. Seattle, United States: Institute for Health Metrics and Evaluation (IHME), 2024. Available from https://vizhub.healthdata.org/gbd-results/

**Workflow integration**: Automatically processed via ``workflow/scripts/prepare_gbd_mortality.py``

**Manual download steps**:

1. Visit https://vizhub.healthdata.org/gbd-results/ and sign in with your IHME account.
2. Reproduce the query parameters above by following this permanent link: https://vizhub.healthdata.org/gbd-results?params=gbd-api-2021-permalink/90f3c59133738e4b70b91072b6fd0db4
3. Export the results as CSV (allow some time for the IHME to process the query) and save to ``data/manually_downloaded``. Rename the file to ``IHME-GBD_2021-dealth-rates.csv`` to match the name expected by the Snakemake workflow.

.. _ihme-relative-risks:

IHME GBD 2019 — Relative Risk Curves
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Provider**: Institute for Health Metrics and Evaluation (IHME)

**Description**: Appendix Table 7a from the Global Burden of Disease Study 2019, listing relative risks by dietary risk factor, outcome, age, and exposure level.

**License**: Free for non-commercial use with attribution (IHME Free-of-Charge Non-commercial User Agreement)

**Citation**: Global Burden of Disease Collaborative Network. Global Burden of Disease Study 2019 (GBD 2019) Results. Seattle, United States of America: Institute for Health Metrics and Evaluation (IHME), 2020.

**Workflow integration**: Automatically processed via ``workflow/scripts/prepare_relative_risks.py``

**Manual download steps**:

1. Navigate to https://ghdx.healthdata.org/record/ihme-data/gbd-2019-relative-risks.
2. Under the Files tab, locate and download the "Relative risks: all risk factors except for ambient air pollution, alcohol, smoking, and temperature [XLSX]" file; it will be named ``IHME_GBD_2019_RELATIVE_RISKS_Y2020M10D15.XLSX``. Log in to your IHME account when requested.
3. Place the downloaded file under ``data/manually_downloaded``; no need to rename.

.. _gdd-dietary-intake:

Global Dietary Database (GDD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Provider**: Tufts University Friedman School of Nutrition Science and Policy

**Description**: Country-level estimates of dietary intake for major food groups and dietary risk factors based on systematic review and meta-analysis of national dietary surveys.

**License**: Free for non-commercial research, teaching, and private study with attribution. Data may not be redistributed or used commercially without Tufts permission.

**Citation**: Global Dietary Database. Dietary intake data by country. https://www.globaldietarydatabase.org/ [Accessed YYYY-MM-DD].

**Workflow integration**: Automatically processed via ``workflow/scripts/prepare_gdd_dietary_intake.py``

**Manual download steps**:

1. Create or sign in to a Global Dietary Database account at https://globaldietarydatabase.org/data-download.
2. When you are signed in, navigate back to the download page, accept the terms and proceed to download the GDD dataset, which will be ~1.6GB zip file.
3. Extract the zip file; you will get a directory named ``GDD_FinalEstimates_01102022``
4. Move this directory to ``data/manually_downloaded`` and rename the directory to ``GDD-dietary-intake``.

Water Resources Data
--------------------

Water Footprint Network — Monthly Blue Water Availability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Provider**: Water Footprint Network (Hoekstra & Mekonnen)

**Description**: Monthly blue water availability for 405 GRDC river basins.

**Format**: Shapefile + Excel workbook

**Access**: https://www.waterfootprint.org/resources/appendix/Report53_Appendix.zip

**License**: No explicit license; citation requested (see below)

**Citation**: Hoekstra, A.Y. and Mekonnen, M.M. (2011). *Global water scarcity: monthly blue water footprint compared to blue water availability for the world's major river basins*, Value of Water Research Report Series No. 53, UNESCO-IHE, Delft, Netherlands.

**Usage**: Constraining irrigated crop production by basin-level water availability

Mock and Placeholder Data
--------------------------

Several CSV files in ``data/`` currently contain **mock placeholder values** and must be replaced with sourced data before publication-quality analysis:

data/foods.csv
~~~~~~~~~~~~~~

**Status**: Mock data

**Description**: Food product definitions and processing relationships

**Needs**: Sourced from food composition databases (e.g., USDA FoodData Central)

data/food_groups.csv
~~~~~~~~~~~~~~~~~~~~

**Status**: Mock data

**Description**: Mapping of foods to dietary food groups

**Needs**: Consistent classification scheme (e.g., USDA food groups, WHO recommendations)

data/nutrition.csv
~~~~~~~~~~~~~~~~~~

**Status**: Mock data

**Description**: Nutritional composition of foods (macronutrients, micronutrients)

**Needs**: USDA FoodData Central, FAO INFOODS, or national food composition tables

data/feed_conversion.csv
~~~~~~~~~~~~~~~~~~~~~~~~~

**Status**: Mock data

**Description**: Crop nutrient content for animal feed

**Needs**: Feedipedia (https://www.feedipedia.org/) — comprehensive livestock feed database

data/feed_to_animal_products.csv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status**: Mock data

**Description**: Feed-to-product conversion ratios for livestock

**Needs**: FAO livestock production data, academic livestock science literature

Data License Summary
--------------------

Most datasets used in this project require attribution. Some disallow redistribution, meaning that food-opt cannot be distributed together with these datasets. Some furthermore prohibit commercial use without prior agreement or a paid-for license.

* **CC BY 4.0** (GAEZ, CROPGRIDS, FAOSTAT): Requires attribution
* **CC BY 3.0 IGO** (UN WPP): Requires attribution to UN
* **Academic use only** (GADM, GBD, GDD): Commercial use requires permission or paid licensed.
