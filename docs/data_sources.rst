.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Data Sources
============

Overview
--------

The model integrates multiple global datasets covering agricultural production, climate, population, health, and water resources. This page documents the key datasets, their licenses, and how to obtain them.

For comprehensive documentation of all datasets, see ``data/DATASETS.md`` in the repository.

Agricultural Production Data
----------------------------

GAEZ (Global Agro-Ecological Zones) v5
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Provider**: FAO/IIASA

**Description**: Global crop suitability and attainable yield estimates under various climate and management scenarios.

**Resolution**: 0.083° × 0.083° (~9 km) for most variables

**Variables used**:
  * RES05-YCX: Attainable yield on current cropland
  * RES05-SX1: Suitability index (fraction suitable)
  * RES05-WDC: Net irrigation water requirement
  * RES02: Growing season start and length

**Access**: https://gaez.fao.org/

**License**: Creative Commons Attribution 4.0 International (CC BY 4.0) + FAO database terms

**Citation**: FAO/IIASA (2023). Global Agro-Ecological Zones v5 (GAEZ v5). Rome: FAO.

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

**Access**: ISIMIP data portal (requires registration)

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

**File used**: ``WPP2024_TotalPopulationBySex.csv.gz``

**Usage**:
  * Scaling per-capita dietary requirements to total demand
  * Age-structured population for health burden calculations

Health and Epidemiology Data
-----------------------------

DIA (Diet Impact Assessment) Model Inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Provider**: WHO-DIA project (Marco Springmann et al.)

**Description**: Epidemiological inputs for translating dietary risk factors to health burdens (DALYs).

**Source repository**: https://github.com/marco-spr/WHO-DIA

**License**: GPL-3.0

**Files used** (snapshots dated 2021):
  * ``diet_05282021.csv``: Baseline dietary intake by country
  * ``RR_int_05282021.csv``: Relative risk breakpoints
  * ``RR_max_05282021.csv``: Maximum relative risk
  * ``dr_05282021.csv``: Dose-response schedules
  * ``lftable_05282021.csv``: Life tables
  * ``VSL_reg_10182021.csv``: Value of statistical life by region

**Citation**: Springmann, M., et al. (2018). Health and nutritional aspects of sustainable diet strategies and their association with environmental impacts. *Nature Sustainability*, 1(11), 624-632.

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

**Recommended source**: USDA FoodData Central (https://fdc.nal.usda.gov/) — comprehensive, regularly updated, public domain

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

Data Retrieval Workflow
------------------------

Most datasets are downloaded automatically by Snakemake rules in ``workflow/rules/retrieve.smk``:

**GADM**::

    rule retrieve_gadm:
        output: "data/downloads/gadm.gpkg"
        # Downloads via HTTP

**GAEZ**::

    rule retrieve_gaez_yield:
        output: "data/downloads/gaez_yield_{params}.tif"
        # Constructs URLs based on config, downloads via HTTP storage plugin

**UN WPP**::

    rule retrieve_un_population:
        output: "data/downloads/WPP_population.csv.gz"
        # Downloads from UN population portal

**FAOSTAT**::

    rule retrieve_faostat_prices:
        output: "processing/{name}/faostat_prices.csv"
        # Uses faostat Python API

**Manual downloads** (not automated):

* **DIA health data**: Copy files from WHO-DIA repository to ``data/health/raw/``
* **Water Footprint Network**: Download and extract ``Report53_Appendix.zip`` to ``data/downloads/``
* **Grassland yields**: Register with ISIMIP, download LPJmL historical run

Data Storage and Caching
-------------------------

**data/downloads/**
  Raw downloaded datasets (excluded from Git via ``.gitignore``)

**processing/{name}/**
  Processed intermediate files (scenario-specific)

**results/{name}/**
  Model outputs and visualizations

**Caching**: Snakemake tracks file timestamps and checksums, avoiding redundant downloads.

Data License Summary
--------------------

When publishing results, ensure compliance with data licenses:

* **CC BY 4.0** (GAEZ, CROPGRIDS, FAOSTAT): Requires attribution
* **CC BY 3.0 IGO** (UN WPP): Requires attribution to UN
* **GPL-3.0** (DIA health data): Derivative works must be GPL-licensed
* **Academic use only** (GADM): Commercial use requires permission
* **Citation requested** (Water Footprint Network): No explicit license but citation expected

Always check the original data provider websites for the most current terms.

Future Data Enhancements
------------------------

Planned dataset additions:

* **Soil carbon maps**: For more accurate land-use change emissions
* **Livestock feed requirement data**: Replace mock conversion ratios
* **Food processing loss factors**: Industry-specific mass balance data
* **Micronutrient databases**: Iron, zinc, vitamin A content
* **Trade flow data**: Historical bilateral trade for calibration

