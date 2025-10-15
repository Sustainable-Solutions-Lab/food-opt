.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Livestock & Grazing
===================

Overview
--------

The livestock module models animal product production (meat, dairy, eggs) through two distinct production systems:

* **Grazing-based**: Animals feed on managed grasslands
* **Feed-based**: Animals consume crops as concentrated feed

Animal Products
---------------

The model includes five major animal product categories configured in ``config/default.yaml``:

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: animal_products ---
   :end-before: # --- section: food_groups ---

Each product can be produced via either production system, with different feed requirements and efficiencies.

Production Systems
------------------

Grazing-Based Production
~~~~~~~~~~~~~~~~~~~~~~~~

**Concept**: Animals graze on managed grasslands, converting grass biomass to animal products.

**Inputs**:
  * Land (per region and resource class, similar to cropland)
  * Managed grassland yields from ISIMIP LPJmL model

**Process**:
  1. Grassland yields (t dry matter/ha/year) are computed per region and resource class
  2. Feed conversion ratios translate grass biomass → animal products
  3. Land allocation to grazing competes with cropland expansion

**Configuration**: Enable/disable grazing with ``grazing.enabled: true``

Feed-Based Production
~~~~~~~~~~~~~~~~~~~~~

**Concept**: Animals consume crops (grains, soybeans, etc.) as concentrated feed.

**Inputs**:
  * Crops from crop production buses
  * Feed conversion ratios (kg crop → kg animal product)

**Process**:
  1. Crops are allocated to animal feed (competing with direct human consumption)
  2. Feed conversion links transform crop inputs to animal products
  3. Multiple crops can contribute (e.g., maize + soybean for poultry)

.. _grassland-yields:

Grassland Yields
----------------

Grazing supply is determined by managed grassland yields from the ISIMIP LPJmL historical simulation.

Data Source
~~~~~~~~~~~

**Dataset**: ISIMIP2b managed grassland yields (historical)

**Resolution**: 0.5° × 0.5° gridded annual yields

**Variable**: Above-ground dry matter production (t/ha/year)

**Processing**: ``workflow/scripts/build_grassland_yields.py``

Aggregation follows the same resource class structure as crops:

1. Load grassland yield NetCDF
2. Aggregate by (region, resource_class) using area-weighted means
3. Output CSV with yields in t/ha/year

.. figure:: _static/figures/grassland_yield.svg
   :alt: Managed grassland yield potential
   :width: 100%
   :align: center

   Global distribution of managed grassland yield potential (tonnes dry matter per hectare per year) from ISIMIP LPJmL historical simulations

Feed Conversion
---------------

The model uses feed conversion ratios to link feed inputs to animal outputs. These are stored in CSV files with mock placeholder data.

data/feed_conversion.csv
~~~~~~~~~~~~~~~~~~~~~~~~~

Maps crops to feed energy/protein content for direct crop-to-feed conversion. Columns:

* ``crop``: Crop name (e.g., "maize", "soybean")
* ``feed_type``: Feed category ("ruminant" or "monogastric")
* ``efficiency``: Digestible fraction of the crop as feed
* ``notes``: Description and source information

This file enables crops to be used directly as animal feed, competing with human consumption and food processing pathways.

.. Note:: Current values are mock data; to be replaced by actual values.

.. _byproduct-feed-conversion:

data/food_feed_conversion.csv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maps food byproducts to feed conversion efficiencies. Columns:

* ``food``: Food byproduct name (e.g., "wheat bran", "sunflower meal")
* ``feed_type``: Feed category ("ruminant" or "monogastric")
* ``efficiency``: Digestible fraction of the byproduct as feed
* ``notes``: Description and source information

This file enables byproducts from food processing (assigned to the ``byproduct`` food group) to be used as animal feed. Byproducts are automatically excluded from human consumption and routed to animal feed systems.

**Typical byproducts**:
  * **Cereal brans**: wheat bran, rice bran, barley bran, oat bran (moderate fiber, good for ruminants)
  * **Oilseed meals**: sunflower meal, rapeseed meal (high protein after oil extraction)
  * **Other**: wheat germ (high protein/fat), buckwheat hulls (high fiber, lower digestibility)

**Feed efficiency differences**: Ruminants generally have higher efficiency for fibrous byproducts due to their multi-stomach digestive system, while monogastrics (pigs, poultry) have lower efficiency for high-fiber materials but may utilize protein-rich meals effectively.

.. Note:: Current values are mock data; to be replaced with actual feed value data from animal nutrition literature.

data/feed_to_animal_products.csv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maps feed requirements to animal product yields. Columns:

* ``animal_product``: Product name (e.g., "cattle meat", "dairy")
* ``feed_type``: Type of feed (e.g., "grass", "grain", "protein")
* ``feed_kg_per_kg_product``: Conversion ratio
* ``production_system``: "grazing" or "feed-based"

.. Note:: Current values are mock data; to be replaced by actual values.

Model Implementation
--------------------

In ``workflow/scripts/build_model.py``, livestock production is represented as multi-bus links:

Grazing Links
~~~~~~~~~~~~~

**Inputs**:
  * ``bus0``: Grassland (land bus for region/class)
  * ``bus2``: Primary resources (water, if constrained)

**Outputs**:
  * ``bus1``: Animal product (to animal product bus)
  * ``bus3``: Emissions (CH₄ from enteric fermentation, N₂O from manure)

**Efficiency**: Grassland yield (t/ha) × feed conversion (t grass → t product)

Feed-Based Links
~~~~~~~~~~~~~~~~

The feed-based production system utilizes two types of feed sources:

**Crop-to-Feed Links**:
  * **Inputs**: Crop buses (e.g., maize, soybean from ``bus0``)
  * **Outputs**: Feed pool buses (``feed_ruminant`` or ``feed_monogastric``)
  * **Efficiency**: Crop digestibility as feed (from ``data/feed_conversion.csv``)
  * Crops compete between human consumption, food processing, and animal feed use

**Byproduct-to-Feed Links**:
  * **Inputs**: Food byproduct buses (e.g., wheat bran, sunflower meal from ``bus0``)
  * **Outputs**: Feed pool buses (``feed_ruminant`` or ``feed_monogastric``)
  * **Efficiency**: Byproduct digestibility as feed (from ``data/food_feed_conversion.csv``)
  * Byproducts from food processing are automatically excluded from human consumption and can only be used as feed

**Feed-to-Animal-Product Links**:
  * **Inputs**: Feed pool buses (``feed_ruminant`` or ``feed_monogastric`` from ``bus0``)
  * **Outputs**: Animal product buses (``bus1``)
  * **Emissions**: CH₄ from enteric fermentation (``bus3``)
  * **Efficiency**: Feed conversion ratios (kg feed → kg product)

Emissions from Livestock
-------------------------

Livestock production generates significant greenhouse gas emissions.

.. Note:: The model currently uses mock data for these emissions.

Enteric Fermentation (CH₄)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ruminants (cattle, sheep) produce methane through digestive fermentation. Emission factors:

* **Cattle**: ~100-300 kg CH₄/head/year (varies by diet and system)
* **Dairy**: Higher per animal but lower per liter milk than meat
* **Pigs/Poultry**: Minimal enteric fermentation

Manure Management (N₂O, CH₄)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Manure storage and application releases:

* **N₂O**: From nitrogen in manure (direct and indirect emissions)
* **CH₄**: From anaerobic manure decomposition (especially in lagoons)

These are incorporated into the production link efficiencies, priced at the configured ``emissions.ghg_price`` (USD/tCO₂-eq).

Configuration Parameters
------------------------

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: animal_products ---
   :end-before: # --- section: trade ---

Disabling grazing (``enabled: false``) forces all animal products to come from feed-based systems or imports, useful for exploring intensification scenarios. All food group minima are zero by default; raise ``food_groups.animal_protein.min_per_person_per_day`` (e.g., to 30 g) to enforce minimum consumption of animal-source foods.

Workflow Rules
--------------

**build_grassland_yields**
  * **Input**: ISIMIP grassland yield NetCDF, resource classes, regions
  * **Output**: ``processing/{name}/grassland_yields.csv``
  * **Script**: ``workflow/scripts/build_grassland_yields.py``

Livestock production is then integrated into the ``build_model`` rule using the grassland yields and feed conversion CSVs.
