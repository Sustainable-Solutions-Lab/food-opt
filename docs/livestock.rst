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

This distinction captures the fundamental trade-off between extensive (land-intensive grazing) and intensive (crop-feed) livestock systems.

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
  * Land (grassland resource classes, similar to cropland)
  * Managed grassland yields from ISIMIP LPJmL model

**Process**:
  1. Grassland yields (t dry matter/ha/year) are computed per region and resource class
  2. Feed conversion ratios translate grass biomass → animal products
  3. Land allocation to grazing competes with cropland expansion

**Advantages**:
  * Can utilize marginal land unsuitable for crops
  * Lower input costs (no crop production needed)

**Disadvantages**:
  * Lower productivity per hectare
  * Higher land use and associated emissions

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

**Advantages**:
  * Higher productivity per unit land (intensive)
  * Can locate near demand centers (doesn't require grazing land)

**Disadvantages**:
  * Competes with human food for crops
  * Depends on crop production infrastructure


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

Since grasslands are perennial, there's no growing season constraint—they're available year-round.

Feed Conversion
---------------

The model uses feed conversion ratios to link feed inputs to animal outputs. These are stored in CSV files with mock placeholder data.

data/feed_conversion.csv
~~~~~~~~~~~~~~~~~~~~~~~~~

Maps crops to feed energy/protein content. Columns:

* ``crop``: Crop name (e.g., "maize", "soybean")
* ``energy_MJ_per_kg``: Energy content
* ``protein_g_per_kg``: Protein content
* Other feed characteristics

**Note**: Current values are mock data; replace with vetted livestock nutrition databases (e.g., Feedipedia).

data/feed_to_animal_products.csv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maps feed requirements to animal product yields. Columns:

* ``animal_product``: Product name (e.g., "cattle meat", "dairy")
* ``feed_type``: Type of feed (e.g., "grass", "grain", "protein")
* ``feed_kg_per_kg_product``: Conversion ratio
* ``production_system``: "grazing" or "feed-based"

**Note**: Current values are mock data; replace with zootechnical conversion factors from FAO or academic literature.

Example Conversion
~~~~~~~~~~~~~~~~~~

For cattle meat via feed-based system:

* 7 kg grain + 2 kg protein feed (e.g., soybean meal) → 1 kg beef
* Feed conversion efficiency: ~11% (9 kg feed → 1 kg beef)

For cattle meat via grazing:

* 40 kg grass dry matter → 1 kg beef
* Extensive system with lower efficiency but uses marginal land

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

**Inputs**:
  * ``bus0``, ``bus2``, ...: Crop buses (e.g., maize, soybean)
  * Multiple crops can feed a single animal product

**Outputs**:
  * ``bus1``: Animal product
  * ``bus3``: Emissions (primarily CH₄ from enteric fermentation)

**Efficiency**: Feed conversion ratios (negative for inputs, positive for output)

Competition and Trade-offs
---------------------------

The model captures several key trade-offs:

Land Use
~~~~~~~~

* **Cropland vs. Grassland**: Expanding crops may reduce grazing land, pushing livestock to feed-based systems or imports
* **Resource Classes**: High-quality land is more valuable for crops, so grazing may concentrate on marginal land

Feed vs. Food
~~~~~~~~~~~~~

* **Crop Allocation**: Crops can go to:

  * Direct human consumption (e.g., wheat → bread)
  * Animal feed (e.g., maize → chicken meat)
  * Non-food uses (e.g., biofuel, waste)

* The model optimizes this allocation based on nutritional needs, environmental costs, and production efficiencies

Environmental Impacts
~~~~~~~~~~~~~~~~~~~~~

* **Grazing**: Higher land use, potential for land-use change emissions, CH₄ from ruminants
* **Feed-based**: Lower land per kg product but higher crop system emissions (fertilizer, machinery)

The optimal mix depends on:

* Carbon price (penalizes ruminant CH₄ and land-use change)
* Land availability (limits grazing expansion)
* Nutritional constraints (protein needs may drive legume production → animal feed)

Emissions from Livestock
-------------------------

Livestock production generates significant greenhouse gas emissions:

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

Visualization
-------------

Livestock production results appear in:

**Food production plots**::

    tools/smk results/{name}/plots/food_production.csv

This CSV includes animal products alongside plant-based foods.

**Objective breakdown**::

    tools/smk results/{name}/plots/objective_breakdown.pdf

Shows contributions of livestock emissions to total environmental costs.

Future Extensions
-----------------

Potential enhancements to the livestock module:

* **Separate grazing types**: Distinguish pasture, rangeland, and mixed systems
* **Manure as fertilizer**: Close nutrient loops by crediting manure N-P-K against synthetic fertilizer
* **Monogastric vs. ruminant**: More detailed emission factor differentiation
* **Dairy vs. beef cattle**: Currently combined; could separate for more realistic herd dynamics
* **Regional suitability**: Some regions better suited for grazing (e.g., arid rangelands) vs. feedlots
