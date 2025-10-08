.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Environmental Impacts
=====================

Overview
--------

The environmental module accounts for greenhouse gas emissions, land use change, water consumption, and nitrogen pollution from food production. These impacts are monetized and included in the objective function via configurable prices/penalties.

Greenhouse Gas Emissions
-------------------------

The model tracks three major greenhouse gases using 100-year global warming potentials (GWP100):

* **CO₂** (GWP = 1): From land use change, fuel combustion
* **CH₄** (GWP = 28): From enteric fermentation (ruminants), rice paddies, manure
* **N₂O** (GWP = 265): From nitrogen fertilizer application, manure

All emissions are aggregated to CO₂-equivalent (tCO₂-eq) for carbon pricing.

Sources of Emissions
~~~~~~~~~~~~~~~~~~~~

**Crop Production**:
  * N₂O from fertilizer (direct and indirect)
  * CH₄ from flooded rice cultivation
  * CO₂ from machinery/fuel (if included)

**Livestock**:
  * CH₄ from enteric fermentation (ruminants)
  * N₂O and CH₄ from manure management
  * CO₂ from feed production (indirect)

**Land Use Change**:
  * CO₂ from clearing vegetation (forest, grassland → cropland)
  * Soil carbon losses

Carbon Pricing
~~~~~~~~~~~~~~

Emissions are priced at a configurable rate:

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: emissions ---
   :end-before: # --- section: crops ---

This creates an economic incentive to reduce emissions, making the model prefer:

* Low-emission crops (legumes fix nitrogen, reducing fertilizer needs)
* Efficient livestock systems (poultry over ruminants)
* Avoiding land use change (using existing cropland)

Typical values:

* **50-100 USD/tCO₂-eq**: Current carbon markets
* **200-300 USD/tCO₂-eq**: Social cost of carbon estimates
* **500+ USD/tCO₂-eq**: Stringent climate mitigation scenarios

Land Use Change
---------------

Expanding agriculture onto non-cropland releases stored carbon.

Carbon Storage by Land Cover
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GAEZ provides carbon density estimates for:

* **Forests**: 100-300 tC/ha (varies by forest type)
* **Grasslands**: 50-150 tC/ha (includes soil carbon)
* **Croplands**: 50-80 tC/ha (lower than natural vegetation)

Conversion emissions:

.. math::

   \text{Emissions (tCO₂)} = (\text{C}_\text{before} - \text{C}_\text{after}) \times 3.67 \times \text{area (ha)}

(Factor 3.67 converts tC → tCO₂)

Model Implementation
~~~~~~~~~~~~~~~~~~~~

Land use change is implicit: when the model allocates land beyond current cropland extent, it incurs:

1. **Opportunity cost**: Lost ecosystem services (not currently monetized)
2. **Carbon emissions**: One-time pulse from clearing
3. **Ongoing emissions**: Reduced soil carbon sequestration

Currently, LUC emissions are approximated using:

* Average carbon density difference (forest/grassland → cropland)
* Regional land use change rates (could use historical GAEZ data)

Future refinements could add:

* Spatially-explicit carbon maps
* Distinction between forest, grassland, wetland conversion
* Afforestation/reforestation credits

Land Constraints
~~~~~~~~~~~~~~~~

To limit unsustainable expansion:

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: downloads ---
   :end-before: # --- section: primary ---

This prevents unrealistic conversion of all suitable land.

Water Use
---------

Irrigated crop production consumes blue water (surface and groundwater), which is limited by basin-level availability.

Water Accounting
~~~~~~~~~~~~~~~~

* **Blue water**: River/lake/aquifer withdrawals for irrigation
* **Green water**: Soil moisture from rainfall (rainfed crops)
* **Grey water**: Pollution dilution capacity (not currently modeled)

The model tracks blue water only, as it's the scarce/contested resource.

Basin-Level Constraints
~~~~~~~~~~~~~~~~~~~~~~~

Water availability data (Water Footprint Network, Hoekstra & Mekonnen 2011) provides monthly blue water availability for 405 GRDC basins.

Constraints:

.. math::

   \sum_{\text{irrigated crops in basin}} \text{water use} \leq \text{blue water available}

This ensures irrigation doesn't exceed sustainable withdrawal limits.

Regional Water Stress
~~~~~~~~~~~~~~~~~~~~~

Regions with high water scarcity (low availability per capita or per cropland area) face:

* Limited irrigated production
* Higher reliance on rainfed crops (lower yields)
* Need to import water-intensive crops (e.g., rice, vegetables)

The model's trade network allows water-rich regions to export to water-scarce regions, effectively trading "virtual water."

Nitrogen Pollution
------------------

Fertilizer application causes nitrogen pollution via:

* **Leaching**: NO₃⁻ contaminating groundwater
* **Runoff**: Eutrophication of rivers/lakes
* **Volatilization**: NH₃ → N₂O emissions

Global Fertilizer Limit
~~~~~~~~~~~~~~~~~~~~~~~

To prevent excessive pollution:

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: primary ---
   :end-before: # --- section: emissions ---

This caps total nitrogen-phosphorus-potassium application globally, forcing efficient use.

Nitrogen Use Efficiency (NUE)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Crop-specific fertilizer requirements (in ``data/crops.csv``) implicitly include NUE. More efficient crops (legumes, which fix nitrogen) require less fertilizer.

Future extensions could:

* Track separate N, P, K limits
* Add regional pollution constraints (e.g., Baltic Sea nitrogen targets)
* Credit manure as organic fertilizer

Emissions Factors
-----------------

Emission factors are stored in ``data/crops.csv`` and embedded in processing scripts:

**Crops** (``data/crops.csv``):
  * ``n2o_kg_per_t``: N₂O from fertilizer per tonne crop
  * ``co2_kg_per_t``: CO₂ from machinery/transport
  * ``ch4_kg_per_t``: CH₄ (e.g., for rice)

**Livestock** (in model code):
  * Enteric CH₄: ~100-300 kg/head/year for cattle
  * Manure N₂O: Function of N excretion rate
  * Manure CH₄: Varies by management system (lagoon vs. pasture)

These are based on IPCC Tier 1/Tier 2 emission factors and can be refined with regional data.

Environmental Cost in Objective
-------------------------------

Total environmental cost:

.. math::

   \text{Environmental cost} = (\text{emissions in tCO₂-eq}) \times \text{carbon price}

Plus potential penalties for:

* Water overuse (shadow price on water constraints)
* Excess nitrogen (if modeled)

Scenario Exploration
--------------------

Environmental parameters enable policy analysis:

**Carbon Pricing Impact**

* 0 USD/tCO₂-eq: Baseline, no emission penalty → high-emission solutions
* 50 USD/tCO₂-eq: Moderate penalty → some emission reduction
* 200 USD/tCO₂-eq: Strong penalty → low-emission crops, less ruminant meat
* 500+ USD/tCO₂-eq: Very stringent → minimal animal products, legume-based diets

**Water Scarcity**

* Tighten basin water limits → shift to rainfed crops, trade virtual water
* Increase irrigation efficiency → allow more irrigated production

**Nitrogen Limits**

* Reduce fertilizer cap → more legumes, extensification, lower yields
* Relax cap → intensification, higher environmental footprint

Visualization
-------------

Environmental results can be visualized:

**Objective breakdown**::

    tools/smk results/{name}/plots/objective_breakdown.pdf

Shows emissions costs vs. production costs vs. health costs.

**Resource usage**::

    tools/smk results/{name}/plots/resource_usage.pdf

Plots water use, fertilizer use, land use by region.

**Water shadow prices**::

    tools/smk results/{name}/plots/water_value_map.pdf

Shows economic value of water in each region (where water constraints bind).
