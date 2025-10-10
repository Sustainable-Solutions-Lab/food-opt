.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Environmental Impacts
=====================

Overview
--------

The environmental module accounts for greenhouse gas emissions, land use change and nitrogen pollution from food production. These impacts are monetized and included in the objective function via configurable prices/penalties.

This is currently a work in progress and not all relevant environmental impacts are implemented and monetized yet.

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

Land Use Change
---------------

Expanding agriculture onto non-cropland releases stored carbon. This is yet to be implemented in the model.

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

Crop-specific fertilizer requirements (in ``data/crops.csv``) implicitly include NUE (currently mock data). More efficient crops (legumes, which fix nitrogen) require less fertilizer.
