.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Configuration
=============

Overview
--------

The food-opt model is configuration-driven: all scenario parameters, crop selections, constraints, and solver options are defined in YAML configuration files under ``config/``. This allows exploring different scenarios without modifying code.

The default configuration is ``config/default.yaml``, structured into thematic sections.

Custom configuration files
~~~~~~~~~~~~~~~~~~~~~~~~~~

Instead of modifying the default configuration file, it is recommended to explore individual scenarios by creating named configuration files, overriding specific parts of the default configuration. Such a named configuration file must contain at the minimum a ``name``. An example could be something like the following: 

.. code-block:: yaml

   # config/my_scenario.yaml
   name: "my_scenario"           # Scenario name → results/my_scenario/
   planning_horizon: 2040        # Override the default 2030 horizon
   primary:
     land:
       regional_limit: 0.6       # Tighten land availability
   emissions:
     ghg_price: 250              # Raise the carbon price above the default

Any keys omitted in your custom file fall back to the defaults shown in the sections below, so you can keep overrides concise.

Results are saved under ``results/{name}/``, allowing multiple scenarios coming from different configuration files to coexist.

To build and solve the model based on the above example configuration, you would run the following::

  tools/smk -j4 --configfile config/my_scenario.yaml

Configuration sections
----------------------

Planning Horizon
~~~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: scenario_metadata ---
   :end-before: # --- section: downloads ---

Matches UN WPP population year and GAEZ climate period.

Download Options
~~~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: downloads ---
   :end-before: # --- section: primary ---

Crop Selection
~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: crops ---
   :end-before: # --- section: macronutrients ---

See :doc:`crop_production` for full list. Add/remove crops to explore specialized vs. diversified production systems.

Country Coverage
~~~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: countries ---
   :end-before: # --- section: data ---

Include countries/territories to model; exclude to reduce problem size. Microstate and countries missing essential data are commented out.

Spatial Aggregation
~~~~~~~~~~~~~~~~~~~

Controls regional resolution and land classification.

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: aggregation ---
   :end-before: # --- section: countries ---

**Trade-offs**:
  * More regions → higher spatial resolution, longer solve time
  * Fewer resource classes → faster solving, less yield heterogeneity

Primary Resource Constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Limits on land, water, and fertilizer availability.

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: primary ---
   :end-before: # --- section: emissions ---

GAEZ Data Parameters
~~~~~~~~~~~~~~~~~~~~

Configures which GAEZ v5 climate scenario and input level to use.

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: data ---
   :end-before: # --- section: irrigation ---

**Scenarios**:
  * SSP126: Strong mitigation (1.5-2°C warming)
  * SSP370: Moderate emissions (~3°C)
  * SSP585: High emissions (~4-5°C)

**Input Levels**:
  * H: Modern agriculture (fertilizer, irrigation, pest control)
  * L: Subsistence farming (minimal external inputs)

Irrigation
~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: irrigation ---
   :end-before: # --- section: solving ---

Restrict irrigation to water-scarce scenarios or explore rainfed-only production.

Macronutrients
~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: macronutrients ---
   :end-before: # --- section: animal_products ---

Use ``min``, ``max``, or ``equal`` constraints.

Food Groups
~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: food_groups ---
   :end-before: # --- section: diet ---

Each food group may specify ``min_per_person_per_day``, ``max_per_person_per_day``,
and ``equal_per_person_per_day``. The defaults leave minima at zero so food group
constraints stay inactive; tighten minima or maxima to guide intakes, or use the
``equal`` field for equality targets.

Diet Controls
~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: diet ---
   :end-before: # --- section: trade ---

Enable ``enforce_gdd_baseline`` to force the optimization to match baseline
consumption from the processed GDD file. Override ``baseline_age`` or
``baseline_reference_year`` if you pre-process alternative cohorts or years.

Animal Products
~~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: animal_products ---
   :end-before: # --- section: food_groups ---

Disable grazing to force intensive feed-based systems.

Trade Configuration
~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: trade ---
   :end-before: # --- section: health ---

Increase trade costs to explore localized food systems; decrease for globalized trade.

Emissions Pricing
~~~~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: emissions ---
   :end-before: # --- section: land use change ---

Land Use Change
~~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: land use change ---
   :end-before: # --- section: crops ---

Controls how land use change emissions and carbon sequestration are modeled over the planning horizon.

**Parameters**:
  * ``horizon_years``: Time horizon (years) for amortizing land use change emissions
  * ``managed_flux_mode``: How to treat emissions from existing managed land (``"zero"`` assumes no net flux from current agricultural land)
  * ``forest_fraction_threshold``: Minimum forest cover fraction (0-1) required for a grid cell to be eligible for regrowth sequestration when land is spared
  * ``spared_land_agb_threshold_tc_per_ha``: Maximum above-ground biomass (tonnes C per hectare) for spared land to be eligible for regrowth sequestration

Health Configuration
~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: health ---
   :end-before: # --- section: aggregation ---

Reduce ``region_clusters`` or ``log_rr_points`` to speed up solving.

Solver Configuration
~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: solving ---
   :end-before: # --- section: plotting ---

**Solver choice**:
  * **HiGHS**: Open-source, fast, good for most problems
  * **Gurobi**: Commercial, often faster for very large problems, requires license (free for academic users)

Plotting Configuration
~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: plotting ---

Customize visualization colors for publication-quality plots. The
``colors.food_groups`` palette is applied consistently across all food-group
charts and maps; extend it if you add new groups to ``data/food_groups.csv``.
