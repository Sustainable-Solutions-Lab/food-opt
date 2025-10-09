.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Configuration System
====================

Overview
--------

The food-opt model is configuration-driven: all scenario parameters, crop selections, constraints, and solver options are defined in YAML configuration files under ``config/``. This allows exploring different scenarios without modifying code.

Main Configuration File
-----------------------

The primary configuration is ``config/default.yaml``, structured into thematic sections.

Scenario Identification
~~~~~~~~~~~~~~~~~~~~~~~

Scenario-specific overrides start by copying ``config/default.yaml`` and editing
only the pieces you need. A minimal custom file might look like this:

.. code-block:: yaml

   # config/my_scenario.yaml (excerpt)
   name: "my_scenario"           # Scenario name → results/my_scenario/
   planning_horizon: 2040        # Override the default 2030 horizon
   primary:
     land:
       regional_limit: 0.6       # Tighten land availability
   emissions:
     ghg_price: 250              # Raise the carbon price above the default

Any keys omitted in your custom file fall back to the defaults shown in the
sections below, so you can keep overrides concise.

Results are saved under ``results/{name}/``, allowing multiple scenarios to coexist.

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

Include countries to model; exclude to reduce problem size. Microstate and countries without level-1 GADM data are commented out.

Spatial Aggregation
-------------------

Controls regional resolution and land classification.

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: aggregation ---
   :end-before: # --- section: countries ---

**Trade-offs**:
  * More regions → higher spatial resolution, longer solve time
  * Fewer resource classes → faster solving, less yield heterogeneity
  * ``land_limit_dataset: "irrigated"`` → uniform land base, simpler but less realistic

Primary Resource Constraints
----------------------------

Limits on land, water, and fertilizer availability.

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: primary ---
   :end-before: # --- section: emissions ---

Tightening these constraints forces more efficient resource use or extensification.

GAEZ Data Parameters
--------------------

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
----------

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: irrigation ---
   :end-before: # --- section: solving ---

Restrict irrigation to water-scarce scenarios or explore rainfed-only production.

Nutritional Requirements
------------------------

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
   :end-before: # --- section: trade ---

Each food group may specify ``min_per_person_per_day``, ``max_per_person_per_day``,
and ``equal_per_person_per_day``. The defaults leave minima at zero so food group
constraints stay inactive; tighten minima or maxima to guide intakes, or use the
``equal`` field for equality targets.

Animal Products
---------------

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: animal_products ---
   :end-before: # --- section: food_groups ---

Disable grazing to force intensive feed-based systems.

Trade Configuration
-------------------

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: trade ---
   :end-before: # --- section: health ---

Increase trade costs to explore localized food systems; decrease for globalized trade.

Emissions Pricing
-----------------

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: emissions ---
   :end-before: # --- section: crops ---

**Values**:
  * 0: No carbon price (baseline)
  * 50-100: Current market prices
  * 200-300: Social cost of carbon
  * 500+: Stringent climate policy

Health Configuration
--------------------

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: health ---
   :end-before: # --- section: diet ---

Reduce ``region_clusters`` or ``log_rr_points`` to speed up solving.

Diet Controls
-------------

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: diet ---
   :end-before: # --- section: aggregation ---

Enable ``enforce_gdd_baseline`` to force the optimization to match baseline
consumption from the processed GDD file. Override ``baseline_age`` or
``baseline_reference_year`` if you pre-process alternative cohorts or years.

Solver Configuration
--------------------

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: solving ---
   :end-before: # --- section: plotting ---

**Solver choice**:
  * **HiGHS**: Open-source, fast, good for most problems
  * **Gurobi**: Commercial, often faster for very large problems, requires license

Plotting Configuration
----------------------

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: plotting ---

Customize visualization colors for publication-quality plots. The
``colors.food_groups`` palette is applied consistently across all food-group
charts and maps; extend it if you add new groups to ``data/food_groups.csv``.

Configuration Workflow
----------------------

Typical workflow for defining a new scenario:

1. **Copy base config**::

       cp config/default.yaml config/my_scenario.yaml

2. **Edit parameters**: Modify crops, constraints, solver options

3. **Set the scenario name** (top of the file)::

       name: "my_scenario"

4. **Run workflow**::

       tools/smk -j4 --configfile config/my_scenario.yaml all

5. **Results appear in**: ``results/my_scenario/``

Multiple scenarios can be run in parallel (different terminal sessions) or sequentially by changing the config file.

Configuration Validation
------------------------

The model performs basic validation:

* Missing required keys → error
* Invalid crop names (not in GAEZ mapping) → error during data retrieval
* Inconsistent constraints (e.g., ``kcal.min > kcal.max``) → solver infeasibility

More sophisticated validation (e.g., checking that food groups sum correctly) is future work.

Configuration Best Practices
----------------------------

**Start small**: Begin with a copy of the default config and relax constraints
  (e.g., fewer regions) before scaling up

**Scale up gradually**: Increase regions/crops/constraints incrementally

**Document changes**: Comment your config file with scenario rationale

**Version control**: Track config files in Git to reproduce results

**Compare scenarios**: Use consistent naming (``baseline``, ``high_carbon_price``, etc.)
