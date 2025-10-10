.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Nutrition
=========

Overview
--------

The nutrition module ensures that the optimized food system meets population dietary requirements. This includes:

* **Macronutrient constraints**: Carbohydrates, protein, fat, and calories per capita
* **Food group constraints**: Consumption of whole grains, fruits, vegetables, etc.
* **Population scaling**: Aggregating per-capita needs to regional/national totals

Macronutrients
--------------

Configuration
~~~~~~~~~~~~~

Macronutrient constraints are specified in ``config/default.yaml``:

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: macronutrients ---
   :end-before: # --- section: animal_products ---

**Constraint types**:

* ``min``: Lower bound (≥)
* ``max``: Upper bound (≤)
* ``equal``: Exact requirement (=)

.. TODO: write in more detail about the implementation
..
   Model Implementation
   ~~~~~~~~~~~~~~~~~~~~

Food Groups
-----------

Beyond macronutrients, the model can also constrains consumption of food groups. Moreover, food groups are used to assess dietary risk factors (see :ref:`health-impacts`).

Configuration
~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: food_groups ---
   :end-before: # --- section: diet ---

Foods are assigned to groups in ``data/food_groups.csv``. Example:

.. TODO: refine this section
..
   Model Implementation
   ~~~~~~~~~~~~~~~~~~~~

   Food group constraints work similarly to macronutrients:

   1. **Food group buses**: Per-country buses (e.g., ``food_group_fruit_USA``)

   2. **Food → Group links**: Foods contribute to their group bus

   3. **Population constraints**: Σ(group consumption) ≥ population × min requirement

   This ensures dietary diversity even if macronutrient needs could be met by a narrow set of crops.

Population Data
---------------

Population projections come from the UN World Population Prospects (WPP) 2024 revision.

Data Processing
~~~~~~~~~~~~~~~

The ``prepare_population`` rule (``workflow/scripts/prepare_population.py``):

1. **Load WPP data**: ``data/downloads/WPP_population.csv.gz``

2. **Filter**:

   * Countries in ``config['countries']``
   * Planning horizon year (``config['planning_horizon']``, e.g., 2030)
   * Medium variant projection

3. **Aggregate**: Sum population by country (converts thousands → persons)

4. **Output**:

   * ``processing/{name}/population.csv``: Total population by country
   * ``processing/{name}/population_age.csv``: Age-structured population for health module

Age Structure
~~~~~~~~~~~~~

Age-structured population is used in the health module to weight dietary risk factors by demographic composition (children vs. adults vs. elderly have different disease burdens).

Nutritional Content Data
-------------------------

The file ``data/nutrition.csv`` (currently mock data) contains nutritional composition for each food product. It is a work in progress to replace this with properly sourced data.

Per-Capita vs. Total Consumption
---------------------------------

The model works with total annual flows (Mt/year) but nutritional requirements are per-capita per-day. Conversion:

.. math::

   \text{Total requirement (Mt/year)} = \frac{\text{per capita (g/day)} \times \text{population} \times 365}{10^{12}}

This is handled internally by ``_per_capita_to_bus_units()`` in ``workflow/scripts/build_model.py``.

From the model's perspective:

* Food buses carry total food availability (Mt)
* Nutrient buses carry total nutrient availability (Mt for mass, Mcal for energy)
* Constraints compare these totals to population-scaled requirements

Dietary Patterns
----------------

The model does not currently prescribe specific dietary patterns (e.g., Mediterranean, vegetarian, EAT-Lancet) but rather:

1. **Lower / upper bounds**: Ensure minimum nutritional adequacy
2. **Cost minimization**: Subject to those bounds, minimize environmental + health costs

Workflow Integration
--------------------

Nutritional constraints are incorporated in the ``build_model`` rule:

1. **Load population**: ``processing/{name}/population.csv``
2. **Load nutrition data**: ``data/nutrition.csv``
3. **Create nutrient buses**: Per-country buses for each nutrient
4. **Create food → nutrient links**: Based on nutritional content
5. **Add global constraints**: Population × requirement bounds

No separate rule needed—nutrition is integrated into the model structure.
