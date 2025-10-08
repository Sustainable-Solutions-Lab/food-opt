.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Nutrition
=========

Overview
--------

The nutrition module ensures that the optimized food system meets population dietary requirements. This includes:

* **Macronutrient constraints**: Minimum carbohydrates, protein, fat, and calories per capita
* **Food group constraints**: Minimum consumption of whole grains, fruits, vegetables, etc.
* **Population scaling**: Aggregating per-capita needs to regional/national totals

Nutritional requirements drive the model to produce not just adequate calories, but a balanced, healthy diet.

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

**Reference values**: Based on WHO/FAO guidelines and EAT-Lancet recommendations for healthy diets.

Model Implementation
~~~~~~~~~~~~~~~~~~~~

Macronutrient constraints are implemented as PyPSA global constraints:

1. **Nutrient buses**: Per-country buses for each nutrient (e.g., ``nutrient_protein_USA``)

2. **Food → Nutrient links**: Each food product contributes to nutrient buses

   * Efficiency: Nutritional content per tonne of food (e.g., 12 kg protein per 100 kg wheat)
   * From ``data/nutrition.csv``

3. **Population constraints**: Σ(nutrient consumption) = population × requirement

Example: Protein constraint for USA with 330M population and 50 g/person/day requirement:

.. math::

   \sum_{\text{foods}} (\text{food consumption}_\text{USA} \times \text{protein content}) \geq 330 \times 10^6 \times 50 \times 365 \times 10^{-9} \text{ Mt}

Unit conversions handled in ``workflow/scripts/build_model.py``:

* g/person/day → Mt/year (for mass nutrients)
* kcal/person/day → Mcal/year (for energy)

Food Groups
-----------

Beyond macronutrients, the model constrains consumption of food groups to promote dietary diversity and quality.

Configuration
~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: food_groups ---
   :end-before: # --- section: trade ---

**Defaults**: All minima are zero in the base config, leaving food group constraints inactive. Raise these values (e.g., 50 g/day for ``whole_grains``) to enforce dietary diversity targets aligned with WHO guidance (≥400 g/day fruit+vegetables, ≥150 g/day whole grains).

Food Group Mapping
~~~~~~~~~~~~~~~~~~

Foods are assigned to groups in ``data/food_groups.csv`` (mock data). Example:

* ``bread`` → ``grain``
* ``whole_wheat_bread`` → ``whole grain``
* ``apple`` → ``fruit``
* ``carrot`` → ``vegetable``
* ``beef`` → ``animal protein``

Some foods may belong to multiple groups or have fractional assignments.

Model Implementation
~~~~~~~~~~~~~~~~~~~~

Food group constraints work similarly to macronutrients:

1. **Food group buses**: Per-country buses (e.g., ``food_group_fruit_USA``)

2. **Food → Group links**: Foods contribute to their group bus

   * Efficiency: 1.0 if fully in group, fractional if partial membership

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

The file ``data/nutrition.csv`` (currently mock data) contains nutritional composition for each food product. Required columns:

* ``food``: Food name (matching ``data/foods.csv``)
* Macronutrients: ``carb_g_per_100g``, ``protein_g_per_100g``, ``fat_g_per_100g``, ``kcal_per_100g``
* Micronutrients (optional): ``calcium_mg_per_100g``, ``iron_mg_per_100g``, etc.

**Units**: Typically per 100g (standard nutrition labeling convention)

**Data sources** (recommended replacements for mock data):

* **USDA FoodData Central**: https://fdc.nal.usda.gov/
* **FAO INFOODS**: https://www.fao.org/infoods/infoods/tables-and-databases/en/
* **National food composition tables**: Many countries maintain detailed databases

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

The model does not prescribe specific dietary patterns (e.g., Mediterranean, vegetarian, EAT-Lancet) but rather:

1. **Lower bounds**: Ensure minimum nutritional adequacy
2. **Cost minimization**: Subject to those bounds, minimize environmental + health costs

This allows the optimizer to discover efficient dietary patterns rather than imposing them a priori.

To explore specific diets, you can:

* **Constrain animal products**: Set ``animal_products.max_per_person_per_day``
* **Require high plant foods**: Increase fruit/vegetable/whole grain minimums
* **Price emissions**: Higher ``emissions.ghg_price`` discourages ruminant meat

Regional Dietary Heterogeneity
-------------------------------

Currently, nutritional constraints are uniform across countries (same g/person/day requirements). Future extensions could add:

* **Country-specific requirements**: Adjust for demographic differences (age structure, body size)
* **Cultural preferences**: Minimum consumption of staple foods (rice in Asia, wheat in Europe)
* **Income-dependent targets**: Lower-income regions may have different nutritional priorities

Micronutrients
--------------

The current model focuses on macronutrients. Micronutrients (vitamins, minerals) can be added by:

1. Extending ``data/nutrition.csv`` with micronutrient columns
2. Adding micronutrient buses and constraints in ``build_model.py``
3. Setting minimum daily intakes (e.g., 10 mg iron, 1000 mg calcium)

This would capture the "hidden hunger" problem where calorie needs are met but micronutrient deficiencies persist.

Validation and Sanity Checks
-----------------------------

After solving, validate nutritional outcomes:

**Total calories**::

    tools/smk results/{name}/plots/food_consumption.pdf

Check that per-capita calorie consumption equals the configured requirement.

**Food group consumption**::

    tools/smk results/{name}/plots/food_consumption.csv

Verify that fruit, vegetable, whole grain consumption meets minimums.

**Macronutrient balance**: Inspect protein/carb/fat ratios (should be reasonable, e.g., 10-35% protein, 20-35% fat, 45-65% carbs by energy)

Workflow Integration
--------------------

Nutritional constraints are incorporated in the ``build_model`` rule:

1. **Load population**: ``processing/{name}/population.csv``
2. **Load nutrition data**: ``data/nutrition.csv``
3. **Create nutrient buses**: Per-country buses for each nutrient
4. **Create food → nutrient links**: Based on nutritional content
5. **Add global constraints**: Population × requirement bounds

No separate rule needed—nutrition is integrated into the model structure.
