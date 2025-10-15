.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Workflow & Execution
====================

Overview
--------

The food-opt model uses `Snakemake <https://snakemake.readthedocs.io/>`__ for workflow orchestration. If you have never used Snakemake before, consider having a look at the official `tutotial <https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html>`__ to get familiar with the basic concepts. The workflow follows these main stages:

1. **Downloads** (GAEZ, GADM, UN WPP, FAOSTAT)
2. **Preprocessing** (regions, resource classes, yields, population, health)
3. **Model Building** (PyPSA network construction)
4. **Solving** (LP optimization with health costs)
5. **Visualization** (plots, maps, CSV exports)

Each stage is defined by Snakemake rules that specify inputs, outputs, and a script or piece of code.

The complete workflow dependency graph is shown below. Each node represents a Snakemake rule, and edges show dependencies between rules.

.. figure:: _static/figures/workflow_rulegraph.svg
   :alt: Workflow dependency graph
   :align: center
   :width: 100%

   Complete workflow dependency graph showing all Snakemake rules and their relationships

Key Snakemake Rules
-------------------

Data Preparation Rules
~~~~~~~~~~~~~~~~~~~~~~

**simplify_gadm**
  * **Input**: ``data/downloads/gadm.gpkg``
  * **Output**: ``processing/shared/gadm-simplified.gpkg``
  * **Script**: ``workflow/scripts/simplify_gadm.py``
  * **Purpose**: Simplify administrative boundaries for faster processing

**build_regions**
  * **Input**: Simplified GADM
  * **Output**: ``processing/{name}/regions.geojson``
  * **Script**: ``workflow/scripts/build_regions.py``
  * **Purpose**: Cluster administrative units into optimization regions

**prepare_population**
  * **Input**: ``data/downloads/WPP_population.csv.gz``
  * **Output**: ``processing/{name}/population.csv``, ``processing/{name}/population_age.csv``
  * **Script**: ``workflow/scripts/prepare_population.py``
  * **Purpose**: Extract population for planning horizon and countries

**compute_resource_classes**
  * **Input**: All GAEZ yield rasters, regions
  * **Output**: ``processing/{name}/resource_classes.nc``
  * **Script**: ``workflow/scripts/compute_resource_classes.py``
  * **Purpose**: Define yield quantile classes within each region

**aggregate_class_areas**
  * **Input**: Resource classes, suitability rasters, regions
  * **Output**: ``processing/{name}/land_area_by_class.csv``
  * **Script**: ``workflow/scripts/aggregate_class_areas.py``
  * **Purpose**: Compute available land area per (region, class, water, crop)

**build_crop_yields**
  * **Wildcards**: ``{crop}`` (crop name), ``{water_supply}`` ("r" or "i")
  * **Input**: Resource classes, GAEZ rasters (yield, suitability, water, growing season)
  * **Output**: ``processing/{name}/crop_yields/{crop}_{water_supply}.csv``
  * **Script**: ``workflow/scripts/build_crop_yields.py``
  * **Purpose**: Aggregate yields by (region, class) for each crop

**build_grassland_yields**
  * **Input**: ISIMIP grassland yield NetCDF, resource classes, regions
  * **Output**: ``processing/{name}/grassland_yields.csv``
  * **Script**: ``workflow/scripts/build_grassland_yields.py``
  * **Purpose**: Aggregate grassland yields for grazing production

**prepare_health_costs**
  * **Input**: Regions, DIA health data, population
  * **Output**: ``processing/{name}/health/*.csv`` (risk breakpoints, dose-response, clusters)
  * **Script**: ``workflow/scripts/prepare_health_costs.py``
  * **Purpose**: Compute health cluster parameters for DALY calculations

Model Building and Solving
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**build_model**
  * **Input**: All crop yields, grassland yields, land areas, population, water availability, static data files (crops.csv, foods.csv with pathway-based processing, etc.)
  * **Output**: ``results/{name}/build/model.nc``
  * **Script**: ``workflow/scripts/build_model.py``
  * **Purpose**: Construct PyPSA network with all components, links, and constraints. Creates multi-output links for crop→food conversion based on processing pathways defined in foods.csv.

**solve_model**
  * **Input**: Built model, health data, food-to-risk mapping
  * **Output**: ``results/{name}/solved/model.nc``
  * **Script**: ``workflow/scripts/solve_model.py``
  * **Purpose**: Add health costs, solve LP, save results

Visualization Rules
~~~~~~~~~~~~~~~~~~~

**Plots and maps** (see ``workflow/rules/plotting.smk``):
  * ``plot_regions_map``: Optimization region boundaries
  * ``plot_resource_classes_map``: Resource class spatial distribution
  * ``plot_crop_production_map``: Crop production by region
  * ``plot_crop_land_use_map``: Land use by crop
  * ``plot_cropland_fraction_map``: Cropland fraction of each region
  * ``plot_water_value_map``: Water shadow prices (economic value)
  * ``plot_health_impacts``: Health risk and baseline maps
  * ``plot_results``: Production, resource usage, objective breakdown
  * ``plot_food_consumption``: Dietary composition
  * ``plot_crop_use_breakdown``: How crops are used (food vs. feed vs. waste)

Execution Commands
------------------

Running the Full Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~

Build, solve, and visualize everything::

    tools/smk -j4 --configfile config/my_scenario.yaml all

* ``-j4``: Use 4 parallel cores (adjust to your CPU count)
* ``--configfile config/my_scenario.yaml``: Specify which configuration file to use
* ``all``: Target rule that depends on all major outputs (strictly speaking optional)

This will:

1. Download datasets (if not already cached)
2. Process data for configured scenario
3. Build and solve the model
4. Generate all plots and exports

Running Specific Stages
~~~~~~~~~~~~~~~~~~~~~~~~

**Build model only** (no solving)::

    tools/smk -j4 --configfile config/my_scenario.yaml results/my_scenario/build/model.nc

**Solve model**::

    tools/smk -j4 --configfile config/my_scenario.yaml results/my_scenario/solved/model.nc

**Regenerate specific plot**::

    tools/smk --configfile config/my_scenario.yaml results/my_scenario/plots/crop_production.pdf

**Prepare data without building model**::

    tools/smk -j4 --configfile config/my_scenario.yaml processing/my_scenario/regions.geojson processing/my_scenario/resource_classes.nc

For any of the above targets, Snakemake will first run any other previous rules in order to generate the necessary inputs for the specified target output/rule.

Checking Workflow Status
~~~~~~~~~~~~~~~~~~~~~~~~~

**Dry-run** (show what would be executed without running)::

    tools/smk --configfile config/my_scenario.yaml -n

**Dependency graph**: See the workflow dependency graph figure at the top of this page. To generate a detailed job-level DAG for a specific configuration (requires Graphviz)::

    tools/smk --dag all | dot -Tpdf > dag.pdf

**List all rules**::

    tools/smk --list

Memory Management
-----------------

The ``tools/smk`` Wrapper
~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to run the workflow directly with the ``snakemake`` command. Food-opt, however, provides a simple shell script, ``tools/smk``, which:

1. Runs Snakemake in a systemd cgroup with hard memory limit (default 10 GB), killing the process group if memory limit is exceeded
2. Disables swap to prevent system instability
3. Sets the ``-j1`` argument (running only one job at a time) by default unless the user sets the ``-j<n>`` option explicitly.

**Default memory limit**: 10 GB (configurable via ``SMK_MEM_MAX`` environment variable)

**Override memory limit**::

    SMK_MEM_MAX=12G tools/smk -j4 all

Parallelization
---------------

Snakemake automatically runs rules concurrently (e.g., downloading multiple GAEZ files, processing yields for different crops), depending on the configured number of parallel rules allowed. This is set with the ``-j<n>`` option, where ``n`` is the number of parallel rules. Note that individual rules (such as the model solving rule) may use more than one processor core. 

Snakemake automatically detects dependencies and runs tasks in correct order.

Incremental Development
-----------------------

**Workflow philosophy**: Snakemake tracks file modification times and only reruns rules whose inputs changed. This includes rule input files, the script associated with the rule as well as rule parameters (relevant configuration sections).

**Example workflow**:

1. Run full workflow: ``tools/smk -j4 all``
2. Modify crop list in config → only crop yield rules rerun
3. Modify solver options → only ``solve_model`` reruns (build model reused)
4. Modify visualization script → only plotting rules rerun

**Rerun specific rule**::

    tools/smk -j4 --configfile config/my_scenario.yaml results/my_scenario/solved/model.nc --forcerun solve_model

**Mark all existing outputs as up to date** (preventing rules from being run due to modification times, etc.)::

    tools/smk --configfile config/my_scenario.yaml --touch
