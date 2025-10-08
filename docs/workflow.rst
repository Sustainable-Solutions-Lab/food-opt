.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Workflow & Execution
====================

Overview
--------

The food-opt model uses Snakemake for workflow orchestration, automating:

1. Data retrieval and preprocessing
2. Model building (PyPSA network construction)
3. Solving (linear program optimization)
4. Postprocessing and visualization

This page describes the workflow stages, key rules, and execution commands.

Workflow Stages
---------------

The workflow follows a dependency graph:

.. code-block:: text

   Downloads (GAEZ, GADM, UN WPP, FAOSTAT)
           ↓
   Preprocessing (regions, resource classes, yields, population, health)
           ↓
   Model Building (PyPSA network construction)
           ↓
   Solving (LP optimization with health costs)
           ↓
   Visualization (plots, maps, CSV exports)

Each stage is defined by Snakemake rules that specify inputs, outputs, and scripts.

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
  * **Input**: All crop yields, grassland yields, land areas, population, water availability, static data files (crops.csv, foods.csv, etc.)
  * **Output**: ``results/{name}/build/model.nc``
  * **Script**: ``workflow/scripts/build_model.py``
  * **Purpose**: Construct PyPSA network with all components, links, and constraints

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

Each visualization rule has:
  * **Input**: ``results/{name}/solved/model.nc`` plus auxiliary data
  * **Output**: ``.pdf`` plots and ``.csv`` data tables
  * **Script**: Corresponding script in ``workflow/scripts/``

Execution Commands
------------------

Running the Full Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~

Build, solve, and visualize everything::

    tools/smk -j4 all

* ``-j4``: Use 4 parallel cores (adjust to your CPU count)
* ``all``: Target rule that depends on all major outputs

This will:

1. Download datasets (if not already cached)
2. Process data for configured scenario
3. Build and solve the model
4. Generate all plots and exports

Running Specific Stages
~~~~~~~~~~~~~~~~~~~~~~~~

**Build model only** (no solving)::

    tools/smk -j4 --configfile config/my_scenario.yaml results/my_scenario/build/model.nc

**Solve existing built model**::

    tools/smk -j4 --configfile config/my_scenario.yaml results/my_scenario/solved/model.nc

**Regenerate specific plot** (assuming model solved)::

    tools/smk --configfile config/my_scenario.yaml results/my_scenario/plots/crop_production.pdf

**Prepare data without building model**::

    tools/smk -j4 --configfile config/my_scenario.yaml processing/my_scenario/regions.geojson processing/my_scenario/resource_classes.nc

Checking Workflow Status
~~~~~~~~~~~~~~~~~~~~~~~~~

**Dry-run** (show what would be executed without running)::

    tools/smk -j4 all --dry-run

**Dependency graph** (requires Graphviz)::

    tools/smk --dag all | dot -Tpdf > dag.pdf

This generates a visual workflow diagram.

**List all rules**::

    tools/smk --list

Memory Management
-----------------

The ``tools/smk`` Wrapper
~~~~~~~~~~~~~~~~~~~~~~~~~

**Never run** ``snakemake`` directly for this project. Always use ``tools/smk``, which:

1. Runs Snakemake in a systemd cgroup with hard memory limit (default 10 GB)
2. Disables swap to prevent system instability
3. Kills the process group if memory limit is exceeded

**Default memory limit**: 10 GB (configurable via ``SMK_MEM_MAX`` environment variable)

**Override memory limit**::

    SMK_MEM_MAX=12G tools/smk -j4 all

**Why this matters**: The model can consume significant memory (especially with many regions/crops), and exceeding system RAM causes thrashing or OOM kills. The wrapper provides graceful failure.

Parallelization
---------------

Snakemake parallelizes independent tasks:

* **Rule-level parallelism**: Different rules run concurrently (e.g., downloading multiple GAEZ files, processing yields for different crops)
* **Within-rule parallelism**: Not used by default (scripts are single-threaded)

**Optimal core count**:
  * **Data prep**: ``-j`` = CPU cores (many independent rules)
  * **Model solving**: ``-j1`` (solver uses all cores internally)

**Example** — 8-core machine::

    # Data prep with parallelism
    tools/smk -j8 --configfile config/my_scenario.yaml processing/my_scenario/resource_classes.nc

    # Solving (solver will use multiple cores)
    tools/smk -j1 --configfile config/my_scenario.yaml results/my_scenario/solved/model.nc

Snakemake automatically detects dependencies and runs tasks in correct order.

Handling Failures
-----------------

**Network downloads fail**:
  * **Cause**: Timeout, connection issues
  * **Solution**: Rerun; Snakemake resumes from where it failed

**Script errors**:
  * **Cause**: Missing data, invalid config, code bug
  * **Solution**: Check error message, fix config/code, rerun

**Memory limit exceeded**:
  * **Cause**: Too many regions, insufficient system RAM
  * **Solution**: Increase ``SMK_MEM_MAX`` or reduce ``aggregation.regions.target_count`` in config

**Solver infeasibility**:
  * **Cause**: Conflicting constraints (e.g., nutritional requirements impossible with available land/crops)
  * **Solution**: Relax constraints, add more crops, increase land availability

Incremental Development
-----------------------

**Workflow philosophy**: Snakemake tracks file modification times and only reruns rules whose inputs changed.

**Example workflow**:

1. Run full workflow: ``tools/smk -j4 all``
2. Modify crop list in config → only crop yield rules rerun
3. Modify solver options → only ``solve_model`` reruns (build model reused)
4. Modify visualization script → only plotting rules rerun

**Force rerun** (ignore timestamps)::

    tools/smk -j4 all --forceall

**Rerun specific rule**::

    tools/smk -j4 --configfile config/my_scenario.yaml results/my_scenario/solved/model.nc --forcerun solve_model

Network Access for Downloads
-----------------------------

Rules that download data (``retrieve_*``) require network access. The ``tools/smk`` wrapper runs in a restricted cgroup by default.

**If downloads fail** due to network restrictions:

1. Confirm you have internet connectivity
2. Check firewall rules
3. Run outside the memory-limited cgroup (not recommended for full workflow)::

       snakemake -j4 data/downloads/gadm.gpkg

   Then use ``tools/smk`` for the rest.

Workflow Customization
----------------------

**Adding a new crop**:

1. Add crop name to ``config.yaml`` ``crops`` list
2. Ensure crop is in ``data/gaez_crop_code_mapping.csv``
3. Rerun: ``tools/smk -j4 all``
4. Snakemake will download new GAEZ files and integrate crop

**Adding a new visualization**:

1. Create script in ``workflow/scripts/plot_*.py``
2. Add rule in ``workflow/rules/plotting.smk``
3. Add output to ``all`` rule dependencies
4. Run: ``tools/smk --configfile config/my_scenario.yaml results/my_scenario/plots/my_new_plot.pdf``

**Changing spatial resolution**:

1. Edit ``aggregation.regions.target_count`` in config
2. Rerun: ``tools/smk -j4 all`` (will rebuild regions and downstream)

Workflow Best Practices
-----------------------

* **Version control**: Commit config changes to track scenario evolution
* **Separate scenarios**: Use different ``name`` values, don't overwrite results
* **Incremental testing**: Test with small region counts (50-100) before full-scale runs
* **Monitor memory**: Watch system resources during first run to gauge memory needs
* **Checkpoint frequently**: For long-running workflows, confirm intermediate outputs (``build/model.nc``) succeed before solving
