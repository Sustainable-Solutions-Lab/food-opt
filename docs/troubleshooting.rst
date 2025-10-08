.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Troubleshooting & FAQ
=====================

Common Issues
-------------

Workflow Execution
~~~~~~~~~~~~~~~~~~

**Problem**: ``ModuleNotFoundError: No module named 'pypsa'`` (or other package)

**Cause**: Dependencies not installed or wrong Python environment

**Solution**::

    uv sync  # Install all dependencies
    uv run snakemake --version  # Verify Snakemake accessible

**Problem**: ``MemoryError`` or system freezes during workflow

**Cause**: Workflow exceeds available RAM

**Solution**:

* Reduce ``aggregation.regions.target_count`` in config (fewer regions = less memory)
* Increase ``SMK_MEM_MAX``::

      SMK_MEM_MAX=16G tools/smk -j4 all

* Close other applications to free RAM

**Problem**: ``Network error`` downloading GAEZ files

**Cause**: Firewall, timeout, or GAEZ server issue

**Solution**:

* Check internet connectivity
* Retry (Snakemake resumes from last successful step)
* Manually download problematic files and place in ``data/downloads/``
* Check GAEZ website for maintenance notices

**Problem**: Workflow hangs at solving stage

**Cause**: Large problem, slow solver convergence, or infeasible model

**Solution**:

* **Check progress**: Look for solver output (HiGHS or Gurobi logs)
* **Reduce problem size**: Fewer regions, crops, or resource classes
* **Increase time limit**: Solvers have default time limits (adjust in config ``solving.options_*``)
* **Switch solvers**: Try Gurobi if using HiGHS (or vice versa)
* **Check feasibility**: If solver reports "infeasible," constraints conflict

Model Infeasibility
~~~~~~~~~~~~~~~~~~~

**Problem**: Solver reports "model is infeasible"

**Cause**: Conflicting constraints (e.g., nutritional requirements exceed production capacity)

**Diagnosis**:

1. **Relax constraints incrementally**:

   * Lower macronutrient minimums
   * Reduce food group requirements
   * Increase land availability (``primary.land.regional_limit``)
   * Add more crops

2. **Check error messages**: Solver may indicate which constraints conflict

3. **Inspect intermediate results**: Verify ``results/{name}/build/model.nc`` created successfully

**Solution**: Adjust config to make problem feasible, or identify which constraint is unrealistic.

**Problem**: Solver reports "unbounded"

**Cause**: Missing upper bounds on decision variables (e.g., infinite production allowed)

**Diagnosis**: Check model for links without capacity limits (``p_nom`` or ``p_nom_max``)

**Solution**: Add realistic upper bounds (e.g., land area limits, maximum fertilizer)

Data Issues
~~~~~~~~~~~

**Problem**: ``KeyError: 'my_crop'`` when building yields

**Cause**: Crop not in GAEZ mapping file

**Solution**:

* Check crop name matches ``data/gaez_crop_code_mapping.csv``
* Add missing crop to mapping (if GAEZ supports it)
* Remove crop from config ``crops`` list if no GAEZ data

**Problem**: ``FileNotFoundError: data/nutrition.csv``

**Cause**: Missing static data file

**Solution**:

* Ensure all required CSV files exist in ``data/``
* Check for typos in file names
* For mock data files, copy from repository (should be committed)

**Problem**: "Mock data warning" in logs

**Cause**: Using placeholder data files (``foods.csv``, ``nutrition.csv``, etc.)

**Solution**:

* This is expected for development
* For analysis, replace with vetted data from USDA FoodData Central, etc.
* See :doc:`data_sources` for recommended replacements

Visualization Errors
~~~~~~~~~~~~~~~~~~~~

**Problem**: ``KeyError`` or ``AttributeError`` in plotting script

**Cause**: Expected component missing from solved model

**Diagnosis**:

* Check if model solved successfully (``results/{name}/solved/model.nc`` exists and is non-empty)
* Inspect model to see if expected links/buses exist::

      import pypsa
      n = pypsa.Network("results/my_scenario/solved/model.nc")
      print(n.links.index)  # List all links

**Solution**: Fix model building script or update plotting script to handle missing components

**Problem**: Plot is empty or shows no data

**Cause**: Optimal solution has zero flows for visualized components

**Diagnosis**: Check if constraints forced zero production (e.g., crop excluded, region has no suitable land)

**Solution**: Verify constraints are reasonable; may be expected behavior

Performance Issues
~~~~~~~~~~~~~~~~~~

**Problem**: Workflow takes hours to complete

**Cause**: Large problem size (many regions, crops, health clusters)

**Solution**:

* **Reduce resolution**: Fewer regions (e.g., 100 instead of 400)
* **Fewer crops**: Start with ~20 key crops, expand later
* **Fewer health clusters**: Reduce ``health.region_clusters`` (e.g., 10 instead of 30)
* **Use Gurobi**: Commercial solver often faster than HiGHS for large problems

**Problem**: High disk I/O during workflow

**Cause**: Reading/writing large raster files

**Solution**:

* Use SSD instead of HDD for ``data/`` and ``processing/`` directories
* Reduce raster resolution if possible (GAEZ lowest resolution: ~0.08°)
* Process fewer crops in parallel (reduce ``-j`` value)

Solver-Specific Issues
~~~~~~~~~~~~~~~~~~~~~~

**Problem**: ``AttributeError: 'Network' object has no attribute 'optimize'`` (PyPSA < 0.20)

**Cause**: Using outdated PyPSA API

**Solution**: Update PyPSA::

    uv sync  # Should install correct version from pyproject.toml

**Problem**: Gurobi license error

**Cause**: Gurobi license not found or expired

**Solution**:

* Obtain academic or commercial license from Gurobi
* Place ``gurobi.lic`` in home directory or set ``GRB_LICENSE_FILE`` environment variable
* Or switch to HiGHS (open-source)::

      # In config.yaml
      solving:
        solver: highs

**Problem**: HiGHS reports "numerical difficulties"

**Cause**: Ill-conditioned problem (large variations in coefficient magnitudes)

**Solution**:

* Check unit consistency (e.g., mixing tonnes and kilograms)
* Rescale variables (e.g., use Mt instead of kg)
* Try barrier method: ``options_highs.solver: "ipm"`` (interior-point)

Frequently Asked Questions
--------------------------

General
~~~~~~~

**Q: How long does the first run take?**

A: First run downloads ~10-20 GB of data (GAEZ, GADM, UN WPP), taking 30-120 minutes depending on internet speed. Subsequent runs reuse cached data and take ~10-60 minutes for a 400-region problem.

**Q: Can I run multiple scenarios in parallel?**

A: Yes, as long as they have different ``name`` values in config. Run in separate terminal sessions::

    # Terminal 1
    # config.yaml: name: "baseline"
    tools/smk -j4 all

    # Terminal 2
    # config.yaml: name: "high_carbon"
    tools/smk -j4 all

Results go to separate directories (``results/baseline/``, ``results/high_carbon/``).

**Q: How do I reduce memory usage?**

A: Key parameters:

* ``aggregation.regions.target_count``: Fewer regions = less memory (try 100-200)
* ``health.region_clusters``: Fewer clusters = less memory (try 10-20)
* Number of crops: Start with ~20 key crops

Also, use ``tools/smk`` wrapper which enforces memory limits.

**Q: Why is solving so slow?**

A: Large linear programs (millions of variables/constraints) are computationally intensive. To speed up:

* Reduce regions/crops as above
* Use Gurobi (faster than HiGHS for large problems)
* Use a machine with more cores (solver parallelizes internally)
* Loosen solver tolerance: ``options_highs.mip_rel_gap: 0.01`` (1% gap instead of 0.1%)

Model Behavior
~~~~~~~~~~~~~~

**Q: Why does the model produce unrealistic diets (e.g., all soy, no variety)?**

A: With only macronutrient constraints, the model finds the cheapest/lowest-emission way to meet requirements, which may be monotonous. Add food group constraints to enforce diversity::

    food_groups:
      whole grain:
        min_per_person_per_day: 100  # Force at least 100g/day

**Q: Why does the model use all available land even with low demand?**

A: Check if there's a constraint forcing production (e.g., minimum food group requirements) or if land use change emissions are too low (set ``emissions.ghg_price`` higher to penalize expansion).

**Q: Why doesn't the model use irrigation even when water is available?**

A: Possible reasons:

* Irrigated crops not enabled: Check ``irrigation.irrigated_crops`` in config
* Water requirements exceed availability: Check ``water_value_map.pdf`` for binding constraints
* Higher cost: Irrigated production may cost more than rainfed + trade

**Q: Why do some regions import crops they can grow locally?**

A: Trade can be cheaper than local production if:

* Local land is more valuable for other crops (opportunity cost)
* Exporting region has comparative advantage (higher yields, lower costs)
* Transport costs are low relative to production cost differences

This is economically efficient but may conflict with food security goals (can constrain trade to explore self-sufficiency).

Technical Details
~~~~~~~~~~~~~~~~~

**Q: What solver should I use?**

A: **HiGHS** (default) for most cases: open-source, fast, no license required.

**Gurobi** for very large problems (>1M variables): faster, but requires license (free academic licenses available).

**Q: How accurate are the results?**

A: Depends on data quality:

* **Yield potentials** (GAEZ): Well-validated, but climate projections uncertain
* **Nutritional content**: Mock data currently — replace with USDA FoodData
* **Feed conversion ratios**: Mock data — needs zootechnical sources
* **Health dose-response**: Based on peer-reviewed GBD studies (high quality)

Treat current results as exploratory; publication-quality analysis requires replacing mock data.

**Q: Can I use this model for a specific country or region?**

A: Yes:

1. Set ``countries`` to your country/countries
2. Reduce ``aggregation.regions.target_count`` to match administrative units
3. Optionally: Replace ``build_regions`` with custom regional definitions

**Q: How do I add a new food product?**

A: Edit ``data/foods.csv`` and ``data/nutrition.csv`` to add the food, then update ``data/food_groups.csv`` to assign it to a group. Rerun workflow.

**Q: How do I export results to GIS software?**

A: Regions with data::

    import geopandas as gpd
    import pypsa

    n = pypsa.Network("results/my_scenario/solved/model.nc")
    regions = gpd.read_file("processing/my_scenario/regions.geojson")

    # Add production/emission columns to regions GeoDataFrame
    # regions["production"] = extract_regional_production(n)

    regions.to_file("results/my_scenario/exports/regions_with_data.geojson")

Then import GeoJSON into QGIS, ArcGIS, etc.

Getting Help
------------

**Documentation**: Read :doc:`index` sections relevant to your issue

**Logs**: Check Snakemake/solver logs for error messages

**Issues**: Open a GitHub issue with:

* Error message (full traceback)
* Configuration file (``config.yaml``)
* Steps to reproduce
* System info (OS, Python version, RAM)

**Discussion**: GitHub Discussions for questions, feature requests, or general Q&A
