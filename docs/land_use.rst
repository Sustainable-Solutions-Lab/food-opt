.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Land Use & Resource Classes
============================

Overview
--------

The land use module translates high-resolution gridded yield potentials (0.05° × 0.05°, approximately 5.6 km at the equator) into aggregated land availability and yield parameters for the optimization model. This aggregation serves two purposes:

1. **Computational efficiency**: Reduces millions of gridcells to hundreds of optimization regions
2. **Yield heterogeneity**: Captures within-region variation in land quality using resource classes

Spatial Aggregation
-------------------

Regional Clustering
~~~~~~~~~~~~~~~~~~~

The model operates at sub-national regional resolution, defined by clustering administrative units (GADM level 1) based on spatial proximity. The clustering process (``workflow/scripts/build_regions.py``) follows these steps:

1. **Simplification**: Simplify GADM geometries to reduce complexity while preserving boundaries

   * Controlled by ``aggregation.simplify_tolerance_km`` and ``aggregation.simplify_min_area_km``
   * Small enclaves are removed to avoid fragmentation

2. **Country selection**: Filter to configured countries (``countries`` list in config)

3. **Clustering**: Aggregate administrative units into target number of regions

   * Method: k-means clustering on region centroids (default)
   * Target count: ``aggregation.regions.target_count`` (e.g., 400)
   * Cross-border clustering: Controlled by ``aggregation.regions.allow_cross_border`` (typically ``false``)

4. **Output**: GeoJSON with region polygons saved to ``processing/{name}/regions.geojson``

The result is a set of contiguous optimization regions that respect national boundaries (unless cross-border clustering is enabled) and balance spatial detail with computational tractability.

Resource Classes
----------------

Concept
~~~~~~~

Within each optimization region, agricultural land is heterogeneous—some areas have high yield potential, others low. To capture this heterogeneity without creating a separate decision variable for each gridcell, the model groups land into **resource classes** based on yield potential quantiles.

For example, with quantiles ``[0.25, 0.5, 0.75]``, each region has 4 resource classes:

* **Class 0**: Bottom 25% of yield potential (lowest quality land)
* **Class 1**: 25th-50th percentile
* **Class 2**: 50th-75th percentile
* **Class 3**: Top 25% of yield potential (highest quality land)

This allows the model to:

* Preferentially allocate high-value crops to high-quality land
* Avoid optimistic bias from averaging yields across heterogeneous land
* Capture marginal land-use decisions

Computation
~~~~~~~~~~~

Resource classes are computed by ``workflow/scripts/compute_resource_classes.py``:

1. **Load yield rasters**: Read GAEZ yield potentials for all crops (both rainfed and irrigated)

2. **Running maximum**: For each gridcell, compute the maximum attainable yield across all crops and water supplies

   * Converts kg/ha → t/ha (tonnes per hectare)
   * This represents the "best possible use" of each piece of land

3. **Regional quantiles**: Within each optimization region, compute yield quantiles

   * Uses ``aggregation.resource_class_quantiles`` from config (e.g., ``[0.25, 0.5, 0.75]``)
   * Only considers gridcells with positive yield potential (deserts/ice don't collapse bins)

4. **Class assignment**: Assign each gridcell to a resource class based on which quantile bin it falls into

5. **Output**: NetCDF file (``processing/{name}/resource_classes.nc``) with:

   * ``resource_class``: Integer class ID for each gridcell
   * ``max_yield``: Maximum yield potential across all crops
   * Coordinate reference system and geotransform for spatial reference

Example
~~~~~~~

For a region with yield distribution [0.1, 0.5, 1.2, 2.0, 3.5, 5.0, 8.0] t/ha with quantiles [0.33, 0.67]:

* Class 0: [0.1, 0.5, 1.2] t/ha (bottom third)
* Class 1: [2.0, 3.5] t/ha (middle third)
* Class 2: [5.0, 8.0] t/ha (top third)

Land Area Aggregation
---------------------

Once resource classes are defined, the next step is to compute how much land is available in each region-class-water combination. This is handled by ``workflow/scripts/aggregate_class_areas.py``.

Suitability-Based Limits
~~~~~~~~~~~~~~~~~~~~~~~~~

Land availability is constrained by **suitability** from GAEZ, which indicates what fraction of each gridcell is suitable for agriculture. GAEZ provides separate suitability rasters for rainfed and irrigated production (SX1 variable: share of cell assessed as "very suitable" or "suitable").

The aggregation process:

1. **Load resource classes**: Read the class assignment raster from previous step

2. **Load suitability**: Read GAEZ suitability rasters for each crop and water supply

3. **Compute suitable area**: For each gridcell, multiply:

   * Cell area (in hectares, accounting for Earth's curvature)
   * Suitability fraction (0-1)
   * Irrigated share (for irrigated crops) or ``1 - irrigated_share`` (for rainfed)

4. **Aggregate by region-class**: Sum suitable area across all gridcells in each (region, resource_class, water_supply, crop) combination

5. **Output**: CSV file (``processing/{name}/land_area_by_class.csv``) with columns:

   * ``region``: Optimization region ID
   * ``resource_class``: Class number (0, 1, 2, ...)
   * ``water_supply``: "r" (rainfed) or "i" (irrigated)
   * ``crop``: Crop name
   * ``area_ha``: Available area in hectares

Land Limit Dataset Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The configuration parameter ``aggregation.land_limit_dataset`` controls how land availability is determined:

* **"suitability"** (default): Use GAEZ suitability rasters per crop and water supply

  * More conservative: only counts land GAEZ deems suitable
  * Crop-specific: different crops have different suitable area

* **"irrigated"**: Use irrigated cropland share for all crops

  * Simpler: all crops have the same suitable area per water supply
  * Based on GAEZ's "land equipped for irrigation" dataset

The "irrigated" option is typically used when you want a uniform land base across all crops, while "suitability" is more realistic but creates more heterogeneity.

Irrigated vs. Rainfed Land
---------------------------

The model distinguishes between irrigated and rainfed production:

**Rainfed**
  * Uses rainfall for water supply
  * Available on all suitable cropland not equipped for irrigation
  * Generally lower yields than irrigated

**Irrigated**
  * Uses irrigation infrastructure
  * Only available on land equipped for irrigation (from GAEZ dataset)
  * Higher yields but requires blue water (see :doc:`crop_production` for water constraints)
  * Controlled by ``irrigation.irrigated_crops`` config (can be "all" or a list)

For each optimization region and resource class, the model maintains separate land variables for rainfed and irrigated production. These compete for the same physical land, so the model includes constraints ensuring that total rainfed + irrigated land in a class doesn't exceed available area.

Land Use Change
---------------

When the model allocates more land to agriculture than is currently used, this represents **land use change** (LUC). The environmental impacts of LUC (carbon emissions from clearing vegetation) are captured in the objective function.

The model does not currently distinguish between:

* Expansion vs. intensification
* Different LUC types (forest → cropland vs. grassland → cropland)
* Historical vs. potential land

These are areas for future refinement. See :doc:`environment` for how LUC emissions are calculated.

Configuration Parameters
------------------------

Key configuration parameters for land use aggregation (in ``config/default.yaml``):

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: aggregation ---
   :end-before: # --- section: countries ---

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: primary ---
   :end-before: # --- section: emissions ---

The ``regional_limit`` parameter applies a global constraint on how much of each region's suitable land can be used (0.7 = 70%). This represents institutional, ecological, or social constraints on agricultural expansion.

Workflow Rules
--------------

The land use aggregation workflow consists of three Snakemake rules:

1. **simplify_gadm**: Simplify administrative boundaries

   * Input: ``data/downloads/gadm.gpkg``
   * Output: ``processing/shared/gadm-simplified.gpkg``
   * Script: ``workflow/scripts/simplify_gadm.py``

2. **build_regions**: Cluster into optimization regions

   * Input: Simplified GADM
   * Output: ``processing/{name}/regions.geojson``
   * Script: ``workflow/scripts/build_regions.py``

3. **compute_resource_classes**: Define resource classes

   * Input: Regions + all GAEZ yield rasters
   * Output: ``processing/{name}/resource_classes.nc``
   * Script: ``workflow/scripts/compute_resource_classes.py``

4. **aggregate_class_areas**: Compute land availability

   * Input: Resource classes + suitability rasters + regions
   * Output: ``processing/{name}/land_area_by_class.csv``
   * Script: ``workflow/scripts/aggregate_class_areas.py``

Visualization
-------------

Regional aggregation can be visualized using the regions map plotting rule::

    tools/smk results/{name}/plots/regions_map.pdf

Resource class distribution can be visualized using::

    tools/smk results/{name}/plots/resource_classes_map.pdf

These maps show the spatial distribution of optimization regions and the quality stratification of land within each region.
