# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
import logging
from exactextract import exact_extract
from exactextract.raster import NumPyRasterSource
from build_crop_yields import calculate_all_cell_areas

logger = logging.getLogger(__name__)


def calculate_region_areas(crop_files, regions_path):
    """Calculate total cropland area per region using fractional coverage.

    Uses exactextract with NumPyRasterSource to handle grid cells that span
    multiple regions by calculating the fractional coverage of each cell.
    """
    logger.info("Processing %d suitability files", len(crop_files))
    logger.info("Using regions from %s", regions_path)

    # Load regions data
    regions_gdf = gpd.read_file(regions_path)
    logger.info("Loaded %d regions", len(regions_gdf))

    # Process first suitability file to get reference parameters
    first_crop = next(iter(crop_files.keys()))
    first_filepath = crop_files[first_crop]

    with rasterio.open(first_filepath) as src:
        height, width = src.shape
        transform = src.transform
        crs = src.crs

        # Ensure regions are in same CRS as raster
        if regions_gdf.crs != crs:
            logger.info("Reprojecting regions from %s to %s", regions_gdf.crs, crs)
            regions_gdf = regions_gdf.to_crs(crs)

        # Calculate cell areas
        cell_areas = calculate_all_cell_areas(src)

        # Initialize max suitability with first crop
        logger.info("Reading first suitability file: %s", first_filepath)
        max_suitability = src.read(1)

        # Handle nodata values and scale
        if src.nodata is not None:
            max_suitability = np.where(
                max_suitability == src.nodata, np.nan, max_suitability
            )
        max_suitability = max_suitability / 10000.0  # Scale to 0-1

    # Process remaining crops to find maximum suitability per pixel
    for crop, filepath in list(crop_files.items())[1:]:
        logger.info("Processing suitability file for %s: %s", crop, filepath)

        with rasterio.open(filepath) as src:
            suitability = src.read(1)

            # Handle nodata values and scale
            if src.nodata is not None:
                suitability = np.where(suitability == src.nodata, np.nan, suitability)
            suitability = suitability / 10000.0  # Scale to 0-1

            # Take maximum across crops
            max_suitability = np.fmax(max_suitability, suitability)

    # Calculate cropland area per pixel: suitability × cell area
    cropland_area = max_suitability * cell_areas

    # Create bounds from transform for NumPyRasterSource
    xmin = transform.c
    ymax = transform.f
    xmax = xmin + width * transform.a
    ymin = ymax + height * transform.e

    # Create NumPyRasterSource directly from the numpy array
    raster_source = NumPyRasterSource(
        cropland_area,
        xmin=xmin,
        ymin=ymin,
        xmax=xmax,
        ymax=ymax,
        nodata=np.nan,
        srs_wkt=crs.to_wkt() if crs else None,
    )

    # Use exactextract to calculate area-weighted sums per region
    logger.info("Calculating fractional coverage cropland areas per region...")
    region_stats = exact_extract(
        raster_source,
        regions_gdf,
        ["sum"],  # Sum of cropland area values
        include_cols=["region"],  # Include region name
        strategy="raster-sequential",
        output="pandas",  # Return pandas DataFrame
    )

    # Convert to pandas Series with region names as index
    region_areas = pd.Series(
        region_stats["sum"].values, index=region_stats["region"].values
    )

    return region_areas


if __name__ == "__main__":
    # Get input files from snakemake
    crop_files = {
        crop: getattr(snakemake.input, crop) for crop in snakemake.config["crops"]
    }

    region_areas = calculate_region_areas(crop_files, snakemake.input.regions)
    region_areas.name = "cropland_area_ha"
    region_areas.index.name = "region"

    logger.info("Total cropland area: %.2f ha", region_areas.sum())

    # Save to CSV
    region_areas.to_csv(snakemake.output[0])
