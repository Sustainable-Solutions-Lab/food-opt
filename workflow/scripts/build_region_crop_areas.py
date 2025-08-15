# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
from rasterio.mask import mask
import logging
from build_crop_yields import calculate_all_cell_areas

logger = logging.getLogger(__name__)


def process_suitability_files(crop_files):
    """Process suitability files to find maximum suitability per pixel across crops."""
    logger.info(
        "Processing %d suitability files to find maximum per pixel", len(crop_files)
    )

    max_suitability = None
    pixel_areas = None
    reference_transform = None
    reference_crs = None

    for i, (crop, filepath) in enumerate(crop_files.items()):
        logger.debug("Processing crop %s: %s", crop, filepath)

        with rasterio.open(filepath) as src:
            # Read suitability data
            suitability = src.read(1)

            # Handle nodata values
            suitability = np.where(
                src.nodata is not None and suitability == src.nodata,
                np.nan,
                suitability,
            )

            # Scale from 1-10000 range to 0-1 range
            suitability = suitability / 10000.0

            if max_suitability is None:
                # Initialize with first crop
                max_suitability = suitability.copy()
                pixel_areas = calculate_all_cell_areas(src, None)
                reference_transform = src.transform
                reference_crs = src.crs
            else:
                # Take maximum across crops, handling NaN values
                max_suitability = np.fmax(max_suitability, suitability)

        # Free memory
        del suitability

    return max_suitability, pixel_areas, reference_transform, reference_crs


def calculate_region_areas(max_suitability, pixel_areas, transform, crs, regions_gdf):
    """Calculate total cropland area per region."""
    logger.info("Calculating cropland areas for %d regions", len(regions_gdf))

    # Ensure regions are in same CRS as raster
    if regions_gdf.crs != crs:
        regions_gdf = regions_gdf.to_crs(crs)

    region_areas = {}

    for _, region_row in regions_gdf.iterrows():
        region_name = region_row["region"]
        region_geom = region_row["geometry"]

        logger.debug("Processing region: %s", region_name)

        # Create a temporary raster dataset for masking
        with rasterio.MemoryFile() as memfile:
            with memfile.open(
                driver="GTiff",
                height=max_suitability.shape[0],
                width=max_suitability.shape[1],
                count=1,
                dtype=max_suitability.dtype,
                crs=crs,
                transform=transform,
            ) as dataset:
                dataset.write(max_suitability, 1)

                try:
                    # Mask raster by region geometry
                    masked_suitability, _ = mask(
                        dataset, [region_geom], crop=True, filled=False
                    )
                    masked_suitability = masked_suitability[0]  # Extract from 3D array

                    # Create a temporary dataset for pixel areas and mask it
                    with rasterio.MemoryFile() as memfile_areas:
                        with memfile_areas.open(
                            driver="GTiff",
                            height=pixel_areas.shape[0],
                            width=pixel_areas.shape[1],
                            count=1,
                            dtype=pixel_areas.dtype,
                            crs=crs,
                            transform=transform,
                        ) as areas_dataset:
                            areas_dataset.write(pixel_areas, 1)
                            region_pixel_areas, _ = mask(
                                areas_dataset, [region_geom], crop=True, filled=False
                            )

                    region_pixel_areas = region_pixel_areas[0]  # Extract from 3D array

                    # Calculate total cropland area: suitability fraction � pixel area
                    # Handle masked values (they are np.ma.masked)
                    valid_mask = ~masked_suitability.mask & ~region_pixel_areas.mask
                    valid_suitability = masked_suitability.data[valid_mask]
                    valid_areas = region_pixel_areas.data[valid_mask]

                    # Remove NaN values
                    finite_mask = np.isfinite(valid_suitability) & np.isfinite(
                        valid_areas
                    )
                    if np.any(finite_mask):
                        total_cropland_area = np.sum(
                            valid_suitability[finite_mask] * valid_areas[finite_mask]
                        )
                    else:
                        total_cropland_area = 0.0

                    region_areas[region_name] = total_cropland_area
                    logger.debug(
                        "Region %s: %.2f ha cropland", region_name, total_cropland_area
                    )

                except Exception as e:
                    logger.warning("Failed to process region %s: %s", region_name, e)
                    region_areas[region_name] = 0.0

    return region_areas


def main():
    """Main function to calculate region crop areas."""
    logger.info("Starting region crop areas calculation")

    # Get input files from snakemake
    crop_files = {
        crop: getattr(snakemake.input, crop) for crop in snakemake.config["crops"]
    }
    regions_file = snakemake.input.regions
    output_file = snakemake.output[0]

    logger.info("Input crop suitability files: %d", len(crop_files))
    logger.info("Regions file: %s", regions_file)
    logger.info("Output file: %s", output_file)

    # Load regions
    regions_gdf = gpd.read_file(regions_file)
    logger.info("Loaded %d regions", len(regions_gdf))

    # Process suitability files to get maximum suitability per pixel
    max_suitability, pixel_areas, transform, crs = process_suitability_files(crop_files)

    # Calculate region areas
    region_areas = calculate_region_areas(
        max_suitability, pixel_areas, transform, crs, regions_gdf
    )

    # Create output DataFrame
    result_df = pd.DataFrame(
        [
            {"region": region, "cropland_area_ha": area}
            for region, area in region_areas.items()
        ]
    )

    # Sort by region name for consistency
    result_df = result_df.sort_values("region").reset_index(drop=True)

    logger.info("Writing results to %s", output_file)
    logger.info("Total regions: %d", len(result_df))
    logger.info("Total cropland area: %.2f ha", result_df["cropland_area_ha"].sum())

    # Save to CSV
    result_df.to_csv(output_file, index=False)

    logger.info("Region crop areas calculation completed")


if __name__ == "__main__":
    main()
