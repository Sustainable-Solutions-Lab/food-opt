# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
from rasterio.mask import mask
from pathlib import Path
from pyproj import Transformer, Geod
import logging

logger = logging.getLogger(__name__)


def calculate_all_cell_areas(src, transformer):
    """Calculate areas for all cells in the raster using geographic coordinates."""
    # Get pixel size in degrees
    pixel_width_deg = abs(src.transform.a)  # degrees longitude
    pixel_height_deg = abs(src.transform.e)  # degrees latitude

    # Create arrays for all pixel centers
    rows, cols = src.shape

    # Get bounds of the raster
    left, bottom, right, top = src.bounds

    # Calculate latitude for each row (pixel centers)
    lats = np.linspace(top - pixel_height_deg / 2, bottom + pixel_height_deg / 2, rows)

    # Use geodesic calculations for accurate area
    geod = Geod(ellps="WGS84")

    # Calculate area for each latitude band
    areas_ha = np.zeros(rows)

    for i, lat in enumerate(lats):
        # Calculate the area of one pixel at this latitude
        # Define corners of a single pixel
        lat_top = lat + pixel_height_deg / 2
        lat_bottom = lat - pixel_height_deg / 2
        lon_left = left  # Any longitude will do for area calculation
        lon_right = left + pixel_width_deg

        # Calculate area using geodesic polygon area
        lons = [lon_left, lon_right, lon_right, lon_left, lon_left]
        lats_poly = [lat_bottom, lat_bottom, lat_top, lat_top, lat_bottom]

        area_m2, _ = geod.polygon_area_perimeter(lons, lats_poly)
        areas_ha[i] = abs(area_m2) / 10000  # Convert to hectares

    # Create 2D array with correct area per latitude
    area_array = np.repeat(areas_ha[:, np.newaxis], cols, axis=1)

    return area_array


def aggregate_yields_by_region(
    yield_path: str,
    suitability_path: str,
    regions_path: str,
    resource_class_quantiles: list,
) -> pd.DataFrame:
    """
    Aggregate crop yields from a GeoTIFF file over regions and resource classes.

    Args:
        yield_path: Path to the GAEZ crop yield GeoTIFF file
        suitability_path: Path to the GAEZ suitability GeoTIFF file
        regions_path: Path to regions GeoJSON file
        resource_class_quantiles: List of quantile thresholds for resource classes

    Returns:
        DataFrame with MultiIndex (regions, resource_classes) containing yield and area
    """
    logger.info("Processing yield data from %s", yield_path)
    logger.info("Processing suitability data from %s", suitability_path)
    logger.info("Using regions from %s", regions_path)
    logger.info("Resource class quantiles: %s", resource_class_quantiles)

    # Load regions (countries) data
    regions_gdf = gpd.read_file(regions_path)
    regions_gdf.set_index("region", inplace=True)
    logger.info("Loaded %d regions", len(regions_gdf))

    # Create transformer for equal-area projection (let it fail if projection not available)
    try:
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:54009", always_xy=True)
        logger.info("Using World Mollweide projection (EPSG:54009)")
    except Exception:
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:3410", always_xy=True)
        logger.info("Using World Cylindrical Equal Area projection (EPSG:3410)")

    # Prepare output data structure
    results_data = []

    # Open both yield and suitability rasters
    with (
        rasterio.open(yield_path) as yield_src,
        rasterio.open(suitability_path) as suit_src,
    ):
        logger.info("Yield raster shape: %s, CRS: %s", yield_src.shape, yield_src.crs)
        logger.info(
            "Suitability raster shape: %s, CRS: %s", suit_src.shape, suit_src.crs
        )

        # Ensure both rasters have the same CRS and shape
        if yield_src.crs != suit_src.crs:
            raise ValueError("Yield and suitability rasters must have the same CRS")
        if yield_src.shape != suit_src.shape:
            raise ValueError("Yield and suitability rasters must have the same shape")

        # Ensure regions are in the same CRS as the rasters
        if regions_gdf.crs != yield_src.crs:
            logger.info(
                "Reprojecting regions from %s to %s", regions_gdf.crs, yield_src.crs
            )
            regions_gdf = regions_gdf.to_crs(yield_src.crs)

        # Calculate cell areas once for the entire raster (this is the expensive operation)
        logger.info("Calculating cell areas for entire raster...")
        cell_areas_global = calculate_all_cell_areas(yield_src, transformer)
        logger.info("Cell area calculation completed")

        # Process each region
        for region, row in regions_gdf.iterrows():
            try:
                # Extract geometry
                geometry = [row.geometry]

                # Mask both rasters with the region geometry
                yield_masked, yield_transform = mask(
                    yield_src, geometry, crop=True, nodata=yield_src.nodata
                )
                suit_masked, suit_transform = mask(
                    suit_src, geometry, crop=True, nodata=suit_src.nodata
                )

                # Get the masked data (first band)
                yield_data = yield_masked[0]
                suit_data = suit_masked[0]

                # Filter out nodata values for both datasets
                valid_mask = np.ones_like(yield_data, dtype=bool)

                if yield_src.nodata is not None:
                    valid_mask &= yield_data != yield_src.nodata
                else:
                    valid_mask &= ~np.isnan(yield_data)
                    valid_mask &= yield_data > 0  # GAEZ yields should be positive

                if suit_src.nodata is not None:
                    valid_mask &= suit_data != suit_src.nodata
                else:
                    valid_mask &= ~np.isnan(suit_data)
                    valid_mask &= suit_data >= 0  # Suitability should be non-negative

                if not np.any(valid_mask):
                    logger.debug("%s: No valid data", region)
                    continue

                # Extract valid data
                valid_yields = yield_data[valid_mask] / 1000  # Convert kg/ha to t/ha
                valid_suitability = suit_data[valid_mask]

                # Create a temporary raster from cell areas and apply the same mask
                with rasterio.MemoryFile() as memfile:
                    with memfile.open(
                        driver="GTiff",
                        height=cell_areas_global.shape[0],
                        width=cell_areas_global.shape[1],
                        count=1,
                        dtype=cell_areas_global.dtype,
                        crs=yield_src.crs,
                        transform=yield_src.transform,
                        nodata=-9999,
                    ) as temp_src:
                        temp_src.write(cell_areas_global, 1)

                    with memfile.open() as temp_src:
                        # Use the same mask operation as for yield data
                        masked_cell_areas, _ = mask(
                            temp_src, geometry, crop=True, nodata=-9999
                        )
                        cell_areas = masked_cell_areas[0]

                valid_cell_areas = cell_areas[valid_mask]
                # Convert suitability from scale 0-10000 to fraction 0-1
                valid_suitability_fraction = valid_suitability / 10000
                valid_suitable_areas = valid_cell_areas * valid_suitability_fraction

                # Create resource classes based on yield quantiles
                if len(valid_yields) > 0:
                    yield_quantiles = np.quantile(
                        valid_yields, [0] + resource_class_quantiles + [1]
                    )

                    for class_idx in range(len(yield_quantiles) - 1):
                        # Find pixels in this yield quantile range
                        class_mask = (valid_yields >= yield_quantiles[class_idx]) & (
                            valid_yields < yield_quantiles[class_idx + 1]
                        )

                        # Handle the last class to include the maximum value
                        if class_idx == len(yield_quantiles) - 2:
                            class_mask = (
                                valid_yields >= yield_quantiles[class_idx]
                            ) & (valid_yields <= yield_quantiles[class_idx + 1])

                        if np.any(class_mask):
                            class_yields = valid_yields[class_mask]
                            class_suitable_areas = valid_suitable_areas[class_mask]

                            mean_yield = np.mean(class_yields)
                            total_suitable_area = np.sum(class_suitable_areas)

                            results_data.append(
                                {
                                    "region": region,
                                    "resource_class": class_idx,
                                    "yield": mean_yield,
                                    "suitable_area": total_suitable_area,
                                }
                            )

                            logger.debug(
                                "%s class %d: yield=%.2f t/ha, area=%.0f ha",
                                region,
                                class_idx,
                                mean_yield,
                                total_suitable_area,
                            )

                # Explicitly clean up large arrays to prevent memory accumulation
                del yield_data, suit_data, valid_mask, valid_yields, valid_suitability
                del cell_areas, valid_cell_areas, valid_suitable_areas

            except Exception as e:
                logger.warning("Error processing %s: %s", region, e)

    # Create DataFrame with MultiIndex
    if results_data:
        results_df = pd.DataFrame(results_data)
        results_df = results_df.set_index(["region", "resource_class"])
        return results_df
    else:
        # Return empty DataFrame with correct structure
        return pd.DataFrame(columns=["yield", "suitable_area"]).set_index(
            ["region", "resource_class"]
        )


if __name__ == "__main__":
    yields_df = aggregate_yields_by_region(
        snakemake.input.yields,
        snakemake.input.suitability,
        snakemake.input.regions,
        snakemake.config["aggregation"]["resource_class_quantiles"],
    )

    # Save the aggregated yields to CSV
    Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)
    yields_df.to_csv(snakemake.output[0])

    logger.info(
        "Saved yields for %d regions and resource classes to %s",
        len(yields_df),
        snakemake.output[0],
    )
    if not yields_df.empty:
        logger.info(
            "Yield statistics: min=%.2f, max=%.2f, mean=%.2f",
            yields_df["yield"].min(),
            yields_df["yield"].max(),
            yields_df["yield"].mean(),
        )
        logger.info(
            "Area statistics: min=%.0f, max=%.0f, total=%.0f ha",
            yields_df["suitable_area"].min(),
            yields_df["suitable_area"].max(),
            yields_df["suitable_area"].sum(),
        )
