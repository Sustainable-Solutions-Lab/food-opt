# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
from pathlib import Path
from pyproj import Geod
import logging
import xarray as xr
from rasterio import features
from exactextract import exact_extract
from exactextract.raster import NumPyRasterSource

logger = logging.getLogger(__name__)


def calculate_all_cell_areas(src):
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
    Aggregate crop yields from a GeoTIFF file over regions and resource classes using vectorized operations.

    This function reads yield and suitability rasters, creates a region mask raster,
    and calculates mean yields and suitable areas for each resource class within each region
    using vectorized operations instead of looping over regions.

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

        cell_areas_global = calculate_all_cell_areas(yield_src)

        # Read full rasters into xarray
        logger.info("Reading raster data into xarray...")
        yield_data = yield_src.read(1)
        suit_data = suit_src.read(1)

        # Create xarray DataArrays with proper coordinates
        height, width = yield_data.shape
        transform = yield_src.transform

        # Create coordinate arrays
        x_coords = np.arange(width) * transform.a + transform.c + transform.a / 2
        y_coords = np.arange(height) * transform.e + transform.f + transform.e / 2

        # Dataset with both yields, suitability and area
        data = xr.Dataset(
            {
                "yield_raw": (["y", "x"], yield_data, {"nodata": yield_src.nodata}),
                "suitability_raw": (["y", "x"], suit_data, {"nodata": suit_src.nodata}),
                "cell_areas": (["y", "x"], cell_areas_global),
            },
            coords={"y": y_coords, "x": x_coords},
        )

        # Create region ID raster using rasterio.features.rasterize
        logger.info("Creating region mask raster...")
        region_shapes = [(geom, idx) for idx, geom in enumerate(regions_gdf.geometry)]
        region_raster = features.rasterize(
            region_shapes,
            out_shape=(height, width),
            transform=transform,
            fill=-1,  # Use -1 for areas outside any region
            dtype=np.int32,
        )

        # Add region raster to dataset
        data = data.assign(region=(("y", "x"), region_raster))

        # Create validity masks and convert units in the dataset
        logger.info("Creating validity masks and converting units...")

        # Check for nodata values
        yield_nodata = data["yield_raw"].attrs.get("nodata")
        suit_nodata = data["suitability_raw"].attrs.get("nodata")

        yield_valid = (
            data["yield_raw"] != yield_nodata
            if yield_nodata is not None
            else ~np.isnan(data["yield_raw"])
        ) & (data["yield_raw"] > 0)

        suit_valid = (
            data["suitability_raw"] != suit_nodata
            if suit_nodata is not None
            else ~np.isnan(data["suitability_raw"])
        ) & (data["suitability_raw"] >= 0)

        # Unit conversions
        data["yield"] = data["yield_raw"] / 1000  # kg/ha to t/ha
        data["suitability"] = data["suitability_raw"] / 10000  # 0-10000 to 0-1
        data["suitable_area"] = data["cell_areas"] * data["suitability"]
        data["valid"] = yield_valid & suit_valid & (data["region"] >= 0)

        # Calculate quantiles by region, then map back to full resolution threshold values
        quantiles = [0] + resource_class_quantiles + [1]
        region_quantiles = data["yield"].groupby(data["region"]).quantile(quantiles)
        thresholds = region_quantiles.sel(region=data.region)  # Has y and x coordinates
        thresholds = thresholds.reset_coords(drop=True)  # Gets rid of region coordinate

        # Resource class assignment using thresholds
        resource_classes = xr.zeros_like(data["yield"], dtype="int8")
        for class_idx in range(len(quantiles) - 1):
            lower_bound = thresholds.isel(quantile=class_idx)
            upper_bound = thresholds.isel(quantile=class_idx + 1)

            # Handle duplicate thresholds by using only the lower bound
            mask = (data["yield"] >= lower_bound) & (data["yield"] < upper_bound)

            # For the highest class, include the upper bound
            if class_idx == len(quantiles) - 2:
                mask = mask | (data["yield"] >= upper_bound)

            resource_classes = resource_classes.where(~mask, class_idx)

        data["resource_class"] = resource_classes

        logger.info("Calculating fractional coverage statistics with exactextract...")

        # Create bounds from transform for NumPyRasterSource
        transform = yield_src.transform
        xmin = transform.c
        ymax = transform.f
        xmax = xmin + width * transform.a
        ymin = ymax + height * transform.e
        crs_wkt = yield_src.crs.to_wkt() if yield_src.crs else None

        results = []

        # Process each resource class separately, but all regions in parallel
        for resource_class_id in range(len(quantiles) - 1):
            # Create masks for this resource class across all regions
            class_mask = (data["resource_class"] == resource_class_id) & data["valid"]

            # Skip if no valid data for this resource class
            if class_mask.sum() == 0:
                logger.info(
                    "No valid data for resource class %d, skipping", resource_class_id
                )
                continue

            # Create composite rasters for this resource class (all regions)
            yield_composite = np.where(class_mask, data["yield"], np.nan)
            area_composite = np.where(class_mask, data["suitable_area"], np.nan)

            # Create NumPyRasterSource objects
            yield_raster = NumPyRasterSource(
                yield_composite,
                xmin=xmin,
                ymin=ymin,
                xmax=xmax,
                ymax=ymax,
                nodata=np.nan,
                srs_wkt=crs_wkt,
            )
            area_raster = NumPyRasterSource(
                area_composite,
                xmin=xmin,
                ymin=ymin,
                xmax=xmax,
                ymax=ymax,
                nodata=np.nan,
                srs_wkt=crs_wkt,
            )

            # Use exactextract to calculate statistics for all regions at once
            # Reset index so 'region' is available as a column for exactextract
            regions_for_extract = regions_gdf.reset_index()

            yield_stats = exact_extract(
                yield_raster,
                regions_for_extract,
                ["mean"],
                include_cols=["region"],
                output="pandas",
            )
            area_stats = exact_extract(
                area_raster,
                regions_for_extract,
                ["sum"],
                include_cols=["region"],
                output="pandas",
            )

            # Combine results for this resource class
            if not yield_stats.empty and not area_stats.empty:
                for idx, row in yield_stats.iterrows():
                    if (
                        not pd.isna(row["mean"])
                        and not pd.isna(row["region"])
                        and row["region"] in area_stats["region"].values
                    ):
                        area_matches = area_stats[
                            area_stats["region"] == row["region"]
                        ]["sum"]
                        if len(area_matches) == 0:
                            logger.warning(
                                f"No area data found for region {row['region']}"
                            )
                            continue
                        area_value = area_matches.iloc[0]
                        results.append(
                            {
                                "region": row["region"],
                                "resource_class": resource_class_id,
                                "yield": row["mean"],
                                "suitable_area": area_value,
                            }
                        )

    results_df = pd.DataFrame(results)
    # results_df = results_df.loc[results_df.region != "NA"]
    results_df = results_df.set_index(["region", "resource_class"])
    return results_df


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
