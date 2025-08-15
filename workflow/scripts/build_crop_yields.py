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
    tif_path: str,
    suitability_path: str,
    regions_path: str,
    resource_class_quantiles: list,
) -> pd.DataFrame:
    """
    Aggregate crop yields from a GeoTIFF file over regions and resource classes.

    Args:
        tif_path: Path to the GAEZ crop yield GeoTIFF file
        suitability_path: Path to the GAEZ suitability GeoTIFF file
        regions_path: Path to regions GeoJSON file (indexed by alpha3 codes)
        resource_class_quantiles: List of quantile thresholds for resource classes

    Returns:
        DataFrame with MultiIndex (regions, resource_classes) containing yield and area
    """
    print(f"Processing yield data from {tif_path}")
    print(f"Processing suitability data from {suitability_path}")
    print(f"Using regions from {regions_path}")
    print(f"Resource class quantiles: {resource_class_quantiles}")

    # Load regions (countries) data
    regions_gdf = gpd.read_file(regions_path)

    # Ensure the index is set to alpha3 codes (should already be from build_regions.py)
    if regions_gdf.index.name != "ISO_A3":
        if "ISO_A3" in regions_gdf.columns:
            regions_gdf = regions_gdf.set_index("ISO_A3")
        else:
            raise ValueError("Regions GeoJSON must have ISO_A3 column or index")

    print(f"Loaded {len(regions_gdf)} regions")

    # Create transformer for equal-area projection (let it fail if projection not available)
    try:
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:54009", always_xy=True)
        print("Using World Mollweide projection (EPSG:54009)")
    except Exception:
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:3410", always_xy=True)
        print("Using World Cylindrical Equal Area projection (EPSG:3410)")

    # Prepare output data structure
    results_data = []

    # Open both yield and suitability rasters
    with (
        rasterio.open(tif_path) as yield_src,
        rasterio.open(suitability_path) as suit_src,
    ):
        print(f"Yield raster shape: {yield_src.shape}, CRS: {yield_src.crs}")
        print(f"Suitability raster shape: {suit_src.shape}, CRS: {suit_src.crs}")

        # Ensure both rasters have the same CRS and shape
        if yield_src.crs != suit_src.crs:
            raise ValueError("Yield and suitability rasters must have the same CRS")
        if yield_src.shape != suit_src.shape:
            raise ValueError("Yield and suitability rasters must have the same shape")

        # Ensure regions are in the same CRS as the rasters
        if regions_gdf.crs != yield_src.crs:
            print(f"Reprojecting regions from {regions_gdf.crs} to {yield_src.crs}")
            regions_gdf = regions_gdf.to_crs(yield_src.crs)

        # Calculate cell areas once for the entire raster (this is the expensive operation)
        print("Calculating cell areas for entire raster...")
        cell_areas_global = calculate_all_cell_areas(yield_src, transformer)
        print("Cell area calculation completed")

        # Process each region
        for alpha3_code, region in regions_gdf.iterrows():
            try:
                # Extract geometry
                geometry = [region.geometry]

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
                    print(f"{alpha3_code}: No valid data")
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
                                    "ISO_A3": alpha3_code,
                                    "resource_class": class_idx,
                                    "yield": mean_yield,
                                    "suitable_area": total_suitable_area,
                                }
                            )

                            print(
                                f"{alpha3_code} class {class_idx}: yield={mean_yield:.2f} t/ha, area={total_suitable_area:.0f} ha"
                            )

                # Explicitly clean up large arrays to prevent memory accumulation
                del yield_data, suit_data, valid_mask, valid_yields, valid_suitability
                del cell_areas, valid_cell_areas, valid_suitable_areas

            except Exception as e:
                print(f"Warning: Error processing {alpha3_code}: {e}")

    # Create DataFrame with MultiIndex
    if results_data:
        results_df = pd.DataFrame(results_data)
        results_df = results_df.set_index(["ISO_A3", "resource_class"])
        return results_df
    else:
        # Return empty DataFrame with correct structure
        return pd.DataFrame(columns=["yield", "suitable_area"]).set_index(
            ["ISO_A3", "resource_class"]
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

    print(
        f"Saved yields for {len(yields_df)} regions and resource classes to {snakemake.output[0]}"
    )
    if not yields_df.empty:
        print(
            f"Yield statistics: min={yields_df['yield'].min():.2f}, "
            f"max={yields_df['yield'].max():.2f}, "
            f"mean={yields_df['yield'].mean():.2f}"
        )
        print(
            f"Area statistics: min={yields_df['suitable_area'].min():.0f}, "
            f"max={yields_df['suitable_area'].max():.0f}, "
            f"total={yields_df['suitable_area'].sum():.0f} ha"
        )
