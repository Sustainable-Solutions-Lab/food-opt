# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import geopandas as gpd
import pandas as pd
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def load_gadm_geojson(output_path: str, country_files: dict, gadm_level: int) -> None:
    """Load GADM data at specified level from downloaded files and save as GeoJSON indexed by region codes.

    Args:
        output_path: Path to save the output GeoJSON file
        country_files: Dictionary mapping country codes to GADM file paths
        gadm_level: GADM level (0 for countries, 1 for states/provinces)
    """

    level_name = "countries" if gadm_level == 0 else "states/provinces"
    logger.info(
        f"Loading {level_name} boundaries from GADM level {gadm_level} files..."
    )

    if not country_files:
        raise ValueError("No country files provided")

    all_regions = []

    for country_code, file_path in country_files.items():
        try:
            logger.debug(f"Loading {country_code} from {file_path}")
            country_data = gpd.read_file(file_path)

            if len(country_data) > 0:
                all_regions.append(country_data)
                logger.debug(f"Loaded {len(country_data)} regions for {country_code}")

        except Exception as e:
            logger.debug(f"Failed to load {country_code}: {e}")
            continue

    if not all_regions:
        raise ValueError(f"No level {gadm_level} data could be loaded from files")

    # Combine all region data
    gdf = pd.concat(all_regions, ignore_index=True)

    # Use appropriate GADM field as region identifier
    if gadm_level == 0:
        # For countries: use ISO_A3 (3-letter country code)
        gdf["region"] = gdf["GID_0"]
    elif gadm_level == 1:
        # For states: use GID_1 (country.state identifier)
        gdf["region"] = gdf["GID_1"]
    else:
        raise ValueError(f"Unsupported GADM level: {gadm_level}")

    # Filter out regions with invalid identifiers
    initial_count = len(gdf)
    valid_mask = (gdf["region"] != "?") & gdf["region"].notna() & (gdf["region"] != "")
    gdf = gdf[valid_mask]

    if len(gdf) < initial_count:
        filtered_count = initial_count - len(gdf)
        logger.warning(
            f"Filtered out {filtered_count} regions with invalid identifiers (?, NaN, or empty)"
        )

    # Make sure there are no duplicates
    region_counts = gdf["region"].value_counts()
    duplicate_ids = region_counts[region_counts > 1].index
    assert len(duplicate_ids) == 0

    # Set index to region identifiers
    gdf = gdf.set_index("region")

    # Create output directory if it doesn't exist
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # Save as GeoJSON indexed by region codes
    gdf.to_file(output_path, driver="GeoJSON")

    logger.info(f"Saved {len(gdf)} {level_name} regions to {output_path}")
    logger.debug(
        f"{level_name.capitalize()} indexed by region codes: %s...",
        ", ".join(sorted(gdf.index)[:10]),
    )


if __name__ == "__main__":
    agg_level = snakemake.config["aggregation"]["regions"]

    # Get input files from Snakemake
    country_files = {c: snakemake.input[c] for c in snakemake.config["countries"]}

    if agg_level == "countries":
        # Use GADM level 0 (country boundaries)
        load_gadm_geojson(snakemake.output[0], country_files, gadm_level=0)
    elif agg_level == "states":
        # Use GADM level 1 (state/province boundaries)
        load_gadm_geojson(snakemake.output[0], country_files, gadm_level=1)
    else:
        raise NotImplementedError(
            f"Aggregation level '{agg_level}' not supported. Options: 'countries', 'states'"
        )
