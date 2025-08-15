# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import geopandas as gpd
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def download_countries_geojson(output_path: str) -> None:
    """Download Natural Earth 1:10m countries data and save as GeoJSON indexed by alpha3 codes."""

    # Natural Earth 1:10m countries dataset (high resolution)
    url = "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip"

    logger.info(
        "Downloading country boundaries from Natural Earth (1:10m resolution)..."
    )

    # Download and read the data directly with geopandas
    gdf = gpd.read_file(url)

    logger.info("Downloaded %d country features", len(gdf))

    # Manually fix some ISO_A3 codes (see https://github.com/geopandas/geopandas/issues/1041)
    gdf.loc[gdf["SOVEREIGNT"] == "France", "ISO_A3"] = "FRA"
    gdf.loc[gdf["SOVEREIGNT"] == "Norway", "ISO_A3"] = "NOR"
    gdf.loc[gdf["SOVEREIGNT"] == "N. Cyprus", "ISO_A3"] = "CYP"
    gdf.loc[gdf["SOVEREIGNT"] == "Somaliland", "ISO_A3"] = "SOM"
    gdf.loc[gdf["SOVEREIGNT"] == "Kosovo", "ISO_A3"] = "XKX"

    # Filter out countries with invalid ISO_A3 codes
    valid_countries = gdf[gdf["ISO_A3"] != "-99"].copy()

    if len(valid_countries) < len(gdf):
        logger.info(
            "Filtered out %d countries with invalid ISO_A3 codes",
            len(gdf) - len(valid_countries),
        )

    # Group by ISO_A3 to ensure unique entries; join geometries
    valid_countries = valid_countries.dissolve(by="ISO_A3", as_index=False)

    # Set index to alpha3 country codes using ISO_A3
    valid_countries = valid_countries.set_index("ISO_A3")

    # Create output directory if it doesn't exist
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # Save as GeoJSON indexed by alpha3 codes
    valid_countries.to_file(output_path, driver="GeoJSON")

    logger.info("Saved %d country regions to %s", len(valid_countries), output_path)
    logger.debug(
        "Countries indexed by ISO alpha3 codes: %s...",
        ", ".join(sorted(valid_countries.index)[:10]),
    )


def main():
    """Main function called by Snakemake."""
    # Get output path from Snakemake
    output_path = snakemake.output[0]

    # Get config to check region type
    config = snakemake.config

    # For now, only support countries as specified
    if config["aggregation"]["regions"] == "countries":
        download_countries_geojson(output_path)
    else:
        raise NotImplementedError(
            f"Region type '{config.get('regions')}' not supported yet"
        )


if __name__ == "__main__":
    main()
