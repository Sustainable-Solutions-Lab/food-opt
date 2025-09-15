# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later


def yield_gap_raster_inputs(wildcards):
    return {
        "potential_yield": _gaez_path(
            "yield", gaez_cfg["water_supply"], gaez_cfg["crops"][wildcards.crop]
        ),
        "actual_yield": _gaez_path(
            "actual_yield", gaez_cfg["water_supply"], gaez_cfg["crops"][wildcards.crop]
        ),
    }


# Average actual/potential by country, using regions dissolved by 'country'.
rule yield_gap_by_country:
    input:
        unpack(yield_gap_raster_inputs),
        regions=f"processing/{name}/regions.geojson",
    params:
        countries=config["countries"],
    output:
        csv=f"processing/{name}/yield_gap_by_country_{{crop}}.csv",
    script:
        "scripts/compute_yield_gap_by_country.py"


def yield_gap_country_csvs(wildcards):
    # Per-crop country CSVs produced by rule yield_gap_by_country
    return [
        f"processing/{name}/yield_gap_by_country_{crop}.csv"
        for crop in config["crops"]
        if gaez_cfg["crops"][crop] in gaez_cfg["actual_yield_crops"]
    ]


# Average actual/potential by country across all crops
rule average_yield_gap_by_country:
    input:
        yield_gap_country_csvs,
    output:
        csv=f"processing/{name}/yield_gap_by_country_all_crops.csv",
    script:
        "scripts/aggregate_yield_gap_all_crops.py"
