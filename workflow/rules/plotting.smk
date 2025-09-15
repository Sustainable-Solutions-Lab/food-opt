# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later


gaez = config["data"]["gaez"]


def yield_gap_inputs(wildcards):
    return {
        "potential_yield": (
            f"data/downloads/gaez_potential_yield_{gaez['climate_model']}_{gaez['time_period']}_{gaez['rcp']}_{gaez['input_management']}_{gaez['water_supply']}_{gaez['co2_fertilization']}_{gaez['crops'][wildcards.crop]}.tif"
        ),
        "actual_yield": (
            f"data/downloads/gaez_actual_yield_{gaez['actual_yield_year']}_{gaez['water_supply']}_{gaez['crops'][wildcards.crop]}.tif"
        ),
    }


rule plot_yield_gap:
    input:
        unpack(yield_gap_inputs),
    output:
        pdf=f"results/{name}/plots/yield_gap_{{crop}}.pdf",
    script:
        "../scripts/plot_yield_gap.py"


rule plot_regions_map:
    input:
        regions=f"processing/{name}/regions.geojson",
    output:
        pdf=f"results/{name}/plots/regions_map.pdf",
    script:
        "../scripts/plot_regions_map.py"


rule plot_results:
    input:
        network="results/{name}/solved/model.nc",
    output:
        crop_pdf="results/{name}/plots/crop_production.pdf",
        food_pdf="results/{name}/plots/food_production.pdf",
        resource_pdf="results/{name}/plots/resource_usage.pdf",
        crop_csv="results/{name}/plots/crop_production.csv",
        food_csv="results/{name}/plots/food_production.csv",
    params:
        output_dir="results/{name}/plots",
    script:
        "../scripts/plot_results.py"


rule plot_crop_production_map:
    input:
        network=f"results/{name}/solved/model.nc",
        regions=f"processing/{name}/regions.geojson",
    output:
        pdf=f"results/{name}/plots/crop_production_map.pdf",
    script:
        "../scripts/plot_crop_production_map.py"


rule plot_cropland_fraction_map:
    input:
        network=f"results/{name}/solved/model.nc",
        regions=f"processing/{name}/regions.geojson",
        land_area_by_class=f"processing/{name}/land_area_by_class.csv",
    output:
        pdf=f"results/{name}/plots/cropland_fraction_map.pdf",
    script:
        "../scripts/plot_cropland_fraction_map.py"


rule plot_average_yield_gap_by_country:
    input:
        regions=f"processing/{name}/regions.geojson",
        csv=f"processing/{name}/yield_gap_by_country_all_crops.csv",
    output:
        pdf=f"results/{name}/plots/yield_gap_by_country_average.pdf",
    script:
        "../scripts/plot_yield_gap_by_country_average.py"
