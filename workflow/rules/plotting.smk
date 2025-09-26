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


rule plot_resource_classes_map:
    input:
        classes=f"processing/{name}/resource_classes.nc",
        regions=f"processing/{name}/regions.geojson",
    output:
        pdf=f"results/{name}/plots/resource_classes_map.pdf",
    script:
        "../scripts/plot_resource_classes_map.py"


rule plot_results:
    input:
        network="results/{name}/solved/model.nc",
    output:
        crop_pdf="results/{name}/plots/crop_production.pdf",
        resource_pdf="results/{name}/plots/resource_usage.pdf",
        crop_csv="results/{name}/plots/crop_production.csv",
        food_csv="results/{name}/plots/food_production.csv",
    params:
        output_dir="results/{name}/plots",
    script:
        "../scripts/plot_results.py"


rule plot_health_impacts:
    input:
        network="results/{name}/solved/model.nc",
        regions=f"processing/{name}/regions.geojson",
        risk_breakpoints=f"processing/{name}/health/risk_breakpoints.csv",
        health_cluster_cause=f"processing/{name}/health/cluster_cause_baseline.csv",
        health_cause_log=f"processing/{name}/health/cause_log_breakpoints.csv",
        health_cluster_summary=f"processing/{name}/health/cluster_summary.csv",
        health_clusters=f"processing/{name}/health/country_clusters.csv",
        food_risk_map="data/health/food_to_risk_factor.csv",
        population=f"processing/{name}/population.csv",
    params:
        ghg_price=config["emissions"]["ghg_price"],
    output:
        breakdown_pdf="results/{name}/plots/objective_breakdown.pdf",
        breakdown_csv="results/{name}/plots/objective_breakdown.csv",
        health_map_pdf="results/{name}/plots/health_risk_map.pdf",
        health_map_csv="results/{name}/plots/health_risk_by_region.csv",
    script:
        "../scripts/plot_health_impacts.py"


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
