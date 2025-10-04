# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later


gaez = config["data"]["gaez"]


def _gaez_actual_yield_raster_path(crop_name: str, water_supply: str) -> str:
    # Wrap helper to provide clearer error message for plotting context.
    try:
        return gaez_path("actual_yield", water_supply, crop_name)
    except ValueError as exc:
        raise ValueError(
            f"Missing RES06 actual yield data for crop '{crop_name}'."
        ) from exc


def yield_gap_raster_inputs(wildcards):
    crop_name = wildcards.crop
    ws = str(gaez["water_supply"]).lower()
    return {
        "potential_yield": gaez_path("yield", ws, crop_name),
        "actual_yield": _gaez_actual_yield_raster_path(crop_name, ws),
    }


rule plot_yield_gap:
    input:
        unpack(yield_gap_raster_inputs),
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
        health_cluster_risk_baseline=f"processing/{name}/health/cluster_risk_baseline.csv",
        food_risk_map="data/health/food_to_risk_factor.csv",
        population=f"processing/{name}/population.csv",
    params:
        ghg_price=config["emissions"]["ghg_price"],
    output:
        breakdown_pdf="results/{name}/plots/objective_breakdown.pdf",
        breakdown_csv="results/{name}/plots/objective_breakdown.csv",
        health_map_pdf="results/{name}/plots/health_risk_map.pdf",
        health_map_csv="results/{name}/plots/health_risk_by_region.csv",
        health_baseline_map_pdf="results/{name}/plots/health_baseline_map.pdf",
        health_baseline_map_csv="results/{name}/plots/health_baseline_by_region.csv",
    script:
        "../scripts/plot_health_impacts.py"


rule plot_crop_production_map:
    input:
        network=f"results/{name}/solved/model.nc",
        regions=f"processing/{name}/regions.geojson",
    output:
        production_pdf=f"results/{name}/plots/crop_production_map.pdf",
        land_pdf=f"results/{name}/plots/crop_land_use_map.pdf",
    script:
        "../scripts/plot_crop_production_map.py"


rule plot_crop_use_breakdown:
    input:
        network=f"results/{name}/solved/model.nc",
    output:
        pdf=f"results/{name}/plots/crop_use_breakdown.pdf",
        csv=f"results/{name}/plots/crop_use_breakdown.csv",
    script:
        "../scripts/plot_crop_use_breakdown.py"


def yield_map_inputs(wildcards):
    if wildcards.item == "pasture":
        return {"raster": "data/downloads/grassland_yield_historical.nc4"}
    else:
        return {
            "raster": gaez_path("yield", wildcards.water_supply, wildcards.item),
            "conversions": "data/yield_unit_conversions.csv",
        }


rule plot_yield_map:
    input:
        unpack(yield_map_inputs),
    output:
        pdf=f"results/{name}/plots/yield_map_{{item}}_{{water_supply}}.pdf",
    params:
        gaez=gaez,
        item=lambda wc: wc.item,
        supply=lambda wc: wc.water_supply,
        unit="t/ha",
        cmap="YlGn",
    script:
        "../scripts/plot_yield_map.py"


rule plot_cropland_fraction_map:
    input:
        network=f"results/{name}/solved/model.nc",
        regions=f"processing/{name}/regions.geojson",
        land_area_by_class=f"processing/{name}/land_area_by_class.csv",
        resource_classes=f"processing/{name}/resource_classes.nc",
    output:
        pdf=f"results/{name}/plots/cropland_fraction_map.pdf",
    script:
        "../scripts/plot_cropland_fraction_map.py"


rule plot_irrigated_cropland_fraction_map:
    input:
        network=f"results/{name}/solved/model.nc",
        regions=f"processing/{name}/regions.geojson",
        land_area_by_class=f"processing/{name}/land_area_by_class.csv",
        resource_classes=f"processing/{name}/resource_classes.nc",
    output:
        pdf=f"results/{name}/plots/irrigated_cropland_fraction_map.pdf",
    params:
        water_supply="i",
        title="Irrigated Cropland Fraction by Region and Resource Class",
        colorbar_label="Irrigated cropland / irrigable land area",
        csv_prefix="irrigated_cropland",
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


rule plot_water_value_map:
    input:
        network=f"results/{name}/solved/model.nc",
        regions=f"processing/{name}/regions.geojson",
    output:
        pdf=f"results/{name}/plots/water_value_map.pdf",
    script:
        "../scripts/plot_water_value_map.py"
