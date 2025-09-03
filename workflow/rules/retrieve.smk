# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later


rule download_gadm_data:
    output:
        "data/downloads/gadm41_{country}_1.json.zip",
    shell:
        "wget -O {output} https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_{wildcards.country}_1.json.zip"


rule download_gaez_yield_data:
    output:
        "data/downloads/gaez_yield_{climate_model}_{time_period}_{rcp}_{input_management}_{water_supply}_{co2_fertilization}_{crop}.tif",
    params:
        url=lambda w: f"https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/res05/{w.climate_model}/rcp{w.rcp}/{w.time_period}H/{config['data']['gaez']['yield_var']}{w.input_management}{w.water_supply}{w.co2_fertilization}_{w.crop}.tif",
    shell:
        "wget -O {output} {params.url}"


rule download_gaez_suitability_data:
    output:
        "data/downloads/gaez_suitability_{climate_model}_{time_period}_{rcp}_{input_management}_{water_supply}_{co2_fertilization}_{crop}.tif",
    params:
        url=lambda w: f"https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/res05/{w.climate_model}/rcp{w.rcp}/{w.time_period}H/{config['data']['gaez']['suitability_var']}{w.input_management}{w.water_supply}{w.co2_fertilization}_{w.crop}.tif",
    shell:
        "wget -O {output} {params.url}"
