# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later


rule download_gadm_zip:
    output:
        temp("data/downloads/gadm_410-levels.zip"),
    params:
        url="https://geodata.ucdavis.edu/gadm/gadm4.1/gadm_410-levels.zip",
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        curl -L --fail --progress-bar -o "{output}" "{params.url}"
        """


rule extract_adm1:
    input:
        zip="data/downloads/gadm_410-levels.zip",
    output:
        "data/downloads/gadm.gpkg",
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        ogr2ogr -f GPKG "{output}" "/vsizip/{input.zip}/gadm_410-levels.gpkg" ADM_1
        """


rule download_gaez_potential_yield_data:
    output:
        "data/downloads/gaez_potential_yield_{climate_model}_{time_period}_{rcp}_{input_management}_{water_supply}_{co2_fertilization}_{crop}.tif",
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


rule download_gaez_actual_yield_data:
    output:
        "data/downloads/gaez_actual_yield_{year}_{water_supply}_{crop}.tif",
    params:
        url=lambda w: f"https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/res06/{w.water_supply.upper()}/{w.year}/{w.crop}_{w.year}_yld.tif",
    shell:
        "wget -O {output} {params.url}"


rule download_gaez_irrigated_cropland_data:
    output:
        "data/downloads/gaez_irrigated_cropland_share.tif",
    params:
        url="https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/LR/wat/GLCSv11_12_5m.tif",
    shell:
        "wget -O {output} {params.url}"


rule download_wpp_population:
    output:
        "data/downloads/WPP_population.csv.gz",
    params:
        url=(
            "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/CSV_FILES/WPP2024_TotalPopulationBySex.csv.gz"
        ),
    shell:
        r"""
        wget -O {output} "{params.url}"
        """
