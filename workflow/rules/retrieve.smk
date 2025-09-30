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


rule retrieve_faostat_prices:
    input:
        mapping="data/faostat_item_map.csv",
    params:
        crops=config["crops"],
    output:
        prices=f"processing/{name}/faostat_prices.csv",
    script:
        "../scripts/retrieve_faostat_prices.py"


rule download_gaez_potential_yield_data:
    output:
        "data/downloads/gaez_potential_yield_{climate_model}_{time_period}_{rcp}_{input_management}_{water_supply}_{co2_fertilization}_{crop}.tif",
    params:
        url=lambda w: f"https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/res05/{w.climate_model}/rcp{w.rcp}/{w.time_period}H/{config['data']['gaez']['yield_var']}{w.input_management}{w.water_supply}{w.co2_fertilization}_{w.crop}.tif",
    shell:
        "wget -O {output} {params.url}"


rule download_gaez_water_requirement_data:
    output:
        "data/downloads/gaez_water_requirement_{climate_model}_{time_period}_{rcp}_{input_management}_{water_supply}_{co2_fertilization}_{crop}.tif",
    params:
        url=lambda w: f"https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/res05/{w.climate_model}/rcp{w.rcp}/{w.time_period}H/{config['data']['gaez']['water_requirement_var']}{w.input_management}{w.water_supply}{w.co2_fertilization}_{w.crop}.tif",
    shell:
        "wget -O {output} {params.url}"


rule download_gaez_growing_season_start_data:
    output:
        f"data/downloads/gaez_growing_season_start_{config['data']['gaez']['climate_model']}_{config['data']['gaez']['time_period']}_{config['data']['gaez']['rcp']}_{config['data']['gaez']['growing_season']['year']}.tif",
    params:
        url=(
            f"https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/res01/"
            f"{config['data']['gaez']['climate_model']}/rcp{config['data']['gaez']['rcp']}/TS/"
            f"{config['data']['gaez']['growing_season']['start_var']}"
            f"_{config['data']['gaez']['climate_model']}"
            f"_rcp{config['data']['gaez']['rcp']}"
            f"_{config['data']['gaez']['growing_season']['year']}.tif"
        ),
    shell:
        "wget -O {output} {params.url}"


rule download_gaez_growing_season_length_data:
    output:
        f"data/downloads/gaez_growing_season_length_{config['data']['gaez']['climate_model']}_{config['data']['gaez']['time_period']}_{config['data']['gaez']['rcp']}_{config['data']['gaez']['growing_season']['year']}.tif",
    params:
        url=(
            f"https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/res01/"
            f"{config['data']['gaez']['climate_model']}/rcp{config['data']['gaez']['rcp']}/TS/"
            f"{config['data']['gaez']['growing_season']['length_var']}"
            f"_{config['data']['gaez']['climate_model']}"
            f"_rcp{config['data']['gaez']['rcp']}"
            f"_{config['data']['gaez']['growing_season']['year']}.tif"
        ),
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


# TODO: license. Different variations?

# See https://data.isimip.org/search/crop/mgr/variable/yield/irrigation/noirr/


# The following is a future projection, but not about yields but primary productivity
# See https://data.isimip.org/search/simulation_round/ISIMIP2b/sector/biomes/model/lpjml/pft/mgr-rainfed/
# url="https://files.isimip.org/ISIMIP2b/OutputData/biomes/LPJmL/gfdl-esm2m/future/lpjml_gfdl-esm2m_ewembi_rcp26_2005soc_2005co2_gpp-mgr-irrigated_global_annual_2006_2099.nc4",
rule download_grassland_yield_data:
    output:
        "data/downloads/grassland_yield_historical.nc4",
    params:
        url="https://files.isimip.org/ISIMIP2a/OutputData/agriculture/LPJmL/watch/historical/lpjml_watch_nobc_hist_co2_yield-mgr-noirr-default_global_annual_1971_2001.nc4",
    shell:
        "wget -O {output} {params.url}"


rule download_wpp_population:
    output:
        "data/downloads/WPP_population.csv.gz",
    params:
        url=(
            "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/CSV_FILES/WPP2024_Population1JanuaryByAge5GroupSex_Medium.csv.gz"
        ),
    shell:
        r"""
        wget -O {output} "{params.url}"
        """


rule download_waterfootprint_appendix:
    output:
        "data/downloads/Report53_Appendix.zip",
    params:
        url="https://www.waterfootprint.org/resources/appendix/Report53_Appendix.zip",
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        wget -O "{output}" "{params.url}"
        """
