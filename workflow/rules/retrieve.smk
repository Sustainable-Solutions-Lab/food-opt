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


rule download_gaez_yield_data:
    output:
        "data/downloads/gaez_yield_{climate_model}_{period}_{scenario}_{input_level}_{water_supply}_{crop}.tif",
    params:
        # GAEZ v5 filename: GAEZ-V5.{VARIABLE}.{PERIOD}.{CLIMATE}.{SCENARIO}.{CROP}.{INPUT}.tif
        # INPUT = {input_level}{water_supply}LM (e.g., HILM, HRLM)
        gcs_url=lambda w: (
            f"gs://fao-gismgr-gaez-v5-data/DATA/GAEZ-V5/MAPSET/{config['data']['gaez']['yield_var']}/"
            f"GAEZ-V5.{config['data']['gaez']['yield_var']}."
            f"{w.period}.{w.climate_model}.{w.scenario}."
            f"{get_gaez_code(w.crop, 'res05')}.{w.input_level}{w.water_supply.upper()}LM.tif"
        ),
    shell:
        "uv run gsutil cp {params.gcs_url} {output}"


rule download_gaez_water_requirement_data:
    output:
        "data/downloads/gaez_water_{climate_model}_{period}_{scenario}_{input_level}_{water_supply}_{crop}.tif",
    params:
        gcs_url=lambda w: (
            f"gs://fao-gismgr-gaez-v5-data/DATA/GAEZ-V5/MAPSET/{config['data']['gaez']['water_requirement_var']}/"
            f"GAEZ-V5.{config['data']['gaez']['water_requirement_var']}."
            f"{w.period}.{w.climate_model}.{w.scenario}."
            f"{get_gaez_code(w.crop, 'res05')}.{w.input_level}{w.water_supply.upper()}LM.tif"
        ),
    shell:
        "uv run gsutil cp {params.gcs_url} {output}"


rule download_gaez_suitability_data:
    output:
        "data/downloads/gaez_suitability_{climate_model}_{period}_{scenario}_{input_level}_{water_supply}_{crop}.tif",
    params:
        gcs_url=lambda w: (
            f"gs://fao-gismgr-gaez-v5-data/DATA/GAEZ-V5/MAPSET/{config['data']['gaez']['suitability_var']}/"
            f"GAEZ-V5.{config['data']['gaez']['suitability_var']}."
            f"{w.period}.{w.climate_model}.{w.scenario}."
            f"{get_gaez_code(w.crop, 'res05')}.{w.input_level}{w.water_supply.upper()}LM.tif"
        ),
    shell:
        "uv run gsutil cp {params.gcs_url} {output}"


rule download_gaez_growing_season_start:
    output:
        "data/downloads/gaez_growing_season_start_{climate_model}_{period}_{scenario}_{input_level}_{water_supply}_{crop}.tif",
    params:
        # RES02-CBD: Beginning of crop growth cycle (day)
        gcs_url=lambda w: (
            f"gs://fao-gismgr-gaez-v5-data/DATA/GAEZ-V5/MAPSET/RES02-CBD/"
            f"GAEZ-V5.RES02-CBD."
            f"{w.period}.{w.climate_model}.{w.scenario}."
            f"{get_gaez_code(w.crop, 'res02')}.{w.input_level}{w.water_supply.upper()}LM.tif"
        ),
    shell:
        "uv run gsutil cp {params.gcs_url} {output}"


rule download_gaez_growing_season_length:
    output:
        "data/downloads/gaez_growing_season_length_{climate_model}_{period}_{scenario}_{input_level}_{water_supply}_{crop}.tif",
    params:
        # RES02-CYL: Length of crop growth cycle (days)
        gcs_url=lambda w: (
            f"gs://fao-gismgr-gaez-v5-data/DATA/GAEZ-V5/MAPSET/RES02-CYL/"
            f"GAEZ-V5.RES02-CYL."
            f"{w.period}.{w.climate_model}.{w.scenario}."
            f"{get_gaez_code(w.crop, 'res02')}.{w.input_level}{w.water_supply.upper()}LM.tif"
        ),
    shell:
        "uv run gsutil cp {params.gcs_url} {output}"


rule download_gaez_actual_yield:
    output:
        "data/downloads/gaez_actual_yield_{water_supply}_{crop}.tif",
    params:
        # RES06-YLD: Actual yields (2010-2019 average)
        # INPUT codes: WSI (irrigated), WSR (rainfed), WST (total)
        # Note: Uses different input naming convention than RES05
        gcs_url=lambda w: (
            f"gs://fao-gismgr-gaez-v5-data/DATA/GAEZ-V5/MAPSET/RES06-YLD/"
            f"GAEZ-V5.RES06-YLD.{get_gaez_code(w.crop, 'res06').lower()}."
            f"WS{w.water_supply.upper()}.tif"
        ),
    shell:
        "uv run gsutil cp {params.gcs_url} {output}"


rule download_gaez_irrigated_landshare_map:
    output:
        "data/downloads/gaez_land_equipped_for_irrigation_share.tif",
    params:
        # LR-IRR: Share of land area equipped for irrigation
        gcs_url="gs://fao-gismgr-gaez-v5-data/DATA/GAEZ-V5/MAP/GAEZ-V5.LR-IRR.tif",
    shell:
        "uv run gsutil cp {params.gcs_url} {output}"





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
