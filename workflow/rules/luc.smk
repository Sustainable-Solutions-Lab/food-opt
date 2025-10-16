# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

rule resample_land_cover:
    input:
        classes=f"processing/{name}/resource_classes.nc",
        land_cover="data/downloads/land_cover_lccs_class.nc",
    output:
        fractions=f"processing/{name}/luc/land_cover_resampled.nc",
    script:
        "../scripts/resample_land_cover.py"


rule resample_regrowth:
    input:
        classes=f"processing/{name}/resource_classes.nc",
        regrowth_raw="data/downloads/forest_carbon_accumulation_griscom_1km.tif",
    output:
        regrowth=f"processing/{name}/luc/regrowth_resampled.nc",
    script:
        "../scripts/resample_forest_carbon_accumulation.py"


rule prepare_luc_inputs:
    input:
        classes=f"processing/{name}/resource_classes.nc",
        land_cover=rules.resample_land_cover.output.fractions,
        regrowth=rules.resample_regrowth.output.regrowth,
        agb="data/downloads/esa_biomass_cci_v6_0.nc",
        soc="data/downloads/soilgrids_ocs_0-30cm_mean.tif",
    params:
        forest_fraction_threshold=config["luc"]["forest_fraction_threshold"],
    output:
        lc_masks=f"processing/{name}/luc/lc_masks.nc",
        agb=f"processing/{name}/luc/agb.nc",
        soc=f"processing/{name}/luc/soc.nc",
        regrowth=f"processing/{name}/luc/regrowth.nc",
    script:
        "../scripts/prepare_luc_inputs.py"


rule build_luc_carbon_coefficients:
    input:
        classes=f"processing/{name}/resource_classes.nc",
        regions=f"processing/{name}/regions.geojson",
        lc_masks=rules.prepare_luc_inputs.output.lc_masks,
        agb=rules.prepare_luc_inputs.output.agb,
        soc=rules.prepare_luc_inputs.output.soc,
        regrowth=rules.prepare_luc_inputs.output.regrowth,
        zone_parameters="data/luc_zone_parameters.csv",
    params:
        horizon_years=config["luc"]["horizon_years"],
        managed_flux_mode=config["luc"]["managed_flux_mode"],
    output:
        pulses=f"processing/{name}/luc/pulses.nc",
        annualized=f"processing/{name}/luc/annualized.nc",
        coefficients=f"processing/{name}/luc/luc_carbon_coefficients.csv",
    script:
        "../scripts/build_luc_carbon_coefficients.py"
