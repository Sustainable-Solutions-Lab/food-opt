# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Rules for generating documentation figures.

These figures are git-tracked and used in the Sphinx documentation.
They use a coarser resolution configuration for faster generation
and to showcase global coverage.
"""

# Documentation figures are generated using the doc_figures config
DOC_FIG_NAME = "doc_figures"

# List of all documentation figures to generate
DOC_FIGURES = [
    # Introduction figures
    "intro_global_coverage",
    # Land use figures
    "land_resource_classes",
    "environment_luc_inputs",
    "environment_luc_lef",
    # Crop production figures
    "crop_yield_wheat",
    "crop_yield_wetland-rice",
    "crop_yield_maize",
    "crop_yield_resource_class_wheat",
    # Water availability figures
    "water_basin_availability",
    "water_region_availability",
    "irrigated_land_fraction",
    # Livestock figures
    "grassland_yield",
    # Trade figures
    "trade_network",
    # Workflow figures
    "workflow_rulegraph",
]


rule doc_fig_intro_global_coverage:
    """Generate global coverage map showing all modeled regions."""
    input:
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
    output:
        svg="docs/_static/figures/intro_global_coverage.svg",
        png="docs/_static/figures/intro_global_coverage.png",
    script:
        "../scripts/doc_figures/intro_global_coverage.py"


rule doc_fig_land_resource_classes:
    """Generate map showing resource class stratification."""
    input:
        classes=f"processing/{DOC_FIG_NAME}/resource_classes.nc",
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
    output:
        svg="docs/_static/figures/land_resource_classes.svg",
        png="docs/_static/figures/land_resource_classes.png",
    script:
        "../scripts/doc_figures/land_resource_classes.py"


rule doc_fig_environment_luc_inputs:
    """Visualise LUC carbon input datasets used in the model."""
    input:
        lc_masks=f"processing/{DOC_FIG_NAME}/luc/lc_masks.nc",
        agb=f"processing/{DOC_FIG_NAME}/luc/agb.nc",
        soc=f"processing/{DOC_FIG_NAME}/luc/soc.nc",
        regrowth="processing/shared/luc/regrowth_resampled.nc",
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
    output:
        svg="docs/_static/figures/environment_luc_inputs.svg",
        png="docs/_static/figures/environment_luc_inputs.png",
    script:
        "../scripts/doc_figures/luc_inputs_map.py"


rule doc_fig_environment_luc_lef:
    """Visualise aggregated land-use change emission factors."""
    input:
        annualized=f"processing/{DOC_FIG_NAME}/luc/annualized.nc",
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
    output:
        svg="docs/_static/figures/environment_luc_lef.svg",
        png="docs/_static/figures/environment_luc_lef.png",
    script:
        "../scripts/doc_figures/luc_lef_map.py"


rule doc_fig_crop_yield:
    """Generate crop yield potential maps for selected crops."""
    input:
        yield_raster=lambda w: gaez_path("yield", "r", w.crop),
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
        conversions="data/yield_unit_conversions.csv",
    output:
        svg="docs/_static/figures/crop_yield_{crop}.svg",
        png="docs/_static/figures/crop_yield_{crop}.png",
    script:
        "../scripts/doc_figures/crop_yield_map.py"


rule doc_fig_crop_yield_resource_class:
    """Generate resource class yield comparison maps."""
    input:
        crop_yields=f"processing/{DOC_FIG_NAME}/crop_yields/{{crop}}_r.csv",
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
    output:
        svg="docs/_static/figures/crop_yield_resource_class_{crop}.svg",
        png="docs/_static/figures/crop_yield_resource_class_{crop}.png",
    params:
        resource_class_1=1,
        resource_class_2=2,
    script:
        "../scripts/doc_figures/crop_yield_resource_class.py"


rule doc_fig_water_basin_availability:
    """Generate basin water availability map."""
    input:
        basin_shapefile="data/downloads/Report53_Appendix/Report53-BlueWaterScarcity-ArcGIS-ShapeFile/Monthly_WS_GRDC_405_basins.shp",
        water_data=f"processing/{DOC_FIG_NAME}/water/blue_water_availability.csv",
    output:
        svg="docs/_static/figures/water_basin_availability.svg",
        png="docs/_static/figures/water_basin_availability.png",
    script:
        "../scripts/doc_figures/water_basin_availability.py"


rule doc_fig_water_region_availability:
    """Generate regional water availability map."""
    input:
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
        water_data=f"processing/{DOC_FIG_NAME}/water/region_growing_season_water.csv",
    output:
        svg="docs/_static/figures/water_region_availability.svg",
        png="docs/_static/figures/water_region_availability.png",
    script:
        "../scripts/doc_figures/water_region_availability.py"


rule doc_fig_irrigated_land_fraction:
    """Generate irrigated land fraction map."""
    input:
        irrigated_fraction="data/downloads/gaez_land_equipped_for_irrigation_share.tif",
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
    output:
        svg="docs/_static/figures/irrigated_land_fraction.svg",
        png="docs/_static/figures/irrigated_land_fraction.png",
    script:
        "../scripts/doc_figures/irrigated_land_fraction.py"


rule doc_fig_grassland_yield:
    """Generate managed grassland yield map."""
    input:
        grassland_yield="data/downloads/grassland_yield_historical.nc4",
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
    output:
        svg="docs/_static/figures/grassland_yield.svg",
        png="docs/_static/figures/grassland_yield.png",
    script:
        "../scripts/doc_figures/grassland_yield_map.py"


rule doc_fig_trade_network:
    """Generate trade network map showing hubs and links."""
    input:
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
    output:
        svg="docs/_static/figures/trade_network.svg",
        png="docs/_static/figures/trade_network.png",
    params:
        n_hubs=config["trade"]["crop_hubs"],
    script:
        "../scripts/doc_figures/trade_network_map.py"


rule doc_fig_workflow_rulegraph_dot:
    """Generate workflow dependency graph in DOT format from Snakemake."""
    output:
        dot="docs/_static/figures/workflow_rulegraph_raw.dot",
    shell:
        """
        snakemake --rulegraph --config name=test > {output.dot}
        """


rule doc_fig_workflow_rulegraph_styled:
    """Apply custom styling and text wrapping to DOT graph."""
    input:
        dot="docs/_static/figures/workflow_rulegraph_raw.dot",
    output:
        dot="docs/_static/figures/workflow_rulegraph.dot",
    script:
        "../scripts/doc_figures/style_workflow_graph.py"


rule doc_fig_workflow_rulegraph:
    """Render workflow dependency graph to SVG and PNG using Graphviz."""
    input:
        dot="docs/_static/figures/workflow_rulegraph.dot",
    output:
        svg="docs/_static/figures/workflow_rulegraph.svg",
        png="docs/_static/figures/workflow_rulegraph.png",
    script:
        "../scripts/doc_figures/render_graph.py"


rule build_docs:
    """Build Sphinx documentation including all figures."""
    input:
        # Figures
        expand("docs/_static/figures/{fig}.svg", fig=DOC_FIGURES),
        expand("docs/_static/figures/{fig}.png", fig=DOC_FIGURES),
        # Documentation source files
        "docs/conf.py",
        "docs/index.rst",
    output:
        "docs/_build/html/index.html",
    shell:
        """
        cd docs && make html
        """
