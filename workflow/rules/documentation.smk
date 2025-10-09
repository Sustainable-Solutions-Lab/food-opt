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
    # Crop production figures
    "crop_yield_wheat",
    "crop_yield_wetland-rice",
    "crop_yield_maize",
]


rule doc_fig_intro_global_coverage:
    """Generate global coverage map showing all modeled regions."""
    input:
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
    output:
        svg="docs/_static/figures/intro_global_coverage.svg",
    script:
        "../scripts/doc_figures/intro_global_coverage.py"


rule doc_fig_land_resource_classes:
    """Generate map showing resource class stratification."""
    input:
        classes=f"processing/{DOC_FIG_NAME}/resource_classes.nc",
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
    output:
        svg="docs/_static/figures/land_resource_classes.svg",
    script:
        "../scripts/doc_figures/land_resource_classes.py"


rule doc_fig_crop_yield:
    """Generate crop yield potential maps for selected crops."""
    input:
        yield_raster=lambda w: gaez_path("yield", "r", w.crop),
        regions=f"processing/{DOC_FIG_NAME}/regions.geojson",
        conversions="data/yield_unit_conversions.csv",
    output:
        svg="docs/_static/figures/crop_yield_{crop}.svg",
    script:
        "../scripts/doc_figures/crop_yield_map.py"


rule generate_all_doc_figures:
    """Generate all documentation figures."""
    input:
        expand("docs/_static/figures/{fig}.svg", fig=DOC_FIGURES),
    output:
        touch("docs/_static/figures/.doc_figures_complete"),


rule build_docs:
    """Build Sphinx documentation including all figures."""
    input:
        # Ensure all figures are generated first
        "docs/_static/figures/.doc_figures_complete",
        # Documentation source files
        "docs/conf.py",
        "docs/index.rst",
    output:
        "docs/_build/html/index.html",
    shell:
        """
        cd docs && make html
        """
