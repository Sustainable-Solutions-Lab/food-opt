.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Introduction
============

Overview
--------

**food-opt** is a global food systems optimization model that addresses the challenge of feeding a growing global population while minimizing environmental impacts and maximizing nutritional outcomes. The model uses a resource flow-based structure implemented with PyPSA/linopy to jointly optimize food production, processing, and consumption patterns.

Key Objectives
~~~~~~~~~~~~~~

The model balances multiple objectives:

* **Environmental sustainability**: Minimize greenhouse gas emissions (CO₂, CH₄, N₂O), land use change, nitrogen pollution, and water use
* **Nutritional adequacy**: Meet population dietary requirements for macronutrients and micronutrients
* **Health outcomes**: Minimize disease burden from dietary risk factors
* **Production constraints**: Respect biophysical limits on crop yields, land availability, and irrigation capacity

Key Features
------------

Comprehensive Food System Coverage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Crop production**: ~70 different crops with spatially-explicit yield potentials
* **Livestock systems**: Multiple production systems (grazing vs. feed-based) for meat and dairy
* **Food processing**: Conversion of raw agricultural products to final food products
* **Nutritional assessment**: Mapping to dietary risk factors and health outcomes

Environmental Impact Assessment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Greenhouse gas emissions from production, land use change, and nitrogen fertilization
* Land use change impacts with spatially-explicit carbon storage estimates
* Water use constraints based on irrigation infrastructure and basin-level availability
* Nitrogen pollution from fertilizer application

Health and Nutrition
~~~~~~~~~~~~~~~~~~~~~

* Integration with Global Burden of Disease dietary risk factors
* Macronutrient and micronutrient constraints
* Population-level health impact assessment in DALYs (Disability-Adjusted Life Years)
* Value of statistical life calculations for health cost valuation

Flexible Spatial Resolution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Input data at high-resolution gridcell level (0.05° × 0.05°)
* Optimization at configurable sub-national regional scale
* Global coverage with detailed country and regional analysis

Getting Started
---------------

Prerequisites
~~~~~~~~~~~~~

* Python >= 3.12
* `uv <https://docs.astral.sh/uv/>`_ for dependency management
* `Snakemake <https://snakemake.readthedocs.io/>`_ workflow management system
* Linear programming solver (HiGHS included, Gurobi optional)

Installation
~~~~~~~~~~~~

1. Clone the repository::

    git clone <repository-url>
    cd food-opt

2. Install dependencies::

    uv sync

3. The workflow will automatically download required datasets when first run.

Quick Start
-----------

Running Your First Model
~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to get started is to run the included toy configuration::

    tools/smk -j4 all

This command will:

1. Download required global datasets (GAEZ, GADM, UN population, etc.)
2. Process and harmonize spatial data for the configured countries
3. Build the linear programming model
4. Solve the optimization problem
5. Generate summary statistics and visualizations

Results will be saved under ``results/toy/``.

Understanding the Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Snakemake workflow is organized into stages:

* **Data preparation**: Population, regions, resource classes, crop yields
* **Model building**: Assemble PyPSA network with all constraints
* **Solving**: Run the linear program with configured solver
* **Visualization**: Generate maps, plots, and CSV exports

You can target individual stages by specifying the output file. For example, to only build the model without solving::

    tools/smk -j4 results/toy/build/model.nc

Or to just prepare regional aggregation::

    tools/smk -j4 processing/toy/regions.geojson

See :doc:`workflow` for detailed information on the workflow stages.

Configuring Your First Scenario
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The toy configuration (``config/config.yaml``) provides a starting point. Key parameters to adjust:

* ``countries``: List of ISO 3166-1 alpha-3 country codes to include
* ``aggregation.regions.target_count``: Number of optimization regions (trade-off between detail and solve time)
* ``crops``: Which crops to include in the model
* ``emissions.ghg_price``: Carbon price in USD/tCO2-eq
* ``macronutrients``: Minimum dietary requirements

After editing the configuration, create a new named scenario by changing the ``name`` field at the top of the file, then run::

    tools/smk -j4 all

Results will be saved under ``results/<your-name>/``.

Project Structure
-----------------

The repository is organized as follows::

    food-opt/
    ├── config/              # Configuration files for scenarios and parameters
    │   └── config.yaml      # Main configuration file
    ├── data/                # Input data (downloaded and processed)
    │   ├── downloads/       # Raw downloaded datasets
    │   ├── crops.csv        # Crop definitions
    │   ├── foods.csv        # Food product definitions (mock data)
    │   └── nutrition.csv    # Nutritional content (mock data)
    ├── processing/          # Intermediate processed datasets
    │   └── {config_name}/   # Processing outputs per scenario
    ├── results/             # Model outputs and analysis
    │   └── {config_name}/   # Results per scenario
    │       ├── build/       # Built model before solving
    │       ├── solved/      # Solved model with optimal values
    │       └── plots/       # Visualizations and CSV exports
    ├── workflow/            # Snakemake workflow
    │   ├── Snakefile        # Main workflow definition
    │   ├── rules/           # Modular rule definitions
    │   └── scripts/         # Data processing and modeling scripts
    ├── tools/               # Utility wrappers
    │   └── smk              # Memory-capped Snakemake wrapper
    ├── notebooks/           # Exploratory analyses
    └── vendor/              # Bundled third-party dependencies

Important Notes
~~~~~~~~~~~~~~~

* The ``results/`` directory contains auto-generated files—never edit these manually
* Several CSV files (``data/foods.csv``, ``data/nutrition.csv``, ``data/feed_conversion.csv``) contain mock placeholder data
* Always use the ``tools/smk`` wrapper to run Snakemake, as it enforces memory limits to prevent system instability
* The first run will take significant time to download global datasets (~several GB)

