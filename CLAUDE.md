<!--
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: CC-BY-4.0
-->

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is food-opt, a global food systems optimization model using linear programming to explore trade-offs between environmental sustainability and nutritional outcomes. The model is built using PyPSA (Python for Power System Analysis) as the optimization framework and Snakemake for workflow management.

## Architecture

### Core Components
- **PyPSA Network Model**: The optimization model is built using PyPSA's network structure with carriers (resources), buses, stores, and links
- **Snakemake Workflow**: Data processing and model execution pipeline in `workflow/Snakefile`
- **Multi-Stage Process**: Data preparation (multiple scripts), model building (`build_model.py`) followed by solving (`solve_model.py`) and finally plotting.

### Key Files
- `workflow/scripts/build_model.py`: Creates the PyPSA network with carriers (land, water, fertilizer, crops, nutrients), buses, and stores
- `workflow/scripts/solve_model.py`: Solves the optimization problem
- `config/config.yaml`: Model configuration including resource limits and crop selection
- `data/`: Input datasets including foods.csv with crop-to-food conversion factors

### Data Flow
1. Configuration defines scenarios, crops, and resource limits
2. Data preparation scripts download and process input data
3. Build script creates network structure from input data
4. Solve script optimizes the network and outputs results
5. Results stored in `results/{config_name}/` directory structure

## Development Commands

### Environment Setup
Prepend the following to any bash command needing snakemake or devel dependencies: `source ~/.bashrc && conda activate food-opt &&`

### Running the Model
In general, it's not necessary to clean up the results directory after modifying code; snakemake will detect when code was changed and rerun the necessary rules.
```bash
# Full workflow (data processing, build, solve)
snakemake --use-conda -j4 all

# Build model only
snakemake --use-conda -j4 results/{config_name}/build/model.nc

# Solve model only
snakemake --use-conda -j4 results/{config_name}/solved/model.nc
```

### Code Quality
Linting and formatting is done using ruff, but this is applied automatically after every code change through a claude code hook and doesn't have to be done manually.

## Configuration System

The model uses a hierarchical configuration system:
- `config/config.yaml` defines scenarios with name, resource limits, and crop selections
- Results organized by config name in `results/{name}/` directories
- Snakemake accesses config via `config["name"]` and nested parameters

## Key Dependencies

- **PyPSA**: Core optimization framework (≥0.35).
- **Gurobi/HiGHS**: Linear programming solvers (gurobipy, highspy)
- **Scientific Stack**: pandas (≥2.1), numpy, matplotlib
- **Workflow**: Snakemake (≥9.0) with conda integration

### PyPSA

For code involving pypsa components, optimization, consult the official documentation using the context7 tool; library ID is pypsa/pypsa.
If you need to look up linopy documentation, search the web.

Links are tricky components. When they have more than 2 buses, the following conventions reign:
- `bus0` is the (first) input bus
- `bus1` is the (first) output bus; `efficiency` controls the efficiency from bus0 to bus1
- `bus2`, etc. are additional inputs or outputs; the corresponding parameter `efficiency2`, etc. controls whether it's an input or output. Positive for output, negative for input. Efficiencies are relative to bus0.
