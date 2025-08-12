<!--
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: CC-BY-4.0
-->

# food-opt

A global food systems optimization model that explores trade-offs between environmental sustainability and nutritional outcomes using linear programming. The model optimizes food production, processing, and consumption patterns while accounting for greenhouse gas emissions, land use change, water use, and dietary health impacts.

## Overview

This model addresses the challenge of feeding a growing global population while minimizing environmental impacts and maximizing nutritional outcomes. It uses a resource flow-based structure implemented with PyPSA/linopy to jointly optimize:

- **Environmental objectives**: Minimize greenhouse gas emissions (CO₂, CH₄, N₂O), land use change, nitrogen pollution, and water use
- **Nutritional objectives**: Meet population dietary requirements while minimizing disease burden from dietary risk factors
- **Production constraints**: Respect biophysical limits on crop yields, land availability, and irrigation capacity

The model operates at flexible sub-national spatial resolution globally, allowing for detailed analysis of regional trade-offs and the spatial distribution of environmental and health impacts.

## Key Features

### Comprehensive Food System Coverage
- **Crop production**: ~170 different crops with spatially-explicit yield potentials
- **Livestock systems**: Multiple production systems (grazing vs. feed-based) for meat and dairy
- **Food processing**: Conversion of raw agricultural products to final food products
- **Nutritional assessment**: Mapping to dietary risk factors and health outcomes

### Environmental Impact Assessment
- Greenhouse gas emissions from production, land use change, and nitrogen fertilization
- Land use change impacts with spatially-explicit carbon storage estimates
- Water use constraints based on irrigation infrastructure
- Nitrogen pollution from fertilizer application

### Health and Nutrition
- Integration with Global Burden of Disease dietary risk factors
- Macronutrient and micronutrient constraints
- Population-level health impact assessment in DALYs (Disability-Adjusted Life Years)

### Flexible Spatial Resolution
- Input data at high-resolution gridcell level (0.05° × 0.05°)
- Optimization at configurable sub-national regional scale
- Global coverage with detailed country and regional analysis

## Project Structure

```
food-model/
├── config/           # Configuration files for scenarios and parameters
├── data/             # Input data (downloaded and processed)
├── results/          # Model outputs and analysis
│   ├── foo/          # Outputs from the config with name 'foo'
│   └── bar/          # Outputs from the config with name 'bar'
└── workflow/         # Snakemake workflow for data processing and model execution
    ├── envs/         # Conda environment specifications
    ├── scripts/      # Data processing and model building scripts
    └── Snakefile     # Main workflow definition
```

## Data Sources

The model integrates multiple global datasets:

- **GAEZ (Global Agro-Ecological Zones)**: Crop suitability and yield potentials
- **FAOSTAT**: Agricultural production, land use, and fertilizer application data
- **CROPGRIDS**: High-resolution crop distribution maps
- **AQUASTAT**: Irrigation infrastructure and water use
- **Global Burden of Disease**: Dietary risk factors and health impacts
- **Global Dietary Database**: Population dietary intake estimates

## Getting Started

### Prerequisites

- Python
- [Snakemake](https://snakemake.readthedocs.io/) workflow management system
- [Conda](https://github.com/conda-forge/miniforge/) for environment management

### Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd food-model
```

2. Install Snakemake (if not already installed):
```bash
conda install -c conda-forge snakemake
```

### Running the Model

1. **Configure the model**: Edit `config/config.yaml` to specify scenarios, spatial resolution, and optimization parameters.

2. **Process data, build and optimize model**:
```bash
snakemake --use-conda -j4 all
```

The workflow will automatically:
- Download required datasets
- Process and harmonize spatial data
- Build the linear programming model
- Solve the optimization problem
- Generate summary statistics and visualizations

Individual stages of the workflow can also be targeted by specifying the name of rule outputs; have a look a `workflow/Snakefile` for an overview of the rules and outputs of the workflow. For example, run `snakemake --use-conda -j4 results/model.nc`

## Model Structure

### Decision Variables
- Crop production by region and crop type
- Livestock production by region and production system
- Food processing and conversion activities
- Land allocation between crops and between agricultural and non-agricultural uses
- Fertilizer application rates

### Constraints
- **Land availability**: Total agricultural land limits by region
- **Yield relationships**: Crop production as function of land and inputs
- **Nutritional requirements**: Minimum intake levels for macro/micronutrients
- **Processing relationships**: Mass balance in food processing chains
- **Environmental limits**: Optional caps on emissions or resource use

### Objective Function
Configurable multi-objective optimization balancing:
- Environmental costs (emissions, land use change, pollution)
- Health costs (dietary risk factors, nutritional deficiencies)
- Economic costs (production costs, opportunity costs)

## Configuration

The `config/config.yaml` file controls key model parameters:

- **Spatial resolution**: Regional aggregation level
- **Scenario definitions**: Population projections, dietary preferences, technology assumptions
- **Environmental constraints**: Emission limits, land use restrictions
- **Objective weights**: Trade-offs between environmental and health objectives

## Results and Analysis

Model outputs include:

- **Optimal food production patterns**: Crop and livestock production by region
- **Land use allocation**: Agricultural land use and land use change
- **Environmental impacts**: Emissions, water use, and pollution by source
- **Nutritional outcomes**: Population diet composition and health impacts
- **Trade flows**: Inter-regional food trade patterns

Results are provided in both tabular format (CSV) and interactive visualizations showing spatial patterns and trade-offs.

## License

All code is licensed under GPL-3.0.

## Contact

[Contact information to be added]

---

*This model is under active development. Please check the repository for the latest updates and documentation.*
