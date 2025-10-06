<!--
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: CC-BY-4.0
-->

# food-opt Documentation

This directory contains the Sphinx documentation for the food-opt global food systems optimization model.

## Building Documentation Locally

### Prerequisites

Ensure documentation dependencies are installed:

```bash
cd ..  # Return to project root
uv sync --dev
```

### Build HTML Documentation

From this directory:

```bash
make html
```

Or directly with sphinx-build:

```bash
uv run sphinx-build -b html . _build/html
```

The HTML documentation will be in `_build/html/`. Open `_build/html/index.html` in your web browser.

### Other Build Formats

- **PDF**: `make latexpdf` (requires LaTeX)
- **EPUB**: `make epub`
- **Clean build**: `make clean` then `make html`

## Documentation Structure

- `index.rst` - Main table of contents
- `introduction.rst` - Getting started guide
- `model_framework.rst` - Mathematical formulation and PyPSA structure
- `land_use.rst` - Spatial aggregation and resource classes
- `crop_production.rst` - Yield modeling and GAEZ data
- `livestock.rst` - Animal production systems
- `food_processing.rst` - Processing chains and trade networks
- `nutrition.rst` - Dietary constraints
- `health.rst` - Health impact assessment (DALYs)
- `environment.rst` - Emissions, water, and nitrogen
- `configuration.rst` - Configuration file reference
- `data_sources.rst` - Dataset documentation
- `workflow.rst` - Snakemake workflow execution
- `results.rst` - Analyzing and visualizing outputs
- `development.rst` - Contributing guidelines
- `troubleshooting.rst` - Common issues and FAQ
- `api/index.rst` - Auto-generated API reference

## Editing Documentation

1. Edit `.rst` files in this directory
2. Rebuild: `make html`
3. Check for warnings in build output
4. Preview in browser: `open _build/html/index.html` (macOS) or `xdg-open _build/html/index.html` (Linux)
5. Commit changes to Git

## Publishing to ReadTheDocs

(Configuration to be added when ReadTheDocs hosting is set up)

1. Create `.readthedocs.yaml` in project root
2. Connect GitHub repository to ReadTheDocs
3. Builds will automatically trigger on push

## Sphinx Configuration

Documentation settings are in `conf.py`:

- Theme: sphinx_rtd_theme (Read the Docs theme)
- Extensions: autodoc, napoleon, viewcode, intersphinx, mathjax
- Intersphinx links to: Python, NumPy, pandas, xarray, PyPSA docs

## License

Documentation content is licensed under CC-BY-4.0 (see SPDX headers in `.rst` files).
