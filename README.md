<!--
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: CC-BY-4.0
-->

<h1>
  <img src="docs/_static/logo.svg" alt="food-opt logo" height="40" style="vertical-align: middle;"> food-opt
</h1>

[![Docs](https://github.com/Sustainable-Solutions-Lab/food-opt/actions/workflows/docs.yml/badge.svg)](https://github.com/Sustainable-Solutions-Lab/food-opt/actions/workflows/docs.yml)

food-opt is a global food-systems optimization model built on [PyPSA](https://pypsa.org/) and [Snakemake](https://snakemake.readthedocs.io). It explores environmental, nutritional, and economic trade-offs through a configuration-driven mixed integer linear program build by reproducible workflow.

## Documentation

Documentation (model design, configuration reference, data provenance, API) lives under `docs/`; the documentation is built automatically by a Github Action and, for now, can be accessed by clicking the documentation badge at the top of this page.

## Quickstart

To install, clone the git repository and make sure you have [uv](https://docs.astral.sh/uv/) installed.

In order to build the model, you will first have to download a few datasets manually. See the documentation for details. Then you are ready to run the model:

```bash
uv sync
tools/smk -j4 --config name=test
```

- `tools/smk` wraps Snakemake with the repository’s resource limits and environment pins.
- Snakemake targets land under `results/{name}/`.

## Repository Layout

- `workflow/` – Snakemake rules and scripts, including the top-level `workflow/Snakefile`.
- `config/` – Scenario YAMLs and shared fragments that parameterize the workflow.
- `docs/` – Sphinx documentation sources (see `docs/README.md` for dev tips).
- `tools/` – Helper wrappers such as `tools/smk` for consistent CLI entry points.
- `results/` – Generated artifacts grouped by configuration (never hand-edit).

Additional contribution guidance can be found in the documentation; dataset provenance is tracked in `data/DATASETS.md`.

## License

food-opt is licensed under GPL-3.0-or-later; documentation content follows CC-BY-4.0. See `LICENSES/` for details.
