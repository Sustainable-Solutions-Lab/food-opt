.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Development & Contributing
===========================

Overview
--------

This page provides guidance for developers contributing to the food-opt project, including code conventions, testing procedures, and best practices.

For AI coding agents, see ``CLAUDE.md`` in the repository root for specific instructions.

Development Setup
-----------------

Prerequisites
~~~~~~~~~~~~~

* Python >= 3.12
* Git
* uv (dependency manager)
* pre-commit (for code quality hooks)

Installation
~~~~~~~~~~~~

1. Clone the repository::

       git clone <repository-url>
       cd food-opt

2. Install dependencies::

       uv sync

3. Install development tools::

       uv sync --dev

4. Set up pre-commit hooks::

       uv run pre-commit install

Code Conventions
----------------

Style Guidelines
~~~~~~~~~~~~~~~~

The project uses **ruff** for linting and formatting, enforcing:

* PEP 8 style (with 88-character line length)
* Import sorting (isort)
* Type hints (where practical)
* Docstrings for public functions

**Run linter**::

    uv run ruff check .

**Auto-format code**::

    uv run ruff format .

**Run from pre-commit** (automatic on ``git commit``)::

    uv run pre-commit run --all-files

Specific Conventions
~~~~~~~~~~~~~~~~~~~~

* **No unused imports**: Ruff removes them automatically
* **No ``from __future__ import annotations``**: Unnecessary with modern Python
* **Fail early**: Validate external inputs; trust internal invariants
* **Concise logic**: Prefer simple control flow; avoid over-engineering
* **Docstrings**: Use NumPy/Google style for functions with non-obvious behavior

Example:

.. code-block:: python

   def aggregate_yields(
       yields: np.ndarray,
       classes: np.ndarray,
       regions: gpd.GeoDataFrame,
   ) -> pd.DataFrame:
       """Aggregate gridded yields to (region, class) combinations.

       Parameters
       ----------
       yields : np.ndarray
           Yield raster (t/ha)
       classes : np.ndarray
           Resource class assignment raster
       regions : gpd.GeoDataFrame
           Region polygons

       Returns
       -------
       pd.DataFrame
           Columns: region, class, mean_yield
       """
       # Implementation...

Licensing
~~~~~~~~~

* **Code**: GPL-3.0-or-later (use SPDX header in ``.py`` files)
* **Documentation**: CC-BY-4.0 (use SPDX header in ``.rst``, ``.md`` files)

SPDX headers (required in all source files):

.. code-block:: python

   # SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
   #
   # SPDX-License-Identifier: GPL-3.0-or-later

Repository Structure
--------------------

::

    food-opt/
    ├── config/              # Scenario configuration files
    ├── data/                # Input data (not committed)
    ├── docs/                # Documentation (Sphinx)
    ├── processing/          # Intermediate outputs (not committed)
    ├── results/             # Model results (not committed)
    ├── workflow/            # Snakemake workflow
    │   ├── Snakefile        # Main workflow definition
    │   ├── rules/           # Modular rule files
    │   └── scripts/         # Python scripts for processing/modeling
    ├── tools/               # Utility wrappers (e.g., smk)
    ├── notebooks/           # Exploratory Jupyter notebooks
    ├── vendor/              # Bundled third-party code (customized PyPSA/linopy)
    ├── .gitignore
    ├── pyproject.toml       # Dependencies and tool config
    ├── README.md
    └── CLAUDE.md            # AI agent guidance

Adding New Features
-------------------

Adding a New Crop
~~~~~~~~~~~~~~~~~

1. **Check GAEZ availability**: Ensure crop is in GAEZ v5 (see ``data/gaez_crop_code_mapping.csv``)

2. **Add to config**:

   .. code-block:: yaml

      crops:
        - wheat
        - maize
        - my_new_crop  # Add here

3. **Add GAEZ mapping** (if needed):

   Edit ``data/gaez_crop_code_mapping.csv``:

   .. code-block:: text

      crop_name,res02_code,res05_code,res06_code
      my_new_crop,CROP_ID,CROP_ID,CROP_ID

4. **Run workflow**::

       tools/smk -j4 results/toy/solved/model.nc

   Snakemake will automatically download new GAEZ files and incorporate the crop.

5. **Add metadata**: Update ``data/crops.csv`` with fertilizer requirements, emission factors, etc.

Adding a New Constraint
~~~~~~~~~~~~~~~~~~~~~~~

1. **Identify constraint type**:

   * Production limit (land, water, fertilizer)
   * Nutritional requirement (macronutrient, food group)
   * Environmental cap (emissions, nitrogen)
   * Policy constraint (minimum animal product, organic share)

2. **Implement in** ``workflow/scripts/build_model.py`` or ``solve_model.py``:

   .. code-block:: python

      # Example: Add minimum legume production constraint
      legume_crops = ["soybean", "chickpea", "lentil"]
      for country in countries:
          legume_production = sum(
              n.links[f"production_crop_{crop}_{country}"]["p0"]
              for crop in legume_crops
          )
          n.add(
              "GlobalConstraint",
              f"min_legume_{country}",
              type=">=",
              carrier_attribute="legume_production",
              constant=min_legume_t,  # From config
          )

3. **Add config parameter**:

   .. code-block:: yaml

      constraints:
        min_legume_production: 1e6  # tonnes globally

4. **Test**: Run with new constraint, verify feasibility

Adding a New Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Create script** ``workflow/scripts/plot_my_metric.py``:

   .. code-block:: python

      import pypsa
      import matplotlib.pyplot as plt

      n = pypsa.Network(snakemake.input.network)

      # Extract and process data
      metric_data = extract_my_metric(n)

      # Plot
      fig, ax = plt.subplots()
      metric_data.plot(kind="bar", ax=ax)
      ax.set_ylabel("My Metric")
      ax.set_title("My Analysis")

      plt.savefig(snakemake.output.plot, bbox_inches="tight")

2. **Add rule** in ``workflow/rules/plotting.smk``:

   .. code-block:: python

      rule plot_my_metric:
          input:
              network="results/{name}/solved/model.nc"
          output:
              plot="results/{name}/plots/my_metric.pdf"
          script:
              "../scripts/plot_my_metric.py"

3. **Add to** ``all`` **rule** (optional):

   .. code-block:: python

      rule all:
          input:
              # ...
              f"results/{name}/plots/my_metric.pdf"

4. **Run**::

       tools/smk results/toy/plots/my_metric.pdf

Testing
-------

Validation with Toy Config
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``toy`` configuration (400 regions, relaxed constraints) serves as an integration test.

**Run full workflow**::

    tools/smk -j4 all

**Expected outcome**: Completes without errors, produces all plots

**Typical runtime**: 30-90 minutes (depending on hardware, whether data is cached)

Unit Testing (Future Work)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently, the project lacks formal unit tests. Future additions should cover:

* **Utility functions**: raster processing, unit conversions
* **Aggregation logic**: Resource class computation, yield averaging
* **Constraint construction**: PyPSA component creation

Use ``pytest`` for unit tests (add to ``pyproject.toml`` dev dependencies).

Workflow Testing
~~~~~~~~~~~~~~~~

Test individual stages::

    # Test region building
    tools/smk processing/toy/regions.geojson

    # Test yield processing for one crop
    tools/smk processing/toy/crop_yields/wheat_r.csv

    # Test model building (no solving)
    tools/smk results/toy/build/model.nc

This allows catching errors early without waiting for full workflow.

Version Control
---------------

Git Workflow
~~~~~~~~~~~~

1. **Branch for features**::

       git checkout -b feature/my-new-feature

2. **Commit frequently** with descriptive messages::

       git commit -m "Add minimum legume production constraint"

3. **Push to remote**::

       git push origin feature/my-new-feature

4. **Create pull request** for review

Commit Messages
~~~~~~~~~~~~~~~

Follow conventional commit style:

* ``feat: Add new crop to GAEZ mapping``
* ``fix: Correct water requirement unit conversion``
* ``docs: Update health module documentation``
* ``refactor: Simplify resource class computation``
* ``test: Add validation for toy config``

What to Commit
~~~~~~~~~~~~~~

**DO commit**:

* Code (``.py``, ``.smk``)
* Configuration (``.yaml``)
* Documentation (``.rst``, ``.md``)
* Static data files (``data/*.csv`` if < 1 MB)

**DO NOT commit**:

* Downloaded datasets (``data/downloads/``)
* Processed intermediate files (``processing/``)
* Results (``results/``)
* Large binary files (> 1 MB)

These are excluded via ``.gitignore``.

Documentation
-------------

Building Documentation Locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    cd docs
    uv run sphinx-build -b html . _build/html
    # Open _build/html/index.html in browser

Or use the Makefile (if created)::

    cd docs
    make html

Updating Documentation
~~~~~~~~~~~~~~~~~~~~~~

1. **Edit** ``.rst`` files in ``docs/``
2. **Rebuild**::

       cd docs && uv run sphinx-build -b html . _build/html

3. **Check** for warnings/errors
4. **Commit** documentation changes

Docstring Guidelines
~~~~~~~~~~~~~~~~~~~~

Use NumPy-style docstrings:

.. code-block:: python

   def my_function(param1: int, param2: str) -> float:
       """One-line summary.

       Longer description if needed, explaining purpose, algorithm, etc.

       Parameters
       ----------
       param1 : int
           Description of param1
       param2 : str
           Description of param2

       Returns
       -------
       float
           Description of return value

       Raises
       ------
       ValueError
           If param1 is negative

       Notes
       -----
       Additional implementation notes, references, etc.
       """

Performance Optimization
------------------------

Profile Before Optimizing
~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``cProfile`` or ``line_profiler`` to identify bottlenecks:

.. code-block:: python

   import cProfile
   cProfile.run("my_function()")

Memory Profiling
~~~~~~~~~~~~~~~~

For memory-intensive operations::

    uv run python -m memory_profiler workflow/scripts/my_script.py

Optimization Tips
~~~~~~~~~~~~~~~~~

* **Vectorize with NumPy**: Avoid Python loops over large arrays
* **Use exactextract for raster aggregation**: Much faster than naive pixel iteration
* **Lazy loading**: Use ``xarray.open_dataset()`` with ``chunks`` for large files
* **Parallelize Snakemake rules**: Independent rules run concurrently with ``-j``

Contributing Guidelines
-----------------------

Before Submitting a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Run linter**: ``uv run ruff check . && uv run ruff format .``
2. **Test workflow**: Verify toy config runs successfully
3. **Update documentation**: If changing user-facing behavior
4. **Add SPDX headers**: To any new files
5. **Write commit messages**: Descriptive and following conventions

Pull Request Process
~~~~~~~~~~~~~~~~~~~~~

1. Fork the repository
2. Create a feature branch
3. Make changes with clear commits
4. Push to your fork
5. Open pull request with description of changes
6. Address review feedback
7. Merge once approved

Community Guidelines
--------------------

* **Be respectful**: Constructive feedback, no harassment
* **Ask questions**: If unsure about approach, open an issue for discussion
* **Document changes**: Help others understand your contributions
* **Iterate**: Expect revisions, embrace code review

