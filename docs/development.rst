.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Development & Contributing
===========================

Overview
--------

This page provides guidance for developers contributing to the food-opt project, including code conventions and best practices.

For AI coding agents, see ``AGENTS.md`` in the repository root for specific instructions.

Development Setup
-----------------

Prerequisites
~~~~~~~~~~~~~

* Python >= 3.12
* Git
* uv (dependency manager)

Installation
~~~~~~~~~~~~

1. Clone the repository::

       git clone <repository-url>
       cd food-opt

2. Install dependencies::

       uv sync

3. Install development tools::

       uv sync --extra dev

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

* **Fail early**: Validate external inputs; trust internal invariants
* **Concise logic**: Prefer simple control flow; avoid over-engineering
* **Docstrings**: Use NumPy style for functions with non-obvious behavior

Licensing
~~~~~~~~~

* **Code**: GPL-3.0-or-later (use SPDX header in ``.py`` files)
* **Documentation**: CC-BY-4.0 (use SPDX header in ``.rst``, ``.md`` files)

SPDX headers (required in all source files):

.. code-block:: python

   # SPDX-FileCopyrightText: 2025 <Author>
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
    └── AGENTS.md            # AI agent guidance

Adding New Features
-------------------

Adding a New Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Create script** ``workflow/scripts/plotting/plot_my_metric.py``:

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
              "../scripts/plotting/plot_my_metric.py"

3. **Add to** ``all`` **rule** (optional):

   .. code-block:: python

      rule all:
          input:
              # ...
              f"results/{name}/plots/my_metric.pdf"

4. **Run**::

       tools/smk --configfile config/my_scenario.yaml results/my_scenario/plots/my_metric.pdf

Version Control
---------------

Git Workflow
~~~~~~~~~~~~

1. **Branch for features**::

       git checkout -b feature/my-new-feature

2. **Commit frequently** with descriptive messages::

       git commit -m "feat: Add minimum legume production constraint"

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
* ``test: Add validation for quickstart config``

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

    tools/smk -j4 --configfile config/doc_figures.yaml -- build_docs
    # Open _build/html/index.html in browser


Updating Documentation
~~~~~~~~~~~~~~~~~~~~~~

1. **Edit** ``.rst`` files in ``docs/``
2. **Rebuild**::

       tools/smk -j4 --configfile config/doc_figures.yaml --forcerun build_docs -- build_docs

   The ``--forcerun build_docs`` is necessary when only documentation text has been updated and no other files/figures.

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

Contributing Guidelines
-----------------------

Before Submitting a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Run linter**: ``uv run ruff check . && uv run ruff format .`` (this is taken care of automatically if you set up ``pre-commit``)
2. **Test workflow**: Verify that the default configuration runs successfully
3. **Update documentation**: If changing user-facing behavior
4. **Write commit messages**: Descriptive and following conventions

Pull Request Process
~~~~~~~~~~~~~~~~~~~~~

1. Fork the repository
2. Create a feature branch
3. Make changes with clear commits
4. Push to your fork
5. Open pull request with description of changes
6. Address review feedback
7. Merge once approved
