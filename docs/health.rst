.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Health Impacts
==============

Overview
--------

The health module converts dietary choices in the optimisation into monetised
health impacts. It combines epidemiological evidence on diet–disease links with
country-level baseline mortality and demographic data, and then represents that
relationship inside the linear programme through carefully constructed
piecewise-linear (SOS2) approximations. The objective therefore weighs
production, environmental and health costs in a consistent monetary unit.

Key ideas:

- Dietary risk factors from the Global Burden of Disease (GBD) study underpin
  the exposure–response curves.
- Countries are grouped into health clusters to keep the optimisation tractable
  while preserving heterogeneity in baseline burden and valuation.
- Relative risks multiply across risk factors, so we work in log space to turn
  the problem into additions that can be linearised.

Data Inputs
-----------

``workflow/scripts/prepare_health_costs.py`` assembles the following datasets:

- **Baseline diet** (``data/health/processed/diet_intake.csv``): average daily
  intake by country and food item.
- **Relative risks** (``data/health/processed/relative_risks.csv``): dose–response
  pairs for each (risk factor, disease cause) combination.
- **Mortality rates** (``data/health/processed/mortality.csv``): cause-specific
  death rates by age, country and year.
- **Population and life tables** (``processing/{name}/population_age.csv`` and
  ``processing/{name}/life_table.csv``): age-structured population counts and
  remaining life expectancy schedules.
- **Value of statistical life (optional)** (``data/health/processed/vsl.csv``):
  used when ``health.value_of_statistical_life`` is set to ``"regional"``.

Preparation Workflow
--------------------

The preprocessing script performs these steps:

1. **Health clustering** – dissolves country geometries, computes equal-area
   centroids and runs K-means to assign each country to one of
   ``health.region_clusters`` clusters. The cluster map is saved as
   ``processing/{name}/health/country_clusters.csv``.
2. **Baseline burden** – combines mortality, population and life expectancy to
   compute years of life lost (YLL) per country and aggregates them to the
   health clusters. The results go into
   ``processing/{name}/health/cluster_cause_baseline.csv`` and
   ``processing/{name}/health/cluster_summary.csv``.
3. **Value per YLL** – either reads the regional VSL dataset or applies the
   configured constant, then converts each cluster’s value of a statistical life
   into a value per YLL using average years lost per death.
4. **Risk-factor breakpoints** – builds dense grids of intake values (including
   observed exposures and configured ``health.intake_grid_step``) and evaluates
   :math:`\log(RR)` for every (risk, cause) pair. These tables are written to
   ``processing/{name}/health/risk_breakpoints.csv``.
5. **Cause-level breakpoints** – as the optimisation needs to recover
   :math:`RR = \exp(\sum_r \log RR_{r})`, the script also constructs breakpoints
   for the aggregated log-relative-risk and its exponential. Stored as
   ``processing/{name}/health/cause_log_breakpoints.csv``.

The generated tables drive the linearisation in
``workflow/scripts/solve_model.py``.

From Diet to Risk Exposure
--------------------------

Per-capita intake
~~~~~~~~~~~~~~~~~

During optimisation, consumption flows are tracked on links named
``consume_<food>_<ISO3>``. For each health cluster :math:`c` and risk factor
:math:`r`, the solver forms a per-capita intake by combining these flows with
shares from ``workflow/scripts/health_food_mapping.py``:

.. math::

   I_{c,r} = \frac{10^{6}}{365\,P_c} \sum_{f \in \mathcal{F}_r} \alpha_{f,r} \; q_{c,f}

where

- :math:`q_{c,f}` is the aggregated flow in million tonnes per year for food
  :math:`f` consumed by cluster :math:`c`;
- :math:`\alpha_{f,r}` is the share of food :math:`f` attributed to risk factor
  :math:`r` (currently 1.0 or 0.0);
- :math:`P_c` is the population represented by the cluster (baseline or updated
  planning population);
- the constant rescales from Mt/year to g/day.

Linearised relative risk curves
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each risk factor :math:`r` affects a subset of causes :math:`g`. The data from
``risk_breakpoints.csv`` provides intake breakpoints
:math:`x_0, \ldots, x_K` and the corresponding
:math:`\log RR_{r,g}(x_k)` values. For every (cluster, risk) pair we introduce
SOS2 “lambda” variables :math:`\lambda_k` that satisfy

.. math::
   \sum_k \lambda_k = 1,\qquad I_{c,r} = \sum_k x_k\,\lambda_k,

and approximate the log-relative-risk as

.. math::
   \log RR_{c,r,g} = \sum_k \lambda_k\, \log RR_{r,g}(x_k).

SOS2 constraints keep only two adjacent :math:`\lambda_k` active, yielding a
piecewise-linear interpolation without binary decision variables when the
solver supports SOS2. When HiGHS is used, the implementation falls back to a
compact binary formulation.

Aggregating across risk factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Epidemiological evidence models the combined effect of multiple risk factors on
one cause as multiplicative:

.. math::
   RR_{c,g} = \prod_{r \in \mathcal{R}_g} RR_{c,r,g}.

Taking logarithms converts this to a sum that remains compatible with linear
programming:

.. math::
   \log RR_{c,g} = \sum_{r \in \mathcal{R}_g} \log RR_{c,r,g}.

The solver accumulates the contributions from each risk factor into
``log_rr_totals`` for every cluster–cause pair.

Recovering total relative risk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optimisation needs :math:`RR_{c,g}` again to price health damages. The
preprocessed ``cause_log_breakpoints.csv`` supplies points
:math:`(z_m, \exp(z_m))` that cover the feasible range of
:math:`z = \log RR_{c,g}`. A second SOS2 interpolation enforces

.. math::
   z = \sum_m z_m \theta_m,\qquad RR_{c,g} = \sum_m e^{z_m} \theta_m,

with :math:`\sum_m \theta_m = 1`. This gives a consistent linearised mapping
from the aggregated log-relative-risk back to the multiplicative relative risk.

Monetising years of life lost
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each cluster–cause pair the preprocessing step stores:

- :math:`\mathrm{YLL}^{\mathrm{base}}_{c,g}` – baseline years of life lost, and
- :math:`V_{c}` – value per YLL derived from the value of a statistical life and
  the average years lost per death.

The solver also records the reference log-relative-risk
:math:`z^{\mathrm{ref}}_{c,g}` (from baseline diets) and its exponential
:math:`RR^{\mathrm{ref}}_{c,g}`. The contribution to the objective is
constructed as

.. math::
   \text{Cost}_{c,g} = V_c\, \mathrm{YLL}^{\mathrm{base}}_{c,g}
   \left( \frac{RR_{c,g}}{RR^{\mathrm{ref}}_{c,g}} - 1 \right).

A constant term subtracts
:math:`V_c\,\mathrm{YLL}^{\mathrm{base}}_{c,g}` so that the baseline diet has
zero health cost and only improvements or deteriorations relative to the
reference affect the optimisation.

Objective Contribution
----------------------

``workflow/scripts/solve_model.py`` adds the summed cost over all clusters and
causes to the PyPSA objective. If the solver exposes SOS2 constraints, the
implementation keeps the formulation linear without integer variables; for
HiGHS a tight binary fallback is activated. The script also records the constant
baseline adjustment in ``network.meta["objective_constant_terms"]["health"]`` to
help interpret objective values ex post.

Configuration Highlights
------------------------

.. code-block:: yaml

   health:
     region_clusters: 30               # Number of geographic health clusters
     reference_year: 2018              # Baseline year for diet and mortality data
     intake_grid_step: 10              # g/day spacing for risk breakpoints
     log_rr_points: 10                 # Points for aggregated log-RR interpolation
     value_of_statistical_life: 3_500_000  # USD; set "regional" to use VSL dataset
     risk_factors:
       - fruits
       - vegetables
       - nuts_seeds
       - legumes
       - fish
       - red_meat
       - prc_meat
       - whole_grains

Lowering ``region_clusters`` or ``log_rr_points`` eases the optimisation at the
cost of coarser health resolution. ``health.intake_grid_step`` controls the
density of the first-stage interpolation grid; smaller values give smoother
curves but produce larger tables.

Outputs
-------

The preprocessing rule saves all intermediate products under
``processing/{name}/health/``. Downstream plotting rules also create quick-look
maps (``results/{name}/plots/health_*.pdf``) and CSV summaries to compare
baseline versus optimised health outcomes.

Limitations and Future Work
---------------------------

- **Mortality focus** – only years of life lost are modelled; years lived with
  disability are currently excluded.
- **Static risk mapping** – all foods mapped to a risk factor contribute in
  fixed proportions; nutrient interactions are not captured.
- **Linearisation error** – SOS2 approximations introduce bounded error that
  depends on the chosen grids. Monitor solver logs if experimenting with coarser
  settings.
- **Valuation assumptions** – constant or regional VSL choices can significantly
  shift policy relevance; document your selections when sharing results.

Future extensions may add morbidity effects, age-dependent optimal intakes, or
multi-period health dynamics that capture delayed impacts of dietary change.
