.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _health-impacts:
   
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

Dietary Risk Factors
--------------------

The model incorporates dietary risk factors as defined by the Global Burden of Disease (GBD) Study 2021 [Brauer2024]_. These risk factors link dietary intake patterns to specific disease outcomes through dose-response relationships.

GBD Risk Factor Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following table reproduces the GBD 2021 dietary risk factor definitions from Brauer et al. (2024, Supplementary Appendix 1, p. 171). All intake quantities are expressed in terms of **fresh (as consumed) weight** unless otherwise specified. The optimal intake levels represent the theoretical minimum risk exposure level (TMREL) used in GBD burden calculations:

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Risk Factor
     - Definition of Exposure
     - Optimal Level or Range
   * - Diet low in fruit
     - Average daily consumption (in grams per day) of fruit including fresh, frozen, cooked, canned, or dried fruit, excluding fruit juices and salted or pickled fruits
     - 340–350 g/day
   * - Diet low in vegetables
     - Average daily consumption (in grams per day) of vegetables, including fresh, frozen, cooked, canned, or dried vegetables and excluding legumes and salted or pickled vegetables, juices, nuts and seeds, and starchy vegetables such as potatoes or corn
     - 306–372 g/day
   * - Diet low in whole grains
     - Average daily consumption (in grams per day) of whole grains (bran, germ, and endosperm in their natural proportion) from breakfast cereals, bread, rice, pasta, biscuits, muffins, tortillas, pancakes, and other sources
     - 160–210 g/day
   * - Diet low in nuts and seeds
     - Average daily consumption (in grams per day) of nuts and seeds, including tree nuts and seeds and peanuts
     - 19–24 g/day
   * - Diet low in fibre
     - Average daily consumption (in grams per day) of fibre from all sources including fruits, vegetables, grains, legumes, and pulses
     - 22–25 g/day
   * - Diet low in seafood omega-3 fatty acids
     - Average daily consumption (in milligrams per day) of eicosapentaenoic acid (EPA) and docosahexaenoic acid (DHA)
     - 470–660 mg/day
   * - Diet low in omega-6 polyunsaturated fatty acids
     - Average daily consumption (in % daily energy) from omega-6 polyunsaturated fatty acids (PUFA) (specifically linoleic acid, γ-linolenic acid, eicosadienoic acid, dihomo-γ-linolenic acid, arachidonic acid)
     - 9–10% of total daily energy
   * - Diet low in calcium
     - Average daily consumption (in grams per day) of calcium from all sources, including milk, yogurt, and cheese
     - 0.72–0.86 g/day (males), 1.1–1.2 g/day (females)
   * - Diet low in milk
     - Average daily consumption (in grams per day) of dairy milk including non-fat, low-fat, and full-fat milk, but excluding plant-based milks, fermented milk products such as buttermilk, and other dairy products such as cheese
     - 280–340 g/day (males), 500–610 g/day (females)
   * - Diet low in legumes
     - Average daily consumption (in grams per day) of legumes and pulses, including fresh, frozen, cooked, canned, or dried legumes
     - 100–110 g/day
   * - Diet high in red meat
     - Average daily consumption (in grams per day) of unprocessed red meat including pork and bovine meats such as beef, pork, lamb, and goat, but excluding all processed meats, poultry, fish, and eggs
     - 0–200 g/day
   * - Diet high in processed meat
     - Average daily consumption (in grams per day) of meat preserved by smoking, curing, salting, or addition of chemical preservatives
     - 0 g/day
   * - Diet high in sugar-sweetened beverages (SSBs)
     - Average daily consumption (in grams per day) of beverages with ≥50 kcal per 226.8 gram serving, including carbonated beverages, sodas, energy drinks, and fruit drinks, but excluding 100% fruit and vegetable juices
     - 0 g/day
   * - Diet high in trans fatty acids
     - Average daily consumption (in percent daily energy) of trans fat from all sources, mainly partially hydrogenated vegetable oils and ruminant products
     - 0–1.1% of total daily energy
   * - Diet high in sodium
     - Average 24-hour urinary sodium excretion (in grams per day)
     - 1–5 g/day

**Notes:**

* All intake quantities are in **fresh (as consumed) weight**, matching the GDD dietary intake data convention (see :doc:`current_diets`)
* **GBD risk factors are evaluated for adult populations (≥25 years)** - the current implementation uses population-weighted "All ages" dietary intake averages, which may underestimate risk for adult-only populations
* The model currently implements a subset of these risk factors based on data availability and model scope
* Risk factor definitions specify both the intake measure (e.g., grams per day) and the threshold or optimal range
* "Diet low in" risk factors specify minimum recommended intakes; "diet high in" risk factors treat any intake as risk-increasing
* Milk/dairy measurements use milk equivalents, where cheese and yogurt are converted to their milk equivalent weight
* See :doc:`current_diets` for detailed mapping between GDD dietary intake data and these risk factors

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
3. **Record cluster totals** – store each cluster’s population for scaling; the
   solver multiplies baseline YLLs by the configured ``health.value_per_yll``
   constant (no external valuation dataset required).
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

For each cluster–cause pair the preprocessing step stores
:math:`\mathrm{YLL}^{\mathrm{base}}_{c,g}` (baseline years of life lost). The
solver also records the reference log-relative-risk
:math:`z^{\mathrm{ref}}_{c,g}` (from baseline diets) and its exponential
:math:`RR^{\mathrm{ref}}_{c,g}`. The contribution to the objective is
constructed as

.. math::
   \text{Cost}_{c,g} = V\, \mathrm{YLL}^{\mathrm{base}}_{c,g}
   \left( \frac{RR_{c,g}}{RR^{\mathrm{ref}}_{c,g}} - 1 \right).

A constant term subtracts
:math:`V\,\mathrm{YLL}^{\mathrm{base}}_{c,g}` so that the baseline diet has zero
health cost and only improvements or deteriorations relative to the reference
affect the optimisation.

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

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: health ---
   :end-before: # --- section: aggregation ---

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

References
----------

.. [Brauer2024] Brauer M, Roth GA, Aravkin AY, et al. Global Burden and Strength of Evidence for 88 Risk Factors in 204 Countries and 811 Subnational Locations, 1990–2021: A Systematic Analysis for the Global Burden of Disease Study 2021. *The Lancet*, 2024;403(10440):2162–203. https://doi.org/10.1016/S0140-6736(24)00933-4
