.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Health Impacts
==============

Overview
--------

The health module translates dietary choices into population health outcomes measured in DALYs (Disability-Adjusted Life Years). This captures the disease burden associated with suboptimal diet composition, based on epidemiological dose-response relationships from the Global Burden of Disease (GBD) study.

Health impacts are incorporated into the objective function, allowing the model to trade off environmental costs against health costs.

Dietary Risk Factors
--------------------

The model tracks eight major dietary risk factors linked to chronic disease:

Configuration
~~~~~~~~~~~~~

.. code-block:: yaml

   health:
     risk_factors:
       - fruits              # Low fruit intake
       - vegetables          # Low vegetable intake
       - nuts_seeds          # Low nuts/seeds intake
       - legumes             # Low legume intake
       - fish                # Low seafood omega-3 intake
       - red_meat            # High red meat intake
       - prc_meat            # High processed meat intake
       - whole_grains        # Low whole grain intake

Each risk factor has:

* **Optimal intake**: Level associated with minimum disease risk
* **Dose-response curve**: How mortality/morbidity changes with intake
* **Attributable diseases**: Which conditions are affected (IHD, stroke, diabetes, colorectal cancer, etc.)

Data Sources
~~~~~~~~~~~~

The health module combines data from multiple sources:

**Mortality rates**: `IHME Global Burden of Disease Study 2021 <https://vizhub.healthdata.org/gbd-results/>`_
  * Cause-specific death rates by country and age
  * Processed via ``workflow/scripts/prepare_gbd_mortality.py``
  * See ``data/DATASETS.md`` for download instructions

**Baseline dietary intake**: `Global Dietary Database <https://www.globaldietarydatabase.org/>`_ (Tufts University)
  * Country-level mean daily intake for major food groups and dietary risk factors
  * Based on systematic review and meta-analysis of national dietary surveys
  * Processed via ``workflow/scripts/prepare_gdd_dietary_intake.py``
  * See ``data/DATASETS.md`` and ``data/manually_downloaded/README.md`` for download instructions

**Risk factor relationships**: `WHO-DIA repository <https://github.com/marco-spr/WHO-DIA>`_ (GPL-3.0 licensed)
  * Relative risk curves linking dietary intake to disease risk
  * Input files (in ``data/health/raw/``):
    - ``RR_int_05282021.csv``: Relative risk breakpoints (intake thresholds)
    - ``RR_max_05282021.csv``: Maximum relative risk at extreme intakes

Food-to-Risk Mapping
---------------------

Foods are mapped to risk factors in ``data/health/food_to_risk_factor.csv``. Example:

* ``apple`` → ``fruits`` (1.0 weight)
* ``beef`` → ``red_meat`` (1.0 weight)
* ``sausage`` → ``prc_meat`` (1.0 weight)
* ``whole_wheat_bread`` → ``whole_grains`` (1.0 weight)

Some foods contribute to multiple risk factors or with fractional weights.

Dose-Response Relationships
----------------------------

The epidemiological model relates intake to relative risk (RR):

.. math::

   RR(\text{intake}) = \exp(\beta \times (\text{intake} - \text{optimal}))

Where:

* :math:`\beta`: Slope parameter from meta-analyses
* **optimal**: Intake level minimizing risk (e.g., 300 g/day fruits)
* **RR**: Relative risk of mortality/disease compared to optimal intake

Piecewise-Linear Approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To keep the optimization linear, the model approximates log(RR) curves with piecewise-linear segments:

1. **Breakpoints**: Divide intake range into bins (configured by ``health.intake_grid_step``)
2. **Linear interpolation**: Between breakpoints, use linear approximation of log(RR)
3. **Binary variables**: Not used—instead, the model uses SOS2 (Special Ordered Set of type 2) constraints or convex combination tricks to maintain linearity

This allows capturing diminishing returns (e.g., going from 0→100 g/day fruits has bigger health benefit than 200→300 g/day) without nonlinear optimization.

Regional Clustering
-------------------

To reduce computational burden, countries with similar health profiles are clustered into health regions.

Configuration
~~~~~~~~~~~~~

.. code-block:: yaml

   health:
     region_clusters: 30         # Number of health clusters
     reference_year: 2018        # Baseline year for health data
     intake_grid_step: 10        # g/day granularity for dose-response
     log_rr_points: 10           # Points for log(RR) linearization

Clustering Process
~~~~~~~~~~~~~~~~~~

The ``prepare_health_costs`` rule (``workflow/scripts/prepare_health_costs.py``):

1. **Load baseline data**: Country-level dietary intake, mortality, demographics
2. **Cluster**: Group countries by similar baseline health burdens (k-means on baseline DALYs)
3. **Compute dose-response**: For each cluster, calculate risk breakpoints and slopes
4. **Valuation**: Apply the configured value of a statistical life (currently a global constant; the code retains optional hooks for future regional datasets).
5. **Output**:

   * ``processing/{name}/health/risk_breakpoints.csv``: Intake thresholds
   * ``processing/{name}/health/cluster_cause_baseline.csv``: Baseline disease burden
   * ``processing/{name}/health/cause_log_breakpoints.csv``: Linearized log(RR) segments
   * ``processing/{name}/health/country_clusters.csv``: Country → cluster mapping

This creates a representative health profile for each cluster, reducing the problem size from ~150 countries to ~30 clusters.

DALY Calculation
----------------

DALYs combine mortality and morbidity:

.. math::

   \text{DALYs} = \text{YLL} + \text{YLD}

* **YLL** (Years of Life Lost): Premature deaths × years lost per death
* **YLD** (Years Lived with Disability): Non-fatal disease burden × disability weights

The model focuses on mortality (YLL) for computational simplicity, as dietary risk factors primarily affect mortality risk.

Calculation Steps
~~~~~~~~~~~~~~~~~

1. **Baseline mortality**: Country-specific death rates by cause (IHD, stroke, diabetes, CRC)
2. **Population Attributable Fraction (PAF)**: % of deaths attributable to suboptimal diet

   .. math::

      PAF = \frac{RR - 1}{RR}

3. **Attributable deaths**: Baseline deaths × PAF
4. **Years of life lost**: Deaths × age-specific life expectancy
5. **Total DALYs**: Σ(attributable deaths × YLL)

Value of Statistical Life
--------------------------

DALYs are monetized using the Value of Statistical Life Year (VSLY) to make health costs commensurable with economic and environmental costs.

Configuration
~~~~~~~~~~~~~

.. code-block:: yaml

   health:
    value_of_statistical_life: 3_500_000  # USD per life (global constant)

Options:

* **Constant**: Single global VSLY (e.g., 3.5M USD, roughly US EPA value)
* **"regional"**: Use region-specific VSL from DIA dataset (higher in high-income countries)

The choice affects optimization priorities:

* **High VSLY**: Model heavily weights health outcomes, may accept higher environmental costs for healthier diets
* **Low VSLY**: Environmental costs dominate, nutrition meets minimums but health optimization is secondary

Health Cost in Objective
-------------------------

Health costs enter the objective function as:

.. math::

   \text{Health cost} = \sum_{\text{clusters}} \text{DALYs}_{\text{cluster}} \times \text{VSLY}

This competes with:

* Production costs
* Trade costs
* Environmental costs (emissions × carbon price)

The optimal solution balances these trade-offs.

Model Integration
-----------------

Health constraints are added during solving (``workflow/scripts/solve_model.py``), not model building, because they require:

1. **Baseline burden**: Loading pre-computed health cluster data
2. **Food consumption variables**: Must be defined first in the model
3. **Risk factor aggregation**: Summing food consumption → risk factor intake
4. **Piecewise-linear constraints**: Linking intake to log(RR) to DALYs

Process:

1. **Load model**: Read built PyPSA network
2. **Load health data**: Risk breakpoints, dose-response, baseline burden
3. **Create risk intake variables**: Σ(food consumption × food-to-risk weights)
4. **Create DALY variables**: Link intake → RR → attributable deaths → YLL
5. **Add to objective**: DALYs × VSLY
6. **Solve**: Optimize with health costs included

Configuration Parameters
------------------------

.. code-block:: yaml

   health:
     region_clusters: 30               # Number of health clusters
     reference_year: 2018              # Baseline year for mortality data
     intake_grid_step: 10              # g/day resolution for breakpoints
     log_rr_points: 10                 # Linearization points for log(RR)
    value_of_statistical_life: 3_500_000  # USD (set "regional" only if dataset provided)
     risk_factors:                     # Which risk factors to include
       - fruits
       - vegetables
       - nuts_seeds
       - legumes
       - fish
       - red_meat
       - prc_meat
       - whole_grains

Reducing ``region_clusters`` or ``log_rr_points`` speeds up solving at the cost of health resolution.

Visualization
-------------

Health impact results can be visualized:

**Health risk map**::

    tools/smk results/{name}/plots/health_risk_map.pdf

Shows spatial distribution of dietary risk-attributable DALYs.

**Health baseline map**::

    tools/smk results/{name}/plots/health_baseline_map.pdf

Shows baseline (pre-optimization) health burden for comparison.

**Regional health breakdown**::

    tools/smk results/{name}/plots/health_risk_by_region.csv
    tools/smk results/{name}/plots/health_baseline_by_region.csv

CSV exports for detailed analysis.

**Objective breakdown**::

    tools/smk results/{name}/plots/objective_breakdown.pdf

Shows contribution of health costs to total objective value.

Scenario Exploration
--------------------

Health module enables exploring:

**Diet Quality vs. Environmental Impact**

* High VSLY → healthier diets (more fruits/vegetables, less red meat) even if higher environmental costs
* Low VSLY → environmentally optimal but potentially lower diet quality

**Trade-offs Between Risk Factors**

* Reducing red meat → lower IHD/CRC risk but may require other protein sources
* Increasing nuts/legumes → health benefits but land use implications

**Effectiveness of Dietary Guidelines**

* Compare optimized diet to EAT-Lancet or national dietary guidelines
* Assess if guidelines balance health, environment, and production constraints

Limitations and Future Work
----------------------------

Current limitations:

* **Mortality focus**: Doesn't capture morbidity (YLD), underestimates full burden
* **Linear approximation**: Piecewise-linear may miss fine-grained nonlinear effects
* **Aggregate risk factors**: Doesn't distinguish subtypes (e.g., processed vs. unprocessed red meat)
* **No nutrient interactions**: Risk factors treated independently

Future enhancements:

* **Morbidity**: Add YLD for diabetes, obesity-related conditions
* **Micronutrient deficiencies**: Iron, vitamin A, zinc deficiency burdens
* **Age-structured**: Different optimal intakes for children vs. adults
* **Dynamic health**: Multi-period model with health transitions
