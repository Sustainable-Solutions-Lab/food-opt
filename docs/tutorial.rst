.. SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Tutorial
========

This tutorial walks you through two complete modelling exercises with
GLADE. It assumes you have finished the :doc:`introduction` (clone,
``pixi install``, credentials, manually-downloaded datasets). You'll leave
each part with solved scenarios, auto-generated plots, and a handful of
hand-rolled comparisons built in a notebook.

Both parts use a reduced spatial resolution of 200 optimisation regions so
that they complete in a few minutes on a laptop (once the one-off raw-data
download — about half an hour, depending on your connection — has run). The
tutorial configs live under ``config/tutorial/``.

.. toctree::
   :hidden:

   tutorials/tutorial_01_analysis
   tutorials/tutorial_02_analysis

Part 1 — GHG prices at a fixed diet
-----------------------------------

In this first exercise, we solve three scenarios that are identical except
for the greenhouse-gas price applied to the objective function, and we hold
consumption at the observed 2020 diet in all three. Because the diet is
fixed, every difference between scenarios comes from how **production** —
which crops are grown where, which livestock systems are used, and where
trade flows — reorganises when emissions become more costly.

Step 1 — Look at the config
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Open ``config/tutorial/01_ghg_prices.yaml``. The file is short — every key
not listed here falls back to ``config/default.yaml``:

.. literalinclude:: ../config/tutorial/01_ghg_prices.yaml
   :language: yaml

A few things to note:

* ``name: "tutorial_01"`` controls the output directory: everything lands
  under ``results/tutorial_01/``.
* ``aggregation.regions.target_count: 200`` keeps the LP small enough to
  solve in minutes. The full-resolution default is 750; values below 200
  fail the per-country clustering step because there are more countries in
  the default list than regions.
* ``planning_horizon`` and ``baseline_year`` both default to 2020, aligning
  the model with the most recent year for which GDD-IA dietary data exist.
* The ``scenarios:`` block defines three scenarios that each set
  ``validation.enforce_baseline_diet: true``. That flag forces consumption
  per food group to equal the observed 2020 diet in every country.
* The health module stays off, which is the default. Health costs are the
  subject of separate documentation (:doc:`health`) — we keep them out of the
  tutorial on purpose, which also means you need none of the manually
  downloaded IHME GBD data to follow along.

If you want to experiment, you can copy this file to a new name (e.g.
``config/tutorial/01_my_variant.yaml``), change the ``name`` field, and edit
any overrides you like.

Step 2 — Dry run
~~~~~~~~~~~~~~~~

Before committing to a full run, it's worth asking Snakemake what it *would*
do:

.. code-block:: bash

   tools/smk -j4 --configfile config/tutorial/01_ghg_prices.yaml -n

The ``-n`` flag prints the planned execution graph without actually running
anything. On a clean checkout you'll see data-preparation rules (downloads,
region clustering, yield aggregation), the model build, three solves (one
per scenario), analysis extraction, and plotting.

Step 3 — Run the workflow
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   tools/smk -j4 --configfile config/tutorial/01_ghg_prices.yaml

The first run is the longest because Snakemake has to download raw datasets
(GAEZ, GADM, UN WPP, FAOSTAT, ESA CCI, …). Subsequent runs of any tutorial
or other configuration reuse the same cached data. The build step itself is
shared across all scenarios; only the three solves and the downstream
analysis/plots are scenario-specific.

When the workflow finishes, you will find:

* ``results/tutorial_01/build/model.nc`` — the PyPSA network before solving.
* ``results/tutorial_01/solved/model_scen-{baseline,ghg_mid,ghg_high}.nc`` —
  the three solved networks.
* ``results/tutorial_01/analysis/scen-{baseline,ghg_mid,ghg_high}/*.parquet`` —
  standardised statistics extracted from each solve (see :doc:`analysis` for
  the full schema).
* ``results/tutorial_01/plots/scen-*/*.pdf`` — auto-generated figures.
* ``results/tutorial_01/plots/comparison/`` — cross-scenario comparison
  plots, produced because we set ``plotting.comparison_scenarios: "all"``.

Step 4 — Analyse in a notebook
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The companion notebook :doc:`tutorials/tutorial_01_analysis` walks through
five quick comparisons across the three scenarios: total agricultural
land, the cropland vs grassland split, net GHG emissions by gas, the
composition of animal feed, and the objective-cost breakdown. Open it in
the docs to browse the rendered outputs, or download it and run it locally
against your own ``results/tutorial_01/`` directory.

To run it yourself:

.. code-block:: bash

   pixi run -e dev jupyter lab docs/tutorials/tutorial_01_analysis.ipynb

Because all three scenarios share the same (baseline) diet, anything that
moves between them reflects production-side reorganisation. Total
agricultural land typically falls sharply as the GHG price rises (because
marginal land is released and the regrowing land sequesters carbon), the
gas composition of net emissions shifts, and the objective's ``ghg_cost``
column becomes strongly negative — at these prices, net emissions are
negative, so ``ghg_price × emissions`` is a revenue term in the
objective.

.. note::

   The notebook opens with a short contextualisation that is worth
   reading: even at ``baseline``, this tutorial's model uses less land
   than the real world and produces net-negative emissions by default.
   Serious studies anchor the model toward observed production using the
   calibrated ``market_response`` curves (see ``config/gsa.yaml``), the
   legacy ``deviation_penalty``, or hard constraints (see
   ``config/validation.yaml``). This tutorial uses the L1 deviation penalty
   because its reduced structural config cannot consume the default
   link-specific market-response artefact.

Part 1 — Summary
~~~~~~~~~~~~~~~~

At this point you've exercised the full end-to-end workflow: config, build,
solve, analysis, and custom post-processing. But because consumption was
held fixed, Tutorial 1 can't tell you whether a different *diet* would
reduce emissions more cheaply — the model had no way to weigh "change what
people eat" against "change how food is produced". Part 2 adds that missing
piece.

Part 2 — Letting the diet respond to policy
--------------------------------------------

In Part 1, we fixed consumption with ``enforce_baseline_diet: true``. That
guarantees realism (nobody is forced to eat something unusual), but it also
rules out dietary shift as a mitigation option. A more interesting model lets
the optimiser decide when giving up some of today's diet is worth the GHG
savings — which requires pricing how much today's diet is *worth*.

GLADE prices it with **calibrated demand curves**, using the same
positive-mathematical-programming machinery that anchors production (the
``market_response`` mechanism, :doc:`calibration`). The idea rests on revealed
preference: today's diet is what consumers actually chose at today's prices,
so treat it as an economic equilibrium rather than a constraint.
Operationally:

1. **Measure** (calibration): solve the model once with every (food, country)
   consumption pinned at its observed intake. The shadow price of each pin is
   the gap between what supplying that food costs and what consumers evidently
   value it at — a per-food willingness to pay.
2. **Reproduce**: attach to each consumption link a marginal-utility curve
   through the observed point, with the fitted willingness to pay as its
   level. Today's diet is now the model's *unconstrained* optimum: solve with
   no policy pressure and you get the observed food system back, with nothing
   pinning it there.
3. **Respond**: away from the observed point the curve falls off at a rate
   set by literature demand elasticities (per food group, e.g. staples ~0.55,
   meat and dairy ~0.72). Under a carbon price, emission-intensive foods
   become dearer and consumption slides down each curve — by a lot where
   demand is elastic, barely at all where it is not.

Step 1 — Calibrate the curves
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Open ``config/tutorial/02_demand_response.yaml``:

.. literalinclude:: ../config/tutorial/02_demand_response.yaml
   :language: yaml

The config inherits the market-response defaults; what it declares is a
dedicated artefact set (``calibration.source: "tutorial-02"``). The curves are
keyed to the exact model structure — at the default ``link`` granularity every
production link carries its own curve — so a 200-region tutorial model cannot
reuse the artefacts fit for the full-resolution default config. Fit the
tutorial's own set (one command, roughly 10–20 minutes; add
``CALIBRATE_PIXI_ENV=default`` to use the open-source solver):

.. code-block:: bash

   tools/calibrate --base config/tutorial/02_demand_response.yaml

The result lands in ``data/curated/calibration/tutorial-02/``. The file to
look at is ``market_response.csv``: one row per curve group, whose
``intercept`` column holds the fitted wedges — for ``demand`` rows, minus the
willingness to pay (bn USD/Mt, numerically USD/kg). Skim the beef and tomato
rows for a feel of what the calibration inferred about consumer valuations.

Step 2 — The reference solve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``reference`` scenario turns GHG pricing off and nothing else:

.. code-block:: bash

   tools/smk -j4 --configfile config/tutorial/02_demand_response.yaml -- \
       results/tutorial_02/solved/model_scen-reference.nc

This solve has **no baseline enforcement, no deviation penalty, and no hard
anchoring constraints** — yet it lands within a fraction of a percent of
observed 2020 production and consumption. That is the exact-calibration
property doing its work, and it is worth pausing on: the observed food system
is now the model's own economic optimum, so every later scenario measures a
*departure caused by policy*, not an artefact of loose calibration.

Step 3 — Price the carbon
~~~~~~~~~~~~~~~~~~~~~~~~~

Solve the remaining scenarios:

.. code-block:: bash

   tools/smk -j4 --configfile config/tutorial/02_demand_response.yaml

Alongside the ``ghg_mid``/``ghg_high`` ladder, two variants make the demand
side itself the experiment:

* ``ghg_high_fixed`` freezes the diet at the observed baseline
  (``components.demand: false`` + ``enforce_baseline_diet``) under the same
  $200 carbon price — production-side mitigation only.
* ``ghg_high_stiff`` halves every demand elasticity via
  ``elasticity_factor: 0.5``. Elasticities enter only at solve time, so
  scanning them needs no recalibration.

Step 4 — Compare in the notebook
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The companion notebook :doc:`tutorials/tutorial_02_analysis` builds three
comparisons:

* **What shifts**: global food-group consumption across the price ladder —
  which foods absorb the carbon price and which barely move.
* **The two sides of the market**: the objective breakdown now carries
  ``supply_response`` and ``demand_response`` columns — the cost of moving
  production and consumption off their observed patterns — next to
  ``ghg_cost``, the force pushing them.
* **The value of flexibility**: net emissions and total cost of ``ghg_high``
  vs ``ghg_high_fixed`` vs ``ghg_high_stiff``. The fixed-vs-flexible gap is a
  direct measure of demand-side mitigation potential; the stiff variant shows
  how much of it hinges on the elasticity assumption.

Gotchas
~~~~~~~

* The intercept artefact is **specific to the model structure it was fit
  against**. Change the region count, crop set, or diet source and the solve
  fails loudly on missing curve groups — rerun
  ``tools/calibrate --base <config>`` for the new structure.
* ``market_response.components.demand`` and
  ``validation.enforce_baseline_diet`` are mutually exclusive: a scenario
  pinning the diet must switch the demand component off (see
  ``ghg_high_fixed``).
* Demand curves value consumption; they do not guarantee adequacy. The
  nutrition constraints (:doc:`nutrition`) still put floors under nutrients,
  so extreme carbon prices reshape diets rather than shrinking them away.

Where to go from here
---------------------

You have now solved two small scenario sets, inspected the output files, and
built a handful of comparisons by hand. Some natural next steps:

* **Scale up the GHG price sweep.** Rather than listing prices by hand,
  generate them programmatically with the
  :doc:`scenario generator DSL <configuration>` -- for example a log-spaced
  sweep from 5 to 500 USD/tCO2-eq:

  .. code-block:: yaml

     scenarios:
       _generators:
         - name: "ghg_{ghg}"
           parameters:
             ghg:
               space: log
               start: 5
               stop: 500
               num: 10
               round: true
           template:
             emissions:
               ghg_price: "{ghg}"
* **Turn on health costs.** :doc:`health` describes the Global Burden of
  Disease integration and how ``health.value_per_yll`` prices diet-related
  disease burden alongside the environmental objectives.
* **Perform a global sensitivity analysis.** :doc:`sensitivity_analysis`
  describes the polynomial-chaos and random-forest surrogate workflows used
  for Sobol-index decomposition.
* **Learn the rule graph.** :doc:`workflow` documents every rule in the
  pipeline; :doc:`results` and :doc:`analysis` document every output file
  and column.
