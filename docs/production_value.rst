.. SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Production value
================

The model tracks each country's *gross agricultural production value*:
every crop and animal product output valued at fixed producer prices,
following the FAOSTAT gross production value convention. An optional
solve-time constraint bounds each country's value from below, protecting
the size of its agricultural sector while leaving the production mix
free -- e.g. to study whether water-scarcity relief can avoid gutting
agriculture in the countries it relocates production away from.

Definition
----------

For country :math:`c`,

.. math::

   V_c = \sum_{\ell \in c} v_\ell \, p_\ell

where :math:`p_\ell` is the dispatch of each crop production link
(single and multi-cropping, in Mha) and animal production link (in Mt
feed DM), and :math:`v_\ell` converts dispatch to value (bnUSD per unit)
by summing ``efficiency * price`` over the link's marketable output
ports. Design choices:

* **Gross, not net**: marketed feed crops (maize, soybean, ...) are
  valued *and* so are the animal products they become (sector size, not
  value added). Grassland output and the fodder crops listed in
  ``production_value.exclude_crops`` (alfalfa, silage maize, biomass
  sorghum) are unmarketed intermediates and are not valued -- fodder
  carries only proxy grain prices, and valuing it lets endogenous fodder
  expansion masquerade as production value (and satisfy the floor).
  This matches the FAOSTAT GPV scope. Animal co-product ports (rendered
  fats) are not valued.
* **Prices**: FAOSTAT PP producer prices per (item, country),
  CPI-deflated and averaged over the same window as production costs
  (see :doc:`costs`), prepared by ``prepare_producer_prices``. Crop
  prices are converted to the dry-matter crop-bus basis via
  ``1/(1 - moisture)``. Meat prices are carcass-basis while model meat
  output is retail weight, so livestock value is conservatively low in
  cross-commodity comparisons; this cancels within the floor constraint.
* **Baseline**: the same coefficients evaluated at the baseline
  quantities (``baseline_area_mha``, ``baseline_feed_use_mt_dm``), so
  floor and realized value share all model conventions and the
  constraint binds on the price-weighted production mix rather than on
  price levels.

Floor constraint
----------------

``production_value.floor`` (solve-time, scenario-overridable)::

   production_value:
     floor:
       enabled: true
       fraction: 0.9          # retain >= 90% of baseline value
       price_basis: observed  # or: uniform (single world price per item)
       min_baseline_bnusd: 0.01

At baseline dispatch :math:`V_c` equals its baseline by construction, so
any fraction <= 1 is feasible whenever the baseline itself is; joint
infeasibility with tight water or emission caps is an informative result
(the country cannot retain that value share within the caps). The
per-country shadow price (bnUSD of system cost per bnUSD of protected
value) is stored on the ``production_value_floor`` global-constraint
rows and reported in the ``production_value`` analysis output.

Interpretation pitfalls
-----------------------

* **Support prices are inside the prices.** Observed producer prices
  include price support (e.g. Indian MSP for paddy), so the floor
  partially protects politically constructed revenue. Use
  ``price_basis: uniform`` to isolate the quantity-mix effect.
* **Calibration-correction steering.** The objective satisfies the
  floor with the cheapest available value, and "cheapest" includes
  cost-calibration corrections; where corrections are large (structural
  data mismatch), the floor may be met by exactly those links. Check the
  value composition in the ``production_value`` output before
  interpreting.
* **Feed-import escape hatch.** Livestock fed on imported feed
  generates gross value with little land or water, so a binding floor
  under a water cap can be met by expanding animal production. Read
  value results jointly with feed trade.
* **Stability-penalty interplay.** The floor relies on reallocating the
  production mix, which the deviation penalties tax; part of a binding
  floor's measured cost is behavioral adjustment friction, not resource
  cost.

Food-energy floor
-----------------

``food_energy.floor`` is the calorie-side twin: a per-country lower
bound on primal production of *human-edible food energy*, with
coefficients built from ``nutrition.csv`` energy densities (crops at
their best-food-use pathway energy; crops with no food pathway carry
zero). This is the production-capacity notion of calorie
self-sufficiency: with fixed diets, keeping domestic calorie production
near baseline caps how much more import-dependent a country can become,
without attributing pooled hub-trade flows. Both floors can be active
simultaneously (value protects the sector economically, energy protects
the calorie base; a country can satisfy one without the other)::

   food_energy:
     floor:
       enabled: true
       fraction: 0.9
       min_baseline_tkcal: 1.0   # Tkcal = 1e12 kcal

Analysis output
---------------

``production_value.parquet`` (per scenario): one row per (country,
carrier, item) with ``realized_value_bnusd``, ``baseline_value_bnusd``,
and per-country ``floor_bnusd`` / ``floor_mu`` when the floor was
active. Coefficients are stamped on every solve, so the output is
available whether or not the floor is enabled.

``food_energy.parquet`` (per scenario): one row per country with primal
energy production (realized and baseline), the energy floor level and
shadow price, gross/net hub-trade calorie imports (crops + foods at
human-edible energy; feed trade reported separately in Mt DM), and
final energy consumption -- the ingredients of an import-dependence
ratio.
