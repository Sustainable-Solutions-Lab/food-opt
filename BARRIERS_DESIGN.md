<!--
SPDX-FileCopyrightText: 2026 Koen van Greevenbroek

SPDX-License-Identifier: CC-BY-4.0
-->

# Barriers to global water-scarcity relief -- design note

Working design document for a proof-of-concept "barriers" study on the
`water-aware` branch. Status: PoC. Not user-facing documentation.

## Thesis

Global agricultural water scarcity is far more relievable than a naive read
suggests: in the fixed-diet model, pricing water cuts accumulated scarcity by
~99% at high water prices. The scientific question is therefore not *whether*
relief is possible but *what it costs on every other axis society cares about*
-- and, empirically, which of those axes actually bind. We reframe the competing
objectives as **guardrails** on the water-relief optimisation and ask, for each:
how much relief survives when the guardrail is honoured, and what is the
marginal cost of honouring it?

Early evidence (the water x carbon x floor grid, 60 solves) already shows the
counter-intuitive shape the paper turns on: relief is *invariant* to carbon
price and cheap to reconcile with value/calorie self-sufficiency, so most
"obvious" barriers are not barriers. The contribution is the **ranking with a
surprise**: identifying the few guardrails that genuinely constrain relief and
quantifying the rest as near-free.

## The uniform framework: guardrail + shadow price

Every barrier is expressed as a **single linear solve-time constraint** on the
existing decision variables (`Link-p`), added after `create_model`, with its
dual persisted to `n.global_constraints` (the `mu` column) so it survives into
analysis. This yields one comparable object per barrier:

- **Relief retained**: accumulated water scarcity with the guardrail enforced
  at its 2020 boundary, relative to the same system with the guardrail absent.
- **Shadow price (cost of protection)**: the constraint dual -- the marginal
  system cost (bnUSD) of honouring one more unit of the guardrail. Directly
  comparable across barriers because they all sit in the same objective.

The primary design samples every guardrail as a binary on/off state while water
price and reallocation flexibility remain continuous. On means "no worse than
the realized frozen-system reference" on that barrier's metric; off means the
constraint is not added. This
avoids mixing mathematical maxima, policy ranges, and optimized reference
outcomes under a nominally common continuous dial.

## The barrier set and why it is not arbitrary

Two established normative frameworks, unioned, fix the set:

- **Planetary boundaries** (Rockstrom 2009; Steffen 2015; food-system
  operationalisations Springmann 2018, Gerten 2020) supply the *environmental*
  guardrails. Water relief eases the freshwater boundary; the guardrails are the
  boundaries agriculture pushes on in exchange -- **climate change** (GHG),
  **land-system change + biosphere integrity** (biodiversity / land conversion),
  and, later, **biogeochemical flows** (nitrogen).
- **FAO's four pillars of food security** supply the *food-system* guardrails:
  **availability** (calorie self-sufficiency), **utilisation** (nutrient / diet
  quality self-sufficiency), **stability** (production concentration), and
  **access** (affordability / diet cost).

Two additions sit outside both frameworks and are labelled as such: **production
value** (rural livelihoods) and **reallocation churn** (feasibility of
adjustment, not a competing objective -- the dominant "barrier" in early runs).

## Operationalisation in GLADE

| Barrier | Framework slot | Instrument | Config | Preliminary verdict |
|---|---|---|---|---|
| Water scarcity (relieved) | freshwater boundary | store price | `water_scarcity.price` | ~99% relievable |
| GHG emissions | climate boundary | store price | `emissions.ghg_price` | free (relief invariant) |
| Biodiversity | land / biosphere | conversion cap (Mha), global + per-region | `biodiversity.cap` | ~free, binds |
| Production concentration | stability | per-crop max-share cap | `production_concentration.cap` | ~free (apple artifact aside) |
| Calorie self-sufficiency | availability | food-energy production floor | `food_energy.floor` | few % |
| Protein self-sufficiency | utilisation | protein production floor | `protein.floor` | few % |
| Production value | livelihoods | GPV floor | `production_value.floor` | ~free |
| Affordability | access | total production-cost cap (global) | `affordability.cost_cap` | to be re-estimated with the 2020 cap |
| Reallocation feasibility | adjustment | per-region deviation cap | `deviation_penalty` | expected dominant |

All nine are solve-time linear guardrails with persisted duals, reported in
`barrier_constraints.parquet`. "Reallocation feasibility" (churn) is labelled
as an *adjustment* axis, not a competing good: it proxies the political/
physical cost of reorganising production, and is the expected dominant barrier.
Affordability is global (a per-country delivered-diet-cost cap is the
distributional "access" variant, deferred).

### New barrier 1 -- biodiversity: natural-land conversion cap

Motivation, from the existing 60-solve grid: cropland expansion (natural land
converted to cropland, carrier `land_conversion`) nearly **doubles** under
water relief, 48.9 Mha at near-zero water price to 95.2 Mha at high water price
(103.7 Mha under joint floors). Relocating production to water-rich regions
requires converting new land there. A biodiversity guardrail limits that.

Constraint (global form):

    sum over links L with carrier in {land_conversion, new_to_pasture} of p[L]
        <= max_conversion_mha

`p[L]` is the converted area in Mha (bus0 flow). Distinct from carbon pricing,
which prices the *LUC CO2* on the same links: an **area** cap also protects
biodiverse-but-low-carbon systems (savanna, steppe, drylands) that a carbon
price undervalues. Report both the conversion carbon pricing alone prevents and
the additional area a hard cap protects -- the wedge is the barrier's
independent content. Global cap is the headline (planetary-boundary "no net
expansion" language); a per-region variant is the leakage-robust extension
(relief can satisfy a global cap by converting in one place and sparing in
another). Dual: bnUSD per Mha of conversion headroom = marginal cost of
protecting natural land.

Config `biodiversity.conversion_cap`: `{enabled, max_conversion_mha}`.

### New barrier 2 -- production concentration: per-crop max-share cap

Stability/fragility guardrail. Relief may concentrate the world's production of
a staple into fewer countries (a shock to a concentrated producer is then
worse). The literal Herfindahl is convex-quadratic; we use the linear
**max-share** form -- no single country may hold more than share `s` of global
production of crop `c`:

    x[c, j]  <=  s * sum_k x[c, k]      for every producing country j

where `x[c, j]` is country j's output of crop c in Mt = sum over that country's
`crop_production` / `crop_production_multi` links and their crop-output ports of
`efficiency_port * p[link]`. Rearranged to `x[c,j] - s * X[c] <= 0` this is
fully linear even with the global total endogenous. One vectorised constraint
block per crop (over its producing countries); the dual per (crop, country) is
the marginal cost of forcing that crop's production to spread -- the cost of
resilience. A richer top-N-share variant (sum-of-k-largest LP epigraph) is a
future extension; max-share is the PoC.

Config `production_concentration.cap`: `{enabled, max_share, min_global_mt}`
(the latter exempts trace crops whose global production is negligible).

Architectural caveat: trade is hub-based (k-means), so *bilateral* import-source
concentration is not representable; production concentration sidesteps the hub
layer entirely and is the cleaner GLADE-native fragility metric.

## The complete study

### The unifying counterfactual

For every competing good, the enforced state is the same normative statement:
the optimized system may not make that outcome worse than its 2020 boundary.
The relaxed state omits the constraint entirely. The headline readout is the
scarcity multiplier from switching the barrier on, averaged over on/off
settings of the other barriers. A negative multiplier is a genuine co-benefit,
not an invalid result: a constraint can redirect the optimum toward a solution
that also uses less scarce water.

Independent 50/50 Bernoulli inputs give the joint GSA an explicit probability
measure. Variance-based indices describe the effect of barrier presence under
that measure; they do not claim that unlike physical units have equal welfare
severity.

### Core experiment

The Sobol design crosses eight binary guardrails with continuous water price
and reallocation flexibility. A 4096-point digital net covers all 256 barrier
combinations exactly 16 times while spreading their continuous coordinates.
The flexible-diet model holds per-country calories at baseline. The frozen 2020
scenario supplies every on-state boundary through its realized guardrail
left-hand sides. This is deliberately not raw build metadata: zero-churn
production requires small feasibility slacks, so only realized values provide
a common feasible witness for every barrier combination.

### Second-order analyses (only for barriers that bite)

- **Onset water-price**: the water price at which each barrier starts eroding
  relief (affordability is expected to bite earliest).
- **Pairwise compounding**: small targeted grids for the biting pairs
  (affordability x reallocation, affordability x nutrient), not the full
  cross-product.
- **Distributional incidence**: per-country/-crop floor and cap duals already
  persisted (desert states, India) -- who bears each barrier.

### Robustness

`nonrenewable_cf` (100 in the primary design), diet source (FBS vs GDD-IA), and
the existing physical-parameter GSA machinery on the headline ranking.

### Analysis & outputs

`barrier_constraints.parquet` (one row per active guardrail: type, crop/country
tag, bound, dual) is the uniform object. `barrier_reference.parquet` records
the frozen per-country, per-crop and global bounds that define "on"; relief and capped quantities come from
`water_metrics`, `luc_breakdown`, `crop_production`, `objective_breakdown`.
Headline figure: barrier-on versus barrier-off scarcity multiplier across the
water-price x reallocation-flexibility plane, plus absolute scarcity contours.

## Resolved design decisions

- Apple/fruit demand-slack degeneracy fixed at the root (demand reconciled with
  expansion-limited production) so concentration and affordability costs are
  clean; `is_slack_pinned` also guards food-price reporting.
- Reallocation feasibility (churn) is IN the barrier set, labelled as the
  adjustment axis.
- Biodiversity: global and per-region/per-country conversion caps.
- Affordability: global production cost may not exceed the frozen 2020 cost
  (per-country diet-cost cap deferred).
- Competing guardrails are binary on/off; reallocation flexibility remains a
  continuous response axis rather than a sampled guardrail.

## Caveats

- The common frozen witness is tested at the all-on corner; any failed sample
  is treated as a model or numerical error rather than silently omitted.
- Fixed diet: food-side barriers are production-side self-sufficiency, not
  consumption adequacy.
- Biodiversity cap overlaps carbon (shared links); report the wedge.
- Global metrics can relieve/convert by relocation; per-region variants are the
  leakage-robust cut.
- Trade is hub-based; concentration is production-side, not bilateral.
