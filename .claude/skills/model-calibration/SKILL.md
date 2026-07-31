---
name: model-calibration
description: Run, refresh, or diagnose the model's calibration pipeline (feed -> food_waste -> food_demand -> cost -> stability) that produces the per-config artefact sets under `data/curated/calibration/<source>/` (the `default` and `gbd-anchored` sets are git-tracked). Covers the dependency order, the `tools/calibrate` wrapper, realistic runtime expectations, when each kind of upstream change forces a re-run, and how to diagnose the most common failure mode: a hidden supply/demand mismatch that inflates the production-stability L1 cost. Use whenever calibration is relevant -- the user touches inputs/build logic that feed the calibration solves, calibration artefacts look off, or a refresh of the artefacts is needed after a model/data change.
---

<!--
SPDX-FileCopyrightText: 2026 Koen van Greevenbroek

SPDX-License-Identifier: CC-BY-4.0
-->

# Model Calibration

The default workflow consumes four calibration artefact groups organized
in per-config *sets* under `data/curated/calibration/<source>/`, selected
by the `calibration.source` config key. Two sets are git-tracked, both
fit against the default GDD-IA diet source: `default` (fit against the
anchoring-off baseline diet of the health-off default config) and
`gbd-anchored` (fit against the GBD-anchored diet; consumed by the
health-enabled configs gsa, gsa_fixed_diet, validation and the doc
configs). No set is shipped for `diet.source: fbs`, which therefore
needs its own calibration run. `tools/calibrate` resolves the base config's
diet.anchor_groups_to_gbd sentinel once and pins it across every
step, and provenance stamps record the *resolved* anchoring. Each
artefact group is produced by a dedicated validation-mode solve and
absorbs a specific class of residual mismatch so that ordinary solves
don't have to. Without these files in place, production-stability,
costs, and food/feed accounting drift from observed 2020 reality.

Every set carries a `provenance.yaml` stamp of the structural config it
was fit against; workflow runs error at DAG time when their config
differs structurally from the consumed set's stamp (see "Artefact sets
and provenance" below).

Authoritative reference: `docs/calibration.rst`. This skill is the operational
companion: when to run, how to run, what to expect, what to watch out for.

## The steps at a glance

Default chain (run by `all`, in order):

| Step | Config | Produces | Absorbs |
|---|---|---|---|
| feed | `config/calibration/feed.yaml` | `grassland_yield.csv`, `fodder_conversion.csv`, `exogenous_forage.csv`, `exogenous_feed.csv` | Per-country forage, and protein- and roughage-feed supply/demand gaps |
| food_waste | `config/calibration/food_waste.yaml` | `food_waste.yaml` | Per-food-group consumer-side waste multiplier (FBS supply vs GDD-IA intake) |
| food_demand | `config/calibration/food_demand.yaml` | `food_demand.csv` | Per-food residual mismatch left over after food_waste |
| market_response | `config/calibration/market_response.yaml` | `market_response.csv` | Per-group PMP curve intercepts (production wedges + demand intercepts with their slope-price basis), fit by sequential pinned solves |

Legacy anchoring steps, outside the default chain, run explicitly (in
this order -- stability depends on cost) for artefact sets consumed by
configs that anchor with them:

| Step | Config | Produces | Absorbs |
|---|---|---|---|
| cost | `config/calibration/cost.yaml` | `crop_cost.csv`, `multi_crop_cost.csv`, `grassland_cost.csv`, `animal_cost.csv` | Additive production-cost corrections from stability-constraint duals |
| stability | `config/calibration/stability.yaml` | `deviation_penalty.yaml` | The L1 penalty triple over the configured components (default: land + feed) that gives ~5% deviation on each axis |

The order is strict (next section explains why). The step configs each
carry `name: "calibration"`, but `tools/calibrate` overrides it per step
with `name: calibration-<source>-<step>`, so each step gets its own
processing tree. That isolation is deliberate -- the steps disable
different calibration blocks, so a shared tree would thrash `build_model`
-- but it means the identical upstream prep is rebuilt once per step
(~100 s and ~250 MB each, so a large share of chain wall-clock and
disk for a full run).

## Strict dependency order

Each step solves against a model whose earlier-step calibrations are
already applied. Re-running out of order leaks residual mismatch into
later artefacts and corrupts them in non-obvious ways:

- `feed` first -- later solves rely on the feed slack already being closed.
- `food_waste` next -- food-bus slack must not be contaminated by feed-side mismatch.
- `food_demand` next -- the per-food residual must not be mis-attributed to waste.
- `market_response` last -- the pinned solves measure wedges against the fully calibrated feed/waste/demand behaviour.
- Legacy: `cost` after `food_demand` -- without it, per-food mismatch leaks into cost duals as spurious sign (olive-oil goes negative, coffee/tea peg at the slack ceiling); `stability` after `cost` -- the L1 Broyden iteration assumes the cost corrections are in place.

If in doubt, run the whole chain. Partial re-runs are only safe when every
earlier artefact is still semantically valid (see "When to re-run").

## Entry point

```bash
tools/calibrate              # the default chain in order
tools/calibrate feed         # one step
tools/calibrate food_waste
tools/calibrate food_demand
tools/calibrate market_response
tools/calibrate cost         # legacy, not in the default chain
tools/calibrate stability    # legacy, not in the default chain
tools/calibrate --check      # per-step staleness + provenance probe (no execution)
tools/calibrate --record     # re-stamp the set as it stands (no execution)
tools/calibrate --base config/<name>.yaml [all|<step>|--check]
                             # calibrate a dedicated set for another config
```

The wrapper defaults to `pixi -e gurobi` -- all calibration configs use
Gurobi. HiGHS is too slow here. Override with `CALIBRATE_PIXI_ENV=<env>`.

With `--base`, the base config must declare its own `calibration.source`
(refusing to overwrite the shared `default` set); a fresh set is
generated in dependency order with no seeding (each step consumes only
its predecessors' artefacts; legacy artefacts appear only when the
`cost`/`stability` steps are run explicitly), and the `all` chain uses
`name: calibration-<source>-<step>` so processing trees don't thrash. After any
successful run the set is (re)stamped with `provenance.yaml`.

Pass extra flags through positionally:

```bash
tools/calibrate cost -j8
```

For local runs, `SMK_MEM_MAX=30G` is a safe ceiling on this workstation;
see [local_machine_smk_defaults.md](../../memory/local_machine_smk_defaults.md).

## When to re-run calibration

### Always probe first

```bash
tools/calibrate --check
```

Takes ~20 s for the `default` set, ~30 s for `gbd-anchored`. It compares
content hashes recorded in `data/curated/calibration/<source>/fingerprint.yaml`
against the step's external DAG leaves (curated/bundled/manually-downloaded
data, rule files, scripts), its config (YAML round-tripped, so comments don't
count), and the artefacts the step wrote. Everything inside the artefact set is
excluded, so the answer is stable: a completed chain reports all-green, and
mtime churn from `git checkout` or `touch` is invisible.

`tools/smk` also prints a one-line mtime-based reminder on every invocation.
That one *is* advisory and does produce false positives after `git pull`;
silence it with `SMK_SKIP_CALIBRATION_HINT=1`.

Editing any rule file or script in a step's DAG marks it stale, including
inert changes (refactors, comments, unrelated rules). When you have
established the artefacts are unaffected, `tools/calibrate --record`
re-stamps the set without solving. It asserts the artefacts are current, so
don't reach for it to make a genuine `[STALE]` go away.

### Triggers for a full re-run (start from `feed`)

Default to a full re-run whenever something upstream of the calibration
solves changes materially. The chain is short and the artefacts are
internally consistent only as a set. Concretely:

- **Feed supply/demand:** GLEAM3 inputs, grassland yields, fodder/forage logic, feed-category routing, feed conversion factors, animal-product baselines.
- **Food-bus mass balance:** food group definitions, food-loss/waste defaults, processing pathways, FBS/QCL preparation, baseline-diet (GDD-IA / FBS) construction, within-group share logic.
- **Crop/grassland/animal cost inputs:** FAOSTAT prices, FADN animal costs, GAEZ/yield calibration touching the cost-baseline solve.
- **Production-stability machinery:** `workflow/scripts/build_model/` changes that touch stability links, constraint sign, slack semantics, or deviation accounting.
- **LUC / land-carbon overhauls** that change cropland or pasture baselines (recent precedent: `4cada27` LUC revert -> `f3f86f3` refresh).
- **Region or resource-class changes** that move the baseline production map.

### Triggers for a partial re-run

You can resume from a later step only if every earlier artefact is still
semantically valid. Rule of thumb:

- Anything upstream of cost in the feed/food chain -> start from `feed`.
- Changes only to cost-extraction logic or to `cost.yaml`/`crop_cost*` inputs -> start from `cost`.
- Changes only to the L1 mechanism or to the `stability.yaml` config -> start from `stability`.

**When in doubt, run the full chain.** A few extra minutes is cheaper than
diagnosing a silently miscalibrated downstream solve.

### What does NOT need a re-run

- Solve-time-only knobs: GHG price, value_per_yll, scenario overrides.
- Analysis, plotting, surrogate fitting, paper notebooks.
- Pure refactors that don't touch model construction, mass balance, or solve-time wiring.

## Realistic runtime (local, gurobi)

End-to-end wall-clock of `tools/calibrate ... all`, measured on this
workstation at `-j8` with the shared `data/` and `processing/shared/`
trees already populated. Each step's time includes rebuilding its own
processing tree (see above), which is why no step comes in under a
minute even though the solves themselves are seconds.

| Step | Solves | `default` set | water study (750 regions, T=4) |
|---|---|---|---|
| feed | 1 validation solve + extract | 96 s | 145 s |
| food_waste | 1 validation solve + extract | 100 s | 145 s |
| food_demand | 1 validation solve + extract | 91 s | 151 s |
| cost | 2 paired solves (hard stability) + extract | 153 s | 198 s |
| stability | Broyden iterations | 193 s (1 iter) | 551 s (5 iters) |
| **full chain** | | **10.6 min** | **19.9 min** |

Of each step's time, only a small part is the solve. On the water study:
`build_model` is 16-18 s and a validation solve 23 s (cost's paired
solves 35 s each); the Broyden loop is the one genuinely long stage at
422 s for 5 iterations. Everything else is the per-step processing-tree
rebuild.

Stability dominates the variance: warm-starting from an existing
`deviation_penalty.yaml` converged in 1 iteration on the `default` set,
while the water set needed 5 from a seeded start. Budget 3-9 minutes for
that step alone.

These are unusually fast because every calibration config disables the
health objective (`value_per_yll: 0`), which makes the problem a pure LP.
A default-config solve carries the health term -- a piecewise-linear
approximation of log/exp products that introduces integer variables --
turning the problem into a MILP and slowing each solve noticeably.
Don't extrapolate these numbers to default-config solves.

They also assume the shared upstream data is already on disk. On a fresh
clone or after a deep model change, add the one-time build of the shared
prep -- 10-30 minutes depending on what has to be rebuilt and whether
retrieval rules fire (those need network access).

HPC offloading isn't worthwhile -- the feed/food/cost steps are single
solves with significant local prereqs, and stability is inherently
sequential (each Broyden iteration depends on the previous solve).

## Output landing zones

- `data/curated/calibration/<source>/*` -- one artefact set per base config, plus its `provenance.yaml` stamp and `fingerprint.yaml` input hashes; the `default` and `gbd-anchored` sets are **git-tracked**. Commit a set together as a refresh; mixed-vintage artefacts are the most common cause of confusing downstream solves.
- `processing/calibration-<source>-<step>/*` -- one upstream prep tree per step (~250 MB each), NOT committed.
- `results/calibration-<source>-<step>/*` -- per-iteration solve logs, NOT committed.
- `results/calibration-<source>-stability/calibration/deviation_penalty_trace.csv` -- per-iter Broyden trace (per-component lambda, achieved deviations, residual norm). Inspect when stability behaves oddly.

## Artefact sets and provenance

- A config selects its set with `calibration.source` (default: `default`); all artefact paths resolve through the `{calibration_source}` placeholder at config-load time.
- Structurally divergent configs must either calibrate their own set (`calibration.source: <name>` + `tools/calibrate --base config/<name>.yaml`), point at a compatible set, or set `calibration.accept_provenance_mismatch: true` (test/tutorial-grade escape hatch: warning instead of error).
- The provenance check covers structural config drift only; input and code staleness is the fingerprint's job. Both run from `tools/calibrate --check`.
- The stamp compares all non-solve-time leaves minus exempt machinery keys (see `PROVENANCE_EXEMPT_PREFIXES` in `workflow/validation/calibration_provenance.py`). Solve-time knobs (GHG price, value_per_yll, deviation_penalty, scenario overrides) never trip it.
- `tests/test_calibration_provenance.py::TestDefaultStampConsistency` fails when `config/default.yaml` changes structurally without a recalibration/restamp -- that is the intended forcing function.

The currently calibrated L1 centre lives in
`data/curated/calibration/default/deviation_penalty.yaml` under
`l1_costs.<component>`. Solves that set
`deviation_penalty.{land,feed,diet}.l1_cost: "calibrated"` resolve the
sentinel from this file at solve time. Per-component
`l1_cost_factor` lets scenarios scan around the calibrated value
without hard-coding absolute numbers.

## Gotchas and pitfalls

### The big one: inflated production-stability L1 cost as a mismatch signal

**Symptom.** After `tools/calibrate stability` finishes,
`l1_costs.land` and/or `l1_costs.feed` in `deviation_penalty.yaml` are
several-fold larger than the previously calibrated centre (current
centre: ~0.10 land, ~0.03 feed). Or Broyden refuses to converge at the
5 % deviation target and oscillates / hits max iterations.

**Mechanism.** A pure cost-minimisation solve is free to reorganise
production arbitrarily; the L1 penalty is what coerces the LP back into
observed allocation. When earlier calibration steps successfully absorb
their respective gaps, the LP can satisfy demand cheaply close to the
baseline pattern, and a modest L1 cost is enough to pin it there. But if
a residual supply/demand mismatch *survives* the feed / food_waste /
food_demand steps -- because of a new bug or a structural change
upstream -- the LP can only match observed production by being physically
coerced through every constraint. The L1 cost is the price of that
coercion, so it inflates.

**Anti-pattern.** Running `tools/calibrate stability` until something
finally converges and committing the inflated L1 cost. That hides the
real bug and silently changes the model's stiffness for every downstream
scenario.

**Diagnosis recipe.**

1. **Solve `config/validation.yaml`** with no stability and inspect food-bus and feed-bus slack:
   ```bash
   tools/smk --configfile config/validation.yaml -- \
       results/validation/solved/model_scen-default.nc
   ```
   The validation scenario has `enforce_baseline_diet: true`,
   `enforce_baseline_feed: true`, no stability penalty, and a modest
   `slack_marginal_cost`. Any slack on the food or feed buses reveals a
   supply/demand gap the calibration chain failed to close.
2. **Compare per-group / per-food / per-feed-category slack** against the
   previously expected post-calibration baseline. The food_waste and
   food_demand intermediates
   (`processing/calibration-<source>-food_waste/food_waste_uncal*`,
   `processing/calibration-<source>-food_demand/food_demand_uncal*`) from
   the previous calibration run are the natural reference for what the gap
   looked like before the change.
3. **A food group with >5 % net slack, or a feed category with a
   shortage that wasn't there before, is the prime suspect.** Trace the
   upstream change before touching the calibration chain.
4. **Other plausible causes** (less common, but real): broken
   food->nutrient routing, a missing trade-hub edge, a mass-balance sign
   flip in a new processing pathway, a unit-conversion bug in a
   feed-conversion factor, an LUC sign error in a pasture flow.

Only after the mismatch is understood and either fixed or accepted should
the calibration chain be re-run. The L1 cost should land back near the
historical centre once the chain absorbs its gaps cleanly.

### Other gotchas

- **Commit an artefact set whole.** The steps write into one
  directory but each is fit against its predecessors, so a set is only
  meaningful as a single vintage. Committing one artefact requires either
  the whole set to be consistent or a careful re-run from the step where
  inconsistency starts. Don't commit them piecemeal.

- **Each step must disable the calibrations of every *later* step.** The
  chain is a single forward pass, not an iteration to a fixed point, so a
  step that consumes a successor's artefact is fit against whatever
  vintage happens to be on disk. That artefact is then an external leaf of
  the step's DAG, so `tools/calibrate --check` goes stale every time the
  successor is regenerated, and the chain never settles. The disable blocks in the step configs form
  a lower triangle (feed disables the other four, food_waste the last
  three, and so on); preserve it when adding a calibration.

  The trap is that a step need not *enable* a successor's calibration to
  depend on it -- inheriting `l1_cost: "calibrated"` from `default.yaml` is
  enough. feed/food_waste/food_demand therefore set
  `deviation_penalty.enabled: false` (production is pinned to actuals, so the
  penalty has nothing to act on) and cost drives stability through hard
  bounds. When adding a calibrated key to `default.yaml`, check what the
  earlier step configs inherit.

- **Step configs must not override structural model keys.** Anything a
  step config pins that isn't calibration machinery, the validation
  regime, or a solver budget silently defeats the point of a per-base
  artefact set, because `provenance.yaml` stamps the *base* config. A
  `water.data.availability` pin in `cost.yaml` did exactly this and moved
  cost corrections by up to 30% on water-resolving configs.

- **Sentinel `"calibrated"` fails loudly when the yaml is stale or
  missing.** `deviation_penalty.{land,feed,diet}.l1_cost: "calibrated"`
  resolves from `deviation_penalty.yaml` at solve time. Scenarios that
  want a fixed numeric value should override the sentinel directly.

- **`penalty_mode` and the sentinel are coupled.** The sentinel resolves
  only under `penalty_mode: "l1"`. If a scenario sets `penalty_mode:
  "hard"`, drop the sentinel and pass an explicit numeric L1 cost (or
  omit -- hard mode doesn't use L1).

- **Calibration components are configurable.** The default subset is
  `deviation_penalty.calibration.components: [land, feed]`, reproducing
  the historical 2D calibration. A separate config with
  `components: [land, feed, diet]` adds diet as a third axis (used for
  specific investigations where the priced optimum reshuffles the diet
  while leaving land use approximately unchanged). The calibrator
  iterates Broyden over the listed components; the output yaml only
  contains keys for the calibrated subset, and the sentinel resolver
  raises if a scenario uses `"calibrated"` on an uncalibrated component.

- **`slack_marginal_cost: 7.5` in `cost.yaml` caps slack-driven duals.**
  Touch with care: raising it masks data errors; lowering it makes the
  high-mismatch foods (buckwheat, plantain, coffee, tea, olive-oil)
  diverge in cost calibration. See `docs/costs.rst`.

- **Step 1 of cost calibration uses hard stability at +/-5 % with
  `enable_slack: true`** -- not a plain validation solve. Step 2
  tightens to +/-1 % and the dual on each tight constraint becomes the
  per-group median cost correction.

- **Warm-starting `stability` is fine and encouraged.** The previously
  calibrated yaml is auto-detected as the Broyden seed. But if upstream
  changed materially, the warm start can mislead the Jacobian. If
  Broyden takes >5 iterations or oscillates, delete (or rename)
  `deviation_penalty.yaml` and start cold with the seeds from
  `config/calibration/stability.yaml`.

- **The canonical `enabled: false, generate: true` pattern.** Every
  calibration block has both flags. Generation always uses
  `enabled: false` so the solve doesn't read the file it is about to
  write. The combination `enabled: true, generate: true` is rejected by
  `workflow/validation/calibration.py`.

## Quick reference

```bash
# Probe staleness without solving
tools/calibrate --check

# Full chain (gurobi env auto-injected)
tools/calibrate

# Single step
tools/calibrate stability

# Diagnose mismatch before re-calibrating
tools/smk --configfile config/validation.yaml -- \
    results/validation/solved/model_scen-default.nc

# Current calibrated L1 centre
cat data/curated/calibration/default/deviation_penalty.yaml

# Per-iter Broyden trace (after a stability run)
cat results/calibration-default-stability/calibration/deviation_penalty_trace.csv
```
