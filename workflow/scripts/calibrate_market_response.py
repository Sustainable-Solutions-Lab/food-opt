# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""PMP phase 1: fit the market-response curve intercepts.

Sequential procedure (see the module docstring of
``workflow.scripts.solve_model.market_response`` for the reasoning):

1. **Production pin** (stage A): every production group pinned to its observed
   activity, demand left to the base config's own regime. The pinning duals
   are the production price wedges.
2. **Slope basis** (stage B0, only with the demand component enabled): the
   demand pin against *zero-intercept* supply curves. Food-bus prices then sit
   at the accounting-cost chain, and the demand duals are realistic reference
   price levels -- the basis for the marginal-utility slopes.
3. **Demand intercepts** (stage B1): the demand pin against the *calibrated*
   supply curves from stage A. The duals are the demand wedges consistent with
   exactly the supply side deployed solves carry.
4. **Refinement sweeps** (``market_response.calibration.sweeps``): the
   production pin repeats with the fitted demand curves active -- stage A
   measured the production wedges under the base regime's demand, while
   deployed solves carry the demand curves -- and the demand intercepts are
   refit against the refined production curves. Each sweep contracts the
   residual cross-side drift by roughly 4x.

Joint pinning of both sides is deliberately avoided: it closes every commodity
chain, leaving the producer/consumer split of each chain's wedge undetermined
(the duals park at the pin slack bound).

The artefact concatenates the production wedges (stage A') with the demand
rows (intercept from B1', ``slope_basis`` from B0), written to the path the
``market_response.intercepts: "calibrated"`` sentinel resolves to.
"""

import logging
from pathlib import Path

import pandas as pd

from workflow.scripts.logging_config import setup_script_logging
from workflow.scripts.solve_model.core import _ShadowPriceLogFilter, run_solve
from workflow.scripts.solve_model.market_response import (
    COMPONENTS,
    DEMAND_COMPONENTS,
)
from workflow.scripts.solve_namespace import (
    build_namespace,
    build_scenario_entry,
    default_path_roots,
)

logger = logging.getLogger(__name__)

PRODUCTION_COMPONENTS = sorted(set(COMPONENTS) - DEMAND_COMPONENTS)


def _pinned_solve(
    base_config: dict, name: str, path_roots: dict, scenario: str, overrides: dict
) -> pd.DataFrame:
    """Run one pinned solve and return its extracted intercept table."""
    entry = build_scenario_entry(
        base_config, scenario, name, path_roots, False, {scenario: overrides}
    )
    Path(entry["log"]).parent.mkdir(parents=True, exist_ok=True)
    n = run_solve(
        build_namespace(entry),
        logger,
        skip_post_processing=True,
        skip_assign_duals=True,
    )
    if n is None:
        raise RuntimeError(f"pinned market-response solve '{scenario}' failed")
    return n._market_response_intercepts


def _log_component_stats(table: pd.DataFrame) -> None:
    for component, sub in table.groupby("component"):
        logger.info(
            "%s: %d groups, wedge median %.4g (p5 %.4g, p95 %.4g), "
            "%d slacked (|slack| sum %.4g)",
            component,
            len(sub),
            sub["intercept"].median(),
            sub["intercept"].quantile(0.05),
            sub["intercept"].quantile(0.95),
            int((sub["slack"].abs() > 1e-6).sum()),
            sub["slack"].abs().sum(),
        )


def main() -> None:
    smk = snakemake  # type: ignore[name-defined]
    name = smk.params.name
    base_config = dict(smk.config)
    base_config.pop("scenarios", None)
    path_roots = default_path_roots(base_config)
    demand_enabled = bool(base_config["market_response"]["components"]["demand"])

    # Stage A: pin production, demand in the base config's own regime.
    production = _pinned_solve(
        base_config,
        name,
        path_roots,
        "_sr_pin",
        {
            "market_response": {
                "pin_baseline": True,
                "intercepts": None,
                "components": {"demand": False},
            }
        },
    )
    table = production

    if demand_enabled:
        stage_dir = Path(path_roots["results"]) / name / "calibration"
        stage_dir.mkdir(parents=True, exist_ok=True)

        # Stage B0: demand pin against zero-intercept supply curves; the duals
        # are the accounting-cost-chain reference prices (slope basis).
        basis = _pinned_solve(
            base_config,
            name,
            path_roots,
            "_sr_basis",
            {"market_response": {"pin_baseline": ["demand"], "intercepts": None}},
        )
        basis = (
            basis[basis["component"] == "demand"]
            .set_index("mr_group")["intercept"]
            .abs()
        )

        def _fit_demand(stage: str, production_table: pd.DataFrame) -> pd.DataFrame:
            """Demand pin against the given production curves (stages B1/B1')."""
            csv = stage_dir / f"market_response_{stage}.csv"
            production_table.to_csv(csv, index=False, float_format="%.6g")
            fitted = _pinned_solve(
                base_config,
                name,
                path_roots,
                f"_sr_{stage}",
                {
                    "market_response": {
                        "pin_baseline": ["demand"],
                        "intercepts": str(csv),
                    }
                },
            )
            fitted = fitted[fitted["component"] == "demand"].copy()
            fitted["slope_basis"] = basis.reindex(fitted["mr_group"]).to_numpy()
            if fitted["slope_basis"].isna().any():
                raise RuntimeError(
                    "slope-basis solve is missing demand groups present in the "
                    "intercept solve; the stages must see the same model"
                )
            return fitted

        # Stage B1: demand intercepts against the stage-A supply curves.
        demand = _fit_demand("b1", production)

        # Refinement sweeps: production pinned again with the fitted demand
        # curves active (measuring its wedges under the deployed demand system
        # instead of the base regime's), then the demand intercepts refit
        # against the refined production curves. Each sweep contracts the
        # residual cross-side drift by roughly 4x.
        sweeps = int(base_config["market_response"]["calibration"]["sweeps"])
        for sweep in range(1, sweeps):
            interim_csv = stage_dir / f"market_response_interim_{sweep}.csv"
            pd.concat([production, demand], ignore_index=True).to_csv(
                interim_csv, index=False, float_format="%.6g"
            )
            production = _pinned_solve(
                base_config,
                name,
                path_roots,
                f"_sr_refine_{sweep}",
                {
                    "market_response": {
                        "pin_baseline": PRODUCTION_COMPONENTS,
                        "intercepts": str(interim_csv),
                    }
                },
            )
            demand = _fit_demand(f"b1_{sweep}", production)
        table = pd.concat([production, demand], ignore_index=True)

    out_path = Path(smk.output.calibrated_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # 6 significant digits: far below any economically meaningful wedge
    # difference, and keeps the git-tracked artefact compact.
    table.to_csv(out_path, index=False, float_format="%.6g")

    _log_component_stats(table)
    logger.info("Wrote %d market-response intercepts to %s", len(table), out_path)


if __name__ == "__main__":
    logger = setup_script_logging(
        log_file=snakemake.log[0] if snakemake.log else None  # type: ignore[name-defined]
    )
    logging.getLogger("pypsa.optimization.optimize").addFilter(_ShadowPriceLogFilter())
    main()
