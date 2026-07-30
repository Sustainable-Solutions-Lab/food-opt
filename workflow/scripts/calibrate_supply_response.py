# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""PMP phase 1: fit the supply-response curve intercepts.

Runs one in-process solve with every supply-curve group pinned to its observed
activity (``supply_response.pin_baseline``, set by the calibration config) and
writes the pinning duals -- the per-group price wedges -- to the artefact the
``supply_response.intercepts: "calibrated"`` sentinel resolves to.
"""

import logging
from pathlib import Path

from workflow.scripts.logging_config import setup_script_logging
from workflow.scripts.solve_model.core import _ShadowPriceLogFilter, run_solve
from workflow.scripts.solve_namespace import (
    build_namespace,
    build_scenario_entry,
    default_path_roots,
)

logger = logging.getLogger(__name__)


def main() -> None:
    smk = snakemake  # type: ignore[name-defined]
    name = smk.params.name
    base_config = dict(smk.config)
    base_config.pop("scenarios", None)
    path_roots = default_path_roots(base_config)

    # The pin regime comes from the calibration config itself; the scenario
    # carries no overrides.
    entry = build_scenario_entry(
        base_config, "_sr_pin", name, path_roots, False, {"_sr_pin": {}}
    )
    Path(entry["log"]).parent.mkdir(parents=True, exist_ok=True)
    n = run_solve(
        build_namespace(entry),
        logger,
        skip_post_processing=True,
        skip_assign_duals=True,
    )
    if n is None:
        raise RuntimeError("pinned supply-response calibration solve failed")

    table = n._supply_response_intercepts
    out_path = Path(smk.output.calibrated_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # 6 significant digits: far below any economically meaningful wedge
    # difference, and keeps the git-tracked artefact compact.
    table.to_csv(out_path, index=False, float_format="%.6g")

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
    logger.info("Wrote %d supply-response intercepts to %s", len(table), out_path)


if __name__ == "__main__":
    logger = setup_script_logging(
        log_file=snakemake.log[0] if snakemake.log else None  # type: ignore[name-defined]
    )
    logging.getLogger("pypsa.optimization.optimize").addFilter(_ShadowPriceLogFilter())
    main()
