# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Solve the model and run analysis in a single process.

Equivalent to running solve_model.py followed by analyze_model.py, but
skips the intermediate .nc write/read cycle.  The in-memory solved
network is passed directly to the analysis extraction functions.

Enabled via ``solving.inline_analysis: true`` in config.
"""

import logging

from workflow.scripts.analysis.analyze_model import (
    run_analysis_from_namespace,
    write_empty_outputs,
)
from workflow.scripts.logging_config import setup_script_logging
from workflow.scripts.solve_model.core import _ShadowPriceLogFilter, run_solve


def main() -> None:
    logger = setup_script_logging(snakemake.log[0])
    logging.getLogger("pypsa.optimization.optimize").addFilter(_ShadowPriceLogFilter())

    # Phase 1: Solve
    n = run_solve(snakemake, logger)

    if n is None:
        logger.warning("Solve failed — writing empty analysis outputs.")
        write_empty_outputs(snakemake.output)
        return

    run_analysis_from_namespace(snakemake, n, logger)

    logger.info("Solve-and-analyze complete.")


if __name__ == "__main__":
    main()
