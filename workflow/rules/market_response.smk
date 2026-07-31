# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""PMP calibration of the market-response curve intercepts.

A sequence of pinned solves run in-process by the script (production pin,
demand slope basis, demand intercepts, then refinement sweeps); the duals of
the pinning constraints are the per-group price wedges written to the
artefact.
"""


_sr_cal_cfg = config["market_response"]["calibration"]

if _sr_cal_cfg["generate"]:

    rule calibrate_market_response:
        input:
            # The pinned solve reads the same calibration artefacts an ordinary
            # solve does (its own output is gated off by generate mode), plus
            # health inputs only when health is enabled.
            unpack(
                lambda w: {
                    **(health_input_paths(name) if health_required() else {}),
                    **calibration_artefact_inputs(config),
                }
            ),
            model=f"<results>/{name}/build/model.nc",
            baseline_diet=f"<processing>/{name}/baseline_diet.csv",
            m49="data/curated/M49-codes.csv",
            food_groups="data/curated/food_groups.csv",
            nutrition="data/curated/nutrition.csv",
            solve_scripts=solve_runtime_code_paths(),
        output:
            calibrated_csv=_sr_cal_cfg["calibrated_csv"],
        params:
            name=name,
        resources:
            runtime=lambda w, attempt: 60 * (1 + attempt),
            mem_mb=lambda w, attempt: 12000 * attempt,
        log:
            f"<logs>/{name}/calibrate_market_response.log",
        benchmark:
            f"<benchmarks>/{name}/calibrate_market_response.tsv"
        script:
            "../scripts/calibrate_market_response.py"
