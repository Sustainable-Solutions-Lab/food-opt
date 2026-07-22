# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Analyze a fixed-price sweep over endpoint-normalized reallocation caps."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

from workflow.scripts.analyze_water_mechanisms import (
    DEMAND_FACTORS,
    _activity_relief,
    _demand_relief,
    _load_state,
    _plot_attribution,
    _plot_detailed_attribution,
    _shapley_relief,
    _source_relief,
)
from workflow.scripts.plot_crop_water_mosaic import (
    _extract_crop_water,
    _plot_intensity_stack,
)


def _cap_from_path(path: Path) -> float:
    scenario = path.stem.removeprefix("model_scen-")
    if not scenario.startswith("wrw"):
        raise ValueError(f"Unexpected reallocation-cap scenario: {scenario}")
    return int(scenario.removeprefix("wrw")) / 10000


def _states(solved_dir: Path, nonrenewable_cf: float) -> list[dict[str, object]]:
    baseline_path = solved_dir / "model_scen-wf00000.nc"
    endpoint_path = solved_dir / "model_scen-wfw10000.nc"
    if not baseline_path.exists() or not endpoint_path.exists():
        raise ValueError("Reallocation analysis requires fresh wf00000 and wfw10000")
    baseline = _load_state(baseline_path, nonrenewable_cf, "wf")
    baseline["price"] = 0.0
    interior = []
    for path in sorted(solved_dir.glob("model_scen-wrw*.nc"), key=_cap_from_path):
        state = _load_state(path, nonrenewable_cf, "wr")
        state["price"] = _cap_from_path(path)
        interior.append(state)
    endpoint = _load_state(endpoint_path, nonrenewable_cf, "wf")
    endpoint["price"] = 1.0
    return [baseline, *interior, endpoint]


def _analyze(states: list[dict[str, object]], output_dir: Path) -> None:
    baseline = states[0]
    baseline_total = float(baseline["burden"]["burden_mm3_eq"].sum())
    baseline_pool = float(np.prod(list(baseline["demand_factors"].values())))
    summary_rows = []
    attribution_rows = []
    demand_rows = []
    source_rows = []
    activity_rows = []
    for state in states:
        cap = float(state["price"])
        burden = state["burden"]
        total = float(burden["burden_mm3_eq"].sum())
        factors = state["demand_factors"]
        land = factors[DEMAND_FACTORS[0]]
        irrigated_share = factors[DEMAND_FACTORS[1]]
        summary_rows.append(
            {
                "reallocation_cap": cap,
                "surface_scarcity_mm3_eq": float(
                    burden["surface_scarcity_mm3_eq"].sum()
                ),
                "renewable_gw_scarcity_mm3_eq": float(
                    burden["renewable_gw_scarcity_mm3_eq"].sum()
                ),
                "scarcity_mm3_eq": float(burden["scarcity_mm3_eq"].sum()),
                "depletion_mm3": float(burden["depletion_mm3"].sum()),
                "priced_burden_mm3_eq": total,
                "relief_mm3_eq": baseline_total - total,
                "pool_demand_mm3": float(state["q"].to_numpy().sum()),
                "demand_balance_error_mm3": state["demand_mismatch_mm3"],
                "crop_production_land_mha": land,
                "irrigated_share": irrigated_share,
                "irrigated_area_mha": land * irrigated_share,
                "nominal_field_requirement_m3_per_ha": factors[DEMAND_FACTORS[2]],
                "pool_consumption_per_field_requirement": factors[DEMAND_FACTORS[3]],
            }
        )
        headline = _shapley_relief(baseline, state)
        demand = _demand_relief(baseline, state)
        source = _source_relief(baseline, state)
        physical_relief = baseline_pool - float(np.prod(list(factors.values())))
        if not np.isclose(sum(headline.values()), baseline_total - total, atol=1e-4):
            raise ValueError(f"Headline attribution does not close at cap {cap}")
        if not np.isclose(sum(demand.values()), physical_relief, atol=1e-4):
            raise ValueError(f"Demand attribution does not close at cap {cap}")
        if not np.isclose(
            sum(source.values()),
            headline["within-region scarcity-burden intensity"],
            atol=1e-4,
        ):
            raise ValueError(f"Supply attribution does not close at cap {cap}")
        attribution_rows.append({"price": cap, **headline})
        demand_rows.append({"price": cap, **demand})
        source_rows.append({"price": cap, **source})
        activity = _activity_relief(baseline, state)
        activity.insert(0, "reallocation_cap", cap)
        activity_rows.append(activity)

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary_rows).to_csv(output_dir / "summary.csv", index=False)
    attribution = pd.DataFrame(attribution_rows).sort_values("price")
    demand = pd.DataFrame(demand_rows).sort_values("price")
    source = pd.DataFrame(source_rows).sort_values("price")
    attribution.rename(columns={"price": "reallocation_cap"}).to_csv(
        output_dir / "attribution.csv", index=False
    )
    demand.rename(columns={"price": "reallocation_cap"}).to_csv(
        output_dir / "demand_attribution.csv", index=False
    )
    source.rename(columns={"price": "reallocation_cap"}).to_csv(
        output_dir / "source_attribution.csv", index=False
    )
    pd.concat(activity_rows, ignore_index=True).to_csv(
        output_dir / "activity_attribution.csv", index=False
    )
    x_label = "Allowed reallocation (fraction of free price-1 response)"
    _plot_attribution(
        attribution,
        output_dir / "attribution.png",
        x_label=x_label,
        linthresh=0.01,
        figure_note="Water scarcity price fixed at 1 USD/m3-eq",
    )
    _plot_detailed_attribution(
        attribution,
        demand,
        source,
        output_dir / "attribution_detailed.png",
        x_label=x_label,
        linthresh=0.01,
        figure_note="Water scarcity price fixed at 1 USD/m3-eq",
    )


def _crop_figure(
    solved_dir: Path, states: list[dict[str, object]], output_dir: Path
) -> None:
    paths = [
        solved_dir / "model_scen-wf00000.nc",
        *sorted(solved_dir.glob("model_scen-wrw*.nc"), key=_cap_from_path),
        solved_dir / "model_scen-wfw10000.nc",
    ]
    caps = [float(state["price"]) for state in states]
    rows = []
    for path, cap in zip(paths, caps, strict=True):
        frame = _extract_crop_water(pypsa.Network(path))
        frame.insert(0, "price", cap)
        rows.append(frame)
    data = pd.concat(rows, ignore_index=True)
    data.rename(columns={"price": "reallocation_cap"}).to_csv(
        output_dir / "crop_water_by_reallocation_cap.csv", index=False
    )
    _plot_intensity_stack(
        data,
        output_dir / "crop_water_intensity_stack.png",
        x_label="Allowed reallocation (fraction of free price-1 response)",
        linthresh=0.01,
        ordering_label="cap-0 intensity",
        title="Irrigation by crop system and intensity (price fixed at 1 USD/m3-eq)",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--solved-dir", type=Path, default=Path("results/water_barriers_gsa/solved")
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/water_barriers_gsa/water_reallocation_cf100"),
    )
    parser.add_argument("--nonrenewable-cf", type=float, default=100.0)
    args = parser.parse_args()
    states = _states(args.solved_dir, args.nonrenewable_cf)
    _analyze(states, args.output_dir)
    _crop_figure(args.solved_dir, states, args.output_dir)


if __name__ == "__main__":
    main()
