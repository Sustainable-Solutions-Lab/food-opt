# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Decompose priced water-burden relief across a water-price sweep.

The script reads solved networks from
``config/water_mechanisms_fixed_diet_cf100.yaml``.
It reports the burden that the water price actually sees: renewable AWARE
scarcity plus the configured characterisation factor times non-renewable
groundwater depletion.

Water sources mix on shared regional buses, so crop burdens are allocated in
proportion to each production link's field-water demand. Region totals remain
exact; crop attribution is an explicit proportional-sharing convention.

The headline decomposition is a four-factor Shapley identity:

    burden = total irrigation demand
             * activity mix
             * geography within activity
             * within-region scarcity-burden intensity

Here an activity is a crop for single-crop links and a crop combination for
multi-cropping links. The four contributions sum exactly to baseline burden
minus scenario burden, including interactions and negative contributions. A
second Shapley identity explains the physical irrigation-volume reduction from
land footprint, irrigated share, nominal field-water requirement, and the
fixed regional delivery calibration. A linear source split partitions the
within-region scarcity-burden contribution by source family.
"""

import argparse
import itertools
from pathlib import Path

from matplotlib.colors import to_rgb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa

NONRENEWABLE = "groundwater_nonrenewable"
RENEWABLE_GW = "groundwater_renewable"
FACTORS = (
    "lower irrigation volume",
    "irrigation activity mix",
    "regional allocation",
    "within-region scarcity-burden intensity",
)
DEMAND_FACTORS = (
    "crop-production land footprint",
    "irrigated share",
    "nominal field-water requirement per irrigated hectare",
    "calibrated pool consumption per unit field requirement",
)
SOURCE_FACTORS = (
    "surface-water AWARE burden",
    "renewable-groundwater AWARE burden",
    "fossil-groundwater depletion burden",
)
HEADLINE_COLORS = ("#4477AA", "#66CCEE", "#228833", "#CCBB44")
DEMAND_COLORS = ("#B8CEE2", "#7FA6CA", "#4477AA", "#294866")
SOURCE_COLORS = ("#F2E8A5", "#CCBB44", "#887A21")
DIRECT_LABELS = {
    "lower irrigation volume": "Irrigation volume",
    "irrigation activity mix": "Activity mix",
    "regional allocation": "Regional allocation",
    "within-region scarcity-burden intensity": (
        "Within-region scarcity-burden intensity"
    ),
    "crop-production land footprint": "Land footprint",
    "irrigated share": "Irrigated share",
    "nominal field-water requirement per irrigated hectare": (
        "Field-water requirement"
    ),
    "calibrated pool consumption per unit field requirement": ("Delivery calibration"),
    "surface-water AWARE burden": "Surface-water AWARE burden",
    "renewable-groundwater AWARE burden": "Renewable-GW AWARE burden",
    "fossil-groundwater depletion burden": "Fossil-GW depletion burden",
}


def _scenario_price(path: Path, prefix: str = "wf") -> float:
    scenario = path.stem.removeprefix("model_scen-")
    if scenario == f"{prefix}00000":
        return 0.0
    priced_prefix = f"{prefix}w"
    if not scenario.startswith(priced_prefix):
        raise ValueError(f"Unexpected water-mechanism scenario: {scenario}")
    return int(scenario.removeprefix(priced_prefix)) / 10000


def _port_number(bus_column: str) -> int:
    return int(bus_column.removeprefix("bus"))


def _efficiency_column(port: int) -> str:
    return "efficiency" if port == 1 else f"efficiency{port}"


def _crop_pool_demand(
    n: pypsa.Network,
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Return pool demand by activity-region and its physical scalar factors."""
    links = n.links.static
    crops = links[links["carrier"].isin(["crop_production", "crop_production_multi"])]
    dispatch = n.links.dynamic.p0.iloc[0].reindex(crops.index).fillna(0.0).clip(0.0)
    bus_carrier = n.buses.static["carrier"]

    field_et = pd.Series(0.0, index=crops.index)
    bus_columns = sorted(
        (column for column in links.columns if column.startswith("bus")),
        key=_port_number,
    )
    for bus_column in bus_columns:
        port = _port_number(bus_column)
        if port == 0:
            continue
        efficiency_column = _efficiency_column(port)
        if efficiency_column not in crops:
            continue
        is_field_water = crops[bus_column].map(bus_carrier).eq("water_field")
        field_et = field_et.add(
            (-dispatch * crops[efficiency_column].fillna(0.0)).where(
                is_field_water, 0.0
            ),
            fill_value=0.0,
        )

    delivery = links[links["carrier"] == "irrigation_delivery"]
    eta = delivery.groupby("region")["efficiency"].first()
    demand = pd.DataFrame(
        {
            "activity": crops["crop"].astype(str),
            "region": crops["region"].astype(str),
            "field_et_mm3": field_et,
            "pool_demand_mm3": field_et / crops["region"].map(eta),
        },
        index=crops.index,
    )
    positive = demand[demand["pool_demand_mm3"] > 1e-9]
    grouped = positive.groupby(["activity", "region"], as_index=False)[
        ["field_et_mm3", "pool_demand_mm3"]
    ].sum()

    cropland_area = float(dispatch.sum())
    irrigated_area = float(dispatch[crops["water_supply"].eq("irrigated")].sum())
    field_et_total = float(field_et.sum())
    pool_total = float(demand["pool_demand_mm3"].sum())
    if min(cropland_area, irrigated_area, field_et_total, pool_total) <= 0:
        raise ValueError("Cannot decompose irrigation demand with a zero factor")
    factors = {
        "crop-production land footprint": cropland_area,
        "irrigated share": irrigated_area / cropland_area,
        "nominal field-water requirement per irrigated hectare": (
            field_et_total / irrigated_area
        ),
        "calibrated pool consumption per unit field requirement": (
            pool_total / field_et_total
        ),
    }
    return grouped, factors


def _regional_burden(n: pypsa.Network, nonrenewable_cf: float) -> pd.DataFrame:
    """Return exact tier draw and priced burden by region."""
    links = n.links.static
    tiers = links[links["carrier"] == "water_supply"]
    draw = n.links.dynamic.p0.iloc[0].reindex(tiers.index).clip(lower=0.0)
    source = tiers["source"].fillna("renewable")
    nonrenewable = source.eq(NONRENEWABLE)
    renewable_gw = source.eq(RENEWABLE_GW)
    surface = ~nonrenewable & ~renewable_gw
    scarcity = draw * tiers["efficiency2"].fillna(0.0)
    frame = pd.DataFrame(
        {
            "region": tiers["region"].astype(str),
            "draw_mm3": draw,
            "surface_scarcity_mm3_eq": scarcity.where(surface, 0.0),
            "renewable_gw_scarcity_mm3_eq": scarcity.where(renewable_gw, 0.0),
            "depletion_mm3": draw.where(nonrenewable, 0.0),
        },
        index=tiers.index,
    )
    out = frame.groupby("region", as_index=False).sum(numeric_only=True)
    out["scarcity_mm3_eq"] = (
        out["surface_scarcity_mm3_eq"] + out["renewable_gw_scarcity_mm3_eq"]
    )
    out["fossil_burden_mm3_eq"] = nonrenewable_cf * out["depletion_mm3"]
    out["burden_mm3_eq"] = out["scarcity_mm3_eq"] + out["fossil_burden_mm3_eq"]
    return out


def _load_state(
    path: Path, nonrenewable_cf: float, scenario_prefix: str
) -> dict[str, object]:
    n = pypsa.Network(path)
    demand, demand_factors = _crop_pool_demand(n)
    burden = _regional_burden(n, nonrenewable_cf)
    q = demand.pivot(index="activity", columns="region", values="pool_demand_mm3")
    q = q.fillna(0.0)
    region_burden = burden.set_index("region")["burden_mm3_eq"]
    q_region = q.sum(axis=0)
    intensity = region_burden.reindex(q.columns, fill_value=0.0).div(q_region)
    intensity = intensity.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    component_columns = {
        "surface-water AWARE burden": "surface_scarcity_mm3_eq",
        "renewable-groundwater AWARE burden": "renewable_gw_scarcity_mm3_eq",
        "fossil-groundwater depletion burden": "fossil_burden_mm3_eq",
    }
    regional = burden.set_index("region")
    intensity_components = {
        name: regional[column].reindex(q.columns, fill_value=0.0).div(q_region)
        for name, column in component_columns.items()
    }
    intensity_components = {
        name: values.replace([np.inf, -np.inf], np.nan).fillna(0.0)
        for name, values in intensity_components.items()
    }

    allocated = q.mul(intensity, axis=1)
    mismatch = float(q_region.sum() - burden["draw_mm3"].sum())
    return {
        "price": _scenario_price(path, scenario_prefix),
        "q": q,
        "intensity": intensity,
        "intensity_components": intensity_components,
        "demand_factors": demand_factors,
        "allocated": allocated,
        "burden": burden,
        "demand_mismatch_mm3": mismatch,
    }


def _align_states(
    baseline: dict[str, object], scenario: dict[str, object]
) -> tuple[list[np.ndarray], list[np.ndarray], pd.Index, pd.Index]:
    q0 = baseline["q"]
    q1 = scenario["q"]
    activities = q0.index.union(q1.index)
    regions = q0.columns.union(q1.columns)
    q0 = q0.reindex(index=activities, columns=regions, fill_value=0.0)
    q1 = q1.reindex(index=activities, columns=regions, fill_value=0.0)
    k0 = baseline["intensity"].reindex(regions, fill_value=0.0)
    k1 = scenario["intensity"].reindex(regions, fill_value=0.0)

    total0 = float(q0.to_numpy().sum())
    total1 = float(q1.to_numpy().sum())
    activity0 = q0.sum(axis=1) / total0
    activity1 = q1.sum(axis=1) / total1

    location0 = q0.div(q0.sum(axis=1).replace(0.0, np.nan), axis=0)
    location1 = q1.div(q1.sum(axis=1).replace(0.0, np.nan), axis=0)
    absent0 = location0.isna().all(axis=1)
    absent1 = location1.isna().all(axis=1)
    location0.loc[absent0] = location1.loc[absent0]
    location1.loc[absent1] = location0.loc[absent1]
    location0 = location0.fillna(0.0)
    location1 = location1.fillna(0.0)

    no_draw0 = q0.sum(axis=0).eq(0.0)
    no_draw1 = q1.sum(axis=0).eq(0.0)
    k0.loc[no_draw0] = k1.loc[no_draw0]
    k1.loc[no_draw1] = k0.loc[no_draw1]

    start = [total0, activity0.to_numpy(), location0.to_numpy(), k0.to_numpy()]
    end = [total1, activity1.to_numpy(), location1.to_numpy(), k1.to_numpy()]
    return start, end, activities, regions


def _identity(factors: list[np.ndarray | float]) -> float:
    total, activity, location, intensity = factors
    return float(total * np.sum(activity * (location @ intensity)))


def _shapley_relief(
    baseline: dict[str, object], scenario: dict[str, object]
) -> dict[str, float]:
    start, end, _, _ = _align_states(baseline, scenario)
    contributions = np.zeros(len(FACTORS))
    permutations = list(itertools.permutations(range(len(FACTORS))))
    for order in permutations:
        state = [value.copy() if hasattr(value, "copy") else value for value in start]
        before = _identity(state)
        for factor in order:
            state[factor] = end[factor]
            after = _identity(state)
            contributions[factor] += before - after
            before = after
    contributions /= len(permutations)
    return dict(zip(FACTORS, contributions, strict=True))


def _demand_relief(
    baseline: dict[str, object], scenario: dict[str, object]
) -> dict[str, float]:
    """Four-factor Shapley decomposition of physical pool-demand relief."""
    start = [baseline["demand_factors"][name] for name in DEMAND_FACTORS]
    end = [scenario["demand_factors"][name] for name in DEMAND_FACTORS]
    contributions = np.zeros(len(DEMAND_FACTORS))
    permutations = list(itertools.permutations(range(len(DEMAND_FACTORS))))
    for order in permutations:
        state = start.copy()
        before = float(np.prod(state))
        for factor in order:
            state[factor] = end[factor]
            after = float(np.prod(state))
            contributions[factor] += before - after
            before = after
    contributions /= len(permutations)
    return dict(zip(DEMAND_FACTORS, contributions, strict=True))


def _source_relief(
    baseline: dict[str, object], scenario: dict[str, object]
) -> dict[str, float]:
    """Split the high-level source-response Shapley value by water source."""
    start, end, _, regions = _align_states(baseline, scenario)
    q0 = baseline["q"].reindex(columns=regions, fill_value=0.0).sum(axis=0)
    q1 = scenario["q"].reindex(columns=regions, fill_value=0.0).sum(axis=0)
    no_draw0 = q0.eq(0.0)
    no_draw1 = q1.eq(0.0)
    components = {}
    for name in SOURCE_FACTORS:
        component0 = baseline["intensity_components"][name].reindex(
            regions, fill_value=0.0
        )
        component1 = scenario["intensity_components"][name].reindex(
            regions, fill_value=0.0
        )
        component0.loc[no_draw0] = component1.loc[no_draw0]
        component1.loc[no_draw1] = component0.loc[no_draw1]
        components[name] = (component0.to_numpy(), component1.to_numpy())

    contributions = dict.fromkeys(SOURCE_FACTORS, 0.0)
    permutations = list(itertools.permutations(range(len(FACTORS))))
    source_factor = FACTORS.index("within-region scarcity-burden intensity")
    for order in permutations:
        state = [value.copy() if hasattr(value, "copy") else value for value in start]
        for factor in order:
            if factor == source_factor:
                total, activity, location, _ = state
                for name, (component0, component1) in components.items():
                    contributions[name] += float(
                        total
                        * np.sum(activity * (location @ (component0 - component1)))
                    )
            state[factor] = end[factor]
    return {name: value / len(permutations) for name, value in contributions.items()}


def _activity_relief(
    baseline: dict[str, object], scenario: dict[str, object]
) -> pd.DataFrame:
    start, end, activities, _ = _align_states(baseline, scenario)
    _, _, _, k0 = start
    _, _, _, k1 = end
    q0 = baseline["q"].reindex(index=activities, fill_value=0.0)
    q1 = scenario["q"].reindex(index=activities, fill_value=0.0)
    regions = q0.columns.union(q1.columns)
    q0 = q0.reindex(columns=regions, fill_value=0.0).to_numpy()
    q1 = q1.reindex(columns=regions, fill_value=0.0).to_numpy()
    k0 = baseline["intensity"].reindex(regions, fill_value=0.0).to_numpy()
    k1 = scenario["intensity"].reindex(regions, fill_value=0.0).to_numpy()

    volume = np.sum((q0 - q1) * (k0 + k1) / 2, axis=1)
    source = np.sum((k0 - k1) * (q0 + q1) / 2, axis=1)
    return pd.DataFrame(
        {
            "activity": activities,
            "demand_relief_mm3_eq": volume,
            "source_relief_mm3_eq": source,
            "total_relief_mm3_eq": volume + source,
        }
    )


def _plot_attribution(
    attribution: pd.DataFrame,
    output: Path,
    x_label: str = "Water price (USD per m3 scarcity-equivalent)",
    linthresh: float = 0.01,
    figure_note: str | None = None,
) -> None:
    data = attribution.sort_values("price")
    x = data["price"].to_numpy()
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    positive_base = np.zeros(len(data))
    negative_base = np.zeros(len(data))
    for factor, color in zip(FACTORS, HEADLINE_COLORS, strict=True):
        values = data[factor].to_numpy() / 1e6
        positive = np.clip(values, 0.0, None)
        negative = np.clip(values, None, 0.0)
        ax.fill_between(
            x,
            positive_base,
            positive_base + positive,
            label=factor,
            color=color,
            alpha=0.9,
        )
        ax.fill_between(
            x,
            negative_base,
            negative_base + negative,
            color=color,
            alpha=0.9,
        )
        positive_base += positive
        negative_base += negative
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xscale("symlog", linthresh=linthresh)
    ax.set_xlim(0.0, x.max())
    ax.set_xlabel(x_label)
    ax.set_ylabel("Relief from unpriced baseline (10^6 m3-equivalent)")
    if figure_note:
        ax.set_title(figure_note, loc="left", fontsize=10.5, color="#555555")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def _stacked_axis(
    ax,
    data: pd.DataFrame,
    factors: tuple[str, ...],
    colors: tuple[str, ...],
    scale: float,
) -> tuple[float, list[tuple[str, str, float]]]:
    """Draw positive and negative stacked areas on an existing axis."""
    x = data["price"].to_numpy()
    positive_base = np.zeros(len(data))
    negative_base = np.zeros(len(data))
    anchors = []
    for factor, color in zip(factors, colors, strict=True):
        values = data[factor].to_numpy() / scale
        positive = np.clip(values, 0.0, None)
        negative = np.clip(values, None, 0.0)
        endpoint = (
            positive_base[-1] + positive[-1] / 2
            if values[-1] >= 0
            else negative_base[-1] + negative[-1] / 2
        )
        anchors.append((factor, color, endpoint))
        ax.fill_between(
            x,
            positive_base,
            positive_base + positive,
            color=color,
            alpha=0.9,
        )
        ax.fill_between(
            x,
            negative_base,
            negative_base + negative,
            color=color,
            alpha=0.9,
        )
        positive_base += positive
        negative_base += negative
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.grid(alpha=0.12)
    return x[-1], anchors


def _direct_area_labels(ax, x: float, anchors: list[tuple[str, str, float]]) -> None:
    """Label stacked bands at the right edge, spreading close labels apart."""
    ymin, ymax = ax.get_ylim()
    span = ymax - ymin
    margin = 0.04 * span
    gap = 0.065 * span
    order = np.argsort([anchor for _, _, anchor in anchors])
    positions = np.array([anchors[index][2] for index in order], dtype=float)
    positions[0] = max(positions[0], ymin + margin)
    for index in range(1, len(positions)):
        positions[index] = max(positions[index], positions[index - 1] + gap)
    overflow = positions[-1] - (ymax - margin)
    if overflow > 0:
        positions -= overflow
    for position, anchor_index in zip(positions, order, strict=True):
        factor, color, anchor = anchors[anchor_index]
        text_color = tuple(0.58 * channel for channel in to_rgb(color))
        ax.annotate(
            DIRECT_LABELS[factor],
            xy=(x, anchor),
            xycoords="data",
            xytext=(1.025, position),
            textcoords=ax.get_yaxis_transform(),
            ha="left",
            va="center",
            color=text_color,
            fontsize=10,
            fontweight="medium",
            arrowprops={"arrowstyle": "-", "color": color, "linewidth": 1.4},
            annotation_clip=False,
        )


def _plot_detailed_attribution(
    attribution: pd.DataFrame,
    demand: pd.DataFrame,
    source: pd.DataFrame,
    output: Path,
    x_label: str = "Water price (USD per m3 scarcity-equivalent)",
    linthresh: float = 0.01,
    figure_note: str | None = None,
) -> None:
    """Plot high-level relief with nested demand and source decompositions."""
    attribution = attribution.sort_values("price")
    demand = demand.sort_values("price")
    source = source.sort_values("price")
    fig, axes = plt.subplots(3, 1, figsize=(11.5, 10.2), sharex=True)
    high_level_labels = _stacked_axis(
        axes[0],
        attribution,
        FACTORS,
        HEADLINE_COLORS,
        1e6,
    )
    demand_labels = _stacked_axis(
        axes[1],
        demand,
        DEMAND_FACTORS,
        DEMAND_COLORS,
        1e3,
    )
    source_labels = _stacked_axis(
        axes[2],
        source,
        SOURCE_FACTORS,
        SOURCE_COLORS,
        1e6,
    )

    axes[0].set_ylabel("Priced burden relief\n(10^6 m3-equivalent)")
    axes[1].set_ylabel("Irrigation consumption reduction\n(km3)")
    axes[2].set_ylabel("Contribution to priced-burden relief\n(10^6 m3-equivalent)")
    axes[0].set_title("a   High-level burden decomposition", loc="left")
    axes[1].set_title("b   Why irrigation volume falls", loc="left")
    axes[2].set_title("c   Sources of within-region scarcity-burden relief", loc="left")
    axes[2].set_xscale("symlog", linthresh=linthresh)
    axes[2].set_xlim(0.0, attribution["price"].max())
    axes[2].set_xlabel(x_label)
    for ax, (x, labels) in zip(
        axes,
        (high_level_labels, demand_labels, source_labels),
        strict=True,
    ):
        _direct_area_labels(ax, x, labels)
    top = 0.93 if figure_note else 0.97
    if figure_note:
        fig.suptitle(
            figure_note,
            x=0.1,
            y=0.985,
            ha="left",
            fontsize=11,
            color="#555555",
        )
    fig.subplots_adjust(left=0.1, right=0.73, bottom=0.075, top=top, hspace=0.35)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def _diagnostic_changes(
    analysis_dir: Path, prices: dict[str, float], baseline_scenario: str
) -> pd.DataFrame:
    """Return selected physical and economic changes relative to the baseline."""
    specifications = {
        "land water system": ("land_use", ["water_supply"], "area_mha"),
        "land crop": ("land_use", ["crop"], "area_mha"),
        "crop production": ("crop_production", ["crop"], "production_mt"),
        "animal production": ("animal_production", ["product"], "production_mt"),
        "feed category": ("feed_by_category", ["category"], "mt_dm"),
        "feed source": ("feed_by_source", ["source_key"], "mt_dm"),
        "food consumption": ("food_consumption", ["food"], "consumption_mt"),
    }

    def aggregate(scenario: str, filename: str, keys: list[str], value: str):
        path = analysis_dir / f"scen-{scenario}" / f"{filename}.parquet"
        frame = pd.read_parquet(path)
        return frame.groupby(keys, dropna=False)[value].sum()

    rows = []
    for domain, (filename, keys, value) in specifications.items():
        baseline = aggregate(baseline_scenario, filename, keys, value)
        for scenario, price in prices.items():
            current = aggregate(scenario, filename, keys, value)
            joined = pd.concat(
                [baseline.rename("baseline"), current.rename("value")], axis=1
            ).fillna(0.0)
            joined["delta"] = joined["value"] - joined["baseline"]
            joined["pct"] = 100 * joined["delta"].div(
                joined["baseline"].replace(0.0, np.nan)
            )
            flat = joined.reset_index()
            flat["item"] = flat[keys].astype(str).agg(" | ".join, axis=1)
            flat.insert(0, "domain", domain)
            flat.insert(0, "price", price)
            rows.append(flat[["price", "domain", "item", *joined.columns]])
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--solved-dir", type=Path, default=Path("results/water_barriers_gsa/solved")
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/water_barriers_gsa/water_mechanisms_cf100"),
    )
    parser.add_argument("--nonrenewable-cf", type=float, default=100.0)
    parser.add_argument("--scenario-prefix", default="wf")
    parser.add_argument("--max-price", type=float, default=1.0)
    parser.add_argument("--skip-diagnostics", action="store_true")
    args = parser.parse_args()

    pattern = f"model_scen-{args.scenario_prefix}*.nc"
    paths = sorted(
        (
            path
            for path in args.solved_dir.glob(pattern)
            if _scenario_price(path, args.scenario_prefix) <= args.max_price
        ),
        key=lambda path: _scenario_price(path, args.scenario_prefix),
    )
    baseline_scenario = f"{args.scenario_prefix}00000"
    if not paths or _scenario_price(paths[0], args.scenario_prefix) != 0.0:
        raise ValueError(
            f"The water-mechanism sweep requires scenario {baseline_scenario}"
        )

    states = [
        _load_state(path, args.nonrenewable_cf, args.scenario_prefix) for path in paths
    ]
    baseline = states[0]
    summary_rows = []
    attribution_rows = []
    demand_rows = []
    source_rows = []
    activity_rows = []
    baseline_total = float(baseline["burden"]["burden_mm3_eq"].sum())
    baseline_pool_demand = float(np.prod(list(baseline["demand_factors"].values())))
    for state in states:
        burden = state["burden"]
        total = float(burden["burden_mm3_eq"].sum())
        price = float(state["price"])
        demand_factors = state["demand_factors"]
        land_footprint = demand_factors["crop-production land footprint"]
        irrigated_share = demand_factors["irrigated share"]
        field_requirement = demand_factors[
            "nominal field-water requirement per irrigated hectare"
        ]
        delivery_multiplier = demand_factors[
            "calibrated pool consumption per unit field requirement"
        ]
        summary_rows.append(
            {
                "price": price,
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
                "crop_production_land_mha": land_footprint,
                "irrigated_share": irrigated_share,
                "irrigated_area_mha": land_footprint * irrigated_share,
                "nominal_field_requirement_m3_per_ha": field_requirement,
                "pool_consumption_per_field_requirement": delivery_multiplier,
            }
        )
        high_level = _shapley_relief(baseline, state)
        demand_detail = _demand_relief(baseline, state)
        source_detail = _source_relief(baseline, state)
        headline_relief = baseline_total - total
        physical_relief = baseline_pool_demand - float(
            np.prod(list(demand_factors.values()))
        )
        if not np.isclose(
            sum(high_level.values()), headline_relief, rtol=1e-10, atol=1e-4
        ):
            raise ValueError(f"High-level attribution does not close at price {price}")
        if not np.isclose(
            sum(demand_detail.values()), physical_relief, rtol=1e-10, atol=1e-4
        ):
            raise ValueError(f"Demand attribution does not close at price {price}")
        if not np.isclose(
            sum(source_detail.values()),
            high_level["within-region scarcity-burden intensity"],
            rtol=1e-10,
            atol=1e-4,
        ):
            raise ValueError(
                f"Supply-exposure attribution does not close at price {price}"
            )
        attribution_rows.append({"price": price, **high_level})
        demand_rows.append({"price": price, **demand_detail})
        source_rows.append({"price": price, **source_detail})
        activity = _activity_relief(baseline, state)
        activity.insert(0, "price", price)
        activity_rows.append(activity)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(summary_rows).sort_values("price")
    attribution = pd.DataFrame(attribution_rows).sort_values("price")
    demand_attribution = pd.DataFrame(demand_rows).sort_values("price")
    source_attribution = pd.DataFrame(source_rows).sort_values("price")
    activity = pd.concat(activity_rows, ignore_index=True)
    summary.to_csv(args.output_dir / "summary.csv", index=False)
    attribution.to_csv(args.output_dir / "attribution.csv", index=False)
    demand_attribution.to_csv(args.output_dir / "demand_attribution.csv", index=False)
    source_attribution.to_csv(args.output_dir / "source_attribution.csv", index=False)
    activity.to_csv(args.output_dir / "activity_attribution.csv", index=False)
    if not args.skip_diagnostics:
        scenario_prices = {
            path.stem.removeprefix("model_scen-"): _scenario_price(
                path, args.scenario_prefix
            )
            for path in paths
        }
        diagnostics = _diagnostic_changes(
            args.solved_dir.parent / "analysis", scenario_prices, baseline_scenario
        )
        diagnostics.to_csv(args.output_dir / "diagnostic_changes.csv", index=False)
    _plot_attribution(attribution, args.output_dir / "attribution.png")
    _plot_detailed_attribution(
        attribution,
        demand_attribution,
        source_attribution,
        args.output_dir / "attribution_detailed.png",
    )


if __name__ == "__main__":
    main()
