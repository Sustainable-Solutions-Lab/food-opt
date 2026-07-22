# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Plot production, irrigation intensity, and water use by crop system.

Single crops are shown by crop type. Multi-cropping links are kept as named
crop systems because a solved link has one combined water input for all crop
cycles. Keeping the bundle intact avoids an arbitrary constituent-crop split
and preserves exact production and water accounting.
"""

import argparse
from pathlib import Path

from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa


def _scenario_price(path: Path) -> float:
    scenario = path.stem.removeprefix("model_scen-")
    if scenario == "wf00000":
        return 0.0
    if not scenario.startswith("wfw"):
        raise ValueError(f"Unexpected CF100 water-mechanism scenario: {scenario}")
    return int(scenario.removeprefix("wfw")) / 10000


def _port_number(column: str) -> int:
    return int(column.removeprefix("bus"))


def _efficiency_column(port: int) -> str:
    return "efficiency" if port == 1 else f"efficiency{port}"


def _crop_output_flows(n: pypsa.Network, crops: pd.DataFrame) -> pd.DataFrame:
    """Return crop-bus output flows by production link and output port."""
    records = []
    bus_crop = n.buses.static["crop"]
    bus_columns = sorted(
        (column for column in crops.columns if column.startswith("bus")),
        key=_port_number,
    )
    for bus_column in bus_columns:
        port = _port_number(bus_column)
        if port == 0:
            continue
        crop = crops[bus_column].map(bus_crop)
        is_crop = crop.notna() & crop.astype(str).ne("")
        if not is_crop.any():
            continue
        dynamic_column = f"p{port}"
        flow = (
            -n.links.dynamic[dynamic_column]
            .iloc[0]
            .reindex(crops.index, fill_value=0.0)
        )
        selected = is_crop & flow.gt(1e-9)
        records.append(
            pd.DataFrame(
                {
                    "link": crops.index[selected],
                    "crop": crop[selected].astype(str).to_numpy(),
                    "production_mt": flow[selected].to_numpy(),
                }
            )
        )
    if not records:
        raise ValueError("No crop production flows found")
    return pd.concat(records, ignore_index=True)


def _field_requirement(n: pypsa.Network, crops: pd.DataFrame) -> pd.Series:
    """Return annual nominal field-water requirement in m3/ha by link."""
    requirement = pd.Series(0.0, index=crops.index)
    bus_carrier = n.buses.static["carrier"]
    for bus_column in sorted(
        (column for column in crops.columns if column.startswith("bus")),
        key=_port_number,
    ):
        port = _port_number(bus_column)
        if port == 0:
            continue
        efficiency_column = _efficiency_column(port)
        if efficiency_column not in crops:
            continue
        is_field_water = crops[bus_column].map(bus_carrier).eq("water_field")
        requirement = requirement.add(
            (-crops[efficiency_column].fillna(0.0)).where(is_field_water, 0.0),
            fill_value=0.0,
        )
    return requirement.clip(lower=0.0)


def _extract_crop_water(n: pypsa.Network) -> pd.DataFrame:
    """Return production and pool-water consumption by exact crop system."""
    links = n.links.static
    crops = links[links["carrier"].isin(["crop_production", "crop_production_multi"])]
    dispatch = n.links.dynamic.p0.iloc[0].reindex(crops.index).fillna(0.0).clip(0.0)
    field_requirement = _field_requirement(n, crops)
    eta = (
        links[links["carrier"] == "irrigation_delivery"]
        .groupby("region")["efficiency"]
        .first()
    )
    pool_water = dispatch * field_requirement / crops["region"].map(eta)
    outputs = _crop_output_flows(n, crops)
    production = outputs.groupby("link")["production_mt"].sum()
    crop_system = crops["crop"].astype(str).copy()
    is_multi = crops["carrier"].eq("crop_production_multi")
    crop_system.loc[is_multi] = "multi:" + crops.loc[is_multi, "combination"].astype(
        str
    )
    allocation = pd.DataFrame(
        {
            "crop_system": crop_system,
            "production_mt": production.reindex(crops.index, fill_value=0.0),
            "water_mm3": pool_water,
        },
        index=crops.index,
    )
    expected = float(pool_water.sum())
    allocated = float(allocation["water_mm3"].sum())
    if not np.isclose(allocated, expected, rtol=1e-10, atol=1e-5):
        raise ValueError(
            f"Crop water allocation does not close: {allocated} != {expected}"
        )
    return allocation.groupby("crop_system", as_index=False)[
        ["production_mt", "water_mm3"]
    ].sum()


def _group_crops(data: pd.DataFrame, top_n: int) -> tuple[pd.DataFrame, list[str]]:
    baseline = data[data["price"].eq(0.0)].set_index("crop_system")["water_mm3"]
    keep = baseline.nlargest(top_n).index.tolist()
    grouped = data.copy()
    grouped["crop_group"] = grouped["crop_system"].where(
        grouped["crop_system"].isin(keep), "other"
    )
    grouped = grouped.groupby(["price", "crop_group"], as_index=False)[
        ["production_mt", "water_mm3"]
    ].sum()
    order = [*keep, "other"]
    return grouped, order


def _display_name(crop_system: str) -> str:
    if crop_system == "other":
        return "Other systems"
    if not crop_system.startswith("multi:"):
        return crop_system.replace("_", " ").title()
    combination = crop_system.removeprefix("multi:")
    if combination.startswith("double_"):
        crop = combination.removeprefix("double_").replace("_", " ")
        return f"Double {crop}"
    if combination.startswith("triple_"):
        crop = combination.removeprefix("triple_").replace("_", " ")
        return f"Triple {crop}"
    return combination.replace("_", " + ").title()


def _baseline_colors(
    data: pd.DataFrame, order: list[str]
) -> dict[str, tuple[float, float, float, float]]:
    baseline = data[data["price"].eq(0.0)].set_index("crop_group").reindex(order)
    intensity = baseline["water_mm3"].div(
        baseline["production_mt"].replace(0.0, np.nan)
    )
    if intensity.isna().any() or intensity.le(0.0).any():
        raise ValueError("Baseline crop-system irrigation intensities must be positive")
    blues = plt.colormaps["Blues"]
    cmap = LinearSegmentedColormap.from_list(
        "crop_intensity_blues", blues(np.linspace(0.28, 0.92, 256))
    )
    norm = Normalize(vmin=float(intensity.min()), vmax=float(intensity.max()))
    return {crop: cmap(norm(value)) for crop, value in intensity.items()}


def _panel_title(ax, letter: str, title: str) -> None:
    ax.set_title(title, loc="left", fontsize=13, pad=7)
    ax.text(
        -0.035,
        1.025,
        letter,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=13,
        fontweight="bold",
    )


def _spread_labels(anchors: np.ndarray, ymin: float, ymax: float) -> np.ndarray:
    positions = anchors.copy()
    gap = 0.06 * (ymax - ymin)
    order = np.argsort(positions)
    ordered = positions[order]
    ordered[0] = max(ordered[0], ymin + 0.03 * (ymax - ymin))
    for index in range(1, len(ordered)):
        ordered[index] = max(ordered[index], ordered[index - 1] + gap)
    overflow = ordered[-1] - 0.96 * ymax
    if overflow > 0:
        ordered -= overflow
    positions[order] = ordered
    return positions


def _spread_x_labels(anchors: np.ndarray, xmax: float) -> np.ndarray:
    """Spread ordered mosaic callouts while keeping leaders uncrossed."""
    xmin = 0.015 * xmax
    label_xmax = 0.78 * xmax
    positions = np.clip(anchors.copy(), xmin, label_xmax)
    gap = 0.065 * xmax
    for index in range(1, len(positions)):
        positions[index] = max(positions[index], positions[index - 1] + gap)
    if positions[-1] > label_xmax:
        positions[-1] = label_xmax
        for index in range(len(positions) - 2, -1, -1):
            positions[index] = min(positions[index], positions[index + 1] - gap)
    if positions[0] < xmin:
        positions = np.linspace(xmin, label_xmax, len(positions))
    return positions


def _plot_water_sweep(
    ax,
    data: pd.DataFrame,
    order: list[str],
    colors: dict[str, tuple[float, float, float, float]],
) -> None:
    pivot = data.pivot(index="price", columns="crop_group", values="water_mm3")
    pivot = pivot.reindex(columns=order, fill_value=0.0).sort_index() / 1000
    x = pivot.index.to_numpy()
    anchors = []
    for crop in order:
        values = pivot[crop].to_numpy()
        ax.plot(x, values, color=colors[crop], linewidth=2.0, solid_capstyle="round")
        anchors.append(values[0])
    ax.set_xscale("symlog", linthresh=0.0001)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=0.0)
    price_ticks = [0.0, 0.0001, 0.001, 0.01, 0.1, 1.0]
    ax.set_xticks(price_ticks, ["0", "0.0001", "0.001", "0.01", "0.1", "1"])
    ax.set_xlabel("Water price (USD/m3-eq)", fontsize=11)
    ax.set_ylabel("Irrigation (km3)", fontsize=11)
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", left=False, labelleft=False, right=True, labelright=True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="both", labelsize=10)
    _panel_title(ax, "a", "Irrigation by crop system")
    ax.grid(alpha=0.12)

    ymin, ymax = ax.get_ylim()
    positions = _spread_labels(np.asarray(anchors), ymin, ymax)
    label_x = -0.03
    label_transform = ax.get_yaxis_transform()
    for crop, anchor, position in zip(order, anchors, positions, strict=True):
        ax.plot(
            [label_x, 0.0],
            [position, anchor],
            transform=label_transform,
            color=colors[crop],
            linewidth=1.4,
            solid_capstyle="butt",
            clip_on=False,
            zorder=3,
        )
        ax.scatter(
            [0.0],
            [anchor],
            s=20,
            color=colors[crop],
            edgecolor="white",
            linewidth=0.5,
            zorder=5,
            clip_on=False,
        )
        ax.text(
            label_x,
            position,
            _display_name(crop),
            transform=label_transform,
            ha="right",
            va="center",
            color="#222222",
            fontsize=10.5,
            fontweight="medium",
            clip_on=False,
            zorder=4,
        )


def _plot_mosaic(
    ax,
    data: pd.DataFrame,
    price: float,
    order: list[str],
    colors: dict[str, tuple[float, float, float, float]],
    xmax: float,
    ymax: float,
    letter: str,
) -> None:
    frame = data[data["price"].eq(price)].set_index("crop_group").reindex(order)
    frame = frame.fillna(0.0)
    frame["intensity_m3_per_t"] = (
        frame["water_mm3"].div(frame["production_mt"].replace(0.0, np.nan)).fillna(0.0)
    )
    frame["tie_order"] = [order.index(crop) for crop in frame.index]
    frame = frame.sort_values(
        ["intensity_m3_per_t", "tie_order"], ascending=[False, True]
    )
    left = 0.0
    label_crops = []
    label_anchors = []
    for crop, row in frame.iterrows():
        width = float(row["production_mt"])
        height = float(row["intensity_m3_per_t"])
        center = left + width / 2
        ax.add_patch(
            Rectangle(
                (left, 0.0),
                width,
                height,
                facecolor=colors[crop],
                edgecolor="white",
                linewidth=0.8,
                alpha=0.9,
            )
        )
        if height >= 0.02 * ymax:
            label_crops.append(crop)
            label_anchors.append((center, height))
        left += width
    label_x = (
        _spread_x_labels(np.asarray([center for center, _ in label_anchors]), xmax)
        if label_anchors
        else np.asarray([])
    )
    for crop, (center, height), text_x in zip(
        label_crops, label_anchors, label_x, strict=True
    ):
        text_y = height + 0.018 * ymax
        is_tall = height > 0.55 * ymax
        if is_tall:
            text_x = max(text_x, center + 0.02 * xmax)
        ax.plot(
            [center, text_x],
            [height, text_y],
            color=colors[crop],
            linewidth=1.0,
            solid_capstyle="butt",
            zorder=3,
        )
        ax.text(
            text_x,
            text_y,
            _display_name(crop),
            ha="left",
            va="center",
            rotation=0 if is_tall else 40,
            rotation_mode="anchor",
            fontsize=8.2,
            fontweight="medium",
            color="#222222",
            clip_on=False,
            zorder=4,
        )
    ax.set_xlim(0.0, xmax)
    ax.set_ylim(0.0, ymax)
    ax.tick_params(axis="x", labelsize=9)
    ax.tick_params(axis="y", labelsize=9.5)
    _panel_title(ax, letter, f"{price:g} USD/m3-eq")
    ax.grid(axis="y", alpha=0.12)


def _plot(
    data: pd.DataFrame,
    order: list[str],
    selected_prices: list[float],
    output: Path,
) -> None:
    colors = _baseline_colors(data, order)
    mosaics = data[data["price"].isin(selected_prices)].copy()
    mosaics["intensity"] = mosaics["water_mm3"].div(
        mosaics["production_mt"].replace(0.0, np.nan)
    )
    xmax = 1.02 * mosaics.groupby("price")["production_mt"].sum().max()
    ymax = 1.42 * mosaics["intensity"].max()

    fig = plt.figure(figsize=(14.5, 9.5))
    grid = fig.add_gridspec(
        3, 2, height_ratios=(1.2, 1.0, 1.0), hspace=0.62, wspace=0.22
    )
    sweep_ax = fig.add_subplot(grid[0, :])
    _plot_water_sweep(sweep_ax, data, order, colors)
    mosaic_axes = [
        fig.add_subplot(grid[row, column]) for row in (1, 2) for column in (0, 1)
    ]
    for index, (ax, price) in enumerate(
        zip(mosaic_axes, selected_prices, strict=True), start=1
    ):
        _plot_mosaic(
            ax,
            data,
            price,
            order,
            colors,
            xmax,
            ymax,
            chr(ord("a") + index),
        )
        if index in (1, 3):
            ax.set_ylabel("Irrigation intensity\n(m3/t DM)", fontsize=10.5)
        else:
            ax.tick_params(axis="y", labelleft=False)
        if index > 2:
            ax.set_xlabel("Crop output (Mt DM)", fontsize=10.5)
    fig.suptitle(
        "Crop irrigation response to water pricing",
        x=0.19,
        y=0.98,
        ha="left",
        fontsize=17,
        fontweight="semibold",
    )
    fig.text(
        0.19,
        0.925,
        "Width = output | height = intensity | area = irrigation | "
        "darker = higher baseline intensity",
        fontsize=10.5,
        color="#444444",
    )
    fig.subplots_adjust(left=0.19, right=0.94, bottom=0.18, top=0.86)
    fig.savefig(output, dpi=200, facecolor="white")
    plt.close(fig)


def _plot_intensity_stack(
    data: pd.DataFrame,
    output: Path,
    x_label: str = "Water scarcity price (USD/m3-eq)",
    linthresh: float = 0.01,
    ordering_label: str = "baseline intensity",
    title: str = "Irrigation by crop system and intensity",
) -> None:
    """Plot exact crop-system water use, colored by irrigation intensity."""
    baseline = data[data["price"].eq(0.0)].set_index("crop_system").copy()
    baseline["intensity"] = baseline["water_mm3"].div(
        baseline["production_mt"].replace(0.0, np.nan)
    )
    order = baseline.sort_values(
        ["intensity", "water_mm3"], ascending=[False, False], na_position="last"
    ).index.tolist()
    water = (
        data.pivot(index="price", columns="crop_system", values="water_mm3")
        .reindex(columns=order, fill_value=0.0)
        .sort_index()
        .fillna(0.0)
        / 1000
    )
    production = (
        data.pivot(index="price", columns="crop_system", values="production_mt")
        .reindex(index=water.index, columns=order, fill_value=0.0)
        .fillna(0.0)
    )
    intensity = water.mul(1000).div(production.replace(0.0, np.nan)).fillna(0.0)
    lower = water.cumsum(axis=1).sub(water)
    upper = lower.add(water)
    x = water.index.to_numpy()

    blues = plt.colormaps["Blues"]
    cmap = LinearSegmentedColormap.from_list(
        "dynamic_intensity_blues", blues(np.linspace(0.18, 0.95, 256))
    )
    vmax = 100 * np.ceil(float(intensity.to_numpy().max()) / 100)
    norm = Normalize(vmin=0.0, vmax=vmax)
    fig, ax = plt.subplots(figsize=(13.5, 7.2))
    polygons = []
    color_values = []
    subdivisions = 12
    for crop in order:
        if water[crop].max() <= 0.0:
            continue
        for index in range(len(x) - 1):
            fractions = np.linspace(0.0, 1.0, subdivisions + 1)
            if x[index] == 0.0:
                segment_x = np.linspace(x[index], x[index + 1], subdivisions + 1)
            else:
                segment_x = np.geomspace(x[index], x[index + 1], subdivisions + 1)
            segment_lower = np.interp(
                fractions,
                [0.0, 1.0],
                [lower[crop].iloc[index], lower[crop].iloc[index + 1]],
            )
            segment_upper = np.interp(
                fractions,
                [0.0, 1.0],
                [upper[crop].iloc[index], upper[crop].iloc[index + 1]],
            )
            segment_intensity = np.interp(
                fractions,
                [0.0, 1.0],
                [intensity[crop].iloc[index], intensity[crop].iloc[index + 1]],
            )
            for step in range(subdivisions):
                polygons.append(
                    [
                        (segment_x[step], segment_lower[step]),
                        (segment_x[step + 1], segment_lower[step + 1]),
                        (segment_x[step + 1], segment_upper[step + 1]),
                        (segment_x[step], segment_upper[step]),
                    ]
                )
                color_values.append(
                    0.5 * (segment_intensity[step] + segment_intensity[step + 1])
                )
    ax.add_collection(
        PolyCollection(
            polygons,
            array=np.asarray(color_values),
            cmap=cmap,
            norm=norm,
            edgecolors="none",
            antialiaseds=False,
            rasterized=True,
        )
    )
    for crop in order[:-1]:
        ax.plot(
            x,
            upper[crop],
            color="white",
            linewidth=0.35,
            alpha=0.8,
            solid_joinstyle="round",
            solid_capstyle="round",
        )
    ax.set_xscale("symlog", linthresh=linthresh)
    ax.set_xlim(0.0, 1.0)
    ymax = 1.02 * float(upper.iloc[0, -1])
    ax.set_ylim(0.0, ymax)
    price_ticks = [0.0, 0.01, 0.1, 1.0]
    ax.set_xticks(price_ticks, ["0", "0.01", "0.1", "1"])
    ax.set_xlabel(x_label, fontsize=12)
    ax.tick_params(axis="both", labelsize=10.5)
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", left=False, labelleft=False, right=True, labelright=True)
    ax.grid(axis="y", alpha=0.15, linewidth=0.7)
    _panel_title(ax, "a", title)

    labelled = baseline["water_mm3"].nlargest(12).index.tolist()
    anchors = np.asarray(
        [0.5 * (lower[crop].iloc[0] + upper[crop].iloc[0]) for crop in labelled]
    )
    positions = _spread_labels(anchors, 0.0, ymax)
    label_x = -0.075
    label_transform = ax.get_yaxis_transform()
    baseline_intensity = intensity.iloc[0]
    for crop, anchor, position in zip(labelled, anchors, positions, strict=True):
        color = cmap(norm(baseline_intensity[crop]))
        ax.plot(
            [label_x, 0.0],
            [position, anchor],
            transform=label_transform,
            color=color,
            linewidth=1.2,
            solid_capstyle="butt",
            clip_on=False,
        )
        ax.text(
            label_x,
            position,
            _display_name(crop),
            transform=label_transform,
            ha="right",
            va="center",
            fontsize=9.5,
            color="#222222",
            clip_on=False,
        )

    colorbar = fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.065, fraction=0.035
    )
    colorbar.set_label("Irrigation intensity (m3/t DM)", fontsize=10.5)
    colorbar.ax.tick_params(labelsize=9.5)
    fig.text(
        0.23,
        0.935,
        f"Global irrigation water (km3); {len(order)} systems ordered by "
        f"{ordering_label}; 12 largest users labelled",
        fontsize=10,
        color="#555555",
    )
    fig.subplots_adjust(left=0.23, right=0.91, bottom=0.14, top=0.88)
    fig.savefig(output, dpi=220, facecolor="white")
    plt.close(fig)


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
    parser.add_argument("--top-crops", type=int, default=8)
    parser.add_argument(
        "--mosaic-prices", type=float, nargs=4, default=[0.0, 0.01, 0.1, 1.0]
    )
    args = parser.parse_args()

    paths = sorted(
        (
            path
            for path in args.solved_dir.glob("model_scen-wf*.nc")
            if _scenario_price(path) <= 1.0
        ),
        key=_scenario_price,
    )
    if not paths or _scenario_price(paths[0]) != 0.0:
        raise ValueError("The crop-water analysis requires scenario wf00000")
    rows = []
    for path in paths:
        frame = _extract_crop_water(pypsa.Network(path))
        frame.insert(0, "price", _scenario_price(path))
        rows.append(frame)
    crop_data = pd.concat(rows, ignore_index=True)
    grouped, order = _group_crops(crop_data, args.top_crops)
    missing_prices = set(args.mosaic_prices) - set(grouped["price"])
    if missing_prices:
        raise ValueError(f"Missing requested mosaic prices: {sorted(missing_prices)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    crop_data.to_csv(args.output_dir / "crop_water_by_price.csv", index=False)
    grouped.to_csv(args.output_dir / "crop_water_grouped_by_price.csv", index=False)
    _plot(
        grouped,
        order,
        list(args.mosaic_prices),
        args.output_dir / "crop_water_mosaic.png",
    )
    _plot_intensity_stack(crop_data, args.output_dir / "crop_water_intensity_stack.png")


if __name__ == "__main__":
    main()
