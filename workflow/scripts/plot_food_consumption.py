# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Plot global food consumption by food group per person per day."""

import logging
from pathlib import Path

import matplotlib
import pandas as pd
import pypsa

matplotlib.use("pdf")
import matplotlib.pyplot as plt


logger = logging.getLogger(__name__)

GRAMS_PER_MEGATONNE = 1e12
KCAL_PER_MCAL = 1e6
DAYS_PER_YEAR = 365


def _select_snapshot(network: pypsa.Network) -> pd.Index | str:
    if "now" in network.snapshots:
        return "now"
    if len(network.snapshots) == 1:
        return network.snapshots[0]
    raise ValueError("Expected snapshot 'now' or single snapshot in solved network")


def _group_from_bus(bus: str) -> str:
    remainder = bus[len("group_") :]
    if "_" in remainder:
        return remainder.rsplit("_", 1)[0]
    return remainder


def _aggregate_group_mass(network: pypsa.Network, snapshot) -> pd.Series:
    loads = network.loads
    if loads.empty or not hasattr(network.loads_t, "p"):
        return pd.Series(dtype=float)

    load_dispatch = network.loads_t.p
    if snapshot not in load_dispatch.index:
        raise ValueError(f"Snapshot '{snapshot}' not present in loads time series")

    snapshot_values = load_dispatch.loc[snapshot]
    totals: dict[str, float] = {}

    for load_name in loads.index:
        bus = str(loads.at[load_name, "bus"])
        if not bus.startswith("group_"):
            continue
        value = abs(float(snapshot_values.get(load_name, 0.0)))
        if value <= 0.0:
            continue
        group = _group_from_bus(bus)
        totals[group] = totals.get(group, 0.0) + value

    return pd.Series(totals, dtype=float)


def _available_legs(links: pd.DataFrame) -> list[int]:
    legs: set[int] = set()
    for column in links.columns:
        if not column.startswith("bus"):
            continue
        if column == "bus0":
            continue
        suffix = column[3:]
        if not suffix:
            continue
        try:
            legs.add(int(suffix))
        except ValueError:
            continue
    return sorted(legs)


def _aggregate_group_calories(network: pypsa.Network, snapshot) -> pd.Series:
    links = network.links
    if links.empty:
        return pd.Series(dtype=float)

    legs = _available_legs(links)
    if not legs:
        return pd.Series(dtype=float)

    time_series_lookup: dict[int, pd.Series] = {}
    for leg in legs:
        attr = f"p{leg}"
        series = getattr(network.links_t, attr, None)
        if series is None or snapshot not in series.index:
            continue
        time_series_lookup[leg] = series.loc[snapshot]

    if not time_series_lookup:
        return pd.Series(dtype=float)

    totals: dict[str, float] = {}

    for link_name in links.index:
        if not link_name.startswith("consume_"):
            continue

        group_name: str | None = None
        kcal_leg: int | None = None

        for leg in legs:
            column = f"bus{leg}"
            if column not in links.columns:
                continue
            bus_value = links.at[link_name, column]
            if pd.isna(bus_value):
                continue
            bus_str = str(bus_value)
            if bus_str.startswith("group_"):
                group_name = _group_from_bus(bus_str)
            if bus_str.startswith("kcal_"):
                kcal_leg = leg

        if group_name is None or kcal_leg is None:
            continue

        series = time_series_lookup.get(kcal_leg)
        if series is None:
            continue

        value = abs(float(series.get(link_name, 0.0)))
        if value <= 0.0:
            continue

        totals[group_name] = totals.get(group_name, 0.0) + value

    return pd.Series(totals, dtype=float)


def _assign_colors(groups: list[str]) -> dict[str, tuple[float, float, float, float]]:
    cmap = plt.get_cmap("tab20")
    return {group: cmap(i % cmap.N) for i, group in enumerate(groups)}


def _plot(
    mass_g_per_person_day: pd.Series,
    calories_kcal_per_person_day: pd.Series,
    output_pdf: Path,
) -> None:
    mass_g_per_person_day = mass_g_per_person_day[mass_g_per_person_day > 0]
    calories_kcal_per_person_day = calories_kcal_per_person_day[
        calories_kcal_per_person_day > 0
    ]

    ordered_groups: list[str] = []
    ordered_groups.extend(
        mass_g_per_person_day.sort_values(ascending=False).index.tolist()
    )
    for group in calories_kcal_per_person_day.sort_values(ascending=False).index:
        if group not in ordered_groups:
            ordered_groups.append(group)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    if not ordered_groups:
        for ax in axes:
            ax.text(
                0.5, 0.5, "No food group consumption data", ha="center", va="center"
            )
            ax.axis("off")
        fig.tight_layout()
        fig.savefig(output_pdf, bbox_inches="tight", dpi=300)
        plt.close(fig)
        return

    colors = _assign_colors(ordered_groups)

    # Mass subplot
    ax_mass = axes[0]
    bottom = 0.0
    for group in ordered_groups:
        value = float(mass_g_per_person_day.get(group, 0.0))
        if value <= 0.0:
            continue
        ax_mass.bar(0, value, bottom=bottom, color=colors[group], label=group)
        bottom += value

    ax_mass.set_xticks([0])
    ax_mass.set_xticklabels(["Mass"])
    ax_mass.set_ylabel("g/person/day")
    ax_mass.set_title("Global Food Consumption (Mass)")
    ax_mass.grid(axis="y", alpha=0.3)

    # Calories subplot
    ax_cal = axes[1]
    bottom = 0.0
    for group in ordered_groups:
        value = float(calories_kcal_per_person_day.get(group, 0.0))
        if value <= 0.0:
            continue
        ax_cal.bar(0, value, bottom=bottom, color=colors[group])
        bottom += value

    ax_cal.set_xticks([0])
    ax_cal.set_xticklabels(["Calories"])
    ax_cal.set_ylabel("kcal/person/day")
    ax_cal.set_title("Global Food Consumption (Calories)")
    ax_cal.grid(axis="y", alpha=0.3)

    handles, labels = ax_mass.get_legend_handles_labels()
    if handles:
        fig.legend(
            handles[::-1],
            labels[::-1],
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
        )
        fig.tight_layout(rect=(0, 0, 0.85, 1))
    else:
        fig.tight_layout()

    fig.savefig(output_pdf, bbox_inches="tight", dpi=300)
    plt.close(fig)


def main() -> None:
    try:
        snakemake  # type: ignore[name-defined]
    except NameError as exc:  # pragma: no cover - Snakemake injects the variable
        raise RuntimeError("This script must be run from Snakemake") from exc

    network_path = snakemake.input.network  # type: ignore[attr-defined]
    population_path = snakemake.input.population  # type: ignore[attr-defined]
    output_pdf = Path(snakemake.output.pdf)  # type: ignore[attr-defined]
    output_pdf.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Loading solved network from %s", network_path)
    network = pypsa.Network(network_path)

    snapshot = _select_snapshot(network)
    logger.info("Using snapshot '%s' for aggregation", snapshot)

    mass = _aggregate_group_mass(network, snapshot)
    calories = _aggregate_group_calories(network, snapshot)

    population_df = pd.read_csv(population_path)
    if "population" not in population_df.columns:
        raise ValueError("Population file must contain a 'population' column")
    total_population = float(population_df["population"].sum())
    if total_population <= 0.0:
        raise ValueError("Total population must be positive for per-capita conversion")

    mass_per_capita = mass * GRAMS_PER_MEGATONNE / (total_population * DAYS_PER_YEAR)
    calories_per_capita = calories * KCAL_PER_MCAL / (total_population * DAYS_PER_YEAR)

    logger.info(
        "Found %d food groups with mass data and %d with calorie data",
        mass_per_capita[mass_per_capita > 0].shape[0],
        calories_per_capita[calories_per_capita > 0].shape[0],
    )

    _plot(mass_per_capita, calories_per_capita, output_pdf)

    logger.info("Food consumption plot saved to %s", output_pdf)


if __name__ == "__main__":
    main()
