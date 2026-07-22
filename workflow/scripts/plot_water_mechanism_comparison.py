# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Compare physical water responses under alternative fossil-water weights."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cf30", type=Path, required=True)
    parser.add_argument("--cf100", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cases = {
        "Fossil weight 30": pd.read_csv(args.cf30),
        "Fossil weight 100": pd.read_csv(args.cf100),
    }
    colors = {"Fossil weight 30": "#4477AA", "Fossil weight 100": "#CC6677"}
    metrics = [
        ("scarcity_mm3_eq", 1e6, "Renewable AWARE scarcity\n(10^6 m3-eq)"),
        ("depletion_mm3", 1e3, "Fossil groundwater\n(km3)"),
        ("pool_demand_mm3", 1e3, "Irrigation consumption\n(km3)"),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(8.2, 8.0), sharex=True)
    for label, frame in cases.items():
        priced = frame[frame["price"] > 0].sort_values("price")
        for ax, (column, scale, ylabel) in zip(axes, metrics, strict=True):
            ax.plot(
                priced["price"],
                priced[column] / scale,
                marker="o",
                markersize=3.2,
                linewidth=1.7,
                color=colors[label],
                label=label,
            )
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.2)

    axes[0].legend(frameon=False)
    axes[-1].set_xscale("log")
    axes[-1].set_xlabel("Water price (USD per m3 scarcity-equivalent)")
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
