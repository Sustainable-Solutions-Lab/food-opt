# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import pypsa
import matplotlib

matplotlib.use("pdf")  # Use PDF backend
import matplotlib.pyplot as plt
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def extract_crop_production(n: pypsa.Network) -> pd.Series:
    """Extract crop production from solved network."""
    crop_production = {}

    # Get all crop production links (those starting with "produce_")
    production_links = [link for link in n.links.index if link.startswith("produce_")]

    for link in production_links:
        # Extract crop name from link name (remove "produce_" prefix)
        crop = link.replace("produce_", "")
        # Get the flow through this link (production level)
        flow = n.links_t.p1.loc["now", link]
        # Convert to absolute value (flows might be negative depending on convention)
        production = abs(flow)
        crop_production[crop] = production

    return pd.Series(crop_production)


def extract_food_production(n: pypsa.Network) -> pd.Series:
    """Extract food production from solved network."""
    food_production = {}

    # Get all food conversion links (those starting with "convert_")
    conversion_links = [link for link in n.links.index if link.startswith("convert_")]

    for link in conversion_links:
        # Extract food name from link name
        # Link format is "convert_{crop}_to_{food}"
        parts = link.split("_to_")
        if len(parts) == 2:
            food = parts[1].replace("_", " ")
            # Get the flow through this link
            flow = n.links_t.p0.loc["now", link]
            # The output flow is the food production
            production = abs(flow * n.links.loc[link, "efficiency"])

            if food in food_production:
                food_production[food] += production
            else:
                food_production[food] = production

    return pd.Series(food_production)


def plot_crop_production(crop_production: pd.Series, output_dir: Path) -> None:
    """Create bar plot for crop production."""
    # Sort by production value for better visualization
    crop_production_sorted = crop_production.sort_values(ascending=False)

    # Filter out crops with zero production
    crop_production_sorted = crop_production_sorted[crop_production_sorted > 1e-6]

    if len(crop_production_sorted) == 0:
        logger.warning("No crop production found")
        return

    plt.figure(figsize=(12, 8))
    plt.bar(range(len(crop_production_sorted)), crop_production_sorted.values)

    # Add value labels on bars
    max_value = crop_production_sorted.max()
    for i, (crop, value) in enumerate(crop_production_sorted.items()):
        plt.text(
            i,
            value + max_value * 0.01,
            f"{value:.1e}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    plt.xlabel("Crops")
    plt.ylabel("Production (tonnes)")
    plt.title("Crop Production by Type")
    plt.xticks(
        range(len(crop_production_sorted)),
        crop_production_sorted.index,
        rotation=45,
        ha="right",
    )
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_dir / "crop_production.pdf", bbox_inches="tight", dpi=300)
    plt.close()

    logger.info("Crop production plot saved to %s", output_dir / "crop_production.pdf")


def plot_food_production(food_production: pd.Series, output_dir: Path) -> None:
    """Create bar plot for food production."""
    # Sort by production value for better visualization
    food_production_sorted = food_production.sort_values(ascending=False)

    # Filter out foods with zero production
    food_production_sorted = food_production_sorted[food_production_sorted > 1e-6]

    if len(food_production_sorted) == 0:
        logger.warning("No food production found")
        return

    plt.figure(figsize=(14, 8))
    plt.bar(range(len(food_production_sorted)), food_production_sorted.values)

    # Add value labels on bars
    max_value = food_production_sorted.max()
    for i, (food, value) in enumerate(food_production_sorted.items()):
        plt.text(
            i,
            value + max_value * 0.01,
            f"{value:.1e}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    plt.xlabel("Food Products")
    plt.ylabel("Production (tonnes)")
    plt.title("Food Production by Type")
    plt.xticks(
        range(len(food_production_sorted)),
        food_production_sorted.index,
        rotation=45,
        ha="right",
    )
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_dir / "food_production.pdf", bbox_inches="tight", dpi=300)
    plt.close()

    logger.info("Food production plot saved to %s", output_dir / "food_production.pdf")


def plot_resource_usage(n: pypsa.Network, output_dir: Path) -> None:
    """Create bar plot for resource usage."""
    resources = ["land", "water", "fertilizer"]
    resource_usage = {}

    for resource in resources:
        # Calculate total resource consumption from links
        total_flow = 0
        for link_name in n.links.index:
            link = n.links.loc[link_name]
            # Check if this link consumes this resource (resource is input bus0)
            if resource == link.get("bus0", ""):
                flow = n.links_t.p0.loc["now", link_name]
                total_flow += abs(flow)
            # Check resource inputs via bus2 (negative efficiency means input)
            elif resource == link.get("bus2", "") and link.get("efficiency2", 0) < 0:
                if "p2" in n.links_t:
                    flow = n.links_t.p2.loc["now", link_name]
                    total_flow += abs(flow)
            # Check resource inputs via bus3 (negative efficiency means input)
            elif resource == link.get("bus3", "") and link.get("efficiency3", 0) < 0:
                if "p3" in n.links_t:
                    flow = n.links_t.p3.loc["now", link_name]
                    total_flow += abs(flow)

        resource_usage[resource] = total_flow

    if not resource_usage or all(v == 0 for v in resource_usage.values()):
        logger.warning("No resource usage data found")
        return

    resource_series = pd.Series(resource_usage)

    plt.figure(figsize=(10, 6))
    plt.bar(resource_series.index, resource_series.values)

    # Add value labels on bars with appropriate units
    units = {"land": "ha", "water": "m³", "fertilizer": "kg"}
    max_value = resource_series.max()
    for i, (resource, value) in enumerate(resource_series.items()):
        unit = units.get(resource, "units")
        plt.text(
            i,
            value + max_value * 0.01,
            f"{value:.1e} {unit}",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    plt.xlabel("Resources")
    plt.ylabel("Usage")
    plt.title("Primary Resource Usage")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_dir / "resource_usage.pdf", bbox_inches="tight", dpi=300)
    plt.close()

    logger.info("Resource usage plot saved to %s", output_dir / "resource_usage.pdf")


if __name__ == "__main__":
    # Load the solved network
    logger.info("Loading solved network...")
    n = pypsa.Network(snakemake.input.network)

    # Create output directory
    output_dir = Path(snakemake.output.plots_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Creating plots in %s", output_dir)

    # Extract data
    logger.info("Extracting crop production data...")
    crop_production = extract_crop_production(n)
    logger.info("Found %d crops with production data", len(crop_production))

    logger.info("Extracting food production data...")
    food_production = extract_food_production(n)
    logger.info("Found %d foods with production data", len(food_production))

    # Create plots
    plot_crop_production(crop_production, output_dir)
    plot_food_production(food_production, output_dir)
    plot_resource_usage(n, output_dir)

    # Save summary data as CSV for reference
    if len(crop_production) > 0:
        crop_production.to_csv(
            output_dir / "crop_production.csv", header=["production_tonnes"]
        )

    if len(food_production) > 0:
        food_production.to_csv(
            output_dir / "food_production.csv", header=["production_tonnes"]
        )

    logger.info("Plotting completed successfully!")
