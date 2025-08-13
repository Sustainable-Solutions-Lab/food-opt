# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import pypsa
import matplotlib

matplotlib.use("pdf")  # Use PDF backend
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path


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
        print("No crop production found")
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

    print(f"Crop production plot saved to {output_dir / 'crop_production.pdf'}")


def plot_food_production(food_production: pd.Series, output_dir: Path) -> None:
    """Create bar plot for food production."""
    # Sort by production value for better visualization
    food_production_sorted = food_production.sort_values(ascending=False)

    # Filter out foods with zero production
    food_production_sorted = food_production_sorted[food_production_sorted > 1e-6]

    if len(food_production_sorted) == 0:
        print("No food production found")
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

    print(f"Food production plot saved to {output_dir / 'food_production.pdf'}")


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
        print("No resource usage data found")
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

    print(f"Resource usage plot saved to {output_dir / 'resource_usage.pdf'}")


def plot_network_topology(n: pypsa.Network, output_dir: Path) -> None:
    """Create network topology visualization with nodes and edges."""
    # Create NetworkX graph
    G = nx.Graph()

    # Define node colors and sizes by type
    node_colors = {}
    node_sizes = {}

    # Color scheme for different types of buses
    color_map = {
        "land": "#8B4513",  # Brown
        "water": "#1E90FF",  # Blue
        "fertilizer": "#FFD700",  # Gold
        "co2": "#FF6347",  # Red
        "ch4": "#FF69B4",  # Pink
        "crop_": "#32CD32",  # Lime green for crops
        "food_": "#FF8C00",  # Orange for foods
        "carb": "#DDA0DD",  # Plum for nutrients
        "protein": "#DDA0DD",
        "fat": "#DDA0DD",
        "grain": "#F4A460",  # Sandy brown for food groups
        "whole grain": "#F4A460",
        "legume": "#F4A460",
        "fruit": "#F4A460",
        "vegetable": "#F4A460",
        "oil": "#F4A460",
        "nuts and seeds": "#F4A460",
        "starchy vegetable": "#F4A460",
    }

    # Add all buses as nodes
    for bus in n.buses.index:
        # Skip empty or invalid bus names
        if not bus or pd.isna(bus):
            continue

        G.add_node(bus)

        # Determine color
        color = "#CCCCCC"  # Default gray
        for key, col in color_map.items():
            if bus.startswith(key) or bus == key:
                color = col
                break

        node_colors[bus] = color

        # Determine size based on type
        if bus in ["land", "water", "fertilizer", "co2", "ch4"]:
            node_sizes[bus] = 800  # Large for resources
        elif bus.startswith("crop_"):
            node_sizes[bus] = 500  # Medium for crops
        elif bus.startswith("food_"):
            node_sizes[bus] = 400  # Medium for foods
        else:
            node_sizes[bus] = 300  # Small for others

    # Add edges from links
    link_colors = {}
    link_widths = {}

    for link_name, link in n.links.iterrows():
        # Get all connected buses for this link
        buses = []
        for bus_col in ["bus0", "bus1", "bus2", "bus3", "bus4"]:
            if bus_col in link.index and pd.notna(link[bus_col]) and link[bus_col]:
                buses.append(link[bus_col])

        # Create edges between all pairs of buses in this link
        for i in range(len(buses)):
            for j in range(i + 1, len(buses)):
                bus1, bus2 = buses[i], buses[j]

                # Skip if either bus is not in the graph
                if not G.has_node(bus1) or not G.has_node(bus2):
                    continue

                # Add edge if not already exists
                if not G.has_edge(bus1, bus2):
                    G.add_edge(bus1, bus2, links=[])

                # Track which links use this edge
                G.edges[bus1, bus2]["links"].append(link_name)

                # Set edge properties
                if link_name.startswith("produce_"):
                    link_colors[(bus1, bus2)] = "#228B22"  # Forest green for production
                    link_widths[(bus1, bus2)] = 2
                elif link_name.startswith("convert_"):
                    link_colors[(bus1, bus2)] = "#4169E1"  # Royal blue for conversion
                    link_widths[(bus1, bus2)] = 2
                elif link_name.startswith("consume_"):
                    link_colors[(bus1, bus2)] = "#DC143C"  # Crimson for consumption
                    link_widths[(bus1, bus2)] = 2
                else:
                    link_colors[(bus1, bus2)] = "#696969"  # Gray for other links
                    link_widths[(bus1, bus2)] = 1

    # Create layout
    print(
        f"Creating network layout with {len(G.nodes)} nodes and {len(G.edges)} edges..."
    )

    # Create custom layered layout
    pos = {}

    # Define layers (y-positions from top to bottom)
    layers = {
        "resources": 3.0,  # Top row: Primary resources
        "crops": 2.0,  # Second row: Crops
        "foods": 1.0,  # Third row: Foods
        "nutrients": 0.0,  # Bottom row: Food groups and nutrients
    }

    # Categorize nodes into layers
    layer_nodes = {"resources": [], "crops": [], "foods": [], "nutrients": []}

    for node in G.nodes():
        if node in ["land", "water", "fertilizer", "co2", "ch4"]:
            layer_nodes["resources"].append(node)
        elif node.startswith("crop_"):
            layer_nodes["crops"].append(node)
        elif node.startswith("food_"):
            layer_nodes["foods"].append(node)
        else:  # Food groups and nutrients
            layer_nodes["nutrients"].append(node)

    # Position nodes in each layer
    for layer_name, y_pos in layers.items():
        nodes_in_layer = layer_nodes[layer_name]
        if nodes_in_layer:
            # Distribute nodes horizontally across the layer
            n_nodes = len(nodes_in_layer)
            for i, node in enumerate(sorted(nodes_in_layer)):
                # Spread nodes evenly across x-axis, centered
                if n_nodes == 1:
                    x_pos = 0.0
                else:
                    x_pos = (i / (n_nodes - 1) - 0.5) * (
                        n_nodes * 0.8
                    )  # Scale to spread nicely
                pos[node] = (x_pos, y_pos)

    # Create plot
    plt.figure(figsize=(20, 16))

    # Draw nodes
    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=[node_colors[node] for node in G.nodes()],
        node_size=[node_sizes[node] for node in G.nodes()],
        alpha=0.8,
    )

    # Draw edges
    for edge in G.edges():
        color = link_colors.get(edge, "#696969")
        width = link_widths.get(edge, 1)
        nx.draw_networkx_edges(G, pos, [edge], edge_color=color, width=width, alpha=0.6)

    # Draw labels
    # Only label important nodes to avoid clutter
    important_nodes = {
        node
        for node in G.nodes()
        if node
        in ["land", "water", "fertilizer", "co2", "ch4", "carb", "protein", "fat"]
        or node.startswith("grain")
        or node.startswith("fruit")
    }

    nx.draw_networkx_labels(
        G,
        pos,
        labels={node: node for node in important_nodes},
        font_size=8,
        font_weight="bold",
    )

    plt.title("Food System Network Topology", size=16, weight="bold")

    # Add categorical labels for each row
    ax = plt.gca()

    # Find the leftmost x position across all nodes to position labels
    all_x_positions = [pos[node][0] for node in pos.keys()]
    leftmost_x = min(all_x_positions) - 1.5  # Position labels to the left of all nodes

    # Add row labels
    layer_labels = {
        "resources": "Primary\nResources",
        "crops": "Crops",
        "foods": "Food\nProducts",
        "nutrients": "Nutrients &\nFood Groups",
    }

    for layer_name, y_pos in layers.items():
        if layer_nodes[layer_name]:  # Only add label if layer has nodes
            label_text = layer_labels[layer_name]
            ax.text(
                leftmost_x,
                y_pos,
                label_text,
                fontsize=12,
                fontweight="bold",
                ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.7),
            )

    plt.axis("off")

    # Create legend
    legend_elements = [
        plt.scatter([], [], c="#8B4513", s=100, label="Resources (land, water, etc.)"),
        plt.scatter([], [], c="#32CD32", s=100, label="Crops"),
        plt.scatter([], [], c="#FF8C00", s=100, label="Foods"),
        plt.scatter([], [], c="#DDA0DD", s=100, label="Nutrients"),
        plt.scatter([], [], c="#F4A460", s=100, label="Food Groups"),
        plt.Line2D([0], [0], color="#228B22", lw=2, label="Production Links"),
        plt.Line2D([0], [0], color="#4169E1", lw=2, label="Conversion Links"),
        plt.Line2D([0], [0], color="#DC143C", lw=2, label="Consumption Links"),
    ]

    plt.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(0, 1))
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_dir / "network_topology.pdf", bbox_inches="tight", dpi=300)
    plt.close()

    print(f"Network topology plot saved to {output_dir / 'network_topology.pdf'}")

    # Also save network statistics
    stats = {
        "nodes": len(G.nodes),
        "edges": len(G.edges),
        "density": nx.density(G),
        "connected_components": nx.number_connected_components(G),
    }

    print(f"Network statistics: {stats}")


if __name__ == "__main__":
    # Load the solved network
    print("Loading solved network...")
    n = pypsa.Network(snakemake.input.network)

    # Create output directory
    output_dir = Path(snakemake.output.plots_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Creating plots in {output_dir}")

    # Extract data
    print("Extracting crop production data...")
    crop_production = extract_crop_production(n)
    print(f"Found {len(crop_production)} crops with production data")

    print("Extracting food production data...")
    food_production = extract_food_production(n)
    print(f"Found {len(food_production)} foods with production data")

    # Create plots
    plot_crop_production(crop_production, output_dir)
    plot_food_production(food_production, output_dir)
    plot_resource_usage(n, output_dir)
    plot_network_topology(n, output_dir)

    # Save summary data as CSV for reference
    if len(crop_production) > 0:
        crop_production.to_csv(
            output_dir / "crop_production.csv", header=["production_tonnes"]
        )

    if len(food_production) > 0:
        food_production.to_csv(
            output_dir / "food_production.csv", header=["production_tonnes"]
        )

    print("Plotting completed successfully!")
