#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Render DOT graph to SVG using Graphviz via pydot."""

import pydot


def main(dot_path: str, svg_path: str):
    """Render DOT file to SVG.

    Args:
        dot_path: Path to input DOT file
        svg_path: Path for output SVG file
    """
    # Read DOT file
    graphs = pydot.graph_from_dot_file(dot_path)
    if not graphs:
        raise ValueError(f"Could not parse DOT file: {dot_path}")
    graph = graphs[0]

    # Render to SVG using Graphviz
    svg_data = graph.create_svg()

    # Write SVG file
    with open(svg_path, "wb") as f:
        f.write(svg_data)


if __name__ == "__main__":
    # Snakemake integration
    main(
        dot_path=snakemake.input.dot,
        svg_path=snakemake.output.svg,
    )
