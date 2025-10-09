# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Shared styling configuration for documentation figures.

This module provides consistent colors, sizes, and styling parameters
for all documentation figures. It complements the matplotlib style sheet
with domain-specific configuration.
"""

from pathlib import Path

# Color palette (matching Furo theme from docs/conf.py)
COLORS = {
    "primary": "#3b745f",  # Brand green (dark)
    "secondary": "#5fa285",  # Brand green (light)
    "accent": "#e07a5f",  # Coral/red
    "neutral": "#f2cc8f",  # Beige
    "highlight": "#81b29a",  # Teal
    "warning": "#f4a261",  # Orange
    "info": "#2a9d8f",  # Cyan
    "danger": "#e76f51",  # Red-orange
}

# Map styling defaults
MAP_STYLE = {
    "edgecolor": "white",
    "linewidth": 0.5,
    "alpha": 0.85,
}

# Common colormaps for different data types
COLORMAPS = {
    "yield": "YlGn",  # Yellow-green for crop yields
    "water": "Blues",  # Blue for water-related
    "emissions": "YlOrRd",  # Yellow-orange-red for emissions/impacts
    "health": "RdYlGn_r",  # Red-yellow-green reversed (red=bad)
    "diverging": "RdYlGn",  # Red-yellow-green for change/difference
    "sequential": "viridis",  # General sequential data
}

# Figure dimensions (width, height in inches)
FIGURE_SIZES = {
    "map_wide": (12, 6),  # Wide maps (global coverage)
    "map_square": (8, 8),  # Square maps (regional focus)
    "map_tall": (6, 9),  # Tall maps (specific regions)
    "chart": (8, 6),  # Standard charts
    "chart_wide": (10, 4),  # Wide charts (timeseries, comparisons)
    "chart_small": (6, 4),  # Small charts (insets, simple plots)
}

# Typography helpers
FONT_SIZES = {
    "title": 12,
    "label": 10,
    "legend": 9,
    "tick": 9,
    "annotation": 8,
}

# Standard DPI for output
OUTPUT_DPI = 100  # For SVG this doesn't matter much, but good for consistency


def get_crop_colors_from_config(config: dict) -> dict:
    """Extract crop color mapping from Snakemake config.

    Args:
        config: Snakemake configuration dictionary

    Returns:
        Dictionary mapping crop names to hex colors
    """
    return config.get("plotting", {}).get("colors", {}).get("crops", {})


def get_crop_fallback_cmap(config: dict) -> str:
    """Get fallback colormap for crops not in config.

    Args:
        config: Snakemake configuration dictionary

    Returns:
        Matplotlib colormap name
    """
    return config.get("plotting", {}).get("fallback_cmaps", {}).get("crops", "Set3")


def apply_doc_style():
    """Apply documentation figure styling to matplotlib.

    This loads the custom style sheet and applies additional
    programmatic styling that can't be set via mplstyle files.
    """
    import matplotlib.pyplot as plt

    # Get path to style sheet (relative to this file)
    style_path = Path(__file__).parent / "doc_figures_style.mplstyle"

    if style_path.exists():
        plt.style.use(str(style_path))
    else:
        # Fallback to minimal styling if file not found
        plt.rcParams.update(
            {
                "font.family": "sans-serif",
                "font.size": 10,
                "axes.grid": True,
                "axes.axisbelow": True,
                "grid.alpha": 0.3,
                "axes.spines.top": False,
                "axes.spines.right": False,
                "figure.facecolor": "white",
                "axes.facecolor": "white",
                "savefig.format": "svg",
                "savefig.bbox": "tight",
            }
        )


def save_doc_figure(fig, output_path: str, format: str = "svg", **kwargs):
    """Save figure with consistent settings.

    Args:
        fig: Matplotlib figure object
        output_path: Output file path (extension will be replaced if needed)
        format: Output format ('svg' or 'png')
        **kwargs: Additional arguments passed to fig.savefig()
    """

    output_path = Path(output_path)

    # Ensure correct extension
    if format == "svg":
        output_path = output_path.with_suffix(".svg")
    elif format == "png":
        output_path = output_path.with_suffix(".png")
    else:
        raise ValueError(f"Unsupported format: {format}")

    # Set defaults
    save_kwargs = {
        "format": format,
        "bbox_inches": "tight",
        "pad_inches": 0.1,
    }

    if format == "png":
        save_kwargs["dpi"] = kwargs.pop("dpi", OUTPUT_DPI)

    # Override with user kwargs
    save_kwargs.update(kwargs)

    # Ensure directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save
    fig.savefig(output_path, **save_kwargs)

    return output_path
