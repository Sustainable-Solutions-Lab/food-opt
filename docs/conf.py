# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: CC-BY-4.0

"""Sphinx configuration for food-opt documentation."""

import os
import sys

# Add project root to path for autodoc
sys.path.insert(0, os.path.abspath(".."))

# Project information
project = "food-opt"
copyright = "2025, Koen van Greevenbroek"
author = "Koen van Greevenbroek"
release = "0.1.0"

# General configuration
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx_autodoc_typehints",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    ".uv-cache",
    "*/.uv-cache/*",
]

# HTML output options
html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "_static/logo.svg"
html_favicon = "_static/logo.svg"
html_title = "food-opt"
html_theme_options = {
    "navigation_with_keys": True,
    "light_css_variables": {
        "color-brand-primary": "#3b745f",
        "color-brand-content": "#2f5e49",
    },
    "dark_css_variables": {
        "color-brand-primary": "#5fa285",
        "color-brand-content": "#7db79e",
    },
}

# Autodoc options
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__",
}

# Napoleon settings for NumPy docstrings
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True

# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
    "pypsa": ("https://pypsa.readthedocs.io/en/latest/", None),
}

# Type hints configuration
typehints_fully_qualified = False
always_document_param_types = True
typehints_document_rtype = True
# Autodoc tweaks
autodoc_typehints = "none"
autodoc_mock_imports = [
    "linopy",
    "pypsa",
]
