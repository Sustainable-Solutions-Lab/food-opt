"""
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: GPL-3.0-or-later
"""

import numpy as np
import rasterio
from pyproj import Geod


def calculate_all_cell_areas(src: rasterio.DatasetReader) -> np.ndarray:
    """Return per-pixel area in hectares for a geographic (lon/lat) raster."""
    pixel_width_deg = abs(src.transform.a)
    pixel_height_deg = abs(src.transform.e)
    rows, cols = src.shape
    left, bottom, right, top = src.bounds
    lats = np.linspace(top - pixel_height_deg / 2, bottom + pixel_height_deg / 2, rows)
    geod = Geod(ellps="WGS84")
    areas_ha = np.zeros(rows)
    for i, lat in enumerate(lats):
        lat_top = lat + pixel_height_deg / 2
        lat_bottom = lat - pixel_height_deg / 2
        lon_left = left
        lon_right = left + pixel_width_deg
        lons = [lon_left, lon_right, lon_right, lon_left, lon_left]
        lats_poly = [lat_bottom, lat_bottom, lat_top, lat_top, lat_bottom]
        area_m2, _ = geod.polygon_area_perimeter(lons, lats_poly)
        areas_ha[i] = abs(area_m2) / 10000.0
    return np.repeat(areas_ha[:, np.newaxis], cols, axis=1)


def scale_fraction(arr: np.ndarray) -> np.ndarray:
    """Scale array to 0..1 if stored as 0..100 or 0..10000; clip to [0,1]."""
    a = arr.astype(float)
    a[~np.isfinite(a)] = np.nan
    vmax = np.nanmax(a)
    if vmax > 1.5:
        a = a / (100.0 if vmax <= 100 else 10000.0)
    return np.clip(a, 0.0, 1.0)


def raster_bounds(transform, width: int, height: int):
    xmin = transform.c
    ymax = transform.f
    xmax = xmin + width * transform.a
    ymin = ymax + height * transform.e
    return xmin, ymin, xmax, ymax
