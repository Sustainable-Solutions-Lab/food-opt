"""
SPDX-FileCopyrightText: 2026 Koen van Greevenbroek

SPDX-License-Identifier: GPL-3.0-or-later

Pack sparse MIRCA-OS monthly growing-area grids for repeated calendar builds.

The source NetCDF files are compact on disk but decode to dense global
``(12, 2160, 4320)`` arrays. This shared artefact stores only cells with a
nonzero value in at least one month, retaining the original float32 values and
subcrop ordering so downstream float64 sums are unchanged.
"""

from pathlib import Path

import numpy as np
import xarray as xr


def prepare_calendar(inputs: dict[str, str], output: str | Path) -> None:
    """Write a sparse, grid-validated NPZ for the supplied MIRCA subcrops."""
    labels = sorted(inputs)
    arrays: dict[str, np.ndarray] = {}
    reference: tuple[np.ndarray, np.ndarray, str] | None = None

    for i, label in enumerate(labels):
        path = inputs[label]
        with xr.open_dataset(path) as ds:
            if "harvested_area" not in ds:
                raise ValueError(f"{path} is missing 'harvested_area'")
            area = ds["harvested_area"]
            if area.dims != ("month", "latitude", "longitude"):
                raise ValueError(
                    f"{path} has dimensions {area.dims}, expected "
                    "('month', 'latitude', 'longitude')"
                )
            months = ds["month"].values
            latitude = ds["latitude"].values.astype(np.float64)
            longitude = ds["longitude"].values.astype(np.float64)
            try:
                crs_wkt = str(ds["spatial_ref"].attrs["crs_wkt"])
            except KeyError as exc:
                raise ValueError(f"{path} is missing spatial_ref.crs_wkt") from exc
            values = area.values

        if not np.array_equal(months, np.arange(1, 13)):
            raise ValueError(f"{path} does not contain numeric months 1..12")
        if latitude[0] <= latitude[-1]:
            raise ValueError(f"{path} does not have a north-up latitude axis")
        if longitude[0] >= longitude[-1]:
            raise ValueError(f"{path} does not have an eastward longitude axis")
        if values.dtype != np.float32:
            raise ValueError(f"{path} has dtype {values.dtype}, expected float32")
        if np.isinf(values).any():
            raise ValueError(f"{path} contains infinite harvested-area values")

        grid = (latitude, longitude, crs_wkt)
        if reference is None:
            reference = grid
        elif not (
            np.array_equal(latitude, reference[0])
            and np.array_equal(longitude, reference[1])
            and crs_wkt == reference[2]
        ):
            raise ValueError(f"{path} does not match the reference MIRCA grid")

        active = np.any(np.isfinite(values) & (values != 0.0), axis=0)
        cell_ids = np.flatnonzero(active).astype(np.int32)
        monthly_area = values[:, active].copy()
        monthly_area[np.isnan(monthly_area)] = 0.0
        arrays[f"cell_ids_{i}"] = cell_ids
        arrays[f"monthly_area_{i}"] = monthly_area

    if reference is None:
        raise ValueError("Expected at least one MIRCA monthly grid")
    latitude, longitude, crs_wkt = reference
    resolution = float(abs(longitude[1] - longitude[0]))
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        labels=np.asarray(labels),
        shape=np.asarray((len(latitude), len(longitude)), dtype=np.int32),
        bounds=np.asarray(
            (
                longitude.min() - resolution / 2,
                latitude.min() - resolution / 2,
                longitude.max() + resolution / 2,
                latitude.max() + resolution / 2,
            ),
            dtype=np.float64,
        ),
        crs_wkt=np.asarray(crs_wkt),
        **arrays,
    )


if __name__ == "__main__":
    raw_inputs = {
        key.removeprefix("nc_"): str(path)
        for key, path in snakemake.input.items()  # type: ignore[name-defined]
        if key.startswith("nc_")
    }
    prepare_calendar(raw_inputs, snakemake.output[0])  # type: ignore[name-defined]
