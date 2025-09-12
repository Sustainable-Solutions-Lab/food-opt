"""
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: GPL-3.0-or-later
"""

import csv
import subprocess
from pathlib import Path


def url_ok(url: str) -> bool:
    """Return True if URL appears downloadable (HTTP 2xx).

    Use wget in spider mode to avoid saving bodies and to ensure non-2xx
    responses (e.g. 403/404 from S3) yield a non-zero exit code.
    """
    res = subprocess.run(
        [
            "wget",
            "-q",
            "--spider",
            "--tries=1",
            "--timeout=10",
            url,
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return res.returncode == 0


if __name__ == "__main__":
    gaez = snakemake.params.gaez  # type: ignore[name-defined]
    crops = snakemake.params.crops  # type: ignore[name-defined]
    require_suit = snakemake.params.require_suitability  # type: ignore[name-defined]

    rows: list[dict] = []
    types = ["g", "s", "d"]  # gravity, sprinkler, drip (as available)
    for crop in crops:
        code = gaez["crops"][crop]
        avail_types: list[str] = []
        for t in types:
            url_y = (
                f"https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/res05/"
                f"{gaez['climate_model']}/rcp{gaez['rcp']}/{gaez['time_period']}H/"
                f"{gaez['yield_var']}{gaez['input_management']}{t}{gaez['co2_fertilization']}_{code}.tif"
            )
            url_s = (
                f"https://s3.eu-west-1.amazonaws.com/data.gaezdev.aws.fao.org/res05/"
                f"{gaez['climate_model']}/rcp{gaez['rcp']}/{gaez['time_period']}H/"
                f"{gaez['suitability_var']}{gaez['input_management']}{t}{gaez['co2_fertilization']}_{code}.tif"
            )
            y_ok = url_ok(url_y)
            s_ok = url_ok(url_s)
            if y_ok and (s_ok or not require_suit):
                avail_types.append(t)

        rows.append({"crop": crop, "code": code, "irrig_types": ";".join(avail_types)})

    out_path = Path(snakemake.output[0])  # type: ignore[name-defined]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["crop", "code", "irrig_types"])
        writer.writeheader()
        writer.writerows(rows)
