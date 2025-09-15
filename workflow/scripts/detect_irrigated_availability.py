"""
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: GPL-3.0-or-later
"""

import pandas as pd
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

    rows: list[dict] = []
    types = ["i", "g", "s", "d"]  # irrigation, gravity, sprinkler, drip (as available)
    for crop in crops:
        code = gaez["crops"][crop]
        available: dict[str, int] = {t: 0 for t in types}
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
            if y_ok and s_ok:
                available[t] = 1

        # First available in the order defined above, or "none" if none
        first = next((t for t in types if available[t]), "none")
        rows.append(
            {
                "crop": crop,
                "code": code,
                **available,
                "first_available": first,
            }
        )

    out_path = Path(snakemake.output[0])  # type: ignore[name-defined]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False)
