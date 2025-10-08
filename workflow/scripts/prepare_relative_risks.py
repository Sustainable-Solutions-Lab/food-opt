# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Convert manually downloaded IHME GBD RR tables into tidy dietary risk curves.

The omega-3 risk factor measures EPA+DHA intake in g/day, but the model tracks fish
consumption. The conversion uses: g_fish = (100 / omega3_per_100g_fish) * g_omega3.
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


# Map IHME dietary risk names to model risk_factor identifiers and exposure conversion factors
RISK_CONFIG = {
    "Diet low in fruits": {"risk_factor": "fruits", "unit": "g/day", "conversion": 1.0},
    "Diet low in vegetables": {
        "risk_factor": "vegetables",
        "unit": "g/day",
        "conversion": 1.0,
    },
    "Diet low in whole grains": {
        "risk_factor": "whole_grains",
        "unit": "g/day",
        "conversion": 1.0,
    },
    "Diet low in legumes": {
        "risk_factor": "legumes",
        "unit": "g/day",
        "conversion": 1.0,
    },
    "Diet low in nuts and seeds": {
        "risk_factor": "nuts_seeds",
        "unit": "g/day",
        "conversion": 1.0,
    },
    "Diet low in seafood omega-3 fatty acids": {
        "risk_factor": "fish",
        "unit": "g/day",
        "conversion": None,
    },
    "Diet high in red meat": {
        "risk_factor": "red_meat",
        "unit": "g/day",
        "conversion": 1.0,
    },
    "Diet high in processed meat": {
        "risk_factor": "prc_meat",
        "unit": "g/day",
        "conversion": 1.0,
    },
}


# Map IHME outcome names to model causes. Any unmapped outcome is ignored.
CAUSE_MAP = {
    "Ischemic heart disease": "CHD",
    "Ischemic stroke": "Stroke",
    "Intracerebral hemorrhage": "Stroke",
    "Subarachnoid hemorrhage": "Stroke",
    "Diabetes mellitus type 2": "T2DM",
    "Colon and rectum cancer": "CRC",
}


VALUE_REGEX = re.compile(r"[-+]?(?:\d+\.\d+|\d+)")


def _parse_rr_value(cell: object) -> Tuple[float, float | None, float | None]:
    """Return (mean, low, high) RR floats parsed from a string cell."""

    if isinstance(cell, (int, float)):
        value = float(cell)
        return value, None, None

    text = str(cell).strip()
    if not text:
        raise ValueError("Empty RR cell")

    matches = VALUE_REGEX.findall(text)
    if not matches:
        raise ValueError(f"Could not parse RR value from '{text}'")

    numbers = [float(v) for v in matches]
    mean = numbers[0]
    low = numbers[1] if len(numbers) > 1 else None
    high = numbers[2] if len(numbers) > 2 else None
    return mean, low, high


def _normalize_exposure(raw: str, conversion: float | None) -> float:
    """Convert exposure text like '100 g/day' into g/day as float."""

    parts = raw.strip().split()
    if len(parts) < 2:
        raise ValueError(f"Unexpected exposure label '{raw}'")

    value = float(parts[0])
    unit = parts[1].lower()

    if unit not in {"g/day", "%energy/day"}:
        raise ValueError(f"Unsupported exposure unit '{unit}' in '{raw}'")

    if unit == "%energy/day":
        raise ValueError(
            "Energy-based exposures are not supported in the current health module"
        )

    if conversion is None:
        raise ValueError("Missing conversion factor for omega-3 exposure")

    return value * conversion


def _extract_risk_blocks(df: pd.DataFrame) -> Dict[str, Tuple[int, int]]:
    """Return mapping from IHME risk name row index to slice bounds."""

    risk_rows = [
        idx for idx, val in df[0].items() if isinstance(val, str) and val in RISK_CONFIG
    ]
    bounds: Dict[str, Tuple[int, int]] = {}
    for i, idx in enumerate(risk_rows):
        next_idx = risk_rows[i + 1] if i + 1 < len(risk_rows) else len(df)
        bounds[df.at[idx, 0]] = (idx + 1, next_idx)
    return bounds


def _parse_relative_risks(
    df: pd.DataFrame,
    omega3_conversion: float,
) -> pd.DataFrame:
    """Parse the Excel sheet into tidy RR records."""

    records: List[Dict[str, float | str]] = []
    skipped_causes: Dict[str, set[str]] = {}
    skipped_units: set[str] = set()
    blocks = _extract_risk_blocks(df)

    for risk_name, (start, end) in blocks.items():
        config = RISK_CONFIG[risk_name]
        risk_id = config["risk_factor"]
        conversion = config["conversion"]

        if risk_id == "fish":
            if omega3_conversion <= 0:
                raise ValueError("omega3_conversion must be positive")
            conversion = omega3_conversion

        block = df.iloc[start:end]
        block = block[block[0].notna()]

        for _, row in block.iterrows():
            outcome = str(row[0]).strip()
            exposure_raw = row[1]

            if not isinstance(exposure_raw, str):
                # Skip cases without quantitative exposure (e.g., "Exposed")
                continue

            if outcome not in CAUSE_MAP:
                skipped_causes.setdefault(risk_id, set()).add(outcome)
                continue

            cause = CAUSE_MAP[outcome]

            try:
                exposure = _normalize_exposure(exposure_raw, conversion)
            except ValueError:
                skipped_units.add(exposure_raw)
                continue

            # Grab the first non-null RR cell from adult ages (prefer all-age, then 25-29+)
            rr_cell = None
            for col_idx in range(4, len(row)):
                value = row[col_idx]
                if isinstance(value, str) and value.strip():
                    rr_cell = value
                    break

            if rr_cell is None:
                continue

            try:
                rr_mean, rr_low, rr_high = _parse_rr_value(rr_cell)
            except ValueError:
                continue

            records.append(
                {
                    "risk_factor": risk_id,
                    "cause": cause,
                    "exposure_g_per_day": float(exposure),
                    "rr_mean": rr_mean,
                    "rr_low": rr_low,
                    "rr_high": rr_high,
                }
            )

    if skipped_causes:
        logger.info("Skipped unmapped outcomes:")
        for risk_id, causes in sorted(skipped_causes.items()):
            items = ", ".join(sorted(causes))
            logger.info(f"  {risk_id}: {items}")

    if skipped_units:
        logger.info("Skipped exposures with unsupported units:")
        for label in sorted(skipped_units):
            logger.info(f"  {label}")

    if not records:
        raise ValueError("No dietary risk records parsed from GBD relative risk file")

    df_out = pd.DataFrame(records)
    df_out = (
        df_out.groupby(["risk_factor", "cause", "exposure_g_per_day"], as_index=False)
        .agg({"rr_mean": "mean", "rr_low": "mean", "rr_high": "mean"})
        .sort_values(["risk_factor", "cause", "exposure_g_per_day"])
    )
    return df_out


def main() -> None:
    snakemake = globals().get("snakemake")  # type: ignore
    if snakemake is None:
        raise RuntimeError("This script must run via Snakemake")

    input_path = Path(snakemake.input["gbd_rr"])
    output_path = Path(snakemake.output["relative_risks"])
    omega3_per_100g = float(snakemake.params["omega3_per_100g"])

    # Convert g omega-3 per 100 g fish to conversion factor g_fish per g omega-3
    if omega3_per_100g <= 0:
        raise ValueError("omega3_per_100g must be positive")
    omega3_conversion = 100.0 / omega3_per_100g

    logger.info(f"Reading {input_path}")
    df = pd.read_excel(input_path, header=None)

    relative_risks = _parse_relative_risks(df, omega3_conversion)

    # Validate that we have all required risk factors and causes
    required_risk_factors = set(snakemake.params["risk_factors"])
    required_causes = set(snakemake.params["causes"])
    output_risk_factors = set(relative_risks["risk_factor"].unique())
    output_causes = set(relative_risks["cause"].unique())

    missing_risk_factors = required_risk_factors - output_risk_factors
    if missing_risk_factors:
        raise ValueError(
            f"[prepare_relative_risks] ERROR: Relative risks data is missing {len(missing_risk_factors)} required risk factors: "
            f"{sorted(missing_risk_factors)}. Available: {sorted(output_risk_factors)}. "
            f"Please ensure the IHME GBD relative risks file includes all risk factors listed in config.health.risk_factors."
        )

    missing_causes = required_causes - output_causes
    if missing_causes:
        raise ValueError(
            f"[prepare_relative_risks] ERROR: Relative risks data is missing {len(missing_causes)} required causes: "
            f"{sorted(missing_causes)}. Available: {sorted(output_causes)}. "
            f"Please ensure the IHME GBD relative risks file includes all causes listed in config.health.causes."
        )

    logger.info(
        "[prepare_relative_risks] ✓ Validation passed: all required risk factors and causes present"
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    relative_risks.to_csv(output_path, index=False)
    logger.info(f"Wrote {len(relative_risks)} records to {output_path}")


if __name__ == "__main__":
    main()
