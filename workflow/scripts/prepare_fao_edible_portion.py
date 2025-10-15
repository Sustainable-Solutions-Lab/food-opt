# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Extract edible portion coefficients from the FAO nutrient table.

The script reads sheet "03" of the FAO Nutrient Conversion Table workbook and
pulls the *edible portion coefficient* for each FAOSTAT item mapped to the
model's configured crops. The resulting CSV lists the coefficient alongside the
FAOSTAT item code and edible portion type flag supplied by FAO.

Special handling is applied to certain crops where FAO's coefficient does not
match the model's yield units:
- Grains (rice, barley, oat, buckwheat): FAO gives milled/hulled conversion;
  we force to 1.0 for whole grain, handling milling separately.
- Oil crops (rapeseed, olive): GAEZ yields are already in kg oil/ha, so we
  force to 1.0 (no further conversion needed).
- Sugar crops (sugarcane, sugarbeet): GAEZ yields are already in kg sugar/ha,
  so we force to 1.0 (no further conversion needed).
"""

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from openpyxl import load_workbook


logger = logging.getLogger(__name__)

ALTERNATE_ITEM_NAMES: Dict[str, str] = {
    "rape or colza seed": "rapeseed or colza seed",
    "oil palm fruit": "palm oil",
}

# Crops for which the edible portion coefficient should be set to 1.0
# despite FAO listing a lower value. Reasons vary by crop type:
# - Grains (rice, barley, oat, buckwheat): FAO coefficient represents milled/
#   hulled grain; we track whole grain and handle milling separately.
# - Oil crops (rapeseed, olive): FAO coefficient represents fruit/seed→oil
#   extraction; GAEZ yields are already in kg oil/ha, so no further conversion.
# - Sugar crops (sugarcane, sugarbeet): GAEZ yields are already in kg sugar
#   per hectare, so no further conversion needed.
EDIBLE_PORTION_EXCEPTIONS: set[str] = {
    "dryland-rice",
    "wetland-rice",
    "barley",
    "oat",
    "buckwheat",
    "rapeseed",
    "olive",
    "sugarcane",
    "sugarbeet",
}


@dataclass
class ComponentRow:
    code: str
    description: str
    edible_coefficient: Optional[float]
    edible_type: Optional[int]
    water_content_g_per_100g: Optional[float]


def _coerce_float(value) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        logger.warning("Could not parse edible coefficient value %r", value)
        return None


def _coerce_int(value) -> Optional[int]:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        logger.warning("Could not parse edible portion type value %r", value)
        return None


def _load_component_values(xlsx_path: Path) -> Dict[str, ComponentRow]:
    wb = load_workbook(filename=xlsx_path, read_only=True, data_only=True)
    try:
        ws = wb["03"]
    except KeyError as exc:
        raise ValueError("Worksheet '03' with component values not found") from exc

    records: Dict[str, ComponentRow] = {}
    for row in ws.iter_rows(min_row=6, values_only=True):
        code = row[0]
        description = row[1]
        if not code or not description:
            continue

        edible_coefficient = _coerce_float(row[4])
        edible_type = _coerce_int(row[5])
        # Column index 8 corresponds to the WATER (g/100g) field in sheet "03".
        water_content = _coerce_float(row[8])

        key = str(description).strip().lower()
        records[key] = ComponentRow(
            code=str(code).strip(),
            description=str(description).strip(),
            edible_coefficient=edible_coefficient,
            edible_type=edible_type,
            water_content_g_per_100g=water_content,
        )

    return records


def _read_crop_mapping(path: Path) -> Dict[str, Optional[str]]:
    mapping: Dict[str, Optional[str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            crop = row.get("crop")
            item = row.get("faostat_item")
            if crop is None:
                continue
            item = item.strip() if item else None
            mapping[crop] = item if item else None
    return mapping


def main() -> None:
    crops: List[str] = list(snakemake.params.crops)  # type: ignore[name-defined]
    xlsx_path = Path(snakemake.input.table)  # type: ignore[name-defined]
    mapping_path = Path(snakemake.input.mapping)  # type: ignore[name-defined]
    output_path = Path(snakemake.output.edible_portion)  # type: ignore[name-defined]

    records_by_name = _load_component_values(xlsx_path)
    crop_to_item = _read_crop_mapping(mapping_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    missing_items: List[str] = []
    missing_components: List[str] = []

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "crop",
                "faostat_item",
                "fao_code",
                "edible_portion_coefficient",
                "edible_portion_type",
                "water_content_g_per_100g",
            ],
        )
        writer.writeheader()

        for crop in crops:
            item = crop_to_item.get(crop)
            if item is None:
                missing_items.append(crop)
                writer.writerow(
                    {
                        "crop": crop,
                        "faostat_item": "",
                        "fao_code": "",
                        "edible_portion_coefficient": "",
                        "edible_portion_type": "",
                        "water_content_g_per_100g": "",
                    }
                )
                continue

            key = item.strip().lower()
            record = records_by_name.get(key)
            if record is None and key in ALTERNATE_ITEM_NAMES:
                record = records_by_name.get(ALTERNATE_ITEM_NAMES[key])
            if record is None:
                missing_components.append(item)
                writer.writerow(
                    {
                        "crop": crop,
                        "faostat_item": item,
                        "fao_code": "",
                        "edible_portion_coefficient": "",
                        "edible_portion_type": "",
                        "water_content_g_per_100g": "",
                    }
                )
                continue

            # Apply exception: force edible_coefficient to 1.0 for certain crops
            edible_coeff = (
                1.0
                if crop in EDIBLE_PORTION_EXCEPTIONS
                else (
                    ""
                    if record.edible_coefficient is None
                    else record.edible_coefficient
                )
            )

            writer.writerow(
                {
                    "crop": crop,
                    "faostat_item": record.description,
                    "fao_code": record.code,
                    "edible_portion_coefficient": edible_coeff,
                    "edible_portion_type": (
                        "" if record.edible_type is None else record.edible_type
                    ),
                    "water_content_g_per_100g": (
                        ""
                        if record.water_content_g_per_100g is None
                        else record.water_content_g_per_100g
                    ),
                }
            )

    if missing_items:
        logger.warning(
            "No FAOSTAT item mapping for crops: %s", ", ".join(sorted(missing_items))
        )
    if missing_components:
        logger.warning(
            "No component entry found in FAO table for items: %s",
            ", ".join(sorted(set(missing_components))),
        )


if __name__ == "__main__":
    main()
