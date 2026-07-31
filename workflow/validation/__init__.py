# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Validation entry points for configuration and data consistency checks."""

from pathlib import Path
from typing import Callable

from snakemake.logging import logger

from .calibration import validate_calibration
from .calibration_provenance import validate_calibration_provenance
from .commodities import validate_commodities
from .config_schema import validate_config_schema
from .country_regions import validate_country_regions
from .crop_food_pathways import validate_crop_food_pathways
from .crop_groups import validate_crop_groups
from .crop_moisture_content import validate_crop_moisture_content
from .cropgrids_crops import validate_cropgrids_crops
from .diet_basis import validate_diet_basis
from .faostat_maps import validate_faostat_maps
from .food_basis import validate_food_basis
from .food_groups import validate_food_groups
from .gaez_crop_mapping import validate_gaez_crop_mapping
from .health_map import validate_health_map
from .m49_codes import validate_m49_codes
from .multi_cropping import validate_multi_cropping
from .nutrition import validate_nutrition
from .secrets import load_secrets_with_env_fallback
from .seed_rates import validate_seed_rates
from .sensitivity_generator import validate_sensitivity_generator
from .yield_unit_conversions import validate_yield_unit_conversions

Validator = Callable[[dict, Path], None]

_SEMANTIC_CHECKS: tuple[tuple[str, Validator], ...] = (
    ("calibration", validate_calibration),
    ("calibration_provenance", validate_calibration_provenance),
    ("commodities", validate_commodities),
    ("country_regions", validate_country_regions),
    ("food_groups", validate_food_groups),
    ("food_basis", validate_food_basis),
    ("diet_basis", validate_diet_basis),
    ("faostat_maps", validate_faostat_maps),
    ("crop_food_pathways", validate_crop_food_pathways),
    ("crop_groups", validate_crop_groups),
    ("crop_moisture_content", validate_crop_moisture_content),
    ("cropgrids_crops", validate_cropgrids_crops),
    ("gaez_crop_mapping", validate_gaez_crop_mapping),
    ("seed_rates", validate_seed_rates),
    ("health_map", validate_health_map),
    ("m49_codes", validate_m49_codes),
    ("multi_cropping", validate_multi_cropping),
    ("nutrition", validate_nutrition),
    ("sensitivity_generator", validate_sensitivity_generator),
    ("yield_unit_conversions", validate_yield_unit_conversions),
)


def validate(
    config: dict,
    project_root: Path | None = None,
) -> None:
    """Run configured validation checks against the active config and data.

    Parameters
    ----------
    config:
        The merged Snakemake configuration dictionary.
    project_root:
        Root directory of the repository. Defaults to the current working directory.
    """
    logger.info("Validating configuration and input datasets")

    root = Path(project_root) if project_root else Path.cwd()
    validate_config_schema(config, root)

    errors: list[str] = []
    for name, check in _SEMANTIC_CHECKS:
        try:
            check(config, root)
        except Exception as exc:
            errors.append(f"{name}: {exc}")

    if errors:
        bullet_list = "\n".join(f" - {msg}" for msg in errors)
        raise RuntimeError(f"Validation failed:\n{bullet_list}")


__all__ = ["load_secrets_with_env_fallback", "validate"]
