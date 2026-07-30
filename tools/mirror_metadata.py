# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Deposition metadata for the combined GLADE input-data mirror on Zenodo.

The mirror is a single Zenodo record (pinned by
``config['data']['mirror']['zenodo_record']``) holding every third-party
input that GLADE redistributes for reproducible, key-free builds:

* the Copernicus ESA CCI land-cover extract (refreshed by
  ``tools/mirror_land_cover.py``), and
* the FAOSTAT bulk zips plus their ``faostat_vintages.yaml`` manifest
  (refreshed by ``tools/mirror_faostat.py``).

Both datasets are CC-BY-4.0, which permits redistribution provided the
attributions below are retained. Each mirror tool refreshes only its own
files (``keep_existing=True``) but regenerates the full record metadata
from this module, so the description stays identical no matter which tool
published last: everything release-specific about FAOSTAT lives in the
``faostat_vintages.yaml`` file inside the record, not in the metadata.
"""

# Source DOI of the Copernicus CDS land-cover dataset (all versions/years).
LAND_COVER_SOURCE_DOI = "10.24381/cds.006f2c9a"

FAOSTAT_URL = "https://www.fao.org/faostat/en/"

# Attribution wording required by the Copernicus product licence (CC-BY-4.0).
COPERNICUS_ATTRIBUTION = (
    "Generated using Copernicus Climate Change Service information {year}. "
    "Neither the European Commission nor ECMWF is responsible for any use that "
    "may be made of the Copernicus information or data it contains."
)

# Attribution for the FAOSTAT bulk data (CC-BY-4.0). Per-domain dataset
# names, access dates and upstream modification dates are recorded in the
# faostat_vintages.yaml manifest uploaded alongside the zips.
FAOSTAT_ATTRIBUTION = (
    "FAO. FAOSTAT statistical database. https://www.fao.org/faostat/en/. "
    "Licence: CC-BY-4.0. Access and last-update dates for each mirrored "
    "domain are recorded in faostat_vintages.yaml in this record."
)


def build_metadata(baseline_year: int, land_cover_version: str) -> dict:
    """Zenodo deposition metadata for the combined GLADE input-data mirror."""
    copernicus_attribution = COPERNICUS_ATTRIBUTION.format(year=baseline_year)
    description = (
        f"<p>Mirror of third-party input datasets for the GLADE global food "
        f"systems model, redistributed here so that builds are reproducible "
        f"(pinned data vintages) and need no upstream API keys. All contents "
        f"are provided under the Creative Commons Attribution 4.0 "
        f"International licence (CC-BY-4.0).</p>"
        f"<p><strong>Copernicus ESA CCI land cover</strong> "
        f"(<code>lccs_class</code> variable only) for {baseline_year}, version "
        f"{land_cover_version}, at 300 m global resolution, extracted from the "
        f"Copernicus Climate Data Store <code>satellite-land-cover</code> "
        f"dataset. Attribution: {copernicus_attribution} Source: Copernicus "
        f"Climate Change Service, Climate Data Store (2019): Land cover "
        f"classification gridded maps from 1992 to present derived from "
        f"satellite observation. "
        f"DOI: https://doi.org/{LAND_COVER_SOURCE_DOI}</p>"
        f"<p><strong>FAOSTAT bulk data</strong> "
        f"(<code>faostat_*.zip</code>): unmodified bulk downloads of the "
        f"FAOSTAT domains used by GLADE. FAO revises these continuously and "
        f"serves only the current release, so this mirror pins the exact "
        f"vintage the model was calibrated against; the "
        f"<code>faostat_vintages.yaml</code> manifest records the source URL, "
        f"upstream last-modified date, retrieval date and SHA-256 digest of "
        f"each zip. Attribution: FAO. FAOSTAT statistical database. "
        f"{FAOSTAT_URL} Licence: CC-BY-4.0.</p>"
    )
    return {
        "title": "GLADE input data mirror",
        "upload_type": "dataset",
        "description": description,
        "creators": [
            {"name": "Copernicus Climate Change Service (C3S)"},
            {"name": "UCLouvain"},
            {"name": "Food and Agriculture Organization of the United Nations (FAO)"},
        ],
        "license": "cc-by-4.0",
        "access_right": "open",
        "keywords": ["GLADE", "land cover", "ESA CCI", "Copernicus", "FAOSTAT", "FAO"],
        "related_identifiers": [
            {
                "identifier": LAND_COVER_SOURCE_DOI,
                "relation": "isDerivedFrom",
                "scheme": "doi",
                "resource_type": "dataset",
            },
            {
                "identifier": FAOSTAT_URL,
                "relation": "isDerivedFrom",
                "scheme": "url",
                "resource_type": "dataset",
            },
        ],
        "notes": f"{copernicus_attribution} {FAOSTAT_ATTRIBUTION}",
    }
