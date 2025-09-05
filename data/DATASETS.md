<!--
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: CC-BY-4.0
-->

# Datasets Used

Brief descriptions of key external datasets used by this project, with links and license notes. Always verify the current terms on the official sites before (re)distribution or alternative uses.

## GADM — Global Administrative Areas

- Description: Global administrative boundary data (levels ADM_0–ADM_5). This project uses level‑1 (ADM_1) regions (e.g., states/provinces).
- Website: https://gadm.org/
- Version/format: GADM 4.1; multi‑layer GeoPackage (`gadm_410-levels.gpkg`), with `ADM_1` extracted to a lighter GPKG for convenience.
- License/terms (summary): Free for academic and other non‑commercial use with attribution; redistribution of the data is not allowed; commercial use requires permission from GADM. See the official GADM license page for full terms and any updates.
  - GADM License: https://gadm.org/license.html

## GAEZ — Global Agro‑Ecological Zones (FAO/IIASA)

- Description: Global suitability and attainable yield datasets by crop and scenario. This project uses crop yield and suitability rasters (e.g., variables `yc`, `sx1`) under selected climate and management scenarios.
- Website: https://gaez.fao.org/
- Version: GAEZ v4.
- License/terms (summary): Governed by FAO’s Database Terms of Use; unless otherwise indicated in dataset metadata, licensing is Creative Commons Attribution-4.0 International licence (CC BY 4.0) with a few additional terms of use.
  - FAO Database Terms of Use: https://www.fao.org/contact-us/terms/db-terms-of-use/en/

