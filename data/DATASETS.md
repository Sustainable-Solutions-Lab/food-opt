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
- License/terms (summary): Datasets disseminated through FAO corporate statistical databases are licensed under Creative Commons Attribution 4.0 International (CC BY 4.0), complemented by FAO’s additional Statistical Database Terms of Use.
  - FAO Database Terms of Use: https://www.fao.org/contact-us/terms/db-terms-of-use/en/

## FAOSTAT — FAO Statistics Division

- Description: FAO’s global statistical database covering food and agriculture domains for 245+ countries and territories from 1961 onward. This project currently uses the Producer Prices (PP) domain to obtain 2015–2024 USD/tonne crop producer prices for cost calibration.
- Website: https://www.fao.org/faostat/en/
- Version/format: Producer Prices (PP) domain via FAOSTAT API, retrieved as CSV through the `faostat` Python client.
- License/terms (summary): Datasets disseminated through FAO corporate statistical databases are licensed under Creative Commons Attribution 4.0 International (CC BY 4.0), complemented by FAO’s additional Statistical Database Terms of Use.
  - FAO Statistical Database Terms of Use: https://www.fao.org/contact-us/terms/db-terms-of-use/en/

## Water Footprint Network — Monthly Blue Water Availability

- Description: Monthly blue water availability for 405 GRDC river basins, provided alongside blue-water scarcity indicators as part of the Water Footprint Network’s Appendix to Value of Water Research Report Series No. 53.
- Download: https://www.waterfootprint.org/resources/appendix/Report53_Appendix.zip
- Version/format: Appendix VII of Hoekstra & Mekonnen (2011); data distributed as an ESRI shapefile (`Monthly_WS_GRDC_405_basins.*`) with basin metadata, plus an Excel workbook (`Report53-Appendices-VI-IX.xls`, sheet “Appendix-VII”) containing monthly availability in Mm³/month.
- License/terms (summary): No explicit license accompanies the dataset. The authors request citation as below; users should evaluate whether their use qualifies as fair use (research is probably allowed) and contact the UNESCO-IHE Institute for Water Education for commercial applications.
- Suggested citation (from the dataset readme): Hoekstra, A.Y. and Mekonnen, M.M. (2011) *Global water scarcity: monthly blue water footprint compared to blue water availability for the world’s major river basins*, Value of Water Research Report Series No. 53, UNESCO-IHE, Delft, the Netherlands. http://www.waterfootprint.org/Reports/Report53-GlobalBlueWaterScarcity.pdf

## UN WPP — World Population Prospects 2024 (UN DESA)

- Description: Official United Nations population estimates and projections prepared by UN DESA’s Population Division. This project uses the TotalPopulationBySex table (Medium variant) to obtain planning-horizon population totals by country and convert them from thousands of persons to persons.
- Website: https://population.un.org/wpp/
- Version/format: 2024 Revision; `WPP2024_TotalPopulationBySex.csv.gz` (CSV, Medium variant) downloaded from the WPP Download Center.
- License/terms (summary): UN population data is made available under the Creative Commons Attribution 3.0 IGO license (CC BY 3.0 IGO)
  - See copyright notice at the bottom of https://population.un.org/wpp/downloads

## DIA Health Impact Inputs (Diet Impact Assessment)

- Description: Epidemiological inputs used by the Diet Impact Assessment (DIA) model to translate dietary exposures into health burdens. We copy a minimal subset covering dietary risk relative-risk schedules, baseline consumption, mortality, demographic structure, life tables, and regional values of a statistical life year.
- Source repository: https://github.com/marco-spr/WHO-DIA
- Version/format: CSV snapshots dated 2021-05-28 (diet, risk schedules, demographics) and 2021-10-18 (VSL region table).
- License/terms (summary): Whole repository licensed under the GPL-3.0
