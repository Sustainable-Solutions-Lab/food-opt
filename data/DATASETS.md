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
- Version: GAEZ v5.
- License/terms (summary): Datasets disseminated through FAO corporate statistical databases are licensed under Creative Commons Attribution 4.0 International (CC BY 4.0), complemented by FAO’s additional Statistical Database Terms of Use.
  - FAO Database Terms of Use: https://www.fao.org/contact-us/terms/db-terms-of-use/en/

## CROPGRIDS

- Description: Global geo-referenced harvested and physical crop area maps for 173 crops around 2020 at 0.05° (~5.6 km) resolution; compiled from Monfreda et al. (2008) plus 28 newer gridded sources aligned to 2020 FAOSTAT statistics.
- Website: https://figshare.com/articles/dataset/CROPGRIDS/22491997
- Version/format: v1.08 release (Figshare v9); we download the NetCDF package `CROPGRIDSv1.08_NC_maps.zip` alongside accompanying country tables.
- License/terms (summary): Creative Commons Attribution 4.0 International (CC BY 4.0); cite Tang, Nguyen, Conchedda, Casse, Tubiello & Maggi (2023), *Scientific Data*, https://doi.org/10.6084/m9.figshare.22491997.v9.
  - License: https://creativecommons.org/licenses/by/4.0/

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

- Description: Official United Nations population estimates and projections prepared by UN DESA’s Population Division. This project uses the total population table (medium variant) to obtain planning-horizon population totals by country. Additionally, we use the abridged life table for years-of-life-lost calculations.
- Website: https://population.un.org/wpp/
- Version/format: 2024 Revision; `WPP2024_TotalPopulationBySex.csv.gz` (CSV, medium variant) and `WPP2024_Life_Table_Abridged_Medium_2024-2100` (CSV, medium variant).
- License/terms (summary): UN population data is made available under the Creative Commons Attribution 3.0 IGO license (CC BY 3.0 IGO)
  - See copyright notice at the bottom of https://population.un.org/wpp/downloads

## DIA Health Impact Inputs (Diet Impact Assessment)

- Description: Epidemiological inputs used by the Diet Impact Assessment (DIA) model to translate dietary exposures into health burdens. We copy a minimal subset covering dietary risk relative-risk schedules, baseline consumption, mortality, demographic structure, and regional values of a statistical life year.
- Source repository: https://github.com/marco-spr/WHO-DIA
- Version/format: CSV snapshots dated 2021-05-28 (diet, risk schedules, demographics) and 2021-10-18 (VSL region table).
- License/terms (summary): Whole repository licensed under the GPL-3.0

## IHME GBD — Global Burden of Disease Study 2021

- Description: Cause-specific mortality rates and dietary risk relative-risk parameters from the Global Burden of Disease (GBD) studies. Death rates feed the baseline disease burden calculation; Appendix Table 7a provides dietary relative risk rates that we resample into optimization breakpoints.
- Website: https://vizhub.healthdata.org/gbd-results/ (mortality); https://ghdx.healthdata.org/record/ihme-data/gbd-2019-relative-risks (relative risks)
- Version/format: GBD 2021 death rates (CSV) and GBD 2019 Appendix Table 7a (XLSX).
- Query configuration for mortality data:
  - Measure: Deaths (Rate per 100,000)
  - Causes: Ischemic heart disease, Stroke, Diabetes mellitus, Colon and rectum cancer, Chronic respiratory diseases, All causes
  - Age groups: <1 year, 12-23 months, 2-4 years, 5-9 years, ..., 95+ years (individual bins, not aggregates)
  - Sex: Both
  - Year: 2021 (or closest available to reference year)
- Permalink (mortality): https://vizhub.healthdata.org/gbd-results?params=gbd-api-2021-permalink/8e5d55f174855a4e62a0ac13c52acf9c
- Permalink (dietary relative risks): https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_RELATIVE_RISKS_Y2020M10D15.XLSX
- License/terms (summary): Free for non-commercial use with attribution; GBD data is made available under the IHME Free-of-Charge Non-commercial User Agreement.
  - Terms: https://www.healthdata.org/data-tools-practices/data-practices/ihme-free-charge-non-commercial-user-agreement
- Citation (mortality): Global Burden of Disease Collaborative Network. Global Burden of Disease Study 2021 (GBD 2021) Results. Seattle, United States: Institute for Health Metrics and Evaluation (IHME), 2024. Available from https://vizhub.healthdata.org/gbd-results/.
- Citation (relative risks): Global Burden of Disease Collaborative Network. Global Burden of Disease Study 2019 (GBD 2019) Results. Seattle, United States of America: Institute for Health Metrics and Evaluation (IHME), 2020.

## GDD — Global Dietary Database (Tufts University)

- Description: Country-level estimates of dietary intake for major food groups and dietary risk factors based on systematic review and meta-analysis of national dietary surveys. This project uses GDD data to establish baseline dietary intake patterns by country for health risk assessment.
- Website: https://www.globaldietarydatabase.org/
- Data download: https://globaldietarydatabase.org/data-download
- Version/format: Downloaded as CSV (~1.6 GB); coverage circa 2015-2020 depending on country survey availability.
- Content: Mean daily intake (g/day per capita) for major food groups including vegetables, fruits, whole grains, legumes, nuts & seeds, red meat, processed meat, and seafood, with uncertainty estimates.
- Coverage: 185+ countries
- License/terms (summary): Free for non-commercial research, teaching, and private study with attribution. Requires user registration. Data may not be redistributed, shared with third parties, or used for commercial purposes without written permission from Tufts.
  - Terms and conditions: https://globaldietarydatabase.org/terms-and-conditions-use
- Citation: Global Dietary Database. Dietary intake data by country. https://www.globaldietarydatabase.org/.
