<!--
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: CC-BY-4.0
-->

# Manually Downloaded Data

This directory contains datasets that must be manually downloaded because they:
- Require interactive query interfaces (e.g., IHME GBD Results Tool)
- Have terms-of-service that preclude automated bulk downloads
- Require authentication or registration

## Current Files

### IHME-GBD_2021-dealth-rates.csv

**Source:** IHME Global Burden of Disease Study 2021
**Download:** https://vizhub.healthdata.org/gbd-results/

Viewing and downloading these results requires a user account on the healthdata.org website.

**Query parameters:**
- **Measure:** Deaths (Rate per 100,000)
- **Causes:**
  - Ischemic heart disease
  - Stroke
  - Diabetes mellitus
  - Colon and rectum cancer
  - Chronic respiratory diseases
  - All causes
- **Age groups:** <1 year, 12-23 months, 2-4 years, 5-9 years, ..., 95+ years
- **Sex:** Both
- **Metric:** Rate
- **Year:** 2021 (or latest available)

This specific query can also be found at the following URL: https://vizhub.healthdata.org/gbd-results?params=gbd-api-2021-permalink/90f3c59133738e4b70b91072b6fd0db4

**Processing:** The Snakemake workflow automatically processes this file via `workflow/scripts/prepare_gbd_mortality.py` to:
1. Map country names to ISO3 codes
2. Map IHME causes to model cause codes
3. Aggregate sub-buckets (12-23 months + 2-4 years → 1-4)
4. Convert rates from per 100k to per 1k
5. Output to `processing/{name}/gbd_mortality_rates.csv`

**License:** IHME Free-of-Charge Non-commercial User Agreement

**Citation:**
> Global Burden of Disease Collaborative Network. Global Burden of Disease Study 2021 (GBD 2021) Results. Seattle, United States: Institute for Health Metrics and Evaluation (IHME), 2024. Available from https://vizhub.healthdata.org/gbd-results/.

---

### GDD-dietary-intake (directory)

**Source:** Global Dietary Database (Tufts University)
**Download:** https://globaldietarydatabase.org/data-download

Downloading requires free user registration and acceptance of terms of use.

**Dataset details:**
- **Content:** Country-level mean daily dietary intake (g/day per capita) for major food groups and dietary risk factors
- **Food groups:** Vegetables, fruits (temperate/tropical/starchy), whole grains, legumes, nuts & seeds, red meat (beef/lamb/pork), processed meat, seafood (fish types + shellfish), grains, dairy, eggs, oils, and others
- **Coverage:** 185+ countries with data circa 2015-2020
- **Format:** CSV (~1.6 GB) with columns for country, food item, mean intake, standard error, and uncertainty bounds
- **Use case:** Baseline dietary patterns for health risk assessment

**Processing:** The Snakemake workflow processes this file via `workflow/scripts/prepare_gdd_dietary_intake.py` to:
1. Filter to baseline (BMK) scenario equivalent
2. Map country names to ISO3 codes
3. Map GDD food items to model dietary risk factors
4. Aggregate multiple food items to risk factor categories
5. Output to `processing/{name}/dietary_intake_baseline.csv`

**License:** Free for non-commercial research, teaching, and private study with attribution. May not be redistributed or used commercially without Tufts permission.

**Citation:**
> Global Dietary Database. Dietary intake data by country. https://www.globaldietarydatabase.org/ [Accessed YYYY-MM-DD].

**Attribution format (when publishing results):**
> Data provided by Global Dietary Database. https://www.globaldietarydatabase.org/ [Date accessed].

---

### IHME_GBD_2019_RELATIVE_RISKS_Y2020M10D15.XLSX

**Source:** IHME Global Burden of Disease Study 2019
**Download:** https://ghdx.healthdata.org/record/ihme-data/gbd-2019-relative-risks

Direct file link: https://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_GBD_2019_RELATIVE_RISKS_Y2020M10D15.XLSX

**Query parameters:**
- **Measure:** Relative Risk
- **Risk factors:** Diet-related risks (high red meat, low vegetables, low fruits, etc.)
- **Causes:** Ischemic heart disease, Stroke, Diabetes, Colon and rectum cancer, Chronic respiratory diseases
- **Age groups:** All age groups
- **Sex:** Both
- **Year:** 2019

**Processing:** The Snakemake workflow processes this file via `workflow/scripts/prepare_relative_risks.py` to:
1. Extract relative risk values for dietary risk factors
2. Map age groups to model age buckets
3. Map causes to model cause codes
4. Output to `processing/{name}/relative_risks.csv`

**License:** IHME Free-of-Charge Non-commercial User Agreement

**Citation:**
> Global Burden of Disease Collaborative Network. Global Burden of Disease Study 2019 (GBD 2019) Relative Risks. Seattle, United States: Institute for Health Metrics and Evaluation (IHME), 2020. Available from https://vizhub.healthdata.org/gbd-results/.

---

## Updating Data

### IHME GBD Mortality Data

When new GBD data is released:

1. Visit https://vizhub.healthdata.org/gbd-results/
2. Configure query with parameters above
3. Download as CSV
4. Replace `IHME-GBD_2021-dealth-rates.csv` (or create new file with updated year)
5. Update `workflow/Snakefile` rule `prepare_gbd_mortality` if filename changes
6. Rerun workflow: `tools/smk processing/{name}/gbd_mortality_rates.csv`

### GDD Dietary Data

When updating GDD data:

1. Visit https://globaldietarydatabase.org/data-download
2. Log in with user account
3. Download the complete dataset CSV
4. Replace `GDD-dietary-intake.csv` in this directory
5. Update access date in citations and documentation
6. Rerun workflow: `tools/smk processing/{name}/dietary_intake_baseline.csv`

### IHME GBD Relative Risks

When new GBD relative risks data is released:

1. Visit https://ghdx.healthdata.org/record/ihme-data/gbd-2019-relative-risks
2. Download the XLSX file (Appendix Table 7a)
3. Replace `IHME_GBD_2019_RELATIVE_RISKS_Y2020M10D15.XLSX` (or create new file with updated year)
4. Update `workflow/Snakefile` rule `prepare_relative_risks` if filename changes
5. Rerun workflow: `tools/smk processing/{name}/relative_risks.csv`
