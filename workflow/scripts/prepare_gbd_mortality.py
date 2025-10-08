# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Process IHME GBD mortality data to replace DIA death rates.

Downloads from IHME GBD Results Tool should include:
- Measure: Deaths (Rate)
- Causes: IHD, Stroke, Diabetes, CRC, Chronic respiratory diseases
- Ages: <1 year, 12-23 months, 2-4 years, 5-9 years, ..., 95+ years
- Sex: Both
- Year: 2019 or 2021
- Metric: Rate (per 100,000)
"""

from pathlib import Path

import pandas as pd
import pycountry


# Manual overrides for country names that pycountry can't match
COUNTRY_NAME_OVERRIDES = {
    "Bolivia (Plurinational State of)": "BOL",
    "Bonaire, Saint Eustatius and Saba": "BES",
    "Cabo Verde": "CPV",
    "Côte d'Ivoire": "CIV",
    "Democratic People's Republic of Korea": "PRK",
    "Democratic Republic of the Congo": "COD",
    "French Guiana": "GUF",  # Use French data for French Guiana
    "Iran (Islamic Republic of)": "IRN",
    "Lao People's Democratic Republic": "LAO",
    "Micronesia (Federated States of)": "FSM",
    "Niger": "NER",  # pycountry fuzzy search confuses with Nigeria (NGA)
    "Republic of Korea": "KOR",
    "Republic of Moldova": "MDA",
    "Republic of the Congo": "COG",
    "Saint Barthélemy": "BLM",
    "Saint Martin (French part)": "MAF",
    "Sint Maarten (Dutch part)": "SXM",
    "The former Yugoslav Republic of Macedonia": "MKD",
    "Türkiye": "TUR",
    "United Kingdom of Great Britain and Northern Ireland": "GBR",
    "United Republic of Tanzania": "TZA",
    "United States of America": "USA",
    "Venezuela (Bolivarian Republic of)": "VEN",
    "Viet Nam": "VNM",
}

# Map IHME cause names to model cause codes
CAUSE_MAP = {
    "Ischemic heart disease": "CHD",
    "Stroke": "Stroke",
    "Diabetes mellitus": "T2DM",
    "Colon and rectum cancer": "CRC",
    "Chronic respiratory diseases": "Resp_Dis",
    "All causes": "all-c",
}

# Map IHME age group names to model age codes
AGE_MAP = {
    "<1 year": "<1",
    "12-23 months": "1-4",  # Map to 1-4 bucket
    "2-4 years": "1-4",  # Map to 1-4 bucket
    "5-9 years": "5-9",
    "10-14 years": "10-14",
    "15-19 years": "15-19",
    "20-24 years": "20-24",
    "25-29 years": "25-29",
    "30-34 years": "30-34",
    "35-39 years": "35-39",
    "40-44 years": "40-44",
    "45-49 years": "45-49",
    "50-54 years": "50-54",
    "55-59 years": "55-59",
    "60-64 years": "60-64",
    "65-69 years": "65-69",
    "70-74 years": "70-74",
    "75-79 years": "75-79",
    "80-84 years": "80-84",
    "85-89 years": "85-89",
    "90-94 years": "90-94",
    "95+ years": "95+",
    "95 plus": "95+",
    "All ages": "all-a",
}


def map_country_name_to_iso3(name: str) -> str | None:
    """Map IHME country name to ISO3 code using pycountry + manual overrides."""
    # Check manual overrides first
    if name in COUNTRY_NAME_OVERRIDES:
        return COUNTRY_NAME_OVERRIDES[name]

    # Try pycountry fuzzy search
    try:
        matches = pycountry.countries.search_fuzzy(name)
        if matches:
            return matches[0].alpha_3
    except LookupError:
        pass

    return None


def main() -> None:
    snakemake = globals().get("snakemake")  # type: ignore
    if snakemake is None:
        raise RuntimeError("This script must run via Snakemake")

    input_path = snakemake.input["gbd_mortality"]
    output_path = snakemake.output["mortality"]
    reference_year = int(snakemake.params["reference_year"])
    set(snakemake.params["countries"])
    set(snakemake.params["causes"])

    print(f"[prepare_gbd_mortality] Reading GBD mortality data from {input_path}")
    df = pd.read_csv(input_path)

    # Validate expected columns
    required_cols = {
        "location_name",
        "age_name",
        "cause_name",
        "year",
        "val",
        "metric_name",
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in GBD data: {missing}")

    # Filter to death rate metric
    df = df[df["metric_name"] == "Rate"].copy()
    if df.empty:
        raise ValueError("No 'Rate' metric found in GBD data")

    # Check available years
    available_years = sorted(df["year"].unique())
    print(f"[prepare_gbd_mortality] Available years: {available_years}")

    if reference_year not in available_years:
        closest_year = min(available_years, key=lambda y: abs(y - reference_year))
        print(
            f"[prepare_gbd_mortality] Reference year {reference_year} not found; "
            f"using closest year {closest_year}"
        )
        target_year = closest_year
    else:
        target_year = reference_year

    df = df[df["year"] == target_year].copy()

    # Map country names to ISO3 (cache unique values to avoid repeated fuzzy searches)
    print("[prepare_gbd_mortality] Mapping country names to ISO3 codes...")
    unique_countries = df["location_name"].unique()
    country_map = {name: map_country_name_to_iso3(name) for name in unique_countries}
    df["country_iso3"] = df["location_name"].map(country_map)

    unmapped = df[df["country_iso3"].isna()]["location_name"].unique()
    if len(unmapped) > 0:
        print(
            f"[prepare_gbd_mortality] Warning: {len(unmapped)} countries could not be "
            "mapped to ISO3:"
        )
        for name in sorted(unmapped)[:10]:
            print(f"  - {name}")
        if len(unmapped) > 10:
            print(f"  ... and {len(unmapped) - 10} more")

    # Drop unmapped countries
    df = df[df["country_iso3"].notna()].copy()

    # Map causes
    df["cause_code"] = df["cause_name"].map(CAUSE_MAP)
    unmapped_causes = df[df["cause_code"].isna()]["cause_name"].unique()
    if len(unmapped_causes) > 0:
        print(
            f"[prepare_gbd_mortality] Warning: {len(unmapped_causes)} causes not mapped:"
        )
        for cause in sorted(unmapped_causes):
            print(f"  - {cause}")

    df = df[df["cause_code"].notna()].copy()

    # Map age groups
    df["age_code"] = df["age_name"].map(AGE_MAP)
    unmapped_ages = df[df["age_code"].isna()]["age_name"].unique()
    if len(unmapped_ages) > 0:
        print(
            f"[prepare_gbd_mortality] Warning: {len(unmapped_ages)} age groups not mapped:"
        )
        for age in sorted(unmapped_ages):
            print(f"  - {age}")

    df = df[df["age_code"].notna()].copy()

    # Check for multiple IHME age buckets mapping to same model bucket (e.g., 1-4)
    # We need to aggregate them by taking population-weighted average
    # For death rates, we'll compute a simple average as an approximation
    # (ideally would weight by age-specific population, but we don't have that here)
    print("[prepare_gbd_mortality] Checking for age bucket aggregation needs...")
    duplicate_keys = df.groupby(["country_iso3", "cause_code", "age_code"]).size()
    needs_aggregation = duplicate_keys[duplicate_keys > 1]

    if len(needs_aggregation) > 0:
        print(
            f"[prepare_gbd_mortality] Found {len(needs_aggregation)} age buckets "
            "needing aggregation (multiple IHME ages -> single model age)"
        )
        print(
            "[prepare_gbd_mortality] Aggregating by simple average "
            "(assumes roughly equal population in sub-buckets)"
        )

        # Group and average
        df = (
            df.groupby(["country_iso3", "cause_code", "age_code", "year"])
            .agg({"val": "mean"})
            .reset_index()
        )

    # Convert rate from per 100,000 to per 1,000
    df["death_rate_per_1000"] = df["val"] / 100.0

    # Build output dataframe matching old format: age, cause, country, year, value
    output = df[
        ["age_code", "cause_code", "country_iso3", "death_rate_per_1000"]
    ].copy()
    output.columns = ["age", "cause", "country", "value"]
    # Use reference year in output to ensure consistency across datasets
    output.insert(3, "year", reference_year)

    # Sort for readability
    output = output.sort_values(["country", "cause", "age"]).reset_index(drop=True)

    # Fill in missing countries using proxy data from similar countries
    # This is for territories/dependencies that don't have separate IHME data
    COUNTRY_PROXIES = {
        "ASM": "WSM",  # American Samoa -> Samoa (if needed)
        "GUF": "FRA",  # French Guiana -> France
        "PRI": "USA",  # Puerto Rico -> USA (if needed)
        "SOM": "ETH",  # Somalia -> Ethiopia (if needed)
    }

    required_countries = set(snakemake.params["countries"])
    output_countries = set(output["country"].unique())
    missing_countries = required_countries - output_countries

    if missing_countries:
        filled = []
        still_missing = []
        for missing in sorted(missing_countries):
            if missing in COUNTRY_PROXIES:
                proxy = COUNTRY_PROXIES[missing]
                if proxy in output_countries:
                    # Duplicate proxy country's data for the missing country
                    proxy_data = output[output["country"] == proxy].copy()
                    proxy_data["country"] = missing
                    output = pd.concat([output, proxy_data], ignore_index=True)
                    filled.append(f"{missing} (using {proxy} data)")
                else:
                    still_missing.append(missing)
            else:
                still_missing.append(missing)

        if filled:
            print(
                f"[prepare_gbd_mortality] Filled {len(filled)} missing countries using proxies:"
            )
            for entry in filled:
                print(f"  - {entry}")

        # Update missing list after filling
        output_countries = set(output["country"].unique())
        missing_countries = required_countries - output_countries

    # Validate that we have all required countries and causes
    required_causes = set(snakemake.params["causes"])
    output_causes = set(output["cause"].unique())
    if missing_countries:
        raise ValueError(
            f"[prepare_gbd_mortality] ERROR: Mortality data is missing {len(missing_countries)} required countries: "
            f"{sorted(list(missing_countries))[:20]}{'...' if len(missing_countries) > 20 else ''}. "
            f"Please ensure the IHME GBD download includes all countries listed in config."
        )

    missing_causes = required_causes - output_causes
    if missing_causes:
        raise ValueError(
            f"[prepare_gbd_mortality] ERROR: Mortality data is missing {len(missing_causes)} required causes: "
            f"{sorted(missing_causes)}. Available causes: {sorted(output_causes)}. "
            f"Please ensure the IHME GBD download includes all causes listed in config.health.causes."
        )

    print(
        "[prepare_gbd_mortality] ✓ Validation passed: all required countries and causes present"
    )

    # Ensure output directory exists
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # Save without header to match old format
    output.to_csv(output_path, index=False, header=False)

    print(f"[prepare_gbd_mortality] Wrote {len(output)} rows to {output_path}")
    print(f"[prepare_gbd_mortality] Countries: {output['country'].nunique()}")
    print(f"[prepare_gbd_mortality] Causes: {sorted(output['cause'].unique())}")
    print(f"[prepare_gbd_mortality] Age groups: {sorted(output['age'].unique())}")


if __name__ == "__main__":
    main()
