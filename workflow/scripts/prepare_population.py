# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

# SPDX-License-Identifier: GPL-3.0-or-later

"""
Download and process UN WPP 2024 population totals by sex (CSV.gz).

Uses the official TotalPopulationBySex dataset with ISO3_code included, filters to
Variant == Medium and Sex == Both sexes, selects the planning horizon year from
Snakemake params, and outputs iso3, country, year, population (persons).
"""

import pandas as pd

if __name__ == "__main__":
    # Read the gzipped CSV from Snakemake input
    df = pd.read_csv(snakemake.input.population_gz, compression="gzip")

    # Filter to Medium variant if present
    if "Variant" in df.columns:
        df = df[df["Variant"].astype(str).str.lower() == "medium"]

    # Filter to planning horizon year (PopTotal is total across sexes)
    year = int(snakemake.params.planning_horizon)
    df = df[df["Time"] == year]

    # Keep rows with ISO3_code only (drop aggregates like regions without ISO3_code)
    df = df[df["ISO3_code"].notna()]

    out = df.loc[:, ["ISO3_code", "Location", "Time", "PopTotal"]].rename(
        columns={
            "ISO3_code": "iso3",
            "Location": "country",
            "Time": "year",
            "PopTotal": "population_thousands",
        }
    )

    # Convert from thousands to persons
    out["population"] = out["population_thousands"].astype(float) * 1e3
    out = out.drop(columns="population_thousands")

    out.to_csv(snakemake.output[0], index=False)
