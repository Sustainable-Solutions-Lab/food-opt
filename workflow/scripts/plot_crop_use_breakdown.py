# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Create a stacked bar chart of global crop uses (human vs. animal feed)."""

import logging
from pathlib import Path

import matplotlib
import pandas as pd
import pypsa

matplotlib.use("pdf")
import matplotlib.pyplot as plt


logger = logging.getLogger(__name__)


def _extract_crop_production(n: pypsa.Network) -> pd.Series:
    """Aggregate crop production (tonnes) over all regions/resource classes."""

    if "now" not in n.snapshots:
        raise ValueError("Expected snapshot 'now' in solved network")

    production = {}
    flows = n.links_t.p1.loc["now"] if "p1" in n.links_t else pd.Series(dtype=float)
    links_df = n.links

    for link in n.links.index:
        if not link.startswith("produce_"):
            continue

        bus1 = str(links_df.at[link, "bus1"]) if link in links_df.index else ""
        if not bus1.startswith("crop_"):
            continue

        crop_token = link.split("_", 2)[1] if "_" in link[len("produce_") :] else ""
        if not crop_token:
            continue

        value = abs(float(flows.get(link, 0.0)))
        if value <= 0.0:
            continue

        production[crop_token] = production.get(crop_token, 0.0) + value

    pasture_total = 0.0
    for link in n.links.index:
        if not link.startswith("graze_"):
            continue

        value = abs(float(flows.get(link, 0.0)))
        if value <= 0.0:
            continue

        pasture_total += value

    if pasture_total > 0.0:
        production["pasture"] = pasture_total

    return pd.Series(production, dtype=float)


def _extract_crop_use(n: pypsa.Network) -> tuple[pd.Series, pd.Series]:
    """Return crop use split into human consumption vs. animal feed (tonnes)."""

    if "now" not in n.snapshots:
        raise ValueError("Expected snapshot 'now' in solved network")

    human_use: dict[str, float] = {}
    feed_use: dict[str, float] = {}

    flows_p0 = n.links_t.p0.loc["now"] if "p0" in n.links_t else pd.Series(dtype=float)
    flows_p1 = n.links_t.p1.loc["now"] if "p1" in n.links_t else pd.Series(dtype=float)

    for link in n.links.index:
        if not link.startswith("convert_"):
            continue

        remainder = link[len("convert_") :]
        if "_to_" not in remainder:
            continue

        crop_token, output = remainder.split("_to_", 1)
        value = abs(float(flows_p0.get(link, 0.0)))
        if value <= 0.0:
            continue

        if output.startswith("ruminant_feed_") or output.startswith(
            "monogastric_feed_"
        ):
            feed_use[crop_token] = feed_use.get(crop_token, 0.0) + value
        else:
            human_use[crop_token] = human_use.get(crop_token, 0.0) + value

    pasture_total = 0.0
    for link in n.links.index:
        if not link.startswith("graze_"):
            continue

        value = abs(float(flows_p1.get(link, 0.0)))
        if value <= 0.0:
            continue

        pasture_total += value

    if pasture_total > 0.0:
        feed_use["pasture"] = feed_use.get("pasture", 0.0) + pasture_total

    return pd.Series(human_use, dtype=float), pd.Series(feed_use, dtype=float)


def _build_dataframe(
    production: pd.Series, human_use: pd.Series, feed_use: pd.Series
) -> pd.DataFrame:
    """Combine the individual series into a single dataframe for plotting/export."""

    df = pd.DataFrame(
        {
            "production_tonnes": production,
            "human_consumption_tonnes": human_use,
            "animal_feed_tonnes": feed_use,
        }
    ).fillna(0.0)

    if df.empty:
        return df

    df = df[(df > 0).any(axis=1)]

    df["residual_tonnes"] = df["production_tonnes"] - df[
        ["human_consumption_tonnes", "animal_feed_tonnes"]
    ].sum(axis=1)
    tolerance = 1e-6
    residual_mask = df["residual_tonnes"].abs() > tolerance
    if residual_mask.any():
        logger.warning(
            "Crop use does not sum to production for: %s",
            ", ".join(df.index[residual_mask]),
        )

    df.sort_values("production_tonnes", ascending=False, inplace=True)
    return df


def _plot(df: pd.DataFrame, output_pdf: Path) -> None:
    """Render stacked bar chart for crop uses."""

    plt.figure(figsize=(12, 7))

    if df.empty:
        plt.text(0.5, 0.5, "No crop production data", ha="center", va="center")
        plt.axis("off")
    else:
        x = range(len(df))
        human = df["human_consumption_tonnes"].to_numpy()
        feed = df["animal_feed_tonnes"].to_numpy()

        plt.bar(x, human, label="Human consumption", color="#1f77b4")
        plt.bar(x, feed, bottom=human, label="Animal feed", color="#ff7f0e")

        plt.ylabel("Tonnes")
        plt.xticks(x, df.index, rotation=45, ha="right")
        plt.grid(axis="y", alpha=0.3)
        plt.legend()
        plt.tight_layout()

    plt.title("Global crop use breakdown")
    plt.tight_layout()
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_pdf, bbox_inches="tight", dpi=300)
    plt.close()
    logger.info("Wrote stacked crop use plot to %s", output_pdf)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    network = pypsa.Network(snakemake.input.network)
    production = _extract_crop_production(network)
    human_use, feed_use = _extract_crop_use(network)

    df = _build_dataframe(production, human_use, feed_use)

    csv_path = Path(snakemake.output.csv)
    pdf_path = Path(snakemake.output.pdf)

    df.to_csv(csv_path)
    logger.info("Wrote crop use data to %s", csv_path)

    _plot(df, pdf_path)
