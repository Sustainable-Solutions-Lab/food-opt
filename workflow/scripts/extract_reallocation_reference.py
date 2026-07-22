# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Extract compact production dispatch from two realized endpoint solutions."""

import gc

import pandas as pd
import pypsa

COMPONENT_CARRIERS = {
    "cropland": ("crop_production", "crop_production_multi"),
    "grassland": ("grassland_production",),
    "feed": ("animal_production",),
}


def _dispatch(path: str, column: str) -> pd.DataFrame:
    network = pypsa.Network(path)
    links = network.links.static
    flow = network.links.dynamic.p0.iloc[0]
    frames = []
    for component, carriers in COMPONENT_CARRIERS.items():
        names = links.index[links["carrier"].isin(carriers)]
        frames.append(
            pd.DataFrame(
                {
                    "component": component,
                    "link": names,
                    column: flow.reindex(names).fillna(0.0).to_numpy(),
                }
            )
        )
    del network
    gc.collect()
    return pd.concat(frames, ignore_index=True)


baseline = _dispatch(snakemake.input.baseline, "baseline_dispatch")
endpoint = _dispatch(snakemake.input.endpoint, "endpoint_dispatch")
reference = baseline.merge(
    endpoint, on=["component", "link"], how="outer", validate="one_to_one"
)
if reference[["baseline_dispatch", "endpoint_dispatch"]].isna().any().any():
    raise ValueError(
        "Reallocation reference endpoints contain different production links"
    )
reference.to_parquet(snakemake.output[0], index=False)
