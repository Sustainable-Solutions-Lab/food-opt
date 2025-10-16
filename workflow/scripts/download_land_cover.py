# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Download land cover data using the ECMWF datastores client.

Snakemake passes the ``snakemake`` object into this module; no standalone CLI
usage is supported.
"""

from pathlib import Path

from ecmwf.datastores import Client


def main(dataset: str, request: dict, output: Path) -> None:
    """Download land cover dataset.

    Parameters
    ----------
    dataset : str
        The dataset identifier (e.g., "satellite-land-cover").
    request : dict
        The request parameters including variable, year, and version.
    output : Path
        The output archive path (ZIP containing the NetCDF payload).
    """
    output.parent.mkdir(parents=True, exist_ok=True)

    client = Client()
    client.retrieve(dataset, request, target=str(output))


if __name__ == "__main__":
    main(
        dataset=snakemake.params.dataset,
        request=snakemake.params.request,
        output=Path(snakemake.output[0]),
    )
