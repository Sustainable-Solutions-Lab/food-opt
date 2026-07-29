# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Mirror the Copernicus ESA CCI land-cover map onto Zenodo.

This is a MAINTAINER tool, not part of the build. It:

1. downloads the ``satellite-land-cover`` dataset for a given year/version from
   the Copernicus Climate Data Store (needs an ECMWF/CDS token),
2. extracts the ``lccs_class`` variable (the only one the model uses), and
3. publishes a new version of the GLADE input-data mirror record on Zenodo
   with the extract, keeping the record's other files (the FAOSTAT bulk zips)
   untouched.

Ordinary builds then fetch the file from Zenodo with no API key (see the
``download_land_cover`` rule and ``config['data']['mirror']['zenodo_record']``).

The 2016-onwards C3S land-cover maps are licensed CC-BY-4.0, which permits this
redistribution provided the Copernicus attribution and source DOI are retained;
both are embedded in the Zenodo deposition metadata (see
tools/mirror_metadata.py).

Run inside the project environment, e.g.::

    # Refresh / new data version (new version of the mirror record):
    pixi run -e dev python tools/mirror_land_cover.py

    # Dry-run against the Zenodo sandbox, leaving an unpublished draft:
    pixi run -e dev python tools/mirror_land_cover.py --sandbox --no-publish

Credentials are read from environment variables (``ECMWF_DATASTORES_URL`` /
``ECMWF_DATASTORES_KEY`` / ``ZENODO_TOKEN``) or from ``config/secrets.yaml``
(``credentials.ecmwf`` and ``credentials.zenodo``).

After publishing, set ``config['data']['mirror']['zenodo_record']`` in
``config/default.yaml`` to the printed record id. When the mirrored year or
version changes, the old ``land_cover_lccs_class_*.nc`` is kept in the new
record version (files are replaced by name); delete it in the Zenodo web UI
before publishing if it should be dropped.
"""

import argparse
from pathlib import Path
import sys

import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _default_config_values() -> tuple[int, str, str]:
    """Read baseline year, land-cover version and mirror record from config."""
    with open(PROJECT_ROOT / "config" / "default.yaml") as handle:
        config = yaml.safe_load(handle)
    return (
        config["baseline_year"],
        config["data"]["land_cover"]["version"],
        config["data"]["mirror"]["zenodo_record"],
    )


def main() -> None:
    default_year, default_version, default_record = _default_config_values()

    parser = argparse.ArgumentParser(
        description="Mirror the Copernicus ESA CCI land-cover map onto Zenodo.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--year", type=int, default=default_year, help="Land-cover map year."
    )
    parser.add_argument(
        "--version", default=default_version, help="ESA CCI land-cover version."
    )
    parser.add_argument(
        "--parent-record",
        default=default_record,
        help=(
            "Mirror record id to publish a new version of. Defaults to "
            "config['data']['mirror']['zenodo_record']; pass an empty string "
            "to create a brand-new record."
        ),
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "downloads",
        help="Directory for the downloaded archive and extracted NetCDF.",
    )
    parser.add_argument(
        "--sandbox",
        action="store_true",
        help="Use sandbox.zenodo.org instead of the production service.",
    )
    parser.add_argument(
        "--no-publish",
        dest="publish",
        action="store_false",
        help="Leave the Zenodo deposition as an unpublished draft for review.",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Reuse an existing extracted NetCDF in --work-dir (skip the CDS download).",
    )
    parser.add_argument(
        "--publish-record",
        default=None,
        help=(
            "Publish an existing draft deposition by record id (e.g. after "
            "reviewing a draft created with --no-publish) and exit. Irreversible."
        ),
    )
    args = parser.parse_args()

    # Make the workflow scripts and sibling tools importable.
    sys.path.insert(0, str(PROJECT_ROOT / "workflow" / "scripts"))
    sys.path.insert(0, str(PROJECT_ROOT / "tools"))
    import download_land_cover
    import extract_land_cover_class
    from mirror_metadata import build_metadata
    from zenodo_publish import load_secret, publish_dataset, publish_draft

    if args.publish_record:
        zenodo_token = load_secret("ZENODO_TOKEN", "zenodo", "token")
        if not zenodo_token:
            parser.error(
                "Missing Zenodo token. Set ZENODO_TOKEN, or credentials.zenodo.token "
                "in config/secrets.yaml."
            )
        result = publish_draft(zenodo_token, args.publish_record, sandbox=args.sandbox)
        print("Published.")
        print(f"  record id : {result['record_id']}")
        print(f"  doi       : {result['doi']}")
        print(
            "\nSet this in config/default.yaml under data.mirror:\n"
            f'    zenodo_record: "{result["record_id"]}"'
        )
        return

    args.work_dir.mkdir(parents=True, exist_ok=True)
    target_name = f"land_cover_lccs_class_{args.year}_{args.version}.nc"
    extracted = args.work_dir / target_name

    if args.skip_download:
        if not extracted.exists():
            parser.error(f"--skip-download given but {extracted} does not exist")
        print(f"Reusing existing {extracted}")
    else:
        ecmwf_url = load_secret("ECMWF_DATASTORES_URL", "ecmwf", "url")
        ecmwf_key = load_secret("ECMWF_DATASTORES_KEY", "ecmwf", "key")
        if not ecmwf_url or not ecmwf_key:
            parser.error(
                "Missing Copernicus CDS credentials. Set ECMWF_DATASTORES_URL and "
                "ECMWF_DATASTORES_KEY, or credentials.ecmwf.{url,key} in "
                "config/secrets.yaml."
            )

        archive = args.work_dir / f"land_cover_{args.year}_{args.version}.zip"
        print(f"Downloading satellite-land-cover {args.version} for {args.year} ...")
        download_land_cover.main(
            dataset="satellite-land-cover",
            request={
                "variable": "all",
                "year": [str(args.year)],
                "version": [args.version],
            },
            output=archive,
            url=ecmwf_url,
            key=ecmwf_key,
        )
        print(f"Extracting lccs_class -> {extracted} ...")
        extract_land_cover_class.main(input_path=archive, output_path=extracted)
        archive.unlink(missing_ok=True)

    zenodo_token = load_secret("ZENODO_TOKEN", "zenodo", "token")
    if not zenodo_token:
        parser.error(
            "Missing Zenodo token. Set ZENODO_TOKEN, or credentials.zenodo.token "
            "in config/secrets.yaml."
        )

    print("Publishing to Zenodo ...")
    result = publish_dataset(
        token=zenodo_token,
        files=[extracted],
        metadata=build_metadata(args.year, args.version),
        parent_record=args.parent_record or None,
        keep_existing=True,
        sandbox=args.sandbox,
        publish=args.publish,
    )

    print("\nDone.")
    print(f"  record id : {result['record_id']}")
    print(f"  doi       : {result['doi']}")
    if result["draft"]:
        edit_link = result["links"].get("html") or result["links"].get("latest_draft")
        print(f"  status    : DRAFT (not published) - review at {edit_link}")
    else:
        print(
            "\nSet this in config/default.yaml under data.mirror:\n"
            f'    zenodo_record: "{result["record_id"]}"'
        )


if __name__ == "__main__":
    main()
