# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Mirror the FAOSTAT bulk zips onto the GLADE input-data mirror on Zenodo.

This is a MAINTAINER tool, not part of the build. FAO's bulk endpoints
(``bulks-faostat.fao.org``) serve only the current release and revise
historical values continuously, so builds never fetch them directly.
Instead this tool:

1. downloads the bulk zip of every FAOSTAT domain the model uses (the
   :data:`workflow.scripts.faostat_bulk.FAOSTAT_BULK_URLS` registry),
2. writes a ``faostat_vintages.yaml`` manifest recording each zip's source
   URL, upstream ``Last-Modified`` date, retrieval date and SHA-256 digest,
3. publishes a new version of the GLADE input-data mirror record on Zenodo
   with the zips (as ``faostat_<code>.zip``) and the manifest, keeping the
   record's other files (the land-cover extract) untouched.

Ordinary builds then fetch the zips from Zenodo with no API key (see the
``download_faostat`` rule and ``config['data']['mirror']['zenodo_record']``).

FAOSTAT datasets are licensed CC-BY-4.0, which permits this redistribution
provided the FAO attribution is retained; it is embedded in the Zenodo
deposition metadata (see tools/mirror_metadata.py).

Run inside the project environment, e.g.::

    # Refresh to the current FAO release (new version of the mirror record):
    pixi run -e dev python tools/mirror_faostat.py

    # Dry-run against the Zenodo sandbox, leaving an unpublished draft:
    pixi run -e dev python tools/mirror_faostat.py --sandbox --no-publish

The Zenodo token is read from ``ZENODO_TOKEN`` or ``config/secrets.yaml``
(``credentials.zenodo.token``).

After publishing, set ``config['data']['mirror']['zenodo_record']`` in
``config/default.yaml`` to the printed record id. Note that a new FAOSTAT
vintage invalidates the calibration artefact sets (the record id is part of
the calibration provenance snapshot), so recalibrate with ``tools/calibrate``
afterwards.
"""

import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import sys

import requests
import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent

VINTAGES_FILENAME = "faostat_vintages.yaml"


def _download(url: str, target: Path) -> dict:
    """Stream *url* to *target*; return its vintage manifest entry."""
    with requests.get(url, stream=True, timeout=60) as response:
        response.raise_for_status()
        digest = hashlib.sha256()
        with open(target, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1 << 20):
                handle.write(chunk)
                digest.update(chunk)
        last_modified = response.headers.get("Last-Modified")
    return {
        "source_url": url,
        "upstream_last_modified": last_modified,
        "sha256": digest.hexdigest(),
        "size_bytes": target.stat().st_size,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mirror the FAOSTAT bulk zips onto Zenodo.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--parent-record",
        default=None,
        help=(
            "Mirror record id to publish a new version of. Defaults to "
            "config['data']['mirror']['zenodo_record'] from config/default.yaml."
        ),
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "downloads" / "faostat_mirror",
        help=(
            "Directory for the downloaded zips and the vintage manifest. "
            "Deliberately separate from data/downloads/faostat/ so mirroring "
            "never touches the build's own inputs."
        ),
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
        help="Reuse existing zips and manifest in --work-dir (skip the FAO downloads).",
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
    sys.path.insert(0, str(PROJECT_ROOT))
    sys.path.insert(0, str(PROJECT_ROOT / "tools"))
    from mirror_metadata import build_metadata
    from zenodo_publish import load_secret, publish_dataset, publish_draft

    from workflow.scripts.faostat_bulk import FAOSTAT_BULK_URLS

    zenodo_token = load_secret("ZENODO_TOKEN", "zenodo", "token")
    if not zenodo_token:
        parser.error(
            "Missing Zenodo token. Set ZENODO_TOKEN, or credentials.zenodo.token "
            "in config/secrets.yaml."
        )

    if args.publish_record:
        result = publish_draft(zenodo_token, args.publish_record, sandbox=args.sandbox)
        print("Published.")
        print(f"  record id : {result['record_id']}")
        print(f"  doi       : {result['doi']}")
        print(
            "\nSet this in config/default.yaml under data.mirror:\n"
            f'    zenodo_record: "{result["record_id"]}"'
        )
        return

    with open(PROJECT_ROOT / "config" / "default.yaml") as handle:
        config = yaml.safe_load(handle)
    parent_record = args.parent_record or config["data"]["mirror"]["zenodo_record"]

    args.work_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.work_dir / VINTAGES_FILENAME
    zip_paths = {
        code: args.work_dir / f"faostat_{code}.zip" for code in FAOSTAT_BULK_URLS
    }

    if args.skip_download:
        missing = [
            str(p) for p in [*zip_paths.values(), manifest_path] if not p.exists()
        ]
        if missing:
            parser.error(f"--skip-download given but missing: {', '.join(missing)}")
        print(f"Reusing existing zips and manifest in {args.work_dir}")
    else:
        vintages = {
            "retrieved": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
            "domains": {},
        }
        for code, url in FAOSTAT_BULK_URLS.items():
            print(f"Downloading {code} ({url}) ...")
            vintages["domains"][code] = _download(url, zip_paths[code])
        with open(manifest_path, "w") as handle:
            handle.write(
                "# Vintage manifest for the FAOSTAT bulk zips in this Zenodo\n"
                "# record; written by tools/mirror_faostat.py.\n"
            )
            yaml.safe_dump(vintages, handle, sort_keys=False)

    print("Publishing to Zenodo ...")
    result = publish_dataset(
        token=zenodo_token,
        files=[*zip_paths.values(), manifest_path],
        metadata=build_metadata(
            config["baseline_year"], config["data"]["land_cover"]["version"]
        ),
        parent_record=parent_record,
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
            f'    zenodo_record: "{result["record_id"]}"\n'
            "then recalibrate the affected artefact sets (tools/calibrate)."
        )


if __name__ == "__main__":
    main()
