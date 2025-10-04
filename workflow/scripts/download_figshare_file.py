# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Download a Figshare file using the public API.

Snakemake passes the ``snakemake`` object into this module; no standalone CLI
usage is supported.
"""

from pathlib import Path

import requests
from tqdm.auto import tqdm

REQUEST_TIMEOUT = (10, 120)  # (connect, read) seconds
CHUNK_SIZE = 4 * 1024 * 1024  # 4 MiB


def _resolve_file(article_id: int, file_name: str) -> int:
    base_url = "https://api.figshare.com/v2"
    response = requests.get(
        f"{base_url}/articles/{article_id}/files", timeout=REQUEST_TIMEOUT
    )
    response.raise_for_status()

    for file_info in response.json():
        if file_info.get("name") == file_name:
            return int(file_info["id"])

    raise FileNotFoundError(
        f"No file named '{file_name}' was found for article {article_id}."
    )


def _download_file(file_id: int, output: Path, show_progress: bool) -> None:
    base_url = "https://api.figshare.com/v2"
    with requests.get(
        f"{base_url}/file/download/{file_id}",
        stream=True,
        timeout=REQUEST_TIMEOUT,
    ) as response:
        response.raise_for_status()
        total_size = response.headers.get("Content-Length")
        total_bytes = int(total_size) if total_size is not None else None

        output.parent.mkdir(parents=True, exist_ok=True)
        desc = f"Downloading {output.name}"
        with (
            output.open("wb") as f,
            tqdm(
                total=total_bytes,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc=desc,
                dynamic_ncols=True,
                disable=not show_progress,
            ) as progress,
        ):
            for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                if not chunk:
                    continue
                f.write(chunk)
                progress.update(len(chunk))


def main(article_id: int, file_name: str, output: Path, show_progress: bool) -> None:
    try:
        file_id = _resolve_file(article_id, file_name)
    except FileNotFoundError as exc:  # pragma: no cover - user-facing error
        raise SystemExit(str(exc)) from exc

    _download_file(file_id, output, show_progress)


if __name__ == "__main__":
    main(
        article_id=int(snakemake.params.article_id),
        file_name=snakemake.params.file_name,
        output=Path(snakemake.output[0]),
        show_progress=bool(snakemake.params.show_progress),
    )
