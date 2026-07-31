# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Content fingerprint of a calibration step's external inputs.

A calibration step is stale when re-running it would change its artefacts,
which depends only on inputs from outside the artefact set: curated and bundled
data, the scripts, and the step config. Those are hashed when the set is
generated and compared when it is checked, so the set's own artefacts and the
chain's intermediates never enter into the answer.

The input list comes from the live Snakemake DAG, built through the Python API:
every file the DAG consumes that no rule in the DAG produces is an external
leaf. Building it afresh on each check means a dependency that a rule or code
change adds or drops is reflected straight away.

Config files are hashed after a YAML round-trip, so reformatting and comment
edits do not register. The base config belongs to ``calibration_provenance``,
which compares it structurally.
"""

import argparse
import contextlib
import hashlib
import io
from pathlib import Path
import sys

from snakemake.api import SnakemakeApi
from snakemake.settings.types import (
    ConfigSettings,
    DAGSettings,
    OutputSettings,
    ResourceSettings,
)
import yaml

from workflow.scripts.solve_namespace import load_merged_config

FINGERPRINT_FILE = "fingerprint.yaml"

# Solve-time keys are absent from the structural provenance stamp.  Select the
# subset that changes a calibration fit here, excluding orchestration fields
# such as output paths, the generate switch, and per-stage pin/intercept values
# that the driver overrides itself.
FIT_CONFIG_PATHS = {
    "market_response": (
        "market_response.n_blocks",
        "market_response.expansion_range",
        "market_response.contraction_range",
        "market_response.width_growth",
        "market_response.granularity",
        "market_response.pin_slack_cost",
        "market_response.elasticities",
        "market_response.elasticity_factor",
        "market_response.components",
        "market_response.calibration.sweeps",
    ),
}

HEADER = """# Content fingerprint of each calibration step's external inputs, written by
# workflow/scripts/calibration_fingerprint.py via tools/calibrate. Compared by
# `tools/calibrate --check`; do not edit by hand.
# Licensing: see the annotation in REUSE.toml.
"""


def sha256_file(path: Path) -> str:
    """Return the SHA-256 of a file, streamed so large inputs stay cheap."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_dir(path: Path) -> str:
    """Return a digest of a directory input's file manifest.

    Manually downloaded sources arrive as directories of hundreds of megabytes.
    Hashing the sorted (relative path, size) listing costs nothing and moves
    whenever a file is added, removed or replaced.
    """
    entries = sorted(
        f"{f.relative_to(path)} {f.stat().st_size}"
        for f in path.rglob("*")
        if f.is_file()
    )
    return hashlib.sha256("\n".join(entries).encode()).hexdigest()


def sha256_yaml(path: Path) -> str:
    """Return the SHA-256 of a YAML file's parsed content.

    Hashing the loaded document rather than the raw bytes keeps comment and
    formatting edits from registering as a change.
    """
    data = yaml.safe_load(path.read_text())
    return hashlib.sha256(
        yaml.safe_dump(data, sort_keys=True, default_flow_style=False).encode()
    ).hexdigest()


def effective_fit_config(
    step: str, configfiles: list[str], workdir: Path
) -> dict[str, object]:
    """Return the resolved config values that determine ``step``'s fit."""
    paths = FIT_CONFIG_PATHS.get(step, ())
    if not paths:
        return {}

    default_config = workdir / "config/default.yaml"
    resolved_files = [default_config]
    for raw in configfiles:
        path = Path(raw)
        path = path if path.is_absolute() else workdir / path
        if path != default_config:
            resolved_files.append(path)
    config = load_merged_config(*resolved_files)

    selected = {}
    for dotted in paths:
        value = config
        for key in dotted.split("."):
            value = value[key]
        selected[dotted] = value
    return selected


def dag_files(
    configfiles: list[str], name: str, targets: list[str], workdir: Path
) -> tuple[set[str], set[str], set[str]]:
    """Return (inputs, outputs, code) over every job in the step's DAG.

    ``code`` is the scripts the DAG's rules run plus the rule files defining
    them; Snakemake carries those outside the input lists, so they are gathered
    separately here.

    The DAG is prepared without taking the workflow lock, so a check is safe to
    run while something else is using the directory.
    """
    with SnakemakeApi(OutputSettings(dryrun=True, quiet={"all"})) as api:
        workflow = api.workflow(
            resource_settings=ResourceSettings(cores=1),
            config_settings=ConfigSettings(
                configfiles=[Path(c) for c in configfiles],
                config={"name": name},
            ),
            snakefile=workdir / "workflow" / "Snakefile",
            workdir=workdir,
        )
        workflow.dag(dag_settings=DAGSettings(targets=frozenset(targets)))
        internal = workflow._workflow
        internal._prepare_dag(
            forceall=False,
            ignore_incomplete=True,
            lock_warn_only=True,
            nolock=True,
        )
        internal._build_dag()

        inputs: set[str] = set()
        outputs: set[str] = set()
        code: set[str] = set()
        for job in internal.dag.jobs:
            inputs.update(str(f) for f in job.input)
            outputs.update(str(f) for f in job.output)
            rule_file = getattr(job.rule, "snakefile", None)
            if rule_file:
                code.add(str(rule_file))
                script = getattr(job.rule, "script", None)
                if script:
                    # Scripts are named relative to the rule file that runs them.
                    code.add(str((Path(rule_file).parent / script).resolve()))
    return inputs, outputs, code


def build_fingerprint(
    configfiles: list[str],
    name: str,
    targets: list[str],
    hash_configs: list[str],
    workdir: Path,
    step: str = "",
) -> dict:
    """Hash a step's external leaves, its config, and the artefacts it wrote."""
    inputs, outputs, code = dag_files(configfiles, name, targets, workdir)

    leaf_hashes = {}
    missing = []
    for leaf in sorted(inputs - outputs):
        path = Path(leaf)
        if path.is_file():
            leaf_hashes[leaf] = sha256_file(path)
        elif path.is_dir():
            leaf_hashes[leaf] = sha256_dir(path)
        else:
            missing.append(leaf)

    def relative(path: Path) -> str:
        """Key code by repo-relative path so a fingerprint stays portable."""
        try:
            return str(path.relative_to(workdir))
        except ValueError:
            return str(path)

    return {
        "config_files": {c: sha256_yaml(Path(c)) for c in hash_configs},
        "effective_config": effective_fit_config(step, configfiles, workdir),
        "code": {
            relative(Path(c)): sha256_file(Path(c))
            for c in sorted(code)
            if Path(c).is_file()
        },
        "inputs": leaf_hashes,
        # The step's own artefacts, so that one edited or checked out behind the
        # chain's back is reported rather than silently trusted.
        "artefacts": {
            tgt: sha256_file(Path(tgt))
            for tgt in sorted(targets)
            if Path(tgt).is_file()
        },
        "missing_inputs": missing,
    }


def load_fingerprints(cal_dir: Path) -> dict:
    path = cal_dir / FINGERPRINT_FILE
    if not path.is_file():
        return {}
    return yaml.safe_load(path.read_text()) or {}


def save_fingerprints(cal_dir: Path, data: dict) -> None:
    cal_dir.mkdir(parents=True, exist_ok=True)
    path = cal_dir / FINGERPRINT_FILE
    path.write_text(HEADER + yaml.safe_dump(data, sort_keys=True))


def compare(recorded: dict, current: dict) -> list[str]:
    """Return human-readable reasons the step is stale, empty when it is not."""
    reasons = []

    for kind, key in (
        ("config", "config_files"),
        ("effective config", "effective_config"),
        ("code", "code"),
        ("input", "inputs"),
        ("artefact", "artefacts"),
    ):
        old, new = recorded.get(key, {}), current.get(key, {})
        for path in sorted(set(old) | set(new)):
            if path not in old:
                reasons.append(f"{path} ({kind} added)")
            elif path not in new:
                reasons.append(f"{path} ({kind} no longer used)")
            elif old[path] != new[path]:
                changed = (
                    "modified since it was generated"
                    if kind == "artefact"
                    else "content changed"
                )
                reasons.append(f"{path} ({kind} {changed})")

    for path in current.get("missing_inputs", []):
        reasons.append(f"{path} (input missing from disk)")

    return reasons


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["record", "check"])
    parser.add_argument("--step", required=True, help="calibration step label")
    parser.add_argument("--cal-dir", required=True, help="artefact set directory")
    parser.add_argument("--name", required=True, help="workflow name override")
    parser.add_argument(
        "--configfile", action="append", default=[], help="DAG config (repeatable)"
    )
    parser.add_argument(
        "--hash-config",
        action="append",
        default=[],
        help="config whose content is part of the fingerprint (repeatable)",
    )
    parser.add_argument(
        "--target",
        action="append",
        default=[],
        help="artefact this step produces (repeatable)",
    )
    args = parser.parse_args()

    workdir = Path.cwd()
    cal_dir = Path(args.cal_dir)

    # Snakemake narrates DAG construction; keep that out of the report, but
    # hand it back if the build fails, since then it is the error message.
    noise = io.StringIO()
    try:
        with contextlib.redirect_stdout(noise), contextlib.redirect_stderr(noise):
            current = build_fingerprint(
                args.configfile,
                args.name,
                args.target,
                args.hash_config,
                workdir,
                step=args.step,
            )
    except Exception as exc:
        print(f"  [error]      {args.step:<10} ({type(exc).__name__}: {exc})")
        sys.stderr.write(noise.getvalue())
        return 2

    if args.mode == "record":
        data = load_fingerprints(cal_dir)
        data.setdefault("steps", {})[args.step] = current
        save_fingerprints(cal_dir, data)
        return 0

    recorded = load_fingerprints(cal_dir).get("steps", {}).get(args.step)
    if recorded is None:
        print(f"  [STALE]      {args.step:<10} (no fingerprint recorded)")
        return 1

    reasons = compare(recorded, current)
    if not reasons:
        print(f"  [up-to-date] {args.step}")
        return 0

    print(f"  [STALE]      {args.step:<10} (run: tools/calibrate {args.step})")
    for reason in reasons[:10]:
        print(f"                 {reason}")
    if len(reasons) > 10:
        print(f"                 ... and {len(reasons) - 10} more")
    return 1


if __name__ == "__main__":
    sys.exit(main())
