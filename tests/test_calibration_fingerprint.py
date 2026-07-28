# SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Unit tests for calibration input fingerprinting."""

from workflow.scripts.calibration_fingerprint import (
    compare,
    sha256_dir,
    sha256_file,
    sha256_yaml,
)


def test_yaml_hash_ignores_comments_and_formatting(tmp_path):
    a = tmp_path / "a.yaml"
    b = tmp_path / "b.yaml"
    a.write_text("# a comment\nfoo:\n  bar: 1\n  baz: [1, 2]\n")
    b.write_text("foo: {baz: [1, 2], bar: 1}\n")
    assert sha256_yaml(a) == sha256_yaml(b)
    assert sha256_file(a) != sha256_file(b)


def test_yaml_hash_tracks_values(tmp_path):
    a = tmp_path / "a.yaml"
    b = tmp_path / "b.yaml"
    a.write_text("foo:\n  bar: 1\n")
    b.write_text("foo:\n  bar: 2\n")
    assert sha256_yaml(a) != sha256_yaml(b)


def test_dir_hash_tracks_membership_and_size(tmp_path):
    root = tmp_path / "src"
    (root / "nested").mkdir(parents=True)
    (root / "nested" / "one.csv").write_text("abc")
    before = sha256_dir(root)

    assert sha256_dir(root) == before

    (root / "two.csv").write_text("d")
    added = sha256_dir(root)
    assert added != before

    (root / "two.csv").write_text("de")
    assert sha256_dir(root) != added


def _fingerprint(**overrides):
    base = {
        "config_files": {"config/calibration/feed.yaml": "cfg"},
        "code": {"workflow/scripts/build_model.py": "code"},
        "inputs": {"data/curated/nutrition.csv": "in"},
        "artefacts": {"data/curated/calibration/unit/exogenous_feed.csv": "art"},
        "missing_inputs": [],
    }
    base.update(overrides)
    return base


def test_identical_fingerprints_are_not_stale():
    assert compare(_fingerprint(), _fingerprint()) == []


def test_changed_input_is_reported():
    current = _fingerprint(inputs={"data/curated/nutrition.csv": "other"})
    (reason,) = compare(_fingerprint(), current)
    assert "nutrition.csv" in reason and "input content changed" in reason


def test_changed_artefact_is_reported_as_hand_modification():
    current = _fingerprint(
        artefacts={"data/curated/calibration/unit/exogenous_feed.csv": "other"}
    )
    (reason,) = compare(_fingerprint(), current)
    assert "modified since it was generated" in reason


def test_added_and_dropped_dependencies_are_reported():
    current = _fingerprint(code={"workflow/scripts/build_luc.py": "code"})
    reasons = compare(_fingerprint(), current)
    assert any("build_luc.py" in r and "code added" in r for r in reasons)
    assert any("build_model.py" in r and "code no longer used" in r for r in reasons)


def test_missing_input_is_reported():
    current = _fingerprint(missing_inputs=["data/manually_downloaded/gbd"])
    (reason,) = compare(_fingerprint(), current)
    assert "input missing from disk" in reason
