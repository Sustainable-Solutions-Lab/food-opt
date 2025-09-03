<!--
SPDX-FileCopyrightText: 2025 Koen van Greevenbroek

SPDX-License-Identifier: CC-BY-4.0
-->

# AGENTS.md

Guidance for AI coding agents contributing to this repository.

If you are Claude Code, see also `CLAUDE.md` for Claude-specific tips.

## Purpose

Provide clear expectations and a safe, efficient workflow so agents can make small, correct, and reversible changes that fit the project’s conventions.

## Project Overview (brief)

- Global food systems optimization using linear programming.
- Built on PyPSA for modeling and Snakemake for workflow orchestration.
- Configuration-driven; results materialized under `results/{config_name}/`.

For a longer overview and concrete file references, see `CLAUDE.md`.

## Core Principles

- Minimal, surgical changes: Tackle the user’s request precisely; avoid drive‑by refactors.
- Keep code concise: Prefer simple control flow; fail early on invalid external inputs.
- Consistent style: Follow existing patterns in nearby files; don’t introduce new paradigms ad hoc.
- Reproducibility: Use the Snakemake targets below to validate changes; don’t hand‑run ad hoc pipelines unless necessary.
- No unused imports: The linter removes them automatically; only add imports when adding code that uses them.
- Documentation-first interfaces: If you change a script’s inputs/outputs, update inline docstrings and any referenced docs/config keys.

## Environment & Tooling

- Dependency manager: `uv` (see `pyproject.toml`).
- Lint/format: `ruff` (auto-enforced via hooks; no manual action usually needed).
- Workflow engine: `snakemake`.

Recommended commands (run inside the uv environment):

```bash
# Install and sync dependencies
uv sync

# Run the full workflow (data prep → build → solve)
uv run snakemake -j4 all

# Build model only
uv run snakemake -j4 results/{config_name}/build/model.nc

# Solve model only (after build)
uv run snakemake -j4 results/{config_name}/solved/model.nc
```

Notes:

- Snakemake tracks code changes and will rerun affected rules; manual cleanup of `results/` is usually unnecessary.
- Prefer small, testable edits and validate by running the narrowest target that exercises your change.

## Repository Conventions

- Scripts used by the workflow live in `workflow/scripts/`.
- Configuration lives under `config/` (e.g., `config.yaml`).
- Input data under `data/`; outputs under `results/` (structured by config name).
- Don’t commit large data or generated results; `.gitignore` and the workflow manage these.

## PyPSA Modeling Notes

- Use standard PyPSA components: carriers, buses, stores, links, etc.
- For multi-bus links:
  - `bus0` is the (first) input bus.
  - `bus1` is the (first) output bus; `efficiency` governs bus0→bus1.
  - `bus2`, `bus3`, … are additional legs with `efficiency2`, `efficiency3`, …
    - Positive efficiency ⇒ output; negative ⇒ input (relative to `bus0`).

## When Implementing Changes

- Keep function/module scope tight; avoid broad rewrites.
- Mirror existing error handling: validate external data; trust internal invariants.
- Add or adjust docstrings where behavior or parameters change.
- If you add a new rule or script, integrate it into the `workflow/Snakefile` and ensure targets are reproducible.
- Don’t introduce network calls or external services in core code unless explicitly required by the task.

## Validation Checklist

- Narrow target runs clean via Snakemake for at least one `config_name`.
- No new linter errors; no unused imports.
- Results land under the expected `results/{config_name}/...` path(s).
- README or relevant docs updated when changing user-visible behavior.

## Agent-Specific Notes

- Claude Code: Follow this document and the additional context in `CLAUDE.md`.
- Other agents (e.g., GitHub Copilot/Cursor/Codex-based tools): Adhere to the same principles—small patches, strong validation, and consistency with existing styles and commands.

## Safety & Licensing

- Respect SPDX headers; keep or add them to new files following repository practice.
- Do not introduce secrets, credentials, or hard-coded local paths.
- Use only licensed datasets and dependencies already declared in `pyproject.toml` unless explicitly instructed to add new ones.

