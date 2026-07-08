# ClimaLand Agent Guide

## Ecosystem Guidelines

Please refer to the shared CliMA agent index for ecosystem-wide rules regarding architecture, performance, code quality, infrastructure, and workflows:

- [docs/dev-guides/AGENTS.md](docs/dev-guides/AGENTS.md) — Shared CliMA agent guidelines.

> Shared guides live at `docs/dev-guides/` and are vendored from the canonical source:
> <https://github.com/CliMA/DeveloperGuides>. Edit shared guides there, not here.

## Before You Act: Agent Autonomy

Before making changes that are externally visible or scientifically consequential (`git push`, version bumps, reproducibility-test edits, CI config changes, public API renames), check [docs/dev-guides/workflow/agent_autonomy.md](docs/dev-guides/workflow/agent_autonomy.md). The boundaries listed there require explicit user approval.

## Repo-Specific Guidelines

Always read the ClimaLand-specific structure guide before working in this repository:

- [docs/src/repo_structure.md](docs/src/repo_structure.md) — directory tree, environments, and test groups for *this* repo.

## Local norms

- For runtime validation of experiments, prefer `julia --project=.buildkite ...` as the `.buildkite` environment manages the necessary dependencies for simulations.
- For package tests, prefer `Pkg.test()` or running `julia --project=test test/runtests.jl` over manually including files.
- Keep edits modular and within the appropriate submodel directory (e.g., `src/standalone` or `src/integrated`).
- Match existing style: explicit names, narrow imports, and use multiple dispatch effectively.
- Keep comments sparse and succinct. Write every comment for a first-time reader of the *final* code, not as a diff narration: do not explain what a change fixes or its history ("this fixes Y", "this newest code fixes the issue introduced by..."), and do not restate the docstring inline. Prefer one short line over a multi-line block; comment only non-obvious rationale (the *why*, not the *what*).
- Follow the software design patterns in [docs/dev-guides/code-quality/software_design_patterns.md](docs/dev-guides/code-quality/software_design_patterns.md) for new code and refactor toward them when touching existing code.
- Run `julia -e 'using JuliaFormatter; format(".")'` before committing code, adhering to the rules in `.JuliaFormatter.toml`.

## Self-correction

- If the code map in [docs/src/repo_structure.md](docs/src/repo_structure.md) is discovered to be stale, update it.
- If the user gives a correction about how work should be done in this repo, add it to `Local norms` or another clearly labeled persistent section in this file or in [docs/src/repo_structure.md](docs/src/repo_structure.md) so future sessions inherit it.
