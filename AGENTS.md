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
- If a persistent Julia REPL is exposed over MCP (e.g. [julia-mcp](https://github.com/aplavin/julia-mcp), which provides a `julia_eval` tool), prefer it over shelling out to `julia` for experiments and interactive work: loaded packages, compiled code, and simulation state persist between calls, so iterating on a long-running simulation skips repeated startup and recompilation. Select the environment with the tool's `env_path` argument instead of `Pkg.activate()` — the repo root for library work, `.buildkite` for experiments. Keep running the test suite through the shell, though: `julia-mcp` treats any `env_path` ending in `test` as a `TestEnv.activate()` request, which fails here (`can not merge projects`) because `test/` carries its own `Project.toml`. Fall back to the plain `julia --project=...` commands when no such server is configured.
- Keep edits modular and within the appropriate submodel directory (e.g., `src/standalone` or `src/integrated`).
- Match existing style: explicit names, narrow imports, and use multiple dispatch effectively.
- Reuse CliMA packages before writing a helper: ClimaAnalysis for operations on `OutputVar`s (e.g. `average_lon`, `weighted_average_lonlat`, `integrate_lonlat`), ClimaCore for operations on `Field`s (e.g. `column_integral_definite!`), and the helpers already in `ext/land_sim_vis/`. When a general capability is missing, propose it upstream, or add it to `src/` with unit tests, rather than as a private helper in a plotting function or script.
- Keep scientific assumptions out of plotting and leaderboard functions: compute derived quantities (such as a steady-state estimate) in their own documented function, unit-tested where possible, and keep the plotting function to loading, plotting, and saving.
- Compute depth integrals of model variables in the model, as 2D diagnostics built on ClimaCore's column integrals (as `iwc` and `soc_int` do), rather than by integrating depth-resolved output afterwards.
- Register a diagnostic that exists only in some configurations with an `add_diagnostics!` method for the model that owns the variable it reads, dispatched on the subcomponent that creates it (e.g. `slit` for a `LandModel` with a `PrognosticCarbonModel`, not for the `CanopyModel`).
- Keep comments sparse and succinct. Write every comment for a first-time reader of the *final* code, not as a diff narration: do not explain what a change fixes or its history ("this fixes Y", "this newest code fixes the issue introduced by..."), and do not restate the docstring inline. Prefer one short line over a multi-line block; comment only non-obvious rationale (the *why*, not the *what*).
- Follow the software design patterns in [docs/dev-guides/code-quality/software_design_patterns.md](docs/dev-guides/code-quality/software_design_patterns.md) for new code and refactor toward them when touching existing code.
- Before committing code, format it with the JuliaFormatter version pinned in `.github/workflows/JuliaFormatter.yml` (currently 2.10.1), following the rules in `.JuliaFormatter.toml`; other versions reformat unrelated files. Install it once in a shared environment with `julia --project=@juliaformatter -e 'using Pkg; Pkg.add(name="JuliaFormatter", version="2.10.1")'`, then run `julia --project=@juliaformatter -e 'using JuliaFormatter; format(".")'`.
- In docs Markdown, never start a bullet or a paragraph with inline math (`- $k$ is the ...`): Julia's Markdown parser turns the leading `$...$` into a centered display equation and pushes the rest of the line to a new paragraph. Lead with words, as in `docs/src/standalone/pages/inland_waters/slab_lake.md` (`- the lake depth, $d$ (m)`).
- When you `export` a new symbol from `ClimaLand`, add it to a `@docs` block on the matching `docs/src/APIs/*.md` page. The docs build sets `checkdocs = :exports`, so an exported name with a docstring that is not included in the manual fails the docs CI job.
- Never tag people on GitHub (`@username`) in PR descriptions, issues, or comments you write; an `@` mention notifies a real person on your behalf. Refer to them by name without the `@`, or cite the PR or issue number instead.
- Merge only what the repo will maintain (run in CI or covered by tests): scripts that only produce artifact data, such as spin-ups, go with the artifact's recipe in ClimaArtifacts or on a side branch, and research tools, such as offline emulators, stay out of ClimaLand.
- Keep restarts reproducible: initial guesses of Newton and other iterative solves (e.g., surface temperature solves) depend only on the state `Y`, the drivers, and constants, never on values cached in `p` from an earlier step, because a restart rebuilds `p` from `Y` alone.

## Self-correction

- If the code map in [docs/src/repo_structure.md](docs/src/repo_structure.md) is discovered to be stale, update it.
- If the user gives a correction about how work should be done in this repo, add it to `Local norms` or another clearly labeled persistent section in this file or in [docs/src/repo_structure.md](docs/src/repo_structure.md) so future sessions inherit it.
