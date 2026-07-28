# Unsupervised calibration loop — prognostic-LAI pipeline

An unattended Claude Code loop that runs a three-stage calibration on NCAR
Derecho, ending with a calibrated prognostic-LAI ClimaLand configuration.

| stage | what it does | config | params | ensemble |
|---|---|---|---|---|
| 1 | calibrate GPP + energy fluxes, prescribed MODIS LAI | `configs/sifgpp_lhf_shf_lwu_rosetta.jl` | 10 | 21 |
| 2 | rebuild the optimal-LAI initial-condition artifact from a 10-year run | `experiments/long_runs/snowy_land_pmodel.jl` | — | — |
| 3 | calibrate prognostic LAI against MODIS | `configs/lai.jl` | 5 | 11 |

The stages are strictly sequential. Stage 1's parameters set the potential GPP
that drives LAI, and are committed to `toml/default_parameters.toml` so stages 2
and 3 inherit them. Stage 2 exists so stage 3 can start near equilibrium: with
`optimal_lai_tau_long_term` at 2 years, a stock initial condition would otherwise
force a multi-year spinup on every ensemble member of every iteration.

## Files

- `loop_prompt.md` — the loop's instructions. Read this to understand the task.
- `STATE.md` — where the pipeline is. Rewritten and committed every iteration.
- `LOG.md` — how it got there. One dated entry per iteration.

## Launching

```bash
cp .claude/unsupervised-loop.settings.json .claude/settings.local.json
export GH_TOKEN=...                # PAT with PR write; gh's keyring is not
                                   # reachable from Claude subprocesses
claude --dangerously-skip-permissions
```

Then feed it `loop_prompt.md`.

Before launching, check the prerequisites listed in the `_note` field of
`.claude/unsupervised-loop.settings.json` — in particular the absolute path to
`.claude/hooks/guard.sh`, which must point at the HPC checkout.

## What this branch added to make the pipeline possible

- `use_rosetta` and `prognostic_lai` on `CalibrateConfig`, threaded through
  `setup_model` (ported and adapted from the stale `ar/calibrate_inversion_nee`).
- The `olf0`, `olvpd` and `olgsl` diagnostics. Since `f0`, `vpd_gs` and `GSL`
  became climate-derived rather than read from a map, they had no diagnostics —
  and without them a spun-up run cannot be turned back into an initial condition.
- Both calibration configs, rewritten against current parameter names and
  conventions rather than copied. See the note in
  `configs/sifgpp_lhf_shf_lwu_rosetta.jl` about `pmodel_α`: it changed meaning in
  PR #1817, and a prior copied unconverted from an older branch would silently
  calibrate against a ~1-day acclimation timescale instead of ~36 days.
