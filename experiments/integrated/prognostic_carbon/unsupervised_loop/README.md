# Unsupervised loop — prognostic carbon pools

An unattended Claude Code loop on NCAR Derecho that implements a simple
prognostic live/dead carbon model in ClimaLand and brings global biomass into
reasonable agreement with observations.

| stage | what it does | where it runs |
|---|---|---|
| 0 | port the parametrized single-pixel ERA5 driver + 20-site battery | CPU, 1 node |
| 1 | carbon pools (sugar, leaf, stem, root) in `biomass.jl` | CPU battery |
| 2 | pool-based autotrophic respiration | CPU battery |
| 3 | prognostic SOC — litter input, decomposition sink | CPU battery |
| 4 | offline spinup of the pools to steady state | CPU, offline |
| 5 | global run, comparison to biomass observations, tuning | GPU |

## Files

- `../MODEL.md` — the model specification. Pools, fluxes, no-PFT climate
  proxies, SOC coupling, parameters. Read this first.
- `loop_prompt.md` — the loop's instructions: process, staging, Derecho
  operational rules.
- `STATE.md` — where the work is. Rewritten and committed every iteration.
- `LOG.md` — how it got there. One dated entry per iteration.

## The three constraints

1. The carbon model does not change GPP or LAI — they are inputs. Autotrophic
   respiration is the one deliberate exception.
2. No PFTs. Global constants, optionally split C3/C4; climate proxies (mean
   annual temperature, mean annual precipitation) only if constants fail.
3. Heavy output goes to `/glade/derecho/scratch/arenchon/claude`, never the repo.

## Launching

```bash
cp .claude/unsupervised-loop.settings.json .claude/settings.local.json
export GH_TOKEN=...                # PAT with PR write; gh's keyring is not
                                   # reachable from Claude subprocesses
claude --dangerously-skip-permissions
```

Then feed it `loop_prompt.md`.

Before launching, check the prerequisites in the `_note` field of
`.claude/unsupervised-loop.settings.json` — in particular the absolute path to
`.claude/hooks/guard.sh`, which must point at this HPC checkout.
