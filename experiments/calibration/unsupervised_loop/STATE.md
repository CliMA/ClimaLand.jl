# Pipeline state

Single source of truth for where the three-stage calibration pipeline is.
The loop reads this first every iteration and commits it every iteration.
See `loop_prompt.md` for what each stage means and how it is gated.

Status values: `not_started` | `running` | `done` | `failed`

All heavy output lives under `/glade/derecho/scratch/arenchon/claude/`.

---

## Stage 1 — GPP + energy fluxes

- **status:** not_started
- **config:** `experiments/calibration/configs/sifgpp_lhf_shf_lwu_rosetta.jl`
- **ensemble:** 21 (10 params, TransformUnscented)
- **PBS job ids:** —
- **output dir:** —
- **calibrated parameters:** —
- **committed to `toml/default_parameters.toml`:** no
- **notes:** —

## Stage 2 — rebuild prognostic-LAI initial conditions

- **status:** not_started
- **driver:** `experiments/long_runs/snowy_land_pmodel.jl`, `PROGNOSTIC_LAI` set, 10 years, GPU
- **required diagnostics:** `lai`, `a0a`, `pra`, `olf0`, `olvpd`, `olgsl`
- **PBS job ids:** —
- **run output:** —
- **new artifact path:** —
- **`~/.julia/artifacts/Overrides.toml` verified:** no
- **notes:** —

## Stage 3 — prognostic LAI

- **status:** not_started
- **config:** `experiments/calibration/configs/lai.jl`
- **ensemble:** 11 (5 params, TransformUnscented)
- **validity mask built for this grid:** no — needs stage 2's long-run
  diagnostics as the reference simdir
- **natural-vegetation mask built for this grid:** yes, pre-built 2026-07-31 at
  `experiments/calibration/natural_vegetation_mask.jld2` (gitignored, so it is
  on disk in the checkout only). 10113 of 22420 land cells kept. Rebuild only if
  `nelements` changes.
- **PBS job ids:** —
- **output dir:** —
- **calibrated parameters:** —
- **committed to `toml/default_parameters.toml`:** no
- **notes:** —

---

## Blockers

None recorded.
