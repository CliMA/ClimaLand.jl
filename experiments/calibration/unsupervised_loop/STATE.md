# Pipeline state

Single source of truth for where the three-stage calibration pipeline is.
The loop reads this first every iteration and commits it every iteration.
See `loop_prompt.md` for what each stage means and how it is gated.

Status values: `not_started` | `running` | `done` | `failed`

All heavy output lives under `/glade/derecho/scratch/arenchon/claude/`.

---

## Stage 1 — GPP + energy fluxes

- **status:** running
- **config:** `experiments/calibration/configs/sifgpp_lhf_shf_lwu_rosetta.jl`
- **ensemble:** 21 (10 params, TransformUnscented)
- **PBS job ids:** orchestrator `6967486.desched1` (submitted 2026-07-31, queue
  `cpu`, 12 h walltime). Member job ids: not yet — they appear once the
  orchestrator reaches iteration 1.
- **output dir:** `/glade/derecho/scratch/arenchon/claude/calibration_stage1_gpp_energy`
- **orchestrator log:** `clima_calibration.o6967486` in the repo root
- **calibrated parameters:** —
- **committed to `toml/default_parameters.toml`:** no
- **queue status (2026-07-31 11:57, iteration 2):** still `Q` after 1 h 08 m
  eligible. PBS reason: `Not Running: Job is requesting an exclusive node and
  node is in use`. The orchestrator inherits `place = scatter:exclhost` from the
  `main` queue, so this ~1-core babysitter process is waiting on a whole
  128-core node. Derecho is busy (1447 queued / 338 running server-wide). Not
  escalating yet: `job_sort_formula` includes `300*(eligible_time/86400)`, so
  the job gains priority the longer it waits, and killing it would forfeit that.
- **escalation to `develop`: RULED OUT (2026-07-31 13:01, iteration 3).** The
  earlier plan was to move the orchestrator to the shared `develop` queue if it
  stayed queued. Two measurements killed that idea:
    1. `develop` caps walltime at 21600 s = 6 h exactly. Verified empirically by
       submitting a probe job at 12 h, which PBS rejected with "your declared
       wallclock time (43200 seconds) exceeds your maximum limit of 21600
       seconds". No job was created. The cap is real, not folklore — no
       `resources_max.walltime` is published on the queue or server.
    2. Ensemble members are themselves `main`-queue GPU jobs. `pbs.jl:29` in
       ClimaCalibrate defaults `queue` to `"main"` and `run_calibration.jl` does
       not override it, so each of the 21 members per iteration is an exclusive
       GPU job at 3 h walltime.
  So 5 EKI iterations of 21 contended 3 h GPU jobs will plausibly run well past
  6 h of orchestrator wall time. Moving to `develop` would trade a queue wait
  for near-certain walltime kills, each orphaning that iteration's in-flight
  members. The 12 h window in `main` is worth waiting for. KEEP WAITING.
  Reducing the orchestrator's `ncpus` would not help either: `main` allocates
  exclusive nodes for any request size.
- **notes:** Pre-flight checks before submitting — `pmodel_α` prior is 0.028
  (post-#1817 `1/τ` convention, τ ≈ 36 d), so the convention trap is clear. Both
  `.jld2` masks are no-ops for this stage's targets (`apply_*_mask` returns early
  for any short name other than `lai`), so stage 1 does not need the validity
  mask that stage 2's long run will supply. Orchestrator and member jobs both
  pin `climacommon/2026_04_08`. Iteration 2 additionally verified all 10 prior
  parameter names resolve in `toml/default_parameters.toml` and that every prior
  mean equals the current default, so the prior really is centered on the
  defaults as the config claims. (Note for whoever repeats this check: the TOML
  mixes bare `[name]` and quoted `["name"]` section headers, so a `^\[name\]`
  grep reports false MISSINGs — match both forms.)

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

## Environment gotchas found by the loop

- `qstat -u $USER` returns EMPTY OUTPUT, silently. The sandbox matches its
  `excludedCommands` list against the literal command string, and the `$USER`
  expansion defeats that match, so the command runs sandboxed and prints
  nothing. Use `qstat -u arenchon`. This one is dangerous rather than merely
  annoying: empty output is indistinguishable from "no jobs running", so an
  iteration that trusts it would resubmit a job that is already queued.
  `loop_prompt.md` has been corrected.
