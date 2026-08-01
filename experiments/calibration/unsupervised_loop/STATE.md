# Pipeline state

Single source of truth for where the three-stage calibration pipeline is.
The loop reads this first every iteration and commits it every iteration.
See `loop_prompt.md` for what each stage means and how it is gated.

Status values: `not_started` | `running` | `done` | `failed`

All heavy output lives under `/glade/derecho/scratch/arenchon/claude/`.

---

## Stage 1 — GPP + energy fluxes

- **status:** running (CHUNKED — see below)
- **config:** `experiments/calibration/configs/sifgpp_lhf_shf_lwu_rosetta.jl`
- **ensemble:** 21 (10 params, TransformUnscented)
- **PBS job ids:**
  - `6967486.desched1` — original `main`-queue submission, 12 h walltime.
    **CANCELLED (`qdel`) 2026-07-31 18:16 after 7 h 26 m queued without ever
    starting.** Never ran, so it produced no output file and wrote nothing.
  - `6970028.desched1` — **CURRENT**, chunk 1 (`N_ITERATIONS=1`), queue
    `develop`→`cpudev`, 5 h 30 m walltime, `select=1:ncpus=4:mem=32gb`. Submitted
    18:17, **started within ~2 minutes**, running on `dec2436`.
  - Member job ids: appear once the orchestrator reaches the ensemble stage.
- **chunking scheme:** each orchestrator job advances the run by ONE iteration
  and exits normally, instead of running to 5 and risking a walltime kill.
  Resubmit with `N_ITERATIONS=2`, then 3, 4, 5 — ClimaCalibrate resumes at
  `last_completed_iteration + 1` (`backends.jl:277`) and stops at
  `n_iterations`. Command (note the command-line overrides beat the script's
  `#PBS -q main` / 12 h directives, verified):

      qsub -q develop -l walltime=05:30:00 -l select=1:ncpus=4:mem=32gb \
          -v CALIBRATION_CONFIG=sifgpp_lhf_shf_lwu_rosetta.jl,CAL_OUTPUT_DIR=/glade/derecho/scratch/arenchon/claude/calibration_stage1_gpp_energy,N_ITERATIONS=<n> \
          experiments/calibration/calibration_orchestrator.pbs
- **output dir:** `/glade/derecho/scratch/arenchon/claude/calibration_stage1_gpp_energy`
- **orchestrator log:** `clima_calibration.o<jobid>` in the repo root, written
  LIVE as the job runs
- **chunk 1 progress (2026-07-31 19:02, 44 min elapsed):** still precompiling,
  no output dir yet, no failures (`✗` count is 0). `run_calibration.sh` begins
  with `Pkg.update()`, which upgraded ClimaCalibrate 0.3.1→0.3.2 and ClimaCore
  0.14.54→0.14.55 and so forced a full rebuild including CUDA. The log's "24
  dependencies failed but may be precompilable after restarting julia" is the
  normal parallel-precompile cache race, not an error: `run_calibration.sh` uses
  three separate `julia` invocations, so the stragglers precompile on the next
  one. Only `.buildkite/Manifest-v1.12.toml` is dirtied; `Project.toml` is not.
  NOTE for chunking: this precompile cost should be paid ONCE, by chunk 1 — the
  depot (`~/.julia`) is shared, so later chunks' `Pkg.update()` is a no-op and
  the caches are warm. If chunk 2 also spends ~1 h precompiling, that assumption
  is wrong and chunking becomes expensive; check it.
- **chunk 1 milestone (2026-07-31 19:44, 1 h 27 m elapsed):** precompile cleared,
  observations built (`land_observation_vector_sifgpp_lhf_shf_lwu.jld2`, 110 MB,
  19:12 — gitignored, so it will not be committed), output dir created 19:25 with
  `interface.jld2`, `prior.jld2`, `eki_file.jld2` and all 21 member dirs, and
  **21 member GPU jobs submitted** (`run_1_1`…`run_1_21` = 10*2+1, as expected).
  Confirms end-to-end that a `develop`-queue orchestrator can submit `main`-queue
  GPU members.
- **⚠ WALLTIME RISK ON CHUNK 1.** Only ~6 members run concurrently (GPU queue at
  77/82 nodes, 26 queued). At the reference ~55 min per member that is 4 waves
  ≈ 3 h 40 m, on top of 1 h 27 m already spent ≈ 5 h 10 m against a 5 h 30 m
  walltime. Margin is roughly 20-40 min. If the orchestrator IS killed:
    1. Do NOT resubmit immediately. Check `qstat -u arenchon` for surviving
       `run_1_*` members and let them finish or `qdel` them first — otherwise the
       resubmitted orchestrator races them for the same member directories.
    2. Then resubmit with `N_ITERATIONS=1`; it will redo iteration 1 from
       scratch (21 fresh member jobs), because an incomplete iteration never
       counts as completed.
  NEXT CHUNKS SHOULD REQUEST MORE WALLTIME than 5 h 30 m is worth — but the
  `develop` cap is a hard 6 h, so the real lesson is that one iteration per chunk
  is near the limit of what fits when GPU concurrency is ~6. If chunk 1 does
  finish comfortably, keep one iteration per chunk; if it is killed, the pipeline
  needs `main`-queue GPU concurrency to improve, or a longer-walltime queue for
  the orchestrator.
- **measured member cost (2026-07-31 20:31): ~1 h 10 m per member**, not the 55 m
  read off the reference run's directory mtimes. `member_001` reported 90 %
  complete at 20:29 with `wall_time_spent = 46 min 58 s` and
  `estimated_sypd = 82.8`, i.e. ~52 min of simulation plus ~16 min of
  startup/compile (member dir created 19:25, simulation began ~19:41). The model
  is healthy — a crashed member would not report sypd or a rising
  `percent_complete`.
  Concurrency rose from 6 to 10 once GPU nodes freed. At 10 concurrent, 21
  members is 3 waves ≈ 3 h 30 m from 19:25, so iteration 1 should finish ~22:55
  against the 23:47 walltime — roughly 50 min of margin. Watch it, but it should
  land. `checkpoint.txt` goes `started` → `completed`; that is the cheap way to
  count finished members without parsing logs.
- **⚠ MEMBER 3 CRASHED — iteration 1 of stage 1 (2026-07-31 ~20:35).**
  `member_003` died with **signal 11 (SIGSEGV), core dumped**, on node
  `deg0069`, at 90 % complete after 47 min of healthy running
  (`estimated_sypd` ~83, `percent_complete` rising normally right up to the
  crash). Its `checkpoint.txt` is stuck at `started` while its PBS job is gone —
  that mismatch, job absent but checkpoint not `completed`, is the signature of
  a crashed member.
  **This is almost certainly transient/hardware, not the parameters.** A bad
  parameter set produces NaN or a solver failure, not a segfault after 47
  minutes at normal throughput.
  **It is handled, and iteration 1 will still complete:**
    - `observation_map.jl:65-70` wraps member processing in try/catch and calls
      `fill_g_ens_col!(..., NaN)` on failure, logging "Error processing member
      $m, filling observation map entry with NaNs". So it degrades rather than
      crashing the orchestrator.
    - EKP's `TransformUnscented` defaults to
      `failure_handler_method = SampleSuccGauss()`
      (`EnsembleKalmanProcess.jl:123-129`); `run_calibration.jl` passes only
      `scheduler`, so that default applies and the update is conditioned on the
      successful particles.
    - 1 failure in 21 is ~4.8 %.
  **WHAT TO WATCH:** if further members segfault in later iterations, especially
  on different nodes, it stops being transient and becomes a real defect worth
  reporting rather than absorbing. Count crashed members per iteration via the
  checkpoint/job mismatch above. Losing a sigma point degrades the unscented
  covariance estimate, so several failures in one iteration would make that
  iteration's update untrustworthy even though EKP will not error.
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
  DECISIVE HAZARD (noted iteration 6), beyond the wasted-time argument: member
  jobs are independent PBS jobs, so killing the orchestrator does NOT kill them.
  They keep running as orphans. A resubmitted orchestrator resumes at
  `last_completed_iteration + 1`, and an iteration cut short never completed, so
  it resubmits that iteration's members — while the orphans are still running
  and writing. Two live member jobs would then target the same
  `iteration_XXX/member_YYY` directory. That is a correctness risk (interleaved
  or corrupted output), not merely wasted GPU time, and it applies to ANY
  walltime kill. It is tolerable once at the end of a 12 h window; it is not
  something to court every 6 h by design.
- **timing evidence (iteration 4).** A previous run of essentially this same
  stage-1 config exists at
  `/glade/derecho/scratch/arenchon/calibration_gpp_energy_rosetta_rebased_v1`
  (2026-06-19). Per-iteration wall times from directory mtimes: 001→002 6 h 48 m
  (cold start: precompile + queue), then 2 h 07 m, 2 h 11 m, 1 h 14 m, 1 h 11 m,
  1 h 12 m, 1 h 49 m. A single member took ~55 min (`iteration_005/member_001`
  `model_log.txt` spans 17:14→18:09). So EXPECT STAGE 1 TO TAKE ~6-10 h of
  orchestrator wall time for 5 iterations, dominated by the cold-start first
  iteration. This confirms the queue decision quantitatively: it fits the 12 h
  `main` window and would certainly be killed by `develop`'s 6 h cap. If the run
  is still on iteration 1 more than ~7 h after the orchestrator starts,
  something is wrong — check member `model_log.txt` rather than waiting.
- **precedent worth knowing:** that reference run went to 9 iterations, while
  this config sets `n_iterations = 5`. That is not evidence 5 is wrong — the
  older run may simply have been left to run on — but when this run reaches
  iteration 5, check convergence properly (loss flat over the last two, means
  stopped moving) before declaring done, and be prepared to extend per the
  ITERATIONS section of `loop_prompt.md`. The same directory is also the natural
  reference for sanity-checking the calibrated parameter values.
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
- **pre-flight done 2026-07-31 (iteration 5), while stage 1 was queued —
  machinery verified, no compute used:**
  - All six IC diagnostics are defined: `lai` comes from the base `CanopyModel`
    list (`default_diagnostics.jl:510`); `a0a`, `pra`, `olf0`, `olvpd`, `olgsl`
    are appended by the `add_diagnostics!` method dispatched on
    `ZhouOptimalLAIModel` (`default_diagnostics.jl:443-453`), which is exactly
    what `PROGNOSTIC_LAI` selects.
  - `snowy_land_pmodel.jl:166-176` requests precisely those six, unioned with
    `get_short_diagnostics(model)`, when `PROGNOSTIC_LAI` is set. So the prompt's
    claim that no extra flag is needed is correct.
  - Env vars behave as documented: `PROGNOSTIC_LAI` is `haskey`-tested (so any
    value, including empty, enables it), `RUN_YEARS` defaults to 2, `OUTPUT_ROOT`
    defaults to `.`.
  - **UNIT CONVERSION for `precip_annual`, confirmed on both sides.** The `pra`
    diagnostic is `units = "m yr^-1"` (`define_diagnostics.jl:374`), but the IC
    reader documents `precip_annual` as `mol H2O m^-2 yr^-1`
    (`spatially_varying_parameters.jl:406`). Multiply by ρ_liq / M_H2O =
    1000 kg m^-3 / 0.018015 kg mol^-1 ≈ **5.5509e4 mol m^-3**. Writing `pra`
    through raw would understate precipitation by ~4.7 orders of magnitude and
    make every cell look hyper-arid.
  - Not yet verified (needs the actual run): that the six fields are finite over
    land, and that the artifact override resolves.
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

- **A `develop`-queue job starts in ~2 minutes; a `main`-queue job waited 7 h 26 m
  and never started.** Measured 2026-07-31 with the CPU queue at 1461 jobs
  queued / 2411 of ~2488 nodes assigned. The GPU queue was nearly EMPTY at the
  same moment (4 queued, 26 running), so the CPU orchestrator was the only real
  bottleneck — the GPU member jobs it submits were never the constraint. Check
  BOTH queues before concluding the machine is busy.
- **Command-line `qsub` flags override the script's `#PBS` directives.** Verified
  with a probe carrying `#PBS -q main` and `#PBS -l walltime=12:00:00`, submitted
  as `qsub -q develop -l walltime=05:30:00`: it landed in `cpudev` and ran. So
  the queue and walltime can be changed without editing
  `calibration_orchestrator.pbs`.
- **`qstat <jobid>` can transiently report "Unknown Job Id"** for a few seconds
  after submission, while the job routes from `develop` to `cpudev`. It is not a
  failure. Re-check with `qstat -u arenchon` before concluding anything.
- **The orchestrator's PBS output file is written LIVE, not at job end.**
  `clima_calibration.o<jobid>` appears in the repo root seconds after the job
  starts and grows as it runs. So "wait for the file to exist" does NOT detect
  completion — watch for the `Orchestrator finished` line (and failure markers)
  in its contents instead.

- `qstat -u $USER` returns EMPTY OUTPUT, silently. The sandbox matches its
  `excludedCommands` list against the literal command string, and the `$USER`
  expansion defeats that match, so the command runs sandboxed and prints
  nothing. Use `qstat -u arenchon`. This one is dangerous rather than merely
  annoying: empty output is indistinguishable from "no jobs running", so an
  iteration that trusts it would resubmit a job that is already queued.
  `loop_prompt.md` has been corrected.
