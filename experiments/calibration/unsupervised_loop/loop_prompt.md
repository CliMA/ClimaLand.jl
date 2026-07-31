================================================================================
UNSUPERVISED CALIBRATION LOOP — prognostic-LAI pipeline on Derecho
================================================================================

You are running unattended on NCAR Derecho, working on the branch
`ar/calibrate_lai_pipeline`. Your job is a three-stage calibration pipeline that
ends with a calibrated prognostic-LAI ClimaLand configuration.

The stages are STRICTLY SEQUENTIAL. Stage 2 consumes stage 1's parameters; stage
3 consumes both stage 1's parameters and stage 2's initial-condition artifact.
Do not start a stage before its predecessor has a recorded, verified result.

Each stage is gated by multi-hour PBS jobs. You do NOT sit and wait: submit,
record the job id in the state file, exit the iteration, and check `qstat` at the
start of the next iteration. Never block on a job.

--------------------------------------------------------------------------------
STATE FILE — read this FIRST, every iteration
--------------------------------------------------------------------------------
`experiments/calibration/unsupervised_loop/STATE.md` is the single source of
truth for where the pipeline is. It records, per stage: status
(not_started / running / done / failed), the PBS job ids submitted, the output
paths, and the resulting parameters or artifact path.

At the start of EVERY iteration:
  1. Read STATE.md.
  2. `qstat -u $USER` to see what is still running.
  3. For each job that finished since last iteration, read its log, decide
     success or failure, and update STATE.md.
  4. Only then decide what to do next.

Commit STATE.md every iteration with `[skip ci]`. It is the handoff between
iterations — if it is wrong or stale, the loop loses its place.

--------------------------------------------------------------------------------
ENVIRONMENT — how to run things here
--------------------------------------------------------------------------------
WHERE OUTPUT GOES. Everything heavy — calibration output directories, long-run
diagnostics, rebuilt artifacts, scratch scripts, plots, extracted data — goes
under `/glade/derecho/scratch/arenchon/claude/`. Never write multi-GB output to
`$HOME` (small quota) or into the repo checkout. Name calibration directories so
the stage is obvious, e.g.
  /glade/derecho/scratch/arenchon/claude/calibration_stage1_gpp_energy
  /glade/derecho/scratch/arenchon/claude/calibration_stage3_lai

LAUNCHING A CALIBRATION — AS A PBS JOB, NOT A BACKGROUND PROCESS. The
orchestrator is a long-lived foreground process that submits one PBS job per
ensemble member and waits for them. You CANNOT run it with `nohup ... &`: in
this sandbox a detached process is killed as soon as the shell that started it
exits (verified — it does not survive to the next tool call). So the
orchestrator itself goes into PBS, using the wrapper committed for this:

  qsub -v CALIBRATION_CONFIG=lai.jl,CAL_OUTPUT_DIR=/glade/derecho/scratch/arenchon/claude/calibration_stage3_lai \
      experiments/calibration/calibration_orchestrator.pbs

Run that from the repository root, bare (see the wrapping trap below). The
wrapper loads `climacommon/2025_02_25` — the same module the forward-model
member jobs are pinned to in `run_calibration.jl` — sets
`HDF5_USE_FILE_LOCKING=FALSE`, and calls `run_calibration.sh`. Submitting jobs
from inside a job works on Derecho; that has been verified, so member jobs
appear in `qstat` under your user as usual.

Record the orchestrator's job id AND the member job ids in STATE.md. The
orchestrator's own log is `clima_calibration.o<jobid>` in the submission
directory.

If the orchestrator hits its 12-hour walltime before finishing all iterations,
resubmit exactly the same command. ClimaCalibrate restarts at
`last_completed_iteration + 1`, so finished iterations are not repeated. Treat a
walltime kill as normal progress, not a failure.

Short interactive work — building the masks, inspecting output — runs fine on
the login node with `julia --project=.buildkite`. Generating the LAI observation
vector does NOT: two attempts were killed with zero output, which is what the
login-node resource arbiter looks like (loading MODIS LAI 2000-2020 is heavy).
Run anything of that weight as a `develop`-queue job instead. This does not
affect the calibration itself — `run_calibration.sh` regenerates the
observations inside the orchestrator job, where it takes about 3 minutes.

A killed process on the login node leaves NO error message, because Julia's
buffered output dies with it. Do not read that silence as "the job is still
running" or as a code bug — check whether the process still exists.

PBS — THE THREE FALSE ALARMS. In the previous unattended run, three separate
job-submission failures were each misread as "Derecho is down for maintenance".
None of them were. Before you ever conclude that the machine is unavailable,
rule out all three:
  1. WRONG ACCOUNT. The allocation is `UCIT0011`. `ucit0012` is not authorized
     and PBS rejects it with "Invalid account".
  2. WRAPPED COMMANDS. Run `qsub`, `qstat`, `qdel`, `qhist` BARE. They are in
     the sandbox's `excludedCommands` list, and that list is matched against the
     literal command, so wrapping one in `timeout ...`, `bash -c ...`, `$(...)`
     or a `for` loop makes it sandboxed and it fails (errno=15008).
  3. TARGETING A SPECIFIC NODE. Do not request a named host or a specific GPU
     node. Let PBS schedule; a single bad node is not an outage.
Only after all three are excluded, and `qstat -B` or the NCAR status notice
actually says so, may you record "maintenance" as a blocker.

SANDBOX. You run sandboxed even with --dangerously-skip-permissions. Writes are
allowed to the repo, `~/.julia`, `~/.cache`, `/glade/derecho/scratch/arenchon`
and `/tmp`; network access is limited to github and julialang. You cannot
download datasets. If something seems to need new external data, say so in
STATE.md as a blocker rather than trying to fetch it.

--------------------------------------------------------------------------------
STAGE 1 — calibrate GPP + energy fluxes
--------------------------------------------------------------------------------
Config: `experiments/calibration/configs/sifgpp_lhf_shf_lwu_rosetta.jl`
Targets: sif_gpp (GOSIF GPP), lhf, shf, lwu (ERA5)
Parameters: 10 (5 P-model/moisture-stress, 5 canopy turbulent/radiative)
Ensemble: TransformUnscented, 10*2+1 = 21 members
Runs with PRESCRIBED MODIS LAI — the canopy is observed, so GPP and the fluxes
are fitted without prognostic-LAI error contaminating them.

  qsub -v CALIBRATION_CONFIG=sifgpp_lhf_shf_lwu_rosetta.jl,CAL_OUTPUT_DIR=/glade/derecho/scratch/arenchon/claude/calibration_stage1_gpp_energy \
      experiments/calibration/calibration_orchestrator.pbs

Nothing to set by hand for the ensemble size: on the Derecho backend
ClimaCalibrate submits one single-task PBS job per member, so 21 member jobs
appear in `qstat` per iteration. (The "match the task count in the batch script"
note in `run_calibration.jl` applies only to the Slurm/WorkerBackend test path.)

Calibrated over ALL land cells. The natural-vegetation mask is stage 3 only —
do not apply it here.

DONE WHEN: the EKI has run its 5 iterations, the final-iteration parameters are
extracted, AND they are written into `toml/default_parameters.toml` (3 significant
figures, matching how the previously calibrated values are stored there). Commit
that TOML change with `[skip ci]` and record the values in STATE.md.

Writing them into the TOML is what makes stages 2 and 3 inherit the calibrated
GPP automatically — there is no override file to keep in sync.

Sanity-check before declaring success:
  - GPP is not uniformly zero or NaN (see the failure signatures below);
  - the loss decreased across iterations;
  - no parameter sits pinned at a prior bound (if one does, say so in STATE.md
    and in the PR comment — it means the prior was mis-specified, and that is a
    finding, not something to quietly widen).

--------------------------------------------------------------------------------
STAGE 2 — rebuild the prognostic-LAI initial conditions
--------------------------------------------------------------------------------
The stock `optimal_lai_inputs` artifact holds a climatology that does not match
this model's own equilibrium. Stage 3's spinup is only 1 year, which is shorter
than `optimal_lai_tau_long_term` (2 years), so it relies on starting near
equilibrium. That is what this stage produces.

  2a. Confirm stage 1's calibrated parameters are already in
      `toml/default_parameters.toml` (stage 1 commits them). The long run picks
      them up with no further action. If they are not there, stage 1 is not
      actually done — go back and finish it.

  2b. Run `experiments/long_runs/snowy_land_pmodel.jl` with PROGNOSTIC_LAI set,
      for 10 years, on GPU. The script honours these env vars:
        PROGNOSTIC_LAI=""   select the optimal-LAI model
        RUN_YEARS=10        length in years
        OUTPUT_ROOT=/glade/derecho/scratch/arenchon/claude   keeps the multi-GB
                            diagnostics off home, alongside the earlier global
                            runs in that directory
      With PROGNOSTIC_LAI set the script already writes the six IC diagnostics
      listed in 2c, on top of the short set — no extra flag needed.

  2c. The IC file needs SIX fields, and all six now have diagnostics:
        lai_init      <- `lai`
        a0_annual     <- `a0a`
        precip_annual <- `pra`    (NOTE: `pra` is m yr^-1; the IC file wants
                                   mol H2O m^-2 yr^-1 — convert via the molar
                                   liquid density, do not write it through raw)
        f0            <- `olf0`
        vpd_gs        <- `olvpd`
        gsl           <- `olgsl`
      `olf0`/`olvpd`/`olgsl` were added for exactly this purpose. Take the LAST
      year of the 10, so the trailing totals have equilibrated.

  2d. Write a new `optimal_lai_inputs.nc` with those six variables on the
      (lon, lat) grid, matching the variable NAMES the reader expects:
      "gsl", "a0_annual", "precip_annual", "vpd_gs", "lai_init", "f0"
      (see `optimal_lai_initial_conditions` in
      `src/standalone/Vegetation/spatially_varying_parameters.jl`).

  2e. Point Julia at it via `~/.julia/artifacts/Overrides.toml`, then VERIFY the
      override is picked up before trusting it: start Julia and check that
      `ClimaLand.Artifacts.optimal_lai_initial_conditions_path()` resolves to
      your new file. A silently-ignored override would make stage 3 look like it
      worked while using the stock IC.

  2f. Sanity-check the new IC against the old one: global mean LAI, the fraction
      of land cells with non-finite values, and the aridity index implied by
      f0. A rebuilt IC that is mostly NaN or has collapsed LAI is a stage-2
      failure, not a stage-3 input.

DONE WHEN: the override resolves, the six fields are finite over land, and
STATE.md records the artifact path.

--------------------------------------------------------------------------------
STAGE 3 — calibrate prognostic LAI
--------------------------------------------------------------------------------
Config: `experiments/calibration/configs/lai.jl`
Target: `lai` (MODIS, via `get_modis_lai_obs_var`)
Parameters: 5 — optimal_lai_z, z_c4, sigma, sigma_c4, alpha
Ensemble: TransformUnscented, 5*2+1 = 11 members, so 11 member jobs per
iteration in `qstat`. Nothing to set by hand.

  qsub -v CALIBRATION_CONFIG=lai.jl,CAL_OUTPUT_DIR=/glade/derecho/scratch/arenchon/claude/calibration_stage3_lai \
      experiments/calibration/calibration_orchestrator.pbs

Build BOTH masks before this (see below) — `run_calibration.sh` regenerates the
observation vector on every launch, and it errors if a mask is missing.

Prerequisites, both of which you must confirm from STATE.md before starting:
  - stage 1's parameters are committed in `toml/default_parameters.toml` (they
    set the potential GPP that drives LAI; calibrating LAI against an
    uncalibrated GPP is meaningless);
  - stage 2's artifact override is active.

DONE WHEN: the EKI has converged and the calibrated optimal-LAI parameters are
written into `toml/default_parameters.toml` (3 significant figures) and committed
with `[skip ci]`, with the before/after values recorded in STATE.md.

TWO MASKS, BOTH REQUIRED BEFORE `generate_observations.jl` RUNS. Neither is
committed (`*.jld2` is gitignored) and both are grid-specific — rebuild them
whenever `nelements` changes. `generate_observations.jl` errors if either is
missing, deliberately, so a stale or absent mask cannot quietly change what is
being calibrated.

  1. `experiments/calibration/model_lai_validity_mask.jld2`, from
     `build_model_lai_validity_mask.jl <reference_simdir>`. Aligns the MODIS obs
     with the simulated-LAI footprint. Needs a reference forward-model run:
     build it from STAGE 2's long-run diagnostics, which are on this grid and
     already exist by the time stage 3 starts. (A member directory from a
     stage-3 iteration also works, but that is circular — observations are
     generated before the first iteration runs.)

  2. `experiments/calibration/natural_vegetation_mask.jld2`, from
     `build_natural_vegetation_mask.jl`. Restricts the LAI target to natural
     undisturbed vegetation: CLM `PCT_NATVEG` >= 80% (which drops cropland,
     urban, lake, wetland and glacier in one criterion) AND GFED4.1s mean annual
     burned fraction <= 5 %/yr (which drops habitually-burning savanna). It
     keeps deserts — they are undisturbed, and the model's arid over-prediction
     is a real residual that belongs in the score. It removes ~43% of land area,
     leaving ~10k of 22k land cells, which is ample for 5 parameters. Takes no
     arguments and needs no simulation; build it early.

The natural-vegetation mask applies to the `lai` target ONLY, by construction
(`apply_natural_vegetation_mask` is a no-op for other short names). Stage 1's
flux targets stay global. When you report the stage-3 fit, say explicitly that
the score is over natural vegetation only, and note the two caveats: the CLM
surface dataset is `simyr2000` so post-2000 cropland expansion is not captured,
and GFED4.1s ends in 2016 so the fire criterion is climatological rather than
specific to the calibrated year.

ITERATIONS. Both configs are set to `n_iterations = 5`, which should be enough.
Check convergence before declaring done: the loss should have flattened over the
last two iterations and the parameter means should have stopped moving
materially. If it has NOT converged, raise `n_iterations` in the config, commit
that with `[skip ci]`, and resubmit the orchestrator with the SAME
`CAL_OUTPUT_DIR` — ClimaCalibrate resumes from `last_completed_iteration + 1`,
so it picks up at iteration 6 rather than starting over. Say in STATE.md and the
PR comment that you extended it and why.

Knobs, in the order to try them if the fit is poor:
  - NOISE_SCALARS["lai"] (currently 0.5): raise toward 1.0 if EKI overfits the
    seasonal amplitude at the cost of timing; lower toward 0.25 to pull harder
    on peak magnitude.
  - spinup: raise to Year(2)+ if the fit looks like it is chasing a spin-up
    transient (diagnose by comparing the first and last scored months).
  - the priors: widen only as a last resort, and say why in the PR comment.

--------------------------------------------------------------------------------
FAILURE SIGNATURES — check for these before believing any result
--------------------------------------------------------------------------------
A calibration that "completes" can still be a total failure. Before recording a
stage as done:

  - Frozen parameters + NaN loss across all iterations usually means every
    ensemble member CRASHED, not that the model is insensitive. Read a member's
    `model_log.txt` — a swallowed GPU compile error looks exactly like a
    numerical NaN from the EKI side.
  - All-zero or all-NaN GPP: check the member logs before touching priors.
  - A parameter pinned at a bound: the prior is mis-specified. Report it.
  - `pmodel_α` in the WRONG CONVENTION is the specific trap on this branch. Since
    PR #1817 it means 1/tau, not 1 - 1/tau. Any prior or parameter value you
    copy from an older branch or an old EKI output must be converted with
    alpha_new = 1 - alpha_old. A value of 0.97 would mean tau ~ 1 day instead of
    ~36 days and would quietly ruin stage 1.

--------------------------------------------------------------------------------
REPORTING
--------------------------------------------------------------------------------
Post a PR comment at the end of each iteration that changed something material:
what you did, what the numbers were, what you concluded, what is next. Prune to
the most recent ~10 comments so the PR stays readable.

Write the durable narrative into
`experiments/calibration/unsupervised_loop/LOG.md` — one dated entry per
iteration. STATE.md is where the pipeline is; LOG.md is how it got there.

--------------------------------------------------------------------------------
RULES
--------------------------------------------------------------------------------
  - `[skip ci]` in EVERY commit message. These commits are bookkeeping; CI on
    them wastes runners.
  - DO commit calibrated parameters into `toml/default_parameters.toml` at the
    end of each calibration stage, at 3 significant figures, in their own
    `[skip ci]` commit. Say in the PR comment which stage produced them and what
    the before/after values were — a default change should never be silent.
    Do not hand-edit any other parameter while you are in there.
  - Never force-push, never `git reset --hard`, never merge or close the PR.
  - All heavy output goes under `/glade/derecho/scratch/arenchon/claude/`;
    nothing large in `$HOME` or in the checkout.
  - Submit `qsub`/`qstat`/`qdel` BARE, and to account UCIT0011. See the three
    false alarms in the ENVIRONMENT section before blaming the machine.
  - Global GPU runs can be submitted in parallel; do so when comparing configs.
  - If you are stuck on the same failure for THREE consecutive iterations, stop
    trying variations. Write up what you tried and what you think is wrong, and
    say clearly in the PR comment that you are blocked and why.
  - Report negative results honestly. A stage that did not work is a finding.
    Do not present a failed calibration as a success with caveats.

================================================================================
