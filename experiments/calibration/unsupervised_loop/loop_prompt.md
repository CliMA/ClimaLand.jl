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
STAGE 1 — calibrate GPP + energy fluxes
--------------------------------------------------------------------------------
Config: `experiments/calibration/configs/sifgpp_lhf_shf_lwu_rosetta.jl`
Targets: sif_gpp (GOSIF GPP), lhf, shf, lwu (ERA5)
Parameters: 10 (5 P-model/moisture-stress, 5 canopy turbulent/radiative)
Ensemble: TransformUnscented, 10*2+1 = 21 members
Runs with PRESCRIBED MODIS LAI — the canopy is observed, so GPP and the fluxes
are fitted without prognostic-LAI error contaminating them.

  CALIBRATION_CONFIG=sifgpp_lhf_shf_lwu_rosetta.jl \
    bash experiments/calibration/run_calibration.sh

Set the task count in the batch script to 21.

DONE WHEN: the EKI has run its 10 iterations, the final-iteration parameters are
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
        OUTPUT_ROOT=/glade/derecho/scratch/$USER   keep multi-GB output off home
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
Ensemble: TransformUnscented, 5*2+1 = 11 members. Set tasks to 11.

  CALIBRATION_CONFIG=lai.jl bash experiments/calibration/run_calibration.sh

Prerequisites, both of which you must confirm from STATE.md before starting:
  - stage 1's parameters are committed in `toml/default_parameters.toml` (they
    set the potential GPP that drives LAI; calibrating LAI against an
    uncalibrated GPP is meaningless);
  - stage 2's artifact override is active.

DONE WHEN: the EKI has converged and the calibrated optimal-LAI parameters are
written into `toml/default_parameters.toml` (3 significant figures) and committed
with `[skip ci]`, with the before/after values recorded in STATE.md.

`build_model_lai_validity_mask.jl` builds a grid-specific
`model_lai_validity_mask.jld2` that the LAI observations are masked by. It is NOT
committed. Build it for this grid BEFORE generating observations, or
`generate_observations.jl` will error.

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
  - Submit `qsub` bare, not inside a for-loop (the loop form gets sandboxed).
  - Global GPU runs can be submitted in parallel; do so when comparing configs.
  - If you are stuck on the same failure for THREE consecutive iterations, stop
    trying variations. Write up what you tried and what you think is wrong, and
    say clearly in the PR comment that you are blocked and why.
  - Report negative results honestly. A stage that did not work is a finding.
    Do not present a failed calibration as a success with caveats.

================================================================================
