================================================================================
UNSUPERVISED LOOP — prognostic carbon pools (ClimaLand.jl)
Branch: ar/prognostic_carbon   |   Machine: NCAR Derecho (login node, PBS)
================================================================================

You are running UNATTENDED on NCAR Derecho. No human is watching. Each wake-up
you make ONE small, safe, committed increment of progress, then stop.
Correctness and reversibility matter far more than speed. When in doubt, do
less, write down why, and push.

GOAL: a simple prognostic live/dead carbon model in ClimaLand that produces
biomass (kg C m⁻²) of the right magnitude and spatial pattern compared to
observations, at a steady state reached under recycled pre-climate-change
forcing. It does NOT reach perfection; it reaches "reasonable and defensible".

READ `experiments/integrated/prognostic_carbon/MODEL.md` FIRST, every iteration.
That is the model specification — pools, fluxes, the no-PFT climate proxy, the
SOC coupling, and the parameter list. This file is the *process*; MODEL.md is
the *physics*. If you deviate from MODEL.md, update MODEL.md in the same commit
and say why in the PR comment.

--------------------------------------------------------------------------------
THE THREE RULES THAT DEFINE THIS MODEL
--------------------------------------------------------------------------------
1. DO NOT CHANGE GPP OR LAI. Both are already modelled (P-model; Zhou optimal
   LAI or prescribed MODIS LAI). The carbon model consumes them. If a change you
   are about to make alters simulated GPP or LAI, stop — it is out of scope, and
   phase 1 is explicitly one-way-coupled (MODEL.md §5). The one deliberate
   exception is autotrophic respiration, which you SHOULD rework (MODEL.md §2.1).
2. NO PFTs. No per-PFT maps, no `pfts.jl`. Allocation fractions and turnover
   times are GLOBAL CONSTANTS, optionally split C3/C4. Only if constants
   demonstrably fail may they vary — and then with mean annual temperature and
   mean annual precipitation (MODEL.md §3), never with a plant functional type.
   This is a hard user constraint.
3. HEAVY OUTPUT GOES TO SCRATCH: /glade/derecho/scratch/arenchon/claude/...
   Never into the repo. Commit figures and small summary tables only.

--------------------------------------------------------------------------------
OPERATIONAL — Derecho gotchas (inherited from ar/derecho_loop, all still true)
--------------------------------------------------------------------------------
!! PBS BARE: run qsub/qstat/qdel BARE — no `timeout ...`, no `for`-loop wrapping,
   no `bash -lc "..."`, no `r=$(qsub ...)`. They are in the sandbox's
   excludedCommands (unsandboxed); wrapping forces them back INTO the sandbox and
   fails with `qsub: cannot connect to server desched (errno=15008)`. That is NOT
   a Derecho outage. If a bare qsub still fails persistently, ask the user to run
   `! qsub <script>` in their shell. Account = UCIT0011 (not ucit0012).

!! NEVER BLOCK ON A JOB. Submit, record the job id in STATE.md, end the
   iteration. Check `qstat -u arenchon` at the START of the next iteration.

!! HDF5: Derecho Lustre throws intermittent HDF5 file-locking errors on
   NetCDF/HDF5 read/write. Always `export HDF5_USE_FILE_LOCKING=FALSE`, including
   in any analysis Julia call that reads NetCDF.

!! BAD GPU NODE: a GPU job that aborts at LAUNCH (Exit_status=-2,
   resources_used.walltime=00:00:00, EMPTY .OU, PBS mail "Aborted by PBS Server")
   repeatedly on the SAME exec_host is a silently-bad node — PBS still reports it
   free and keeps assigning it. Diagnose with
   `qstat -x -f <jid> | grep -iE "Exit_status|exec_host|walltime"`. Fix: add
   `#PBS -m n`, and resubmit until it lands on a different node.

!! GPU PARALLEL: global GPU runs can run in parallel. When comparing configs
   (e.g. a τ_wood / f_root sweep), submit them ALL AT ONCE rather than serially.

!! ENVIRONMENT: `julia --project=.buildkite` for experiments and simulations;
   `julia --project=test test/runtests.jl` for package tests.

--------------------------------------------------------------------------------
STEP 0 — ORIENT (every iteration, in this order)
--------------------------------------------------------------------------------
1. Read `experiments/integrated/prognostic_carbon/unsupervised_loop/STATE.md`.
   It is the source of truth for where the work is. Trust it over memory.
2. Read `experiments/integrated/prognostic_carbon/MODEL.md`.
3. `qstat -u arenchon`. Reconcile EVERY job in STATE.md's job table:
   - still queued/running -> leave it, do NOT resubmit the same test;
   - finished             -> read its log and outputs, record pass/fail plus the
                             key numbers in STATE.md, and decide the follow-up.
4. `git status` and `git log --oneline -5`.
5. Only then decide what to do.

--------------------------------------------------------------------------------
HOW TO TEST — cheap and parallel by default, global occasionally
--------------------------------------------------------------------------------
The workhorse is a battery of SINGLE-PIXEL ERA5 columns run in PARALLEL on ONE
CPU node. Each column is ~1 core and ~2–3 GB, they are fully independent, and a
20-site battery spanning tropical forest / savanna / desert / temperate forest /
grassland / mediterranean / boreal / tundra covers the whole biomass gradient
for a few core-hours. Every code change gets tested this way FIRST.

Global GPU runs are for milestones only: after a stage passes the battery, and
for the stage-5 parameter sweep (submitted all at once, in parallel). Do not
burn GPU hours debugging something a single column would have caught.

Both the parametrized driver and the battery runner already exist on
`origin/ar/derecho_loop`:
  experiments/integrated/era5/single_site.jl     (SITE_LON/SITE_LAT/SITE_NAME/
                                                  START/STOP/DT via ENV)
  experiments/integrated/generic_site/run_battery.pbs   (NPAR concurrent sites,
                                                  sidecar battery.conf config)
  experiments/integrated/generic_site/test_sites.csv    (the 20-site battery)
PORT them (stage 0), do not copy blindly: that branch predates this one, several
of its `ONLINE_*` env switches are now permanent behaviour here, and its `REPO`
path points at /glade/u/home/arenchon/GitHub/ClimaLand.jl while this checkout is
/glade/u/home/arenchon/ClimaLand.jl.

--------------------------------------------------------------------------------
STAGES — strictly ordered. Do not start one before its predecessor is verified.
--------------------------------------------------------------------------------

STAGE 0 — test harness
  Port and adapt the single-site driver, `run_battery.pbs` and `test_sites.csv`
  from `origin/ar/derecho_loop` onto this branch, fixing the repo path and
  dropping switches that no longer exist. DONE WHEN: the full 20-site battery
  runs green on unmodified code, with per-site output under scratch, and you
  have recorded the baseline GPP / LAI / Ra / Rh per site in STATE.md. That
  baseline is what rule 1 is checked against later.

STAGE 1 — implement the pools (CPU, site battery)
  Add the carbon model to `src/standalone/Vegetation/biomass.jl`. It composes
  with an existing LAI model rather than replacing it: prescribed MODIS LAI and
  `ZhouOptimalLAIModel` must both work underneath it, and neither may change
  behaviour when the carbon model is absent.
  Four prognostic pools: C_sugar, C_leaf, C_stem, C_root (MODEL.md §1–2).
  Sugar is regulated by the allocation ramp of MODEL.md §2.3, not by a hard
  clamp — it should swing seasonally and stay off zero on its own. If it pins at
  zero, that is a signal to investigate, not to clamp harder.
  DONE WHEN: package tests pass; the site battery completes; a CARBON
  CONSERVATION test passes — d(ΣC)/dt equals GPP − Ra − litter to round-off
  (write that test); and GPP and LAI at every site are UNCHANGED from the
  stage-0 baseline.

STAGE 2 — pool-based autotrophic respiration
  New subtype of `AbstractAutotrophicRespirationModel` computing Rm from the
  pools and Rg = (1−a)·S (MODEL.md §2.1–2.2). Keep the JULES model as the
  default constructor until stage 5 says otherwise; select the new one
  explicitly. DONE WHEN: the battery's annual NPP/GPP lands in 0.3–0.6 at the
  forest and grassland sites, and Ra neither collapses to zero nor exceeds GPP
  in the annual mean at any site.

STAGE 3 — turn SOC on
  `dY.soilco2.SOC = I_litter(z) − Sm`. Debit `Sm` from SOC in `MicrobeProduction`
  (today it is added to CO₂ but never removed from SOC — the C balance is open).
  Add the canopy→soil litter coupling in the integrated model. Standalone
  `SoilCO2Model` keeps zero/prescribed litter, and its existing tests must still
  pass. DONE WHEN: soil C conservation test passes, and the battery shows SOC
  drifting slowly (not collapsing, not exploding) with Rh within a factor of ~2
  of its stage-0 baseline at every site.

STAGE 4 — offline spinup to steady state
  Build the offline pool integrator of MODEL.md §6: dump the carbon model's
  drivers to scratch, integrate the pools alone under recycled forcing for as
  long as stem and SOC need, and write an initial-condition file. Develop and
  validate it on the site battery — where the offline integrator can be checked
  against a genuinely long coupled column run — before applying it globally.
  Verify the artifact override actually resolves before trusting it; a silently
  ignored override makes stage 5 look like it worked.
  DONE WHEN: the offline integrator reproduces a coupled run's pools over the
  overlapping period to within a few percent at every battery site, and the
  global equilibrium pools are finite and non-negative everywhere over land.

STAGE 5 — global run, compare to observations, tune
  Global GPU run from the spun-up IC. Compare cVeg to the biomass reference
  dataset and cSoil to SoilGrids. Tune the 4–6 parameters MODEL.md §8 identifies
  — τ_stem, f_stem, f_root, r_stem, a — with a parallel sweep, not one run at a
  time. Do the coarse sweep on the site battery (cheap, and biomass at 20 sites
  already separates the good parameter sets from the bad), then confirm the best
  few globally. Report RMSE and bias by broad climate zone (tropics / temperate /
  boreal / semi-arid), and attach the maps.
  DONE WHEN: global biomass is in the right ballpark and the pattern is
  defensible, OR you can say precisely which term is wrong and why.

--------------------------------------------------------------------------------
FAILURE SIGNATURES — check before believing any result
--------------------------------------------------------------------------------
  - C_sugar pinned at zero everywhere: a real imbalance, not a numerical one —
    Rm too big, or the GPP unit conversion wrong. Check mol→kg C
    (M_C = 0.012011) before touching parameters, and do not "fix" it with a
    clamp. Sugar should swing seasonally; permanently empty means the plant is
    being asked to pay more maintenance than it earns.
  - C_stem growing without bound: turnover is not being applied, or τ_stem came
    out in the wrong units (seconds vs years).
  - σl_implied = C_leaf/LAI far from ~0.03–0.1 kg C m⁻² leaf: the allocation
    fractions and the LAI model disagree. Report it; it is the main diagnostic
    that the constant-fraction assumption is or is not holding.
  - cVeg two orders of magnitude off: almost always a unit error, not a
    parameter error. Check units first, always.
  - SOC collapsing to zero within a few years: litter input is not reaching the
    soil, or `Sm` is being double-debited.
  - A run that "completes" with frozen parameters and NaN loss usually means
    every member CRASHED. Read a member's log; a swallowed GPU compile error
    looks exactly like a numerical NaN from the outside.
  - Any change in GPP or LAI relative to the same run without the carbon model:
    that is a bug in phase 1 by definition. Check it explicitly at each stage.

--------------------------------------------------------------------------------
REPORTING
--------------------------------------------------------------------------------
Post a PR comment at the end of each iteration that changed something material:
what you did, the numbers, what you concluded, what is next. Show results, not
just descriptions of code. Prune to the most recent ~10 comments.

`STATE.md` is where the work is; `LOG.md` is how it got there — one dated entry
per iteration. Commit both every iteration.

--------------------------------------------------------------------------------
RULES
--------------------------------------------------------------------------------
  - `[skip ci]` in EVERY bookkeeping commit.
  - Format before committing: `julia -e 'using JuliaFormatter; format(".")'`.
  - Never force-push, never `git reset --hard`, never merge or close the PR,
    never push to main.
  - Exported symbols need a `@docs` entry in `docs/src/APIs/*.md` — the docs
    build uses `checkdocs = :exports` and will fail CI otherwise.
  - If you are stuck on the SAME failure for THREE consecutive iterations, stop
    trying variations. Write up what you tried and what you think is wrong, and
    say clearly in the PR comment that you are blocked and why.
  - Report negative results honestly. A stage that did not work is a finding.
    Do not present a failed result as a success with caveats.

================================================================================
