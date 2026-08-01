# Iteration log

One dated entry per loop iteration: what was tried, what the numbers were, what
was concluded. `STATE.md` is where the pipeline is; this is how it got there.

---

## 2026-07-31 — iteration 1: stage 1 submitted

Pipeline was at `not_started` across all three stages, with nothing in the
queue, so this iteration's job was to get stage 1 moving.

**Pre-flight.** Checked the two things most likely to make a "successful"
stage-1 run meaningless:

- `pmodel_α` prior is `constrained_gaussian("pmodel_α", 0.028, 0.012, 0.001,
  0.15)` — the post-#1817 `1/τ` convention (τ ≈ 36 d), not the old `1 - 1/τ`
  form that would read 0.97 and mean τ ≈ 1 d. The config carries a comment
  saying as much. Trap cleared.
- Neither `.jld2` mask is needed for this stage. `apply_model_validity_mask`
  and `apply_natural_vegetation_mask` in `generate_observations.jl` both return
  early for any short name other than `lai`, and stage 1's targets are
  `sif_gpp`/`lhf`/`shf`/`lwu`. So the missing validity mask — which cannot be
  built until stage 2's long run exists — does not block stage 1.

**Submitted.** Orchestrator job `6967486.desched1`, queue `cpu`, 12 h walltime,
output to `/glade/derecho/scratch/arenchon/claude/calibration_stage1_gpp_energy`.
Confirmed queued (state `Q`).

**Two documentation corrections.** Both in `loop_prompt.md`:

1. `qstat -u $USER` returns empty output, silently — and the prompt told every
   iteration to use it as the "what is still running" check. The sandbox matches
   `excludedCommands` against the literal command string, so the `$USER`
   expansion defeats the match and the command runs sandboxed, printing nothing.
   `qstat -u arenchon` works. This is the dangerous kind of wrong: empty output
   reads as "no jobs running", so the next iteration would have resubmitted a
   job that was already queued. Caught it only because a submission I had just
   made failed to appear. Prompt now spells the username out and explains why.
2. The prompt said the orchestrator wrapper loads `climacommon/2025_02_25`. It
   loads `2026_04_08`, which is what `run_calibration.jl` pins for the member
   jobs — so the actual code is self-consistent and only the prose was stale.

**Next.** Wait on `6967486`. At the next check: if it is still `Q`/`R`, do
nothing; if it has finished, read `clima_calibration.o6967486` and the member
logs before believing any parameter values, per the failure signatures.

## 2026-07-31 — iteration 2: still queued; diagnosed why

Woken by the fallback heartbeat, not the log-file waiter, so stage 1 had not
finished. It had not started either: `6967486` still `Q` at 1 h 08 m eligible.

**Why it is waiting.** `qstat -f` gives the reason outright — `Not Running: Job
is requesting an exclusive node and node is in use`, with `place =
scatter:exclhost`. The `main` queue allocates whole nodes, so the orchestrator,
which uses about one core to submit member jobs and poll, is queued behind a
free 128-core node. Derecho is busy: 1447 queued against 338 running.

**Chose not to escalate.** The server's `job_sort_formula` includes
`300*(eligible_time/86400)`, so a waiting job climbs; killing and resubmitting
would forfeit the accrued eligible time and restart that climb. Recorded an
explicit rule in STATE.md instead — if still `Q` at ~3 h eligible, move to the
shared `develop` queue (`cpudev` defaults to `pack:shared` and had 88 running
against 1 queued, so it starts near-immediately), after first confirming its
walltime cap. No `resources_max.walltime` is published on `develop`, `cpudev`,
or the server, so the widely-assumed 6 h is unverified and must be checked
rather than assumed before trading away a 12 h window.

**Worth recording about the member jobs**, since it changes where the real
bottleneck is: `run_calibration.jl` gives each ensemble member
`:gpus_per_task => 1` and `CLIMACOMMS_DEVICE = CUDA`, at 3 h walltime. So the
orchestrator is the only CPU job; the 21 members per iteration are GPU jobs.
Orchestrator queue time is a startup cost, not the throughput limit.

**Used the wait for a pre-flight check.** Verified all 10 prior parameter names
resolve in `toml/default_parameters.toml`, and that every prior mean equals the
current default — so the prior is centered on the defaults, as the config
claims, and `pmodel_α` is 0.028 on both sides. One trap for anyone repeating
this: the TOML mixes bare `[name]` and quoted `["name"]` headers, so grepping
`^\[name\]` reports most of them MISSING. My first pass did exactly that and
looked like a serious finding until I checked the file format.

Two `cpudev` jobs (`6967717`, `6967718`, `pc_base_p*`) are running under this
account but were not submitted by the loop. Left alone.

**Next.** Same as before, plus: apply the escalation rule if `6967486` is still
`Q` at ~3 h eligible.

## 2026-07-31 — iteration 3: escalation plan measured, then abandoned

`6967486` still `Q` at 2 h 11 m eligible, same exclusive-node reason, no
estimated start time offered. Rather than burn the escalation iteration on
investigation, settled the open question early — and the answer reversed the
plan.

**`develop` caps at 6 h, confirmed by experiment rather than assumption.**
Submitted a trivial probe job requesting 12 h; PBS rejected it outright: "your
declared wallclock time (43200 seconds) exceeds your maximum limit of 21600
seconds ... the queue limit is 21600 seconds". Rejection means no job was
created, so nothing needed cleaning up. 21600 s is 6 h exactly.

**The members are GPU jobs in `main`, which is the fact that matters.**
`ClimaCalibrate/src/pbs.jl:29` defaults the member queue to `"main"`, and
`run_calibration.jl` never overrides it — it only sets `gpus_per_task = 1` and
`time = 180`. So every one of the 21 members per iteration is an exclusive-node
GPU job with a 3 h walltime, competing in the same contended queue.

**Therefore: do not move to `develop`.** Five EKI iterations of 21 contended 3 h
GPU jobs will plausibly exceed 6 h of orchestrator wall time. Switching would
trade a bounded queue wait for repeated walltime kills, and each kill orphans
that iteration's in-flight member jobs — ClimaCalibrate resumes at
`last_completed_iteration + 1`, so the work is not lost, but the running members
are wasted. The 12 h window in `main` is worth waiting for. Also checked whether
shrinking the request would schedule sooner: it would not, because `main`
allocates an exclusive node for any request size.

This is the second time on this branch that a plausible-sounding assumption
about Derecho turned out to be checkable in about a minute. Worth continuing to
check rather than reason about.

**Next.** Keep waiting on `6967486`; the waiter fires when it ends. If it is
still queued around 8-10 h eligible, that stops being normal contention and
should be reported as a blocker rather than waited out further.

## 2026-07-31 — iteration 4: put a number on the queue decision

`6967486` still `Q` at ~3 h 15 m eligible. Rather than re-argue the queue choice
from first principles, went looking for evidence, and it was already on disk.

**A previous run of essentially this stage-1 config exists**, at
`/glade/derecho/scratch/arenchon/calibration_gpp_energy_rosetta_rebased_v1`,
from 2026-06-19. Per-iteration wall times, read off directory mtimes:

| iteration | finished | Δ |
|---|---|---|
| 001 | 06:07 | — |
| 002 | 12:55 | 6 h 48 m (cold start) |
| 003 | 15:02 | 2 h 07 m |
| 004 | 17:13 | 2 h 11 m |
| 005 | 18:27 | 1 h 14 m |
| 006 | 19:38 | 1 h 11 m |
| 007 | 20:50 | 1 h 12 m |
| 008 | 22:39 | 1 h 49 m |

A single member takes about 55 min (`iteration_005/member_001/model_log.txt`
spans 17:14→18:09). So five iterations should take roughly **6-10 h** of
orchestrator wall time, dominated by the cold-start first iteration.

That converts the earlier queue argument into a measurement: 6-10 h fits the
12 h `main` window and would certainly be killed by `develop`'s 6 h cap. The
decision to wait stands, now on evidence rather than reasoning.

It also gives a stall detector worth having: if the run is still on iteration 1
more than ~7 h after the orchestrator starts, read a member `model_log.txt`
rather than continuing to wait.

**One flag for later.** That reference run went to **9** iterations; this config
sets `n_iterations = 5`. That is not proof 5 is too few — the older run may just
have been left going — but it is a reason to check convergence carefully at
iteration 5 rather than assuming, and to be ready to extend per the ITERATIONS
section of the prompt. The same directory is the natural reference for sanity-
checking the calibrated parameter values when they arrive.

**Next.** Unchanged: wait for `6967486` to start; the waiter catches completion.

## 2026-07-31 — iteration 5: contention confirmed system-wide; stage 2 pre-flighted

`6967486` still `Q` at ~4 h 18 m eligible.

**The wait is the machine, not the job.** `qstat -Qf cpu` shows 2373 nodes
assigned out of Derecho's ~2488, with 1330 queued against 259 running, and
`backfill_depth = 5` on this queue. The cluster is essentially full. This is
worth stating explicitly because the prompt warns about misreading contention as
an outage — the opposite error is also available, and the queue statistics rule
out anything specific to this job.

**Used the wait to pre-flight stage 2**, since a defect there would otherwise
surface only after a 10-year GPU run had completed and been wasted. All checks
are static; no compute was used.

- All six IC fields have diagnostics. `lai` is in the base `CanopyModel` list;
  `a0a`, `pra`, `olf0`, `olvpd`, `olgsl` are appended by the `add_diagnostics!`
  method dispatched on `ZhouOptimalLAIModel`, which is what `PROGNOSTIC_LAI`
  selects.
- `snowy_land_pmodel.jl:166-176` requests exactly those six, unioned with the
  short set, when `PROGNOSTIC_LAI` is set — confirming the prompt's claim that no
  extra flag is needed.
- `PROGNOSTIC_LAI` is `haskey`-tested, so `PROGNOSTIC_LAI=""` does enable it;
  `RUN_YEARS` defaults to 2 and must be set to 10; `OUTPUT_ROOT` defaults to `.`
  and must be set, or multi-GB diagnostics land in the checkout.

**Pinned down the `precip_annual` unit trap with both sides quoted.** The prompt
flags it; now it is exact. `pra` is `units = "m yr^-1"`
(`define_diagnostics.jl:374`), the reader documents `precip_annual` as
`mol H2O m^-2 yr^-1` (`spatially_varying_parameters.jl:406`). The factor is
ρ_liq / M_H2O = 1000 / 0.018015 ≈ 5.5509e4 mol m^-3. Writing `pra` through raw
would understate precipitation by about 4.7 orders of magnitude — every land
cell would look hyper-arid, and the resulting IC would be quietly wrong rather
than obviously broken. Recorded in STATE.md so stage 2 does not have to
rediscover it.

**Next.** Unchanged. Blocker threshold still ~8-10 h eligible.

## 2026-07-31 — iteration 6: still queued; one hazard worth recording

`6967486` still `Q` at ~5 h 21 m eligible. No change, and no further pre-flight
work available — stage 3's remaining prerequisite, the model-LAI validity mask,
genuinely needs stage 2's long-run diagnostics as its reference simdir.

Revisited the `develop`-queue idea once more, because five hours of waiting is a
real cost and the earlier argument was about wasted time, which is arguable.
Found a stronger reason not to, and it is about correctness rather than waste.

**Member jobs are independent PBS jobs, so killing the orchestrator does not
kill them.** They keep running as orphans. A resubmitted orchestrator resumes at
`last_completed_iteration + 1`; an iteration cut short never completed, so it
resubmits that iteration's members while the orphans are still running. Two live
member jobs would then be writing the same `iteration_XXX/member_YYY` directory.
That is interleaved or corrupted output, not just duplicated GPU hours.

This applies to any walltime kill, including one in `main`. The difference is
frequency: tolerable once at the end of a 12 h window, not something to design
around every 6 h. It also means that if `6967486` ever does get walltime-killed
mid-iteration, the right move is to check for still-running member jobs and let
them finish or `qdel` them BEFORE resubmitting the orchestrator — not to
resubmit immediately.

**Next.** Unchanged. Blocker threshold still ~8-10 h eligible, i.e. ~19:00-21:00.

## 2026-07-31 — iteration 7: reported as blocked; chunked fallback verified

`6967486` still `Q` at 6 h 23 m. The reason this iteration is different from the
last two: the backlog is GROWING, so "keep waiting" stopped being obviously
convergent.

| | it. 5 (15:07) | it. 7 (17:12) |
|---|---|---|
| CPU nodes assigned | 2373 / ~2488 | 2406 |
| queued | 1330 | 1411 |
| running | 259 | 279 |

**Found and verified a third option: chunking.**
`ClimaCalibrate/src/backends.jl:277` computes
`first_iter = last_completed_iteration(output_dir) + 1` and runs
`first_iter:n_iterations`. So raising `n_iterations` one step at a time makes
each job run exactly one more iteration and then exit NORMALLY. Each chunk is
~1 h 15 m, well inside `develop`'s 6 h cap, and because no job is ever
walltime-killed, the orphaned-member hazard from iteration 6 never arises. That
hazard was the decisive argument against the naive "just use `develop`" idea;
chunking sidesteps it rather than accepting it.

**Why I reported instead of just switching.** The two approaches are mutually
exclusive: two orchestrators must not share a `CAL_OUTPUT_DIR`, so adopting
chunking means `qdel`-ing `6967486` and forfeiting 6 h 23 m of accrued
`eligible_time` priority. It also deviates from the committed 12 h `main`
wrapper, and the cold-start first chunk is the least certain to fit in 6 h. That
is a genuine trade-off rather than a clear improvement, and the user is away, so
the honest move is to state it, commit to a trigger, and keep working.

**Trigger set:** implement chunking if `6967486` is still queued at ~9-10 h
eligible (≈20:00-21:00 MDT). Posted to the PR as a blocked-status comment so the
situation is visible rather than buried in another quiet tick.

**Next.** Wait; act on the trigger if it fires.

## 2026-07-31 — iteration 8: switched to chunked `develop`; stage 1 IS RUNNING

Pulled the trigger early — at 7 h 26 m rather than the announced 9-10 h — because
two measurements taken this iteration were exactly the evidence the trigger was
waiting for, not impatience.

**1. A `develop` job starts instantly.** Submitted a probe carrying
`#PBS -q main` and `#PBS -l walltime=12:00:00`, overridden on the command line
with `-q develop -l walltime=05:30:00`. It was accepted (proving the walltime
override took effect, since 12 h would have been rejected), landed in `cpudev`,
ran, and finished within seconds. So the queue and walltime can be changed
WITHOUT editing the committed PBS wrapper.

**2. The GPU queue is nearly empty.** 4 queued, 26 running, against a CPU queue
with 1461 queued. The members are GPU jobs; the CPU orchestrator was the entire
bottleneck. That also collapses the cold-start risk that made me hesitate: with
the GPU queue clear, iteration 1 should be precompile (~30 min) + observations
(~3 min) + one member wave (~55 min) ≈ 2-2.5 h, comfortably inside 6 h.

Waiting two more hours purely to honour a self-imposed number would have been
worse work than acting on new information.

**What was done.**
- Added an `N_ITERATIONS` env override to the stage-1 config
  (`n_iterations = parse(Int, get(ENV, "N_ITERATIONS", "5"))`, default
  unchanged) and passed it through the PBS wrapper. This is the enabler for
  chunking: each job advances the run by one iteration and exits NORMALLY, so no
  walltime kill, so no orphaned members.
- `qdel 6967486` — 7 h 26 m queued, never ran, wrote nothing.
- Submitted chunk 1: `6970028`, `N_ITERATIONS=1`, `develop`, 5 h 30 m,
  `ncpus=4:mem=32gb`. **Started within ~2 minutes** on `dec2436`. Confirmed from
  its log: `n_iterations cap: 1`.

**Three environment facts learned, all recorded in STATE.md.**
- `qstat <jobid>` transiently reports "Unknown Job Id" while a job routes from
  `develop` to `cpudev`. Briefly alarming; not a failure.
- The PBS output file is written LIVE and appears seconds after the job starts,
  not at the end. My iteration-1 waiter design — "wait for the file to exist" —
  would therefore have fired immediately and meant nothing. Replaced with a
  waiter that greps the log for `Orchestrator finished` plus failure markers.
  Worth flagging that the earlier design was only ever correct by accident: it
  was never exercised, because that job never ran.
- `pkill -f <pattern>` self-matches, because the pattern appears in pkill's own
  command line; it kills its own shell first (exit 143). Use a bracket class
  (`o696748[6]`) to avoid it. The stale waiter from iteration 1 survives but is
  inert — it watches for an output file that a cancelled job never produced — so
  it was left alone rather than spending more turns on it.

**Next.** Watch `6970028`. When it exits cleanly, verify iteration 1's output
(members present, no all-NaN/all-zero GPP, loss finite) BEFORE submitting chunk 2
with `N_ITERATIONS=2`. The point of chunking is that each step is checkable.

## 2026-07-31 — iteration 12: first real failure — member 3 segfaulted

Wave 1 finished. **5 of 6 members completed; `member_003` crashed.**

**How it was spotted.** All of `run_1_1`…`run_1_6` had left the queue, but only
five `checkpoint.txt` files read `completed`. Member 3's still read `started`.
Job gone + checkpoint not `completed` is the signature of a crashed member, and
it is cheaper to check than parsing logs.

**What happened.** The last line of `member_003/model_log.txt`:

    deg0069.hsn.de.hpc.ucar.edu: rank 0 died from signal 11 and dumped core

SIGSEGV at 90 % complete, after 47 minutes of entirely healthy running —
`estimated_sypd` ~83, `percent_complete` climbing normally, `wall_time_per_step`
steady at ~29.7 ms right up to the crash. **This is transient or hardware, not
the parameters.** A bad parameter set gives NaN or a solver failure; it does not
segfault after 47 good minutes at full throughput.

**It is handled, and I verified that rather than assuming it.** Two layers:

1. `observation_map.jl:65-70` wraps `process_member_data!` in try/catch and
   calls `fill_g_ens_col!(g_ens_builder, m, NaN)` on failure. So the orchestrator
   degrades instead of dying.
2. EKP's `TransformUnscented` defaults to
   `failure_handler_method = SampleSuccGauss()`
   (`EnsembleKalmanProcess.jl:123-129`). `run_calibration.jl` overrides only
   `scheduler`, so that default is in force and the update is conditioned on the
   successful particles.

One failure in 21 is ~4.8 %, so iteration 1 should produce a valid update.

**What I will watch.** Losing a sigma point degrades the unscented covariance
estimate. EKP will not error on several failures, so "it completed" is not
evidence the update is sound. If more members segfault — particularly on
different nodes — this stops being transient and becomes a defect to report
rather than absorb.

**Timing improved.** GPU concurrency rose from 6 to 15; all remaining members
(7-21) are now running and should finish ~22:15, comfortably inside the 23:47
walltime. The earlier "20-40 min of margin" worry has receded.

## 2026-07-31 — iteration 13: chunk 1 done; I was wrong about member 3

**Chunk 1 finished cleanly at 22:14**, 3 h 57 m after starting, well inside its
5 h 30 m walltime. It wrote `iteration_002/eki_file.jld2` (336 MB) and staged the
next member dirs, so the chunking scheme is proven end-to-end. Chunk 2
(`6970860`, `N_ITERATIONS=2`) went in at 22:16 and started immediately.

**Correction: member 3's data is valid, and iteration 1 used all 21 members.**
Last iteration I reported that `member_003`'s segfault would degrade it to a NaN
column, leaving the update on 20/21 sigma points. That was wrong, and it matters
because it changes how sound iteration 1's update is.

What actually happened: the simulation ran to completion and wrote all its
diagnostics at 20:35 — its own predicted finish time — and then segfaulted
during teardown, before `checkpoint.txt` could be flipped to `completed`. So
ClimaCalibrate flagged `Failed ensemble members: [3]` off the checkpoint, while
`observation_map` read the file tree without complaint.

Verified rather than assumed, three ways:
- `member_003` has the same 10 `.nc` files as a healthy member;
- the same 36 time steps, with an identical time coordinate (first=0,
  last=92102400) to `member_001` and `member_005`;
- the orchestrator log contains ZERO "Error processing member" lines, and the G
  ensemble reports `size=(647824, 21)` with **0 NaN columns and finite fraction
  1.0**.

**The general lesson: `checkpoint.txt` is a process-exit signal, not a science
signal.** It says the process ended cleanly, not that the output is complete —
and the two can disagree in both directions. A member that crashed EARLY would
leave short or missing files, and the NaN path exists for exactly that; it
simply did not trigger here. Check the diagnostics before believing a flagged
member lost data.

**Iteration 1 results.** Misfit `[0.157]` — a single value, so nothing can be
said about convergence yet. Parameter means moved substantially from the prior,
and one deserves attention:

**`canopy_z_0b_coeff` went 0.04422 → 0.09765 against an upper bound of 0.1** —
97.7 % of the way there after ONE iteration. Not yet pinned, but heading for it.
Its prior σ is 0.01 on a mean of 0.0444, so the bound sits only 5.6σ out; the
data appear to want a larger roughness ratio than the prior permits. Per the
rules this is a finding to report, not something to quietly widen, so it goes in
the PR comment now and will be re-checked at iteration 2. If it reaches ~0.0999
it is effectively pinned and the prior is mis-specified.

Several other parameters also moved multiple σ in a single step (`pmodel_β_c4`
5.12 → 18.0, `canopy_z_0m_coeff` 0.352 → 0.158, `pmodel_β_c3` halved). That is
possible with `TransformUnscented` + `DataMisfitController`, but iterations 2-5
need watching for oscillation rather than convergence.

**Next.** Chunk 2 runs iteration 2. On completion: re-check
`canopy_z_0b_coeff` against its bound, confirm the misfit fell, and submit chunk
3.

## 2026-07-31 — iteration 14: warm start confirmed, full GPU concurrency

Chunk 2 is an hour in and both open questions resolved in the pipeline's favour.

**Warm start holds — the assumption chunking depends on.** Chunk 2's setup took
27 min (job start 22:16 → members submitted ~22:43) against chunk 1's 68 min
(18:17 → 19:25); its log is 255 lines against 1355. The shared `~/.julia` depot
meant `Pkg.update()` was a no-op and the caches were warm. So each chunk costs
about 27 min of overhead rather than an hour, and there is no reason to pack
more iterations into a chunk to amortise it. Had this gone the other way,
chunking would have cost ~4 extra hours across stage 1 and I would have had to
reconsider.

**All 21 members are running concurrently**, versus chunk 1's 6-then-15 ramp.
The GPU contention that shaped the earlier timing estimates has cleared.

Together these cut the per-chunk time from 3 h 57 m to a projected ~2 h 15 m.
Revised estimate: chunks 3-5 at ~2 h 15 m each ≈ 6 h 45 m, so stage 1 should
finish around 07:00-07:30 on 2026-08-01 if concurrency holds. Chunk 2's walltime
does not expire until 03:46, so there is no pressure.

**Next.** On chunk 2 completion: check whether the misfit fell from 0.157,
re-check `canopy_z_0b_coeff` against its 0.1 bound (0.09765 after iteration 1),
count crashed members, and submit chunk 3.
