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
