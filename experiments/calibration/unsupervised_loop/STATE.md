# Pipeline state

Single source of truth for where the three-stage calibration pipeline is.
The loop reads this first every iteration and commits it every iteration.
See `loop_prompt.md` for what each stage means and how it is gated.

Status values: `not_started` | `running` | `done` | `failed`

All heavy output lives under `/glade/derecho/scratch/arenchon/claude/`.

---

## Stage 1 — GPP + energy fluxes

- **status:** **done** (2026-08-01 08:47, terminated by EKP at iteration 8)
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
- **MEMBER 3 CRASHED BUT ITS DATA IS VALID — corrected finding (2026-07-31 22:15).**
  `member_003` died with signal 11 (SIGSEGV) on node `deg0069`, and its
  `checkpoint.txt` is stuck at `started`, so ClimaCalibrate logged
  `Failed ensemble members: [3]`. My first reading — that it would degrade to a
  NaN column and iteration 1 would run on 20/21 sigma points — WAS WRONG.
  What actually happened: the simulation ran to completion and wrote all its
  diagnostics at 20:35 (its own predicted finish time), then segfaulted during
  teardown, before `checkpoint.txt` could be updated. Verified directly:
  `member_003` has the same 10 `.nc` files as a healthy member, the same 36
  time steps, and an identical time coordinate (first=0, last=92102400).
  Consequently `observation_map` read it successfully — the log contains ZERO
  "Error processing member" lines — and the G ensemble has **0 NaN columns,
  finite fraction 1.0 across all 21 members**. Iteration 1 used the full
  ensemble; the unscented covariance is NOT degraded.
  **Lesson for future iterations: `checkpoint.txt` is not a reliable success
  signal.** It reports whether the process exited cleanly, not whether the
  science completed. Check the diagnostics (file count, time length) before
  concluding a flagged member actually lost data. Conversely, a member that
  crashes EARLY would leave short/missing files, and that is what the NaN path
  is for — it just did not trigger here.
- **ITERATION 1 COMPLETE (2026-07-31 22:14).** Chunk 1 (`6970028`) exited
  cleanly in 3 h 57 m, well inside its 5 h 30 m walltime, and wrote
  `iteration_002/eki_file.jld2` (336 MB) plus staged member dirs. Chunking works
  end-to-end. **Chunk 2 (`6970860`, `N_ITERATIONS=2`) submitted 22:16, running.**
  G ensemble: `size=(647824, 21)`, **0 members with NaN, finite fraction 1.0**.
  Misfit history so far: `[0.157]` — one value, so convergence cannot be judged
  yet.
  Constrained parameter means, prior (= TOML default) → after iteration 1:

  | parameter | prior mean | after it.1 | bounds | note |
  |---|---|---|---|---|
  | pmodel_cstar | 0.4295 | 0.4617 | 0.2–0.7 | small |
  | pmodel_β_c3 | 87.05 | 41.44 | 10–300 | halved |
  | pmodel_β_c4 | 5.121 | 18.00 | 1–100 | ~3.5×, large |
  | pmodel_α | 0.0268 | 0.02869 | 0.001–0.15 | small; τ ≈ 35 d, convention OK |
  | moisture_stress_c | 0.5947 | 0.7168 | 0.05–1.0 | moderate |
  | leaf_Cd | 0.06922 | 0.04393 | 0–Inf | moderate |
  | canopy_z_0m_coeff | 0.3520 | 0.1578 | 0–0.5 | large, ~3.9σ |
  | canopy_z_0b_coeff | 0.04422 | **0.09765** | 0–**0.1** | ⚠ **97.7 % of upper bound** |
  | canopy_d_coeff | 0.05446 | 0.05019 | 0–1.0 | small |
  | canopy_K_lw | 0.9168 | 1.229 | 0–2.0 | moderate |

  **⚠ `canopy_z_0b_coeff` is at 97.7 % of its upper bound (0.09765 vs 0.1) after
  one iteration.** Not yet pinned, but heading there. Per `loop_prompt.md`, a
  parameter pinned at a bound means the prior is mis-specified, and that is a
  finding to REPORT, not something to quietly widen. Watch at iteration 2: if it
  sits at ~0.0999 it is effectively pinned and must go in the PR comment. Its
  prior σ is only 0.01 against a mean of 0.0444, so the bound is 5.6σ out — a
  tight prior, and the data appear to want a larger roughness ratio than it
  allows.
  Several other parameters moved multiple σ in one step, which is possible with
  `TransformUnscented` + `DataMisfitController` but worth watching for
  oscillation rather than convergence across iterations 2-5.
- **CHUNK 2 / ITERATION 2 (2026-07-31 23:21, 1 h 01 m elapsed): warm start
  CONFIRMED, full concurrency.** The assumption chunking rests on holds:
    - Setup took **27 min** (job start 22:16 → members submitted ~22:43) versus
      **68 min** for chunk 1 (18:17 → 19:25). The shared `~/.julia` depot meant
      `Pkg.update()` was a no-op and the precompile caches were warm. Chunk 2's
      log is 255 lines against chunk 1's 1355. So chunking costs ~27 min of
      overhead per chunk, not ~1 h — it is cheap, and there is no reason to pack
      more iterations per chunk.
    - **All 21 members are running CONCURRENTLY** (`run_2_1`…`run_2_21`), versus
      chunk 1's 6-then-15 ramp. GPU contention has cleared.
  Expected: members finish ~23:53-00:07, iteration 2 completes ~00:15-00:30,
  i.e. ~2 h 15 m for the chunk against chunk 1's 3 h 57 m, and far inside the
  5 h 30 m walltime (which expires 03:46). REVISED ESTIMATE for the remaining
  work: chunks 3-5 at ~2 h 15 m each ≈ 6 h 45 m, so stage 1 should finish around
  07:00-07:30 on 2026-08-01 if concurrency holds.
- **ITERATION 2 COMPLETE (2026-08-01 00:08).** Chunk 2 ran in **1 h 52 m** — the
  warm start and full 21-way concurrency delivered. **All 21 members completed,
  zero failures, zero segfaults**, which supports reading member 3's crash as
  transient/hardware: no recurrence. Chunk 3 (`6971281`, `N_ITERATIONS=3`)
  submitted 00:10 and running.

  | parameter | prior | it.1 | it.2 | bounds |
  |---|---|---|---|---|
  | pmodel_cstar | 0.4295 | 0.4617 | 0.4066 | 0.2–0.7 |
  | pmodel_β_c3 | 87.05 | 41.44 | 32.42 | 10–300 |
  | pmodel_β_c4 | 5.121 | 18.00 | 9.125 | 1–100 |
  | pmodel_α | 0.0268 | 0.02869 | 0.03244 | 0.001–0.15 |
  | moisture_stress_c | 0.5947 | 0.7168 | 0.6197 | 0.05–1.0 |
  | leaf_Cd | 0.06922 | 0.04393 | 0.06503 | 0–Inf |
  | canopy_z_0m_coeff | 0.3520 | 0.1578 | 0.1834 | 0–0.5 |
  | canopy_z_0b_coeff | 0.04422 | 0.09765 | **0.08586** | 0–0.1 |
  | canopy_d_coeff | 0.05446 | 0.05019 | 0.03263 | 0–1.0 |
  | canopy_K_lw | 0.9168 | 1.229 | 1.215 | 0–2.0 |

  **`canopy_z_0b_coeff` backed OFF the bound** (0.09765 → 0.08586). It is not
  pinned; the iteration-1 alarm has receded. Keep watching, but do not treat it
  as a mis-specified prior on current evidence.

- **⚠ THE MISFIT SERIES IS NOT A CONVERGENCE DIAGNOSTIC HERE — important.**
  Misfit went 0.157 → **0.173**, i.e. UP. That is NOT evidence of divergence,
  because `minibatch_size = 2`: `run_calibration.jl:117` builds the
  `ObservationSeries` with `minibatcher_over_samples(16, 2)`, so EACH ITERATION
  SCORES A DIFFERENT PAIR OF YEARS. The two numbers are computed on different
  data and are not comparable. `loop_prompt.md`'s "the loss decreased across
  iterations" check is therefore misleading for this config, and must not be
  applied naively.
  **Use instead:** parameter means ceasing to move materially, which is the
  other criterion the prompt gives. For a like-for-like loss comparison one would
  have to score a fixed batch, which the pipeline does not currently do.
- **⚠ 5 ITERATIONS DOES NOT COMPLETE ONE EPOCH.** 16 samples at minibatch 2 is
  **8 iterations per epoch**, so `n_iterations = 5` sees only 10 of the 16
  sample years — some data is never used, and no parameter has been scored
  against the full record. This is a principled reason the reference run's **9**
  iterations may have been the right choice (one full epoch plus one), rather
  than an arbitrary preference. RECOMMENDATION to weigh at iteration 5: extend
  to at least 8, ideally 9, per the ITERATIONS section of `loop_prompt.md`
  (raise `n_iterations`, commit `[skip ci]`, resubmit with the same
  `CAL_OUTPUT_DIR`). Report the reasoning in the PR rather than extending
  silently.
- **ITERATION 3 COMPLETE (2026-08-01 01:48).** Chunk 3 ran in 1 h 38 m;
  **21/21 members completed, zero segfaults**. Chunk 4 (`6971902`,
  `N_ITERATIONS=4`) submitted 01:50, running. All G blocks remain
  `size=(647824, 21)` with 0 NaN and finite fraction 1.0.
  Misfit: `[0.157, 0.173, 0.154]` — fluctuating, as expected given each
  iteration scores a different minibatch. NOT a convergence signal; see above.

  | parameter | prior | it.1 | it.2 | it.3 | trend |
  |---|---|---|---|---|---|
  | pmodel_cstar | 0.4295 | 0.4617 | 0.4066 | 0.3945 | drifting down |
  | pmodel_β_c3 | 87.05 | 41.44 | 32.42 | 32.28 | **settled** |
  | pmodel_β_c4 | 5.121 | 18.00 | 9.125 | 8.558 | settling |
  | pmodel_α | 0.0268 | 0.02869 | 0.03244 | 0.02896 | oscillating |
  | moisture_stress_c | 0.5947 | 0.7168 | 0.6197 | 0.6090 | settling |
  | leaf_Cd | 0.06922 | 0.04393 | 0.06503 | 0.06658 | settling |
  | canopy_z_0m_coeff | 0.3520 | 0.1578 | 0.1834 | 0.1382 | still moving |
  | canopy_z_0b_coeff | 0.04422 | 0.09765 | 0.08586 | 0.09301 | oscillating high |
  | canopy_d_coeff | 0.05446 | 0.05019 | 0.03263 | 0.02785 | drifting down |
  | canopy_K_lw | 0.9168 | 1.229 | 1.215 | 1.211 | **settled** |

  Read: `pmodel_β_c3` and `canopy_K_lw` have essentially stopped moving;
  `pmodel_β_c4`, `moisture_stress_c` and `leaf_Cd` are settling; but
  `canopy_z_0m_coeff`, `canopy_d_coeff` and `pmodel_cstar` are still drifting
  materially at iteration 3 of 5. **Not converged.** This is direct evidence for
  the epoch argument — stopping at 5 would freeze three parameters mid-drift.
  `canopy_z_0b_coeff` keeps orbiting 0.086-0.098 against its 0.1 bound: never
  pinned, but persistently pressed against it, so it stays on the watch list
  even though the iteration-1 alarm was withdrawn.
- **NOTE on extending:** no config edit is needed to go past 5. `N_ITERATIONS` is
  an env override, so chunks 6-9 are just `N_ITERATIONS=6…9`. The config default
  should still be updated to 9 if that becomes the conclusion, so the file
  reflects reality rather than relying on how the loop happened to invoke it.
- **ITERATION 4 COMPLETE (2026-08-01 03:22).** Chunk 4 ran in 1 h 32 m;
  **21/21 members, zero segfaults** (as in iterations 2 and 3 — member 3's crash
  remains a one-off). Chunk 5 (`6972224`, `N_ITERATIONS=5`) submitted 03:25.
  Misfit: `[0.157, 0.173, 0.154, 0.150]` — still not comparable across
  iterations (minibatching), so still not a convergence signal.

  Parameter means, and the change over the LAST step (it.4 → it.5 column):

  | parameter | it.1 | it.2 | it.3 | it.4 | it.5 | last step |
  |---|---|---|---|---|---|---|
  | pmodel_cstar | 0.4295 | 0.4617 | 0.4066 | 0.3945 | 0.3823 | −3.1 % |
  | pmodel_β_c3 | 87.05 | 41.44 | 32.42 | 32.28 | 29.49 | **−8.6 %** |
  | pmodel_β_c4 | 5.121 | 18.00 | 9.125 | 8.558 | 8.389 | −2.0 % |
  | pmodel_α | 0.0268 | 0.02869 | 0.03244 | 0.02896 | 0.02782 | −3.9 % |
  | moisture_stress_c | 0.5947 | 0.7168 | 0.6197 | 0.6090 | 0.5936 | −2.5 % |
  | leaf_Cd | 0.06922 | 0.04393 | 0.06503 | 0.06658 | 0.07191 | **+8.0 %** |
  | canopy_z_0m_coeff | 0.3520 | 0.1578 | 0.1834 | 0.1382 | 0.1105 | **−20 %** |
  | canopy_z_0b_coeff | 0.04422 | 0.09765 | 0.08586 | 0.09301 | 0.08608 | −7.5 % |
  | canopy_d_coeff | 0.05446 | 0.05019 | 0.03263 | 0.02785 | 0.04199 | **+51 %** |
  | canopy_K_lw | 0.9168 | 1.229 | 1.215 | 1.211 | 1.199 | −1.0 % |

  **NOT CONVERGED, decisively.** Only `canopy_K_lw` (−1 %) and `pmodel_β_c4`
  (−2 %) have settled. `canopy_z_0m_coeff` is still moving 20 % per step,
  `canopy_d_coeff` jumped 51 %, and `pmodel_β_c3` moved −8.6 % after appearing
  flat at iteration 3 — a reminder that two similar consecutive values are not
  convergence. **Stopping at 5 would freeze at least four parameters mid-flight.**
  Combined with the epoch argument (5 iterations sees 10 of 16 sample years),
  extending is the clear call.

- **⚠ NEW: `canopy_z_0m_coeff` is heading for its LOWER bound.** Monotonic across
  every iteration: 0.352 → 0.158 → 0.183 → 0.138 → 0.1105, against bounds
  0–0.5 and a prior of (0.349, σ=0.05). It is now ~4.8σ below the prior mean and
  still falling. If it keeps going it will press against 0, which is the mirror
  image of the `canopy_z_0b_coeff` concern and would likewise mean a
  mis-specified prior. Watch it at iterations 5+.
- **⚠ PHYSICAL PLAUSIBILITY, worth a domain expert's eye.** These two coefficients
  multiply canopy height to give the momentum and scalar roughness lengths, and
  they are moving in OPPOSITE directions: `z_0b/z_0m` has gone from
  0.04422/0.3520 ≈ **0.13** at the prior to 0.08608/0.1105 ≈ **0.78** at
  iteration 5. Scalar roughness is conventionally much smaller than momentum
  roughness (ratio ~0.1), so a ratio approaching unity is atypical. The data may
  genuinely want this — it is compensating for something in the flux fit — but it
  should be reported rather than absorbed, and NOT quietly corrected by changing
  priors mid-run.
- **ITERATION 5 COMPLETE (2026-08-01 04:46), AND THE RUN WAS EXTENDED TO 9.**
  Chunk 5 ran in 1 h 21 m; 21/21 members, zero segfaults. `n_iterations` raised
  5 → 9 in the config and committed; chunk 6 (`6972345`, `N_ITERATIONS=6`)
  submitted 04:50 and running.
  Misfit: `[0.157, 0.173, 0.154, 0.150, 0.151]` — now fluctuating in a narrow
  band, but still not a convergence signal (different minibatch each iteration).

  | parameter | it.4 | it.5 | it.6 | last step |
  |---|---|---|---|---|
  | pmodel_cstar | 0.3945 | 0.3823 | 0.3870 | +1.2 % |
  | pmodel_β_c3 | 32.28 | 29.49 | 30.37 | +3.0 % |
  | pmodel_β_c4 | 8.558 | 8.389 | 9.315 | **+11 %** |
  | pmodel_α | 0.02896 | 0.02782 | 0.02630 | −5.5 % |
  | moisture_stress_c | 0.6090 | 0.5936 | 0.5875 | −1.0 % |
  | leaf_Cd | 0.06658 | 0.07191 | 0.06908 | −3.9 % |
  | canopy_z_0m_coeff | 0.1382 | 0.1105 | 0.1009 | −8.7 % |
  | canopy_z_0b_coeff | 0.09301 | 0.08608 | 0.08427 | −2.1 % |
  | canopy_d_coeff | 0.02785 | 0.04199 | 0.03645 | **−13 %** |
  | canopy_K_lw | 1.211 | 1.199 | 1.195 | −0.3 % |

  **Converging, not yet converged — extending was the right call.** Step sizes
  are broadly smaller than the previous iteration (`canopy_z_0m_coeff` −20 % →
  −8.7 %), and `canopy_K_lw`, `moisture_stress_c` and `pmodel_cstar` are inside
  a few percent. But `pmodel_β_c4` (+11 %) and `canopy_d_coeff` (−13 %) are still
  swinging. Had the run stopped at 5, those two would have been frozen at
  arbitrary points on their trajectories.
- **`canopy_z_0m_coeff` lower-bound concern is RECEDING.** Sequence:
  0.352 → 0.158 → 0.183 → 0.138 → 0.1105 → 0.1009. Still monotonic, but the
  decline is decelerating (−20 %, then −8.7 %), consistent with asymptoting near
  0.09-0.10 rather than running into the bound at 0. Keep watching to iteration
  9, but do not report it as pinned on current evidence.
- **The roughness-ratio oddity persists and is the main thing for a domain
  expert.** `z_0b/z_0m` = 0.08427/0.1009 ≈ **0.84**, up from ≈0.13 at the prior.
  The two coefficients have converged toward each other rather than one being an
  order of magnitude below the other. This is stable across iterations 2-5, so it
  is a real feature of the fit, not a transient. Report it with the final
  numbers; do not tune priors to remove it.
- **ITERATION 6 COMPLETE (2026-08-01 06:06).** Chunk 6 ran in 1 h 16 m; 21/21
  members, zero segfaults. Chunk 7 (`6972557`, `N_ITERATIONS=7`) submitted 06:10.
  Misfit: `[0.157, 0.173, 0.154, 0.150, 0.151, 0.154]`.

  | parameter | it.5 | it.6 | it.7 | last step |
  |---|---|---|---|---|
  | pmodel_cstar | 0.3823 | 0.3870 | 0.3728 | −3.7 % |
  | pmodel_β_c3 | 29.49 | 30.37 | 25.90 | **−14.7 %** |
  | pmodel_β_c4 | 8.389 | 9.315 | 8.352 | −10.3 % |
  | pmodel_α | 0.02782 | 0.02630 | 0.02516 | −4.3 % |
  | moisture_stress_c | 0.5936 | 0.5875 | 0.5886 | +0.2 % |
  | leaf_Cd | 0.07191 | 0.06908 | 0.06781 | −1.8 % |
  | canopy_z_0m_coeff | 0.1105 | 0.1009 | 0.08303 | **−17.7 %** |
  | canopy_z_0b_coeff | 0.08608 | 0.08427 | 0.09546 | **+13.3 %** |
  | canopy_d_coeff | 0.04199 | 0.03645 | 0.04532 | **+24.3 %** |
  | canopy_K_lw | 1.199 | 1.195 | 1.175 | −1.7 % |

- **CORRECTION to the iteration-5 read.** I recorded that
  `canopy_z_0m_coeff`'s decline was decelerating and "consistent with asymptoting
  near 0.09-0.10". **That was premature.** It fell another 17.7 % to 0.08303, so
  the deceleration was a one-iteration pause, not an asymptote. Full sequence:
  0.352 → 0.158 → 0.183 → 0.138 → 0.1105 → 0.1009 → 0.08303 — monotonic decline
  throughout, at a fluctuating rate. It is still heading toward its lower bound
  of 0 and the concern is live again. Lesson, the same one iteration 3 taught
  about `pmodel_β_c3`: a single slower step is not convergence.
- **⚠⚠ THE ROUGHNESS RATIO HAS NOW INVERTED — this is the headline finding.**
  `canopy_z_0b_coeff` (0.09546) has overtaken `canopy_z_0m_coeff` (0.08303), so
  `z_0b/z_0m` ≈ **1.15**. Scalar roughness now EXCEEDS momentum roughness.
  Standard surface-layer theory has z_0b roughly an order of magnitude SMALLER
  than z_0m (ratio ~0.1); the prior encoded ~0.13. A ratio above 1 is not merely
  atypical, it inverts the expected physics. Trajectory of the ratio: 0.13
  (prior) → 0.62 → 0.47 → 0.67 → 0.78 → 0.84 → **1.15**, i.e. monotonic and
  still climbing. This is a stable, developing feature of the fit across seven
  iterations, not noise.
  **Do not tune priors to suppress it.** Report it with the final numbers and
  flag it for a domain expert: either the flux data genuinely require it (the
  pair may be compensating for a missing process in the canopy turbulence
  scheme), or the parameterisation is being pushed outside its valid regime — and
  which of those it is cannot be settled from the EKI output alone.
- **convergence status: NOT converging monotonically.** Iteration 6's steps are
  LARGER than iteration 5's for four parameters (`β_c3` +3.0 → −14.7 %,
  `z_0m` −8.7 → −17.7 %, `d_coeff` −13 → +24.3 %, `z_0b` −2.1 → +13.3 %). With a
  different minibatch each iteration the means chase the current pair of years,
  so per-iteration steps are noisy by construction. This is a further argument
  for running the full 9 (one epoch plus one) and, when judging the final answer,
  looking at the trend over the last few iterations rather than the last step.
- **ITERATION 7 COMPLETE (2026-08-01 07:24).** Chunk 7 ran in 1 h 14 m; 21/21
  members, zero segfaults. Chunk 8 (`6972840`, `N_ITERATIONS=8`) submitted 07:28.
  Misfit: `[…, 0.1539, 0.1535]`.

  | parameter | it.6 | it.7 | it.8 | last step |
  |---|---|---|---|---|
  | pmodel_cstar | 0.3870 | 0.3728 | 0.3655 | −2.0 % |
  | pmodel_β_c3 | 30.37 | 25.90 | 24.94 | −3.7 % |
  | pmodel_β_c4 | 9.315 | 8.352 | 8.891 | +6.5 % |
  | pmodel_α | 0.02630 | 0.02516 | 0.02514 | **−0.08 %** |
  | moisture_stress_c | 0.5875 | 0.5886 | 0.5907 | +0.4 % |
  | leaf_Cd | 0.06908 | 0.06781 | 0.07003 | +3.3 % |
  | canopy_z_0m_coeff | 0.1009 | 0.08303 | 0.08225 | −0.9 % |
  | canopy_z_0b_coeff | 0.08427 | 0.09546 | 0.07297 | **−23.6 %** |
  | canopy_d_coeff | 0.03645 | 0.04532 | 0.04160 | −8.2 % |
  | canopy_K_lw | 1.195 | 1.175 | 1.156 | −1.6 % |

- **CORRECTION: I overstated the roughness-ratio trend last iteration.** I
  described it as "monotonic and still climbing" and called the inversion the
  headline finding. Two things were wrong:
    1. It was never monotonic. The sequence I quoted contains a decrease I
       glossed over (0.62 → 0.47 between iterations 2 and 3).
    2. The inversion did not persist. `canopy_z_0b_coeff` fell 23.6 % this
       iteration, so `z_0b/z_0m` = 0.07297/0.08225 ≈ **0.89**, back below 1.
  Accurate sequence: **0.126 (prior) → 0.62 → 0.47 → 0.67 → 0.78 → 0.84 → 1.15
  → 0.89.** Rising strongly overall but oscillating, and the >1 excursion was a
  single iteration.
  **The substantive finding survives, restated correctly:** the ratio has settled
  around **0.8-0.9 over the last four iterations, versus ~0.13 in the prior and a
  conventional ~0.1** — roughly a factor of 7 higher than surface-layer theory
  expects, with scalar and momentum roughness nearly equal. That is still well
  worth a domain expert's attention. It is NOT a stable inversion, and I should
  not have called it one on a single data point.
- **`canopy_z_0m_coeff` has nearly stopped** (−0.9 %, 0.08303 → 0.08225) after
  −17.7 % the iteration before. Deliberately NOT calling this convergence — I
  made exactly that mistake at iteration 5 and it resumed falling. Judge it at
  iteration 9 on the trend, not on one small step.
- **`pmodel_α` is effectively converged** at 0.02514 (−0.08 %), i.e. τ ≈ 40 d.
  Still firmly in the post-#1817 `1/τ` convention.
- **STAGE 1 IS DONE. The EKI TERMINATED ITSELF AT ITERATION 8 (2026-08-01
  08:47).** Chunk 8's log ends with:

      Warning: Termination condition of scheduler `DataMisfitController` has been
      exceeded, returning `true` from `update_ensemble!` and preventing futher updates

  `run_calibration.jl:140` sets `DataMisfitController(terminate_at = 100)`, and
  `LearningRateSchedulers.jl:284` terminates once `sum_Δt >= T` — the inversion
  has consumed its allotted algorithmic time. **This is a principled stopping
  criterion, not a failure.** Iteration 8's members ran and were processed, but
  the resulting update was refused, which is why `iteration_009/eki_file.jld2`
  contains no new block.
  It stopped at **exactly one epoch (8 iterations at minibatch 2 over 16
  samples)**, which independently corroborates the epoch reasoning used to raise
  `n_iterations` from 5 to 9.
  **Chunk 9 (`6973790`) was cancelled** — with updates blocked it would have
  burned 21 GPU jobs to produce nothing.

- **CALIBRATED PARAMETERS, written to `toml/default_parameters.toml` at 3 s.f.
  and committed.** Final ensemble mean (u[8]); spread is the ensemble std.

  | parameter | old default | **calibrated** | spread | CV |
  |---|---|---|---|---|
  | pmodel_cstar | 0.43 | **0.366** | 0.00071 | 0.2 % |
  | pmodel_β_c3 | 91.1 | **24.9** | 0.20 | 0.8 % |
  | pmodel_β_c4 | 5.36 | **8.89** | 0.27 | 3.0 % |
  | pmodel_α | 0.028 | **0.0251** | 0.00011 | 0.4 % |
  | moisture_stress_c | 0.59 | **0.591** | 0.0015 | 0.3 % |
  | leaf_Cd | 0.0726 | **0.07** | 0.00053 | 0.8 % |
  | canopy_z_0m_coeff | 0.349 | **0.0823** | 0.0028 | 3.4 % |
  | canopy_z_0b_coeff | 0.0444 | **0.073** | 0.0046 | 6.3 % |
  | canopy_d_coeff | 0.0573 | **0.0416** | 0.0127 | **30 %** |
  | canopy_K_lw | 0.919 | **1.16** | 0.0040 | 0.4 % |

- **Sanity checks, all passed.**
  - G ensemble finite over ALL 8 iterations: `size=(647824, 21)`, 0 NaN columns,
    finite fraction 1.0. GPP is not zero or NaN anywhere.
  - **No parameter pinned at a prior bound.** Closest is `canopy_z_0b_coeff` at
    73 % of its 0–0.1 range. `canopy_d_coeff` sits at 4.2 % of its 0–1 range but
    its prior mean (0.0573) was already near the bottom of that range and it is
    3.3σ above zero, so it is poorly CONSTRAINED, not pinned.
  - Ensemble spreads are mostly under 1 % CV, i.e. the ensemble has genuinely
    collapsed. The exception is `canopy_d_coeff` at 30 % CV — the flux targets
    barely constrain it, which is worth knowing before anyone treats 0.0416 as a
    precise result.
  - The prompt's "loss decreased across iterations" check is NOT applicable:
    `minibatch_size = 2` means each iteration scores a different pair of years.
    EKP's own termination is the sound stopping signal here.

- **⚠ FINDING FOR A DOMAIN EXPERT: the roughness ratio.** Final
  `z_0b/z_0m` = 0.073/0.0823 ≈ **0.89**, against ≈0.13 in the prior and a
  conventional ≈0.1 — roughly 7× higher, i.e. scalar and momentum roughness end
  up nearly equal where surface-layer theory expects an order of magnitude
  between them. `canopy_z_0m_coeff` fell from 0.349 to 0.0823 (a factor of 4)
  while `canopy_z_0b_coeff` rose from 0.0444 to 0.073. Both stayed inside their
  bounds, and the ensemble is tight, so this is what the flux data actually
  prefer under this parameterisation. Either the pair is compensating for a
  missing process in the canopy turbulence scheme, or the scheme is being pushed
  outside its valid regime. Not resolvable from EKI output; reported, not tuned
  away.
- **committed to `toml/default_parameters.toml`:** **yes** — commit `ec4774837`.
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

- **status:** running — long run submitted 2026-08-01 09:35
- **PBS job:** `6973822.desched1`, queue `gpu` (`main`), **12 h walltime**,
  `select=1:ncpus=4:ngpus=1`, via the new
  `experiments/long_runs/long_run_gpu.pbs` wrapper:

      qsub -v LONG_RUN_SCRIPT=snowy_land_pmodel.jl,PROGNOSTIC_LAI=1,RUN_YEARS=10,OUTPUT_ROOT=/glade/derecho/scratch/arenchon/claude \
          experiments/long_runs/long_run_gpu.pbs

  **Why `main`/12 h and not the faster `develop` queue:** a long run has NO
  checkpoint to resume from, so a walltime kill loses everything and starts over.
  The estimate is ~4.5 h, which leaves too little margin against `develop`'s hard
  6 h cap. The calibration could safely use `develop` only because ClimaCalibrate
  resumes at `last_completed_iteration + 1`; that argument does not transfer here.
  `PROGNOSTIC_LAI=1` rather than an empty value: the driver tests
  `haskey(ENV, ...)`, and a non-empty value avoids any risk of PBS dropping an
  empty one.
- **log:** `clima_long_run.o6973822` in the repo root
- **expected output:** `/glade/derecho/scratch/arenchon/claude/global_runs/snowy_land_pmodel_opt_lai_longrun_gpu`
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
- **sizing the stage-2 run (2026-08-01 01:12, from prior runs on scratch).**
  `/glade/derecho/scratch/arenchon/claude/global_runs/` holds many earlier
  `PROGNOSTIC_LAI` long runs. Taking
  `snowy_land_pmodel_opt_lai_zhouparams_longrun_gpu` as representative: file
  mtimes span 02:46 → 04:18, i.e. **1 h 32 m end-to-end for a 3-year run**
  (36 monthly outputs), including the trailing plotting/ILAMB step. Scaling the
  simulation part to `RUN_YEARS=10` gives roughly **3.5-4 h**, so stage 2 fits
  inside `develop`'s 6 h cap as well as `main`. Prefer `develop` for the fast
  start, same as the calibration chunks.
  Those older runs carry `pra` and `a0a` but NOT `olf0`/`olvpd`/`olgsl`, which
  confirms those three diagnostics are new on this branch and were added for the
  IC rebuild — so no earlier run can substitute for stage 2's, even though the
  directory is full of superficially similar output.
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
