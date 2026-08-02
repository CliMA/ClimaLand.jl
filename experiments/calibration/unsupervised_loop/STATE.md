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

- **status:** **done** (2026-08-01 13:20)
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
- **progress (2026-08-01 12:02, 3 h 06 m in): STILL PRECOMPILING, and it is
  genuine work, not a loop.** 847 `.ji` files written since the job started,
  across **707 DISTINCT PACKAGES** — repeats are negligible (max 6, for one small
  extension). So this is a real first-time build of the whole environment, not
  cache thrashing. Cause: the calibration's `Pkg.update()` bumped ClimaCore, a
  deep dependency, invalidating most of the tree, and `snowy_land_pmodel.jl`
  additionally loads the CairoMakie/GeoMakie visualization stack for its trailing
  ILAMB/plot step, which the calibration never touched. 707 packages at ~10 s
  each is ~2 h, matching what is observed. Memory 11.1 GB and climbing.
  **HOW TO TEST FOR A PRECOMPILE LOOP** (worth keeping, it is not obvious): count
  `.ji` files by their PARENT DIRECTORY, which is the package name —

      find ~/.julia/compiled -name '*.ji' -newermt '<job start>' -printf '%h\n' \
        | sed 's|.*/||' | sort | uniq -c | sort -rn | head

  Many files but few distinct directories means thrashing; many distinct
  directories means honest progress. Do NOT key on the `.ji` FILENAME — its
  suffix encodes the build configuration, not the package, so grouping by it
  shows a handful of "slugs" and looks exactly like a loop when nothing is wrong.
- **progress (2026-08-01 11:01, 2 h 06 m in): PRECOMPILING, not stuck.**
  `find ~/.julia/compiled -name '*.ji' -newermt '2026-08-01 08:55' | wc -l`
  gives **274** caches written since the job started, with writes continuing past
  10:53. This is the CairoMakie/GeoMakie visualization stack that
  `snowy_land_pmodel.jl` pulls in for its trailing ILAMB/plot step and that the
  calibration path never loaded, so it was not already warm. Memory is climbing
  (8.6 → 9.9 GB) and `cpupercent = 153`.
  **This vindicates asking for 12 h in `main`.** A 5 h 30 m `develop` job would
  have spent 2 h precompiling and had only ~3.5 h left for a ~4 h simulation —
  it would very likely have been killed, and a long run has no checkpoint.
  ⚠ **`find -newermt "-20 minutes"` GAVE A FALSE NEGATIVE** here — it reported
  0 recent caches while absolute-timestamp queries showed writes minutes earlier,
  which briefly made a healthy job look hung. Use absolute timestamps
  (`-newermt "2026-08-01 08:55"`) for this check, not relative ones.
- **progress (2026-08-01 09:57, 1 h 02 m in):** running on `deg0019`. The PBS log
  shows only the wrapper's echo lines and no output dir yet, which looks alarming
  but is not: Julia block-buffers stdout to a file. `qstat -f` confirms real work
  — `cpupercent = 153`, `cput = 00:40:42` over 1 h walltime, 8.5 GB resident. It
  is compiling. Do NOT read the quiet log as a hang; check `resources_used`.
- **AN OVERRIDE IS ALREADY IN PLACE — do not assume stage 2e starts from
  nothing.** `~/.julia/artifacts/Overrides.toml` already maps the
  `optimal_lai_inputs` artifact (hash `2ce5e5e0…`) to
  `/glade/derecho/scratch/arenchon/artifacts_local/optimal_lai_inputs`. So the
  mechanism is live and the *stock* artifact is already shadowed by a local copy.
  Write the rebuilt IC to a NEW directory and repoint the override, rather than
  overwriting that file — otherwise the current IC is destroyed with no way back.
- **IC file template, read off the existing file** (match it): dims `lon = 360`,
  `lat = 180`; `lon`/`lat` as `double` with `degrees_east`/`degrees_north`; six
  `double` variables named `gsl`, `a0_annual`, `precip_annual`, `vpd_gs`,
  `lai_init`, `f0`.
  Its own `precip_annual` description independently confirms the conversion
  factor — "1 mm = 55.51 mol H2O m^-2", i.e. 5.551e4 mol m^-3, matching the
  ρ_liq/M_H2O = 5.5509e4 derived earlier.
- **builder written and ready:**
  `/glade/derecho/scratch/arenchon/claude/build_optimal_lai_ic.jl`, takes
  `<simdir> <out.nc>`, takes the LAST time slice of each of the six diagnostics,
  applies the `pra` conversion, and reports per-field finite counts and ranges so
  a mostly-NaN or collapsed field is visible immediately rather than after
  stage 3 starts.
- **driver:** `experiments/long_runs/snowy_land_pmodel.jl`, `PROGNOSTIC_LAI` set, 10 years, GPU
- **required diagnostics:** `lai`, `a0a`, `pra`, `olf0`, `olvpd`, `olgsl`
- **PBS job ids:** —
- **run output:** —
- **new artifact path:** `/glade/derecho/scratch/arenchon/artifacts_local/optimal_lai_inputs_calibrated/optimal_lai_inputs.nc`
  (the previous local IC at `.../optimal_lai_inputs/` is UNTOUCHED, and the
  pre-edit `Overrides.toml` is backed up under
  `/glade/derecho/scratch/arenchon/claude/Overrides.toml.bak-*`, so this is
  revertible by editing one line back)
- **`~/.julia/artifacts/Overrides.toml` verified:** **yes** — `ClimaLand.Artifacts.optimal_lai_initial_conditions_path()` resolves to the rebuilt file and it exists
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
- **STAGE 2 RESULT (2026-08-01 13:20).** The 10-year run produced all six IC
  diagnostics with 120 monthly steps each. IC built by
  `/glade/derecho/scratch/arenchon/claude/build_optimal_lai_ic.jl`.

  | field | rebuilt | existing IC | note |
  |---|---|---|---|
  | lai_init | mean 0.578, max 5.34 | mean 0.502, max 5.83 | comparable |
  | a0_annual | mean 276, max 628 | mean 126, max 479 | higher — calibrated GPP |
  | precip_annual | mean 3.302e4 | mean 3.281e4 | **within 0.6 %** |
  | vpd_gs | mean 1455 | mean 630 | higher |
  | gsl | mean 182, min 0.21 | mean 193, min 30 | see floor note |
  | f0 | mean 0.457 | mean 0.428 | comparable |

  **The `precip_annual` means agree to 0.6 %.** That is an independent check on
  the unit conversion — a raw write would have been ~4.7 orders of magnitude out,
  so this confirms the ρ_liq/M_H2O factor end to end.

- **Two judgement calls, both deliberate:**
  1. **NaN fill over zero-vegetation cells.** `a0a` and `olvpd` were NaN on 6921
     land cells — 6843 of them in lat −90..−60, i.e. Antarctica, and EVERY one
     with `lai < 0.01`. Potential GPP is undefined with no vegetation and
     `vpd_gs = VPDA0/A0` is 0/0. The reader ingests values directly, so a NaN
     would propagate into the model state; the runtime `Ao_annual_safe` guards
     protect the computation, not the IC. Filled `a0_annual = 0` (physically
     right, and the stock file's minimum is also 0), and `vpd_gs`/`f0` with the
     median of their finite values — with A0 = 0 those cells are inert, so the
     values cannot affect the solution and only need to be finite and safe in a
     denominator. All six fields are now finite over all 22420 land cells.
  2. **`gsl` is NOT floored at 30 days**, even though the stock file's own
     description advertises "Minimum GSL = 30.0 days" and ours reaches 0.21.
     `optimal_lai.jl:224-231` puts `GSL` in the NUMERATOR of Eq. 20 and guards
     the denominator with `Ao_annual_safe`/`fAPAR_max_safe`, commenting that
     "m ~ 0 naturally" as A0 → 0. So the floor is a generation-time choice in the
     stock file, not a model requirement — and imposing it would overwrite the
     model's own equilibrium, which is the entire point of rebuilding the IC.
- **notes:** —

## Stage 3 — prognostic LAI

- **status:** **complete, but NOT converged — parameters deliberately NOT written to the TOML.** See the verdict below.
- **config:** `experiments/calibration/configs/lai.jl`
- **ensemble:** 11 (5 params, TransformUnscented)
- **validity mask built for this grid:** **yes** — built 2026-08-01 13:30 from
  stage 2's long-run diagnostics
  (`.../snowy_land_pmodel_opt_lai_longrun_gpu/global_diagnostics/output_0000`).
  **Verified it is a valid reference despite the domains differing.** The long run
  uses `global_box_domain` defaults while the calibration sets
  `nelements = (180, 360, 15)`, so the footprints could in principle disagree.
  Compared the NaN footprint of the long run's `lai` against a stage-1 member's
  `gpp`: both 22420 finite cells, **identical footprint, zero disagreement in
  either direction**. So the prompt's claim that the long run is "on this grid"
  holds empirically, not just by assertion.
  (Note when reading these files directly: the arrays are `(time, lon, lat)`,
  not `(lon, lat, time)`.)
- **natural-vegetation mask built for this grid:** yes, pre-built 2026-07-31 at
  `experiments/calibration/natural_vegetation_mask.jld2` (gitignored, so it is
  on disk in the checkout only). 10113 of 22420 land cells kept. Rebuild only if
  `nelements` changes.
- **PBS job ids:** `6976711` (chunk 1, `N_ITERATIONS=1`). An earlier submission
  `6976710` was cancelled within a minute — see the config note below.
- **output dir:** `/glade/derecho/scratch/arenchon/claude/calibration_stage3_lai`
- **⚠ `lai.jl` NEEDED THE SAME `N_ITERATIONS` OVERRIDE, and did not have it.**
  The chunking override was only ever added to the stage-1 config, so the first
  stage-3 submission silently ignored `N_ITERATIONS=1` and would have attempted
  all 5 iterations in one 5 h 30 m job. At ~1 h 15 m per iteration that overruns
  the walltime, and a kill mid-iteration orphans member jobs that then race the
  resubmitted orchestrator for the same member directories. Cancelled it,
  added the override to `lai.jl` (default still 5), and resubmitted.
- **DIFFERENT FROM STAGE 1: the misfit series here IS a valid convergence
  signal.** `lai.jl` sets `minibatch_size = 1` with a single
  `sample_date_ranges` entry, so every iteration scores the SAME data. The
  reasoning that invalidated the "loss decreased across iterations" check for
  stage 1 (16 samples at minibatch 2, a different pair of years each iteration)
  does not apply here. Use the loss trend for stage 3.
- **ITERATION 1 COMPLETE (2026-08-01 15:12).** Chunk 1 ran in 1 h 39 m; **11/11
  members, zero segfaults, zero failed members.** Chunk 2 (`6977375`,
  `N_ITERATIONS=2`) submitted 15:20.
  G ensemble `(40424, 11)`, **0 NaN columns, finite fraction 1.0** — both masks
  are doing their job (40424 ≈ 10113 natural-vegetation cells × 4 scored months).
  Misfit: `[0.2229]`.

- **THE REBUILT IC WORKS — this was the real test of stage 2.** `member_001`'s
  simulated LAI shows a clean two-cycle seasonal signal: global mean oscillating
  0.55 → 1.36, maxima ~5.4-5.6, finite throughout. No collapse, no NaN.
  Crucially the winter minimum (~0.55) sits essentially at the IC's `lai_init`
  mean (0.578), i.e. **the model starts near its own equilibrium**, which is
  exactly what stage 2 existed to achieve — the 1-year spinup is adequate despite
  `optimal_lai_tau_long_term` being 2 years.

  | parameter | prior mean | after it.1 | bounds | % of range |
  |---|---|---|---|---|
  | optimal_lai_z | 14.92 | 9.326 | 1–40 | 21 % |
  | optimal_lai_z_c4 | 14.92 | **36.45** | 1–40 | **91 %** ⚠ |
  | optimal_lai_sigma | 0.9293 | 0.4933 | 0.1–3.0 | 14 % |
  | optimal_lai_sigma_c4 | 0.9293 | **2.734** | 0.1–3.0 | **91 %** ⚠ |
  | optimal_lai_alpha | 0.0685 | 0.1869 | 0.01–0.3 | 61 % |

- **⚠ WATCH: both C4 parameters are at ~91 % of their range after ONE step.**
  `optimal_lai_z_c4` → 36.45 against a bound of 40, and `optimal_lai_sigma_c4` →
  2.734 against 3.0. Deliberately NOT calling this pinned yet: at stage 1
  `canopy_z_0b_coeff` reached 97.7 % of its range at iteration 1 and then backed
  off, so one iteration proves nothing.
  There is a plausible reason it may be real, though. The C3/C4 split is NEW in
  this PR, and `lai.jl:82-85` gives the C4 parameters priors IDENTICAL to the C3
  ones (`z_c4` mean 15.0 like `z`; `sigma_c4` mean 0.939 like `sigma`) — the C4
  priors are borrowed C3 defaults rather than independently chosen. If C4
  vegetation genuinely wants a much larger z and sigma, that borrowed range is
  simply too narrow, and the parameters will pin. **If they are still ≥95 % at
  iteration 3-5, report it as a mis-specified prior — do NOT widen the bounds
  mid-run.**
- **ITERATION 2 COMPLETE (2026-08-01 16:36).** Chunk 2 ran in 1 h 16 m; 11/11
  members, zero segfaults. Chunk 3 (`6978126`) submitted 16:40. G ensembles still
  `(40424, 11)` with 0 NaN and finite fraction 1.0.

  | parameter | it.1 | it.2 | it.3 | bounds | % of range (it.3) |
  |---|---|---|---|---|---|
  | optimal_lai_z | 14.92 | 9.326 | 20.28 | 1–40 | 49 % |
  | optimal_lai_z_c4 | 14.92 | 36.45 | 31.04 | 1–40 | 77 % |
  | optimal_lai_sigma | 0.9293 | 0.4933 | 0.8209 | 0.1–3.0 | 25 % |
  | optimal_lai_sigma_c4 | 0.9293 | 2.734 | 2.349 | 0.1–3.0 | 78 % |
  | optimal_lai_alpha | 0.0685 | 0.1869 | 0.1661 | 0.01–0.3 | 54 % |

- **The C4 bound concern has RECEDED.** `z_c4` 91 % → 77 %, `sigma_c4` 91 % →
  78 %. Same pattern as `canopy_z_0b_coeff` at stage 1, which is why it was worth
  not calling it pinned after one step. Keep watching, but on current evidence
  the borrowed C3 priors are not obviously too narrow.

- **⚠ THE MISFIT INCREASED: 0.2229 → 0.3154 (+41 %).** Unlike stage 1, this
  number IS meaningful — `lai.jl` uses `minibatch_size = 1` over a single sample
  range, so both values score the SAME data. So this is a real rise in
  misfit, not a change of scoring set.
  Not yet acting on it. Two points are not a trend, EKI misfit can rise
  transiently while `DataMisfitController` adapts its step, and the parameters
  are still swinging hard (`z` 9.3 → 20.3, `sigma` 0.49 → 0.82), which is
  consistent with a large early step rather than divergence.
  **Decision rule for the remaining iterations:** if the misfit is still above
  the iteration-1 value of 0.223 at iteration 5, that is a genuine negative
  result and must be reported as such, NOT papered over. The knobs the prompt
  offers, in order, would be `NOISE_SCALARS["lai"]` (currently 0.5; raise toward
  1.0 if the fit is overfitting seasonal amplitude at the cost of timing) and
  then the spinup. Any such change gets its own commit and PR note explaining
  why — and is a next-run decision, not a mid-run one.
- **ITERATION 3 COMPLETE (2026-08-01 17:50).** Chunk 3 ran in 1 h 10 m; 11/11
  members reached `completed` and there were no segfaults — but see the crash
  below. Chunk 4 (`6978373`) submitted 17:53.
  Misfit: `[0.2229, 0.3154, 0.2233]` — **the rise reversed**. Iteration 3 is back
  to essentially the iteration-1 value, so the +41 % at iteration 2 was a
  transient large step, as suspected, not divergence. It is now FLAT rather than
  decreasing, which is its own question for iterations 4-5.

  | parameter | it.2 | it.3 | it.4 | bounds |
  |---|---|---|---|---|
  | optimal_lai_z | 9.326 | 20.28 | 27.55 | 1–40 |
  | optimal_lai_z_c4 | 36.45 | 31.04 | **8.17** | 1–40 |
  | optimal_lai_sigma | 0.4933 | 0.8209 | 1.062 | 0.1–3.0 |
  | optimal_lai_sigma_c4 | 2.734 | 2.349 | 1.570 | 0.1–3.0 |
  | optimal_lai_alpha | 0.1869 | 0.1661 | 0.1676 | 0.01–0.3 |

  **The C4 parameters are oscillating violently**, not converging: `z_c4` has gone
  14.9 → 36.5 → 31.0 → 8.2 (a 74 % drop in the last step), `sigma_c4`
  0.93 → 2.73 → 2.35 → 1.57. Only `optimal_lai_alpha` looks settled (0.1661 →
  0.1676, +1 %). With 5 iterations configured this is very unlikely to converge.

- **⚠ A MEMBER CRASHED WITH A `DomainError`, AND THIS ONE IS NOT TRANSIENT.**
  Iteration 3's G ensemble has **1 NaN column of 11** (finite fraction 0.909);
  the orchestrator logged `Error processing member 6, filling observation map
  entry with NaNs`. Its `model_log.txt` (41 MB of repeated output) contains:

      ERROR: a DomainError was thrown during kernel execution on thread ...
      Error: ClimaLand simulation crashed.

  Parameters were `z = 20.28, z_c4 = 31.20, sigma = 0.8208, sigma_c4 = 2.351,
  alpha = 0.1895` — **all inside the prior bounds.**
  This is CATEGORICALLY DIFFERENT from stage 1's `member_003`, which segfaulted
  once on a bad node and whose data turned out valid. This is deterministic: a
  parameter combination the optimal-LAI model cannot evaluate, reached by the
  EKI while exploring its own prior. A tell worth remembering — **member 6
  finished conspicuously EARLY** (it was the only one complete while ten still
  ran), because it crashed rather than finished.
  Handled for now (NaN column + `SampleSuccGauss`), but it is a real robustness
  finding to REPORT: with only 11 sigma points, losing one degrades the unscented
  covariance, and the violent C4 oscillation above may partly be a consequence.
  **Watch iterations 4-5 for recurrence, and record which parameter sets crash.**
- **ITERATION 4 COMPLETE (2026-08-01 19:32).** 11/11 members, **no DomainError
  crash** (0 in iteration 4), G ensemble clean again. Chunk 5 — the final
  configured iteration — submitted 19:35. So the iteration-3 crash was a one-off
  parameter combination, not a systematic failure.
  **Misfit: `[0.2229, 0.3154, 0.2233, 0.2032]` — now the LOWEST of the run**, and
  below the iteration-1 value. The stage-3 fit is genuinely improving. (Recall
  this metric IS valid here: `minibatch_size = 1` over a single sample, so every
  iteration scores the same data.)

  | parameter | it.1 | it.2 | it.3 | it.4 | it.5 |
  |---|---|---|---|---|---|
  | optimal_lai_z | 14.92 | 9.326 | 20.28 | 27.55 | 21.02 |
  | optimal_lai_z_c4 | 14.92 | 36.45 | 31.04 | 8.17 | **34.85** |
  | optimal_lai_sigma | 0.9293 | 0.4933 | 0.8209 | 1.062 | 0.9361 |
  | optimal_lai_sigma_c4 | 0.9293 | 2.734 | 2.349 | 1.570 | 1.992 |
  | optimal_lai_alpha | 0.0685 | 0.1869 | 0.1661 | 0.1676 | 0.1260 |

- **⚠ THE C4 PARAMETERS ARE NOT CONVERGING, AND THE LOSS IMPROVES ANYWAY.**
  `optimal_lai_z_c4` has gone 14.9 → 36.5 → 31.0 → 8.2 → 34.9, traversing almost
  its entire 1–40 prior range while the misfit falls. A parameter that can be
  moved across its whole range without hurting the fit is not being identified by
  the data.

- **Tested a hypothesis for why, and REFUTED it.** I expected the
  natural-vegetation mask to be the cause, since it drops cropland and much C4 is
  maize. Measured instead: majority-C4 cells are **20.8 % of scored cells versus
  18.6 % of all vegetated land**, and mean `fc3` over scored cells is 0.731
  (C4 fraction 0.269). **The mask slightly ENRICHES C4 rather than removing it.**
  Cropland removal is not the explanation, and it would have been an easy and
  wrong thing to assert.

- **Better-supported explanation: a degeneracy between `z_c4` and `sigma_c4`.**
  The two move TOGETHER across iterations — (36.5, 2.73), (31.0, 2.35),
  (8.2, 1.57) — which is what a compensating ridge looks like: raising one and
  raising the other produce similar LAI, so the EKI slides along the ridge
  instead of converging to a point. Both enter `compute_m` (Eq. 20) in the same
  direction, which is a plausible mechanism.
  **This is a finding to REPORT, not to tune away.** If it is a genuine
  degeneracy, more iterations will not fix it; the C3/C4 split may need either an
  independent constraint on C4 (a C4-specific target) or the two parameters tied
  rather than free.
- **STAGE 3 FINISHED ALL 5 ITERATIONS (2026-08-01 20:40).** Every iteration
  11/11 members except iteration 3 (one `DomainError` crash, absorbed as a NaN
  column). Final misfit series: **`[0.2229, 0.3154, 0.2233, 0.2032, 0.1984]`** —
  falling steadily over the last three iterations and ending at the best value of
  the run, **11 % better than iteration 1**. This metric is valid here
  (`minibatch_size = 1`, one sample, same data every iteration).

  | parameter | it.1 | it.2 | it.3 | it.4 | it.5 | **final** | bounds |
  |---|---|---|---|---|---|---|---|
  | optimal_lai_z | 14.92 | 9.326 | 20.28 | 27.55 | 21.02 | **24.29** | 1–40 |
  | optimal_lai_z_c4 | 14.92 | 36.45 | 31.04 | 8.17 | 34.85 | **14.61** | 1–40 |
  | optimal_lai_sigma | 0.9293 | 0.4933 | 0.8209 | 1.062 | 0.9361 | **1.053** | 0.1–3.0 |
  | optimal_lai_sigma_c4 | 0.9293 | 2.734 | 2.349 | 1.570 | 1.992 | **1.398** | 0.1–3.0 |
  | optimal_lai_alpha | 0.0685 | 0.1869 | 0.1661 | 0.1676 | 0.1260 | **0.1295** | 0.01–0.3 |

- **VERDICT: the loss improved but the parameters are NOT identified. They are
  NOT being written to `toml/default_parameters.toml`.**
  `optimal_lai_z_c4` visited 14.9 → 36.5 → 31.0 → 8.2 → 34.9 → 14.6, i.e. most of
  its 1–40 prior range, while the misfit fell monotonically over the last three
  steps. A parameter that can be moved across its whole range without hurting the
  fit has not been determined by the data, and the final value is an arbitrary
  point on a ridge rather than an optimum. Writing it as a model default would
  present a partial failure as a success, which the rules explicitly forbid. The
  existing defaults are left in place; the numbers above are recorded here and in
  the PR for a human to act on.

- **QUANTIFIED THE DEGENERACY (this is the substantive finding).**
  Within the FINAL ensemble (11 members):

  | pair | Pearson r |
  |---|---|
  | `optimal_lai_z` ↔ `optimal_lai_sigma` | **+0.860** |
  | `optimal_lai_z_c4` ↔ `optimal_lai_sigma_c4` | +0.452 |
  | `optimal_lai_z` ↔ `optimal_lai_z_c4` | −0.438 |

  Across-iteration correlation of the ensemble MEANS: `z`↔`sigma` **+0.875**,
  `z_c4`↔`sigma_c4` **+0.822**, `sigma`↔`sigma_c4` −0.774, `sigma_c4`↔`alpha`
  +0.822.
  Final ensemble spreads: `z` CV **1.1 %**, `sigma` CV **0.7 %**,
  `sigma_c4` 5.3 %, `alpha` 7.1 %, `z_c4` **19.8 %**.

  **Read:** the ensemble is very tight (CV ~1 %) yet its mean wanders by a factor
  of three between iterations. Tight spread + wandering mean + r ≈ +0.87 between
  `z` and `sigma` is a RIDGE: the data constrains a combination of the two, not
  each separately. Both enter `compute_m` (Eq. 20) in the same direction —
  `m = (sigma * GSL * LAI_max) / (A0_annual * fAPAR_max)` with `z` setting
  `LAI_max` — so raising either produces nearly the same LAI. The C4 pair shows
  the same structure and is worse identified because C4 is a minority of the
  scored cells.
  **More iterations will not fix this.** It is structural, not a convergence
  problem.
- **calibrated parameters:** recorded above; **deliberately NOT committed**
- **committed to `toml/default_parameters.toml`:** **no — and this is the correct
  outcome, not an omission.**
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

---

## Stage 3 follow-up — experiments to improve LAI RMSE and bias

**Measured baselines** (leaderboard method: mean over months of area-weighted
`ClimaAnalysis.global_rmse` / `global_bias`; matched years 2017-18). Validated
against ClimaLand's own leaderboard, whose LAI RMSE curve runs 0.806-1.059 and
reproduces the unmasked numbers exactly.

| configuration | RMSE unmasked | RMSE masked | bias unmasked | bias masked |
|---|---|---|---|---|
| default optimal-LAI params | 1.046 | 0.941 | +0.225 | +0.081 |
| single-year calibrated (u[5]) | 0.926 | **0.885** | +0.081 | **−0.098** |

So the single-year calibration cut RMSE 6 % masked / 11.5 % unmasked and cut the
unmasked bias by ~64 %, but OVERSHOT the masked bias from +0.081 to −0.098. On
natural vegetation the uncalibrated model was already nearly unbiased.

⚠ **Two measurement traps, both of which caught me:**
1. RMSE of the ANNUAL-MEAN field (~0.674) is ~25 % lower than the mean of
   per-month RMSEs (0.885) because monthly errors partly cancel. The leaderboard
   uses the per-month form. Always say which is quoted.
2. `iteration_00N/member_001` is the mean of u[N], i.e. the SECOND-TO-LAST
   ensemble — not the final calibrated parameters u[N+1]. The 0.885 above is for
   u[5] (`z_c4` = 35.0), while the committed values are u[6] (`z_c4` = 14.6).
   The committed set's RMSE is UNMEASURED; it needs one forward run.

### Experiment A — more data (running)

`configs/lai_multiyear.jl`, output
`/glade/derecho/scratch/arenchon/claude/calibration_stage3_lai_multiyear`,
chunked the same way (`N_ITERATIONS=1,2,…`). First chunk `6979038`.
Four annual cycles (2015-2018) at `minibatch_size = 2`, `n_iterations = 8`
(four full epochs). Priors, noise scalar, spinup and masks IDENTICAL to
`lai.jl`, so this isolates the effect of the target's information content.
Confirmed live: "The number of samples is 4", members span ~3 years, so the cost
is ~1.6x per member rather than 4x — `minibatcher_over_samples` batches
CONSECUTIVE samples and `models/snowy_land.jl:171` spans only the minibatch.

**What it tests.** Whether the ridge (`z`↔`sigma` r = +0.86; `z_c4` roaming most
of its 1-40 range while the loss fell) is a DATA problem or a STRUCTURAL one.
- If `z_c4` pins down with four years → data. Extend further.
- If it still roams → structural degeneracy. Next experiment: drop the C3/C4
  split so `z_c4`/`sigma_c4` stop being free, or reparameterise onto
  `log(z·sigma)` / `log(z/sigma)`.

**Judge it on RMSE and bias, not the EKI loss.** The single-year run cut its loss
11 % while masked RMSE fell only 6 % and masked bias got worse — the normalised
misfit and the leaderboard metric are different functionals and decouple on a
ridge. Use `/glade/derecho/scratch/arenchon/claude/lai_compare.jl`, which takes a
simdir and reports both, masked and unmasked.

### Experiment A RESULT — more data does NOT help (2026-08-02)

**Falsified my own top recommendation.** Quadrupling the target from one annual
cycle to four reproduced the single-year trajectory almost exactly, including the
pathology:

| iteration | z single / multi | sigma single / multi | z_c4 single / multi | misfit single / multi |
|---|---|---|---|---|
| 2 | 9.33 / 7.75 | 0.493 / 0.47 | 36.45 / 39.27 | 0.315 / 0.317 |
| 3 | 20.28 / 21.02 | 0.821 / 0.825 | 31.04 / 36.19 | 0.223 / 0.239 |
| 4 | 27.55 / 27.44 | 1.062 / 1.065 | 8.17 / 15.09 | 0.203 / 0.205 |
| 5 | 21.02 / 20.43 | 0.936 / 0.923 | 34.85 / 36.62 | 0.198 / — |

`z` and `sigma` agree to three significant figures at every iteration, and
`z_c4` still swings 15.09 → 36.62. Four times the observations changed neither
the parameter path, the loss, nor the oscillation.

**Therefore the z-sigma ridge is structural, not a data shortage.** "Calibrate
against more years" was the top-ranked suggestion in the report; it is now
falsified and should be struck from the recommendations. The remaining
suggestions that address the ridge itself — reparameterising onto
`log(z·sigma)` / `log(z/sigma)`, fixing one of the pair from theory, or dropping
the C3/C4 split — are unaffected, and are the ones worth pursuing.

**Experiment B** (`configs/lai_noc4split.jl`, output
`.../calibration_stage3_lai_noc4split`) tests the last of those. ⚠ Caveat it
must be read with: `z_c4` and `sigma_c4` are FROZEN at 14.6 / 1.4, values
inherited from the unconverged single-year run. B started with a worse misfit
(0.387 vs A's 0.231) partly for that reason, so a worse loss for B does NOT show
the split was useful — it may only show the frozen values are poor. Judge on
RMSE/bias.

**Experiment A fit result — no improvement either.** Measured on
`iteration_004/member_001` (u[4]), 24 months 2016-12..2018-11:

| configuration | RMSE unmasked | RMSE masked | bias unmasked | bias masked |
|---|---|---|---|---|
| default params | 1.046 | 0.941 | +0.225 | +0.081 |
| single-year calibrated | 0.926 | **0.885** | +0.081 | −0.098 |
| experiment A, 4 years | 0.912 | 0.895 | +0.036 | −0.146 |

Masked RMSE is marginally WORSE than the single-year run and the masked bias is
worse. Everything is within ~1.5 % on RMSE. Parameter trajectory and fit agree:
four years bought nothing.

### A pattern worth acting on: every calibration drives the masked bias negative

+0.081 (default) → −0.098 (single year) → −0.146 (four years). Calibration
systematically pushes LAI DOWN on natural vegetation. `z` rises 15 → 20-27 in
every run, and larger `z` means lower peak LAI (see the prior docstring), so this
is the optimiser trading one large over-prediction against a widespread mild
under-prediction — which pure squared error rewards.

**Implication: the bias is not a tuning failure, it is what the objective asks
for.** No parameter search will remove it while the loss is unpenalised squared
error. Options, if unbiasedness on natural vegetation matters as much as RMSE:
add an explicit bias penalty to the objective; weight cells by observed LAI so
high-LAI regions are not sacrificed; or score seasonal amplitude and phase
separately so amplitude errors cannot be paid for with a mean offset.
This is a DIFFERENT problem from the z-sigma ridge: the ridge governs whether
parameters are identifiable, this governs where the optimum sits.

### ⚠ The `DomainError` crash is reproducible and correlates with high `alpha`

Three occurrences so far, all `DomainError` on the SAME GPU thread/block
(156,1,1)/(25,1,1), i.e. the same grid cell:

| run | z | z_c4 | sigma | sigma_c4 | alpha |
|---|---|---|---|---|---|
| stage 3 single-year, it3 member 6 | 20.28 | 31.20 | 0.8208 | 2.351 | **0.1895** |
| experiment B, it2 member 3 | 13.26 | (fixed) | 0.6452 | (fixed) | **0.2322** |
| experiment B, it2 member 4 | 13.26 | (fixed) | 0.6297 | (fixed) | **0.2416** |

All parameter sets are INSIDE the priors. The common factor is elevated `alpha`
(the LAI acclimation rate, τ = 1 day / alpha, so 0.24 is τ ≈ 4 days — very fast).

**The stack trace does NOT point at the cause.** It surfaces in
`interpolate!` / `_collect_interpolated_values!` inside the ClimaDiagnostics
NetCDF writer, but the `CUDA/exceptions.jl:39` and `synchronization.jl` frames
show this is a DEFERRED kernel exception: the fault happened in an earlier
kernel and is only reported at the next synchronisation, which happens to be the
diagnostics remap. Anyone debugging this should NOT start in the writer. Re-run
the failing parameter set with `CUDA_LAUNCH_BLOCKING=1`, or reproduce on CPU, to
get the true frame.

**This hurt experiment B disproportionately.** B has only 7 members (3 params),
so losing 2 in iteration 2 is 29 % of the ensemble against 9 % for an 11-member
run. Fewer parameters means fewer sigma points means each crash costs more — an
unintended downside of the reduced-parameter design, and a reason to read B's
covariance with caution.

### Experiment C, iteration 1 — and a caution about comparing misfits

C's misfit is 0.4453 against the baseline's 0.2229 — almost exactly 2x. That is
arithmetic, not a worse fit: the misfit is normalised by the covariance, so
halving the variance (0.5 -> 0.25) doubles it for an IDENTICAL residual.
**Misfit values are not comparable across noise settings.** Do not read C's
higher loss as a regression.

The parameter step is nearly unchanged from the baseline:

| | baseline (0.5) | C (0.25) |
|---|---|---|
| z | 9.326 | 8.837 |
| sigma | 0.4933 | 0.4835 |
| alpha | 0.1869 | 0.1877 |
| z_c4 | 36.45 | 37.54 |

**Hypothesis to check at iterations 2-5: the noise scalar may be largely a
no-op here**, because `DataMisfitController` sets its step from the misfit
magnitude, so scaling the covariance scales the misfit and the controller
normalises it back out. If that holds, `NOISE_SCALARS` is not the lever the
`lai.jl` docstring claims ("tighten it -> 0.25 to pull harder on peak-LAI
magnitude"), and that docstring should be corrected. The accumulated-Δt schedule
can still diverge later, so do not conclude from iteration 1 alone.

### Crash rate tracks `alpha` across experiments — independent confirmation

| experiment | crashes / member-runs | rate | alpha range reached |
|---|---|---|---|
| A (5 params, 4 yr) | 1 / 77 | 1.3 % | 0.135 – 0.210 |
| B (3 params, 4 yr) | 3 / 28 | **10.7 %** | **0.233 – 0.242** |
| C (1 yr, noise 0.25) | 0 / 22 | 0 % | 0.188 |

B crashes ~8x more often than A, and B is the run that pushes `alpha` highest.
The correlation therefore holds BETWEEN experiments, not only within individual
failures, which is much stronger evidence than the three-case table above.

**Working conclusion: `optimal_lai_alpha` above roughly 0.19 (τ ≲ 5 days) risks a
`DomainError` in a GPU kernel.** The prior allows up to 0.3, and every
calibration drives alpha upward from 0.0701, so the optimiser reliably walks into
the unsafe region. Two independent fixes, both worth doing:
1. Guard the expression that faults (find it with `CUDA_LAUNCH_BLOCKING=1` — NOT
   by reading the deferred stack trace, which points at the diagnostics writer).
2. Until then, cap the `optimal_lai_alpha` prior at ~0.18 so the ensemble cannot
   reach the unsafe region. This costs little: no run has settled above 0.19
   anyway, and it would have prevented all four crashes seen so far.

### Experiment C RESULT — the noise scalar is very nearly a no-op

| | baseline (0.5) | C (0.25) | ratio |
|---|---|---|---|
| misfit it1 | 0.2229 | 0.4453 | 2.00x |
| misfit it2 | 0.3174 | 0.6329 | 1.99x |
| z at it3 | 20.28 | 19.27 | −5 % |
| sigma at it3 | 0.8209 | 0.8012 | −2 % |
| alpha at it3 | 0.1661 | 0.1713 | +3 % |

The misfit scales by exactly the variance ratio while the parameters move by
under 5 %. `DataMisfitController` sets its timestep from the misfit magnitude, so
scaling the covariance uniformly is normalised back out and the trajectory is
essentially unchanged.

**`NOISE_SCALARS` is therefore not the lever `lai.jl`'s docstring claims.** That
docstring says "tighten it (-> 0.25) to pull harder on peak-LAI magnitude" and
should be corrected: with an adaptive-timestep scheduler, a uniform variance
rescaling changes almost nothing. It would still matter for the RELATIVE
weighting between several targets (as in stage 1's four variables), just not for
a single-target calibration like this one.

### Running score: three config knobs tested, three null results

| lever | experiment | effect on the calibration |
|---|---|---|
| more data (1 yr -> 4 yr) | A | none — trajectory identical to 3 s.f.; masked RMSE 0.895 vs 0.885 |
| tighter noise (0.5 -> 0.25) | C | none — misfit rescales, parameters move <5 % |
| fewer parameters (5 -> 3) | B | inconclusive so far, and hampered by a 10.7 % crash rate |

Nothing reachable from the config moves LAI RMSE. The binding constraints are
the z-sigma ridge (identifiability) and the model's peak-LAI ceiling (bias), and
both need a change to the model or the objective, not to settings.

### Experiment B, iteration 4 — the first non-null result, but partial

| | A (5 params) | B (3 params) |
|---|---|---|
| `z` range visited | 7.75 – 27.44 (span 19.7) | 13.26 – 25.9 (span **12.6**) |
| `sigma` range | 0.47 – 1.065 (span 0.60) | 0.630 – 1.064 (span **0.43**) |
| misfit | 0.231 / 0.317 / 0.239 / 0.205 | 0.387 / 0.484 / 0.358 / 0.337 |

**Dropping the C3/C4 split damps the oscillation by roughly 35 %.** That is the
first lever of the three to show any effect on the calibration's behaviour.

**But it does not converge.** B's `z` still moved −15 % in the last step
(25.9 → 22.01) and `alpha` is now its loosest parameter (CV 12.4 %). This is what
the correlation analysis predicted: the DOMINANT ridge is `z`↔`sigma`
(r = +0.86), so removing `z_c4`/`sigma_c4` eliminates a secondary degeneracy and
leaves the primary one untouched.

**Conclusion so far: worth doing for stability, not sufficient on its own.** The
reparameterisation onto `log(z·sigma)` / `log(z/sigma)` targets the actual ridge
and remains the one untested idea with a mechanism behind it.

B's crashes all fell in iterations 2-3 (2 then 1, none since), tracking alpha's
early peak at 0.233 and its settling to 0.104 — consistent with the alpha
threshold above.

### Experiment A COMPLETE (8/8) — the EKI is in a PERIOD-2 LIMIT CYCLE

| iteration | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|
| z | 20.43 | 24.54 | 20.48 | 24.09 | 20.29 |
| sigma | 0.9229 | 1.069 | 0.9176 | 1.059 | 0.9126 |
| z_c4 | 36.62 | 17.64 | 37.77 | 23.95 | 37.64 |
| misfit | 0.2107 | 0.2012 | 0.2116 | 0.2002 | — |

The parameters alternate between two states (z ≈ 20.3 ↔ 24.3, sigma ≈ 0.91 ↔
1.06) and **the misfit alternates with them** (0.211 ↔ 0.200). Two points on the
ridge with almost identical loss; the filter orbits between them.

**This is a much sharper diagnosis than "it does not converge".** A limit cycle
means MORE ITERATIONS CANNOT HELP — it is a closed orbit, not slow convergence.
Stage 1's `DataMisfitController` terminated on its own; here nothing stops it.

⚠ **It also breaks the obvious convergence test.** "Check the last two iterations
agree" would be FOOLED: iterations 7 and 9 agree to 1 % (20.48, 20.29) while
iteration 8 sits 18 % away. Any convergence check on this problem must compare
iterations TWO apart, or track the cycle amplitude, not adjacent pairs.

A's final crash tally: 1 in 88 member-runs (iteration 8 had one).

### Experiment A FINAL FIT — a complete null result

Measured on `iteration_008/member_001` (u[8], the HIGH state of the limit cycle:
z = 24.09, sigma = 1.059), 24 months 2016-12..2018-11:

| configuration | RMSE masked | bias masked | RMSE unmasked | bias unmasked |
|---|---|---|---|---|
| default params | 0.941 | +0.081 | 1.046 | +0.225 |
| single-year, 5 iterations | **0.885** | −0.098 | 0.926 | +0.081 |
| four-year, 8 iterations | **0.887** | −0.088 | 0.939 | +0.093 |

**Four times the data and eight iterations instead of five land within 0.3 % of
the single-year run.** A is a null result on the parameters (identical
trajectory), on the loss, and now on the fit. Struck from the recommendations.
