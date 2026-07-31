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
