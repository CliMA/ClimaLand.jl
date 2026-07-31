# Prognostic carbon loop — narrative log

One dated entry per iteration. `STATE.md` is where the work is; this is how it
got there. Record what was tried, what the numbers were, and what was concluded
— including negative results.

---

## 2026-07-31 — setup

Branch `ar/prognostic_carbon` created off `ar/climate_responsive_lai_inputs`.
Model specification written to `../MODEL.md`; loop scaffolding
(`loop_prompt.md`, `STATE.md`, `README.md`, this file) and the staged sandbox
settings (`.claude/unsupervised-loop.settings.json`, `.claude/hooks/guard.sh`,
ported from `origin/ar/derecho_loop` with the repo path fixed to this checkout)
in place. No code written yet — design under discussion with the user.

## 2026-07-31 — iteration 1: stage-0 harness ported

Ported the single-column test harness from `origin/ar/derecho_loop` into a
self-contained `experiments/integrated/prognostic_carbon/harness/`:
`site_driver.jl`, `run_battery.pbs`, `test_sites.csv`. Kept it separate from
`experiments/integrated/era5/single_site.jl` — that file is this branch's
hard-coded desert experiment, and nothing in CI or the docs references either,
so a new directory avoids repurposing a working experiment.

What changed relative to the `derecho_loop` original, and why:

- `REPO` now points at `/glade/u/home/arenchon/ClimaLand.jl` (the loop's
  checkout), not `~/GitHub/ClimaLand.jl`.
- Dropped the `BETA_IN_A0`, `ONLINE_F0`, `F0_SCALE`, `ONLINE_VPD_GS` and
  `ONLINE_GSL` switches. Verified against `toml/default_parameters.toml`: those
  five parameters no longer exist on this branch, their behaviour now being
  unconditional. Setting them would have written an override for a parameter
  that is never read — silently inert, which is the worst failure mode. The
  overrides that survive are the ones still in the toml: `optimal_lai_z`,
  `z_c4`, `sigma`, `sigma_c4`, `alpha`, `f0`, `z_a0`, `online_c3c4`.
- New `LAI_MODE` switch: `prescribed` (MODIS, default) or `prognostic` (Zhou).
  The settled configuration wants prescribed LAI for stages 1–3 and prognostic
  from stage 4, and the carbon model must work under both. Both `LandModel`
  constructors exist — passing `LAI` selects `PrescribedBiomassModel`, omitting
  it builds `ZhouOptimalLAIModel` — so the switch is a genuine one-liner.
- Replaced the LAI-vs-MODIS scoring block with the stage-0 baseline the loop
  actually needs: post-spinup means of GPP, Ra, Rh, LAI (plus air temperature
  and soil water content for context) to `carbon_metrics.txt` per site, which
  `run_battery.pbs` collects into `baseline_summary.tsv`. GPP/Ra/Rh are
  converted from mol CO₂ m⁻² s⁻¹ to g C m⁻² day⁻¹ (× M_C × 86400 × 1000).
- Job name `pc_battery`, output under
  `/glade/derecho/scratch/arenchon/claude/prognostic_carbon/`, `#PBS -m n`.

Checked the metric lookup rather than assuming it: `output_short_name` is
`"<short>_<schedule>_<suffix>"` (e.g. `gpp_1d_average`), built by
`ClimaDiagnostics.descriptive_short_name`. The driver therefore matches on
`startswith(name, short * "_")`. The original's `occursin` would have been
loose enough to cross-match in principle; anchoring removes the question.

Submitted job **6967513** — a 4-site smoke test (`amazon_central`, `sahara`,
`ozark_us`, `us_great_plains`), 1-year window, prescribed LAI — rather than the
full 20 straight away. If the driver has a porting bug, four sites surface it
just as well as twenty. The full battery follows once this is green.

Operational finding: `module load gh` must be run **bare**. `module` is a shell
function, so piping it (`module load gh | tail`) runs it in a subshell and the
`PATH` change is discarded — which looks exactly like "gh is not installed".
Bare, it works and is already authenticated (account AlexisRenchon).

## 2026-07-31 — iteration 2: smoke test green, 20-site baselines launched

Job 6967513 finished in 7 minutes, `Exit_status=0`, 4/4 sites PASS. The port
works: per-site output under scratch, `carbon_metrics.txt` written, and
`baseline_summary.tsv` collected. Annual means (g C m⁻² day⁻¹):

| site | GPP | Ra | Rh | LAI | NPP/GPP |
|---|---|---|---|---|---|
| amazon_central | 8.83 | 4.08 | 1.60 | 4.78 | 0.54 |
| ozark_us | 3.75 | 2.04 | 0.49 | 1.13 | 0.46 |
| us_great_plains | 2.67 | 1.69 | 0.45 | 0.72 | 0.37 |
| sahara | 0.00 | 1.00 | 0.06 | 0.00 | — |

Amazon GPP of 8.8 g C m⁻² day⁻¹ is ≈ 3.2 kg C m⁻² yr⁻¹, right for tropical
rainforest, and the three vegetated NPP/GPP ratios sit inside the 0.3–0.6 band
stage 2 has to hit. So the existing JULES Ra is not badly wrong where there is
vegetation.

**The desert is a different story, and it is the finding of this iteration.**
Sahara: GPP = 0, LAI = 0, Ra = 1.00 g C m⁻² day⁻¹ — about 365 g C m⁻² yr⁻¹ of
autotrophic respiration with no photosynthesis anywhere in the year. The cause
is structural, not numerical: the JULES term respires *prescribed constant*
area indices (`Rd_ref·RAI`, `Rd_ref·μs·SAI`), and SAI/RAI do not vanish in a
desert. The model is emitting carbon it never fixed, from a pool that does not
exist. That is precisely the open budget MODEL.md §2.1 sets out to close, and
it gives stage 2 a falsifiable prediction: once Ra is computed from the pools,
Sahara Ra must fall to ≈ 0, because the pools there equilibrate near zero.

Two harness fixes, both from things that actually bit:

- **Per-job configs.** `run_battery.pbs` sourced one shared
  `battery.conf`, so two concurrent batteries would have silently shared a
  configuration — the second run would report the first's settings. It now
  sources `battery_${PBS_JOBNAME}.conf` with `battery.conf` as fallback, and
  `RUNROOT` carries the job name. Deleted the smoke-test `battery.conf` so its
  4-site `SUBSET` could not leak into a full run.
- **PBS output naming.** The iteration-1 monitor watched
  `pc_battery.o6967513` and never fired. With `-j oe` and `-o <dir>`, Derecho
  writes `<jobid>.OU` — `6967513.desched1.OU`. A monitor on the wrong path is
  indistinguishable from a job that is still running, which is the failure mode
  worth remembering.

Also noted: a job id can be missing from `qstat` for a few seconds after `qsub`
returns it (`Unknown Job Id`). Propagation delay, not a failed submission —
6967718 looked lost and was in fact already running. Do not resubmit on that
evidence.

Submitted the real stage-0 baselines, both 20 sites, 2-year window with the
first year excluded from the means, running concurrently:

- **6967717** `pc_base_pres` — prescribed MODIS LAI (the stages 1–3 setup)
- **6967718** `pc_base_prog` — prognostic `ZhouOptimalLAIModel` (stage 4+)

Both baselines matter: stage 1's tests must cover both LAI models, so rule 1
needs a reference table under each.

## 2026-07-31 — iteration 3: prescribed baseline recorded; Q10 is 1.0, not 2.0

Job 6967717 finished **20/20 PASS**. The prescribed-LAI stage-0 baseline is
recorded in STATE.md and committed as `harness/baseline_prescribed_lai.tsv`.
6967718 (prognostic LAI) is slower and still running at 33 minutes.

GPP looks right: 8.7 / 8.2 / 7.4 g C m⁻² day⁻¹ at congo / amazon / borneo
(≈3.0–3.2 kg C m⁻² yr⁻¹, correct for tropical rainforest), falling to 0 at the
true deserts. LAI tracks MODIS as it should, since it *is* MODIS here.

**Ra is where the model is wrong, and the baseline says so quantitatively.**
Sorting the 20 sites by NPP/GPP = 1 − Ra/GPP:

- Ra **exceeds** GPP at `mojave_sw_us` (NPP/GPP = −0.024) and
  `alaska_north_slope` (−0.093) — negative annual NPP.
- Ra is positive with GPP = 0 at both true deserts (0.81 g C m⁻² day⁻¹).
- Three cold forest sites fall well below the 0.3–0.6 target band:
  `canada_boreal` 0.13, `ne_china` 0.17, `central_siberia` 0.21.
- The 12 warmer/wetter sites land inside 0.3–0.6 and look healthy.

So the **current model already fails stage 2's acceptance criterion** ("Ra
neither collapses to zero nor exceeds GPP in the annual mean at any site") at
four of twenty sites, all of them cold or arid. Useful to know now: stage 2 is
not a refinement, it is fixing something measurably broken.

**Root cause found, and it changes what stage 2 must do.** Chasing why Ra was so
insensitive: across `alaska_north_slope` (262 K), `mojave_sw_us` (281 K) and
`sahel` (300 K) — a 38 K span — Ra is 1.068, 1.066, 1.094, constant to ±3%,
where `Q10 = 2` predicts a 14× range. The two zero-LAI deserts agree to 4.5e-8
relative despite a 3.45 K difference.

`autotrophic_respiration_Q10` in `toml/default_parameters.toml` is **1.0**, not
2.0. Its own description says it is "fixed at 1.0 (not calibrated), keeping the
maintenance baseline seasonally flat". So `f_T = 1^x ≡ 1` and the temperature
term is disabled. Confirming the arithmetic: at the deserts
Ra = `Rd_ref` × 1.0000085, exactly `Rd_ref` with `f_T = 1` and `RAI + μs·SAI = 1`.

This matters because MODEL.md §2.1's new `Rm` reuses the same `Q10`, and §8
listed it as "2.0, 298.15 K (existing)" — wrong about the existing value. Had
stage 2 been written against that, it would have shipped a `Rm` with
`f_T(T_canopy)` and `f_T(T_soil)` terms that read correctly and do nothing.
Corrected MODEL.md §8 and added a note to §2.1 recording the measurement.
Setting `Q10 = 2.0` is a deliberate change to respiration behaviour, so it
belongs to stage 2 and gets reported, not slipped in.

**A monitor gave a false alarm, and the reason is worth recording.** The
monitor watching 6967718 tested liveness with `! qstat <jid>` and fired
`JOB_GONE`; a bare `qstat` showed the job Running at 33 minutes. `qstat` is in
the sandbox's excludedCommands, so wrapping it in a script forces it back into
the sandbox, where it fails — and a failed `qstat` is indistinguishable from a
finished job. Monitors now key only on files the job writes; bare `qstat`
liveness checks happen at the start of an iteration, where they work. This is
the second monitor bug in three iterations, both of the same family: the
monitor's own failure looking like a result.

**Stage 0 closed later the same iteration.** 6967718 finished **20/20 PASS**
(the rewritten file-based monitor fired correctly on `BATTERY_DONE`). Both
baselines are recorded and committed: `harness/baseline_prescribed_lai.tsv` and
`harness/baseline_prognostic_lai.tsv`.

The two LAI models disagree a lot, which is worth knowing before writing any
carbon code. Prognostic Zhou LAI runs much higher at grassland, savanna and
temperate sites — `us_great_plains` 0.66 → 2.31, `iberia` 0.74 → 2.03,
`n_australia_savanna` 1.17 → 4.52, `ozark_us` 1.12 → 2.52 — and lower at
`california_vaira` 1.74 → 0.82, `central_siberia` 1.29 → 0.72,
`alaska_north_slope` 0.40 → 0.14. At `sahel` it collapses to 0.00 and takes GPP
with it, turning a savanna into bare ground.

Under rule 1 that is emphatically not mine to fix, and I am not touching it. It
is recorded because the carbon model has to work under both, and the pools will
come out very differently between them — a constant allocation-fraction set
that looks right under MODIS LAI may not under Zhou LAI, which is exactly what
the `σl_implied` diagnostic exists to expose.

The respiration problem is the same under both LAI models, which is
reassuring: it is a property of the respiration scheme, not of the LAI input.
Under prognostic LAI the middle of the distribution is healthier (13 sites in
0.4–0.55), but `alaska_north_slope` is much worse — GPP halves to 0.49 while Ra
stays at 0.93, giving NPP/GPP = −0.915 — and `central_siberia` is marginal at
0.028.

Stage 0 is therefore **done**: harness ported, battery green 20/20 under both
LAI models, per-site output in scratch, baselines committed. Stage 1 is
unblocked and is the next work — the four pools in `biomass.jl`, the carbon
conservation test, and the rule-1 check against these two tables.

## 2026-07-31 — iteration 4: stage 1 model and tests

Wrote the carbon pools. Commit `2419a0c73`; 56 new tests pass under Float32 and
Float64, and the existing 468 canopy tests still pass.

**The design decision that shaped everything: a wrapper, not a replacement.**
`:biomass` is a single component slot, and both `PrescribedBiomassModel` and
`ZhouOptimalLAIModel` occupy it. Since the specification requires the carbon
model to compose with either, `PrognosticCarbonModel` wraps an LAI model and
delegates every area-index and LAI decision to it. That makes rule 1 structural
rather than something to remember to check — the LAI path *is* the wrapped
model's, untouched — and the test confirms LAI, GPP and Ra are bit-identical
with and without the pools.

Four things that had to be got right, all found by running rather than by
reading:

- **`update_fractional_c3!` must be forwarded.** It dispatches on the biomass
  model. Unforwarded, the Zhou C3/C4 competition silently stops being applied
  when the carbon model wraps it — which changes GPP and breaks rule 1 quietly.
- **`height` and `rooting_depth` are read as direct fields** across
  plant hydraulics, root interactions and diagnostics, so the wrapper mirrors
  them. Safe only because phase 1 never mutates either.
- **The `CanopyModel` constructor asserted `biomass.plant_area_index.LAI == LAI`**,
  which no wrapper can satisfy. Replaced with `prescribed_lai_input(biomass)`,
  a small accessor the wrapper forwards. One line in `Canopy.jl`.
- **Accessors must be hoisted out of `@.` blocks.**
  `@. lazy(M_C * get_GPP(p, canopy.photosynthesis))` broadcasts over `p` itself
  and throws "broadcasting over dictionaries and `NamedTuple`s is reserved".

The conservation test is the real gate and it is a genuine test, not a
tautology: it integrates one step from a non-trivial pool state and checks the
summed tendencies against `GPP - Ra - litter` independently. Summing the four
pool equations gives

    d(ΣC)/dt = GPP - Rm - S + a·S - litter = GPP - (Rm + (1-a)S) - litter

so the balance only closes if growth respiration is exactly `(1-a)S` and
allocation is neither double-counted nor dropped. Residual is at round-off.

**A parameter decision worth flagging.** The carbon model gets its own
`carbon_Q10 = 2.0`, deliberately separate from `autotrophic_respiration_Q10`,
which iteration 3 found ships at 1.0. Sharing the parameter would have meant
stage 1's pools inheriting a temperature-blind maintenance term — the exact trap
iteration 3 identified. Keeping them separate also means stage 1 changes no
existing behaviour: the legacy scheme still drives the canopy fluxes untouched,
and the decision about its Q10 belongs to stage 2, where it will be reported.

**Deliberately deferred, and stated rather than hidden.** Phenological leaf
shedding is not implemented; only background `C_leaf/τ_leaf` turnover is.
`dLAI/dt` is directly available for the Zhou model but needs finite-differencing
for prescribed LAI, and shipping it for one LAI mode only would make the two
modes structurally different — worse than not shipping it. Expect deciduous leaf
carbon to lag LAI through autumn in the battery.

Stage 1 is **not** done: the gate also requires the site battery to complete
with GPP and LAI unchanged from the stage-0 baseline at every site. Next
iteration wires the carbon model into `harness/site_driver.jl` behind a
`CARBON=1` switch and runs the battery under both LAI modes.
