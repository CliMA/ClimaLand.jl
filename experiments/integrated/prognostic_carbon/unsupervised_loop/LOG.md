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
