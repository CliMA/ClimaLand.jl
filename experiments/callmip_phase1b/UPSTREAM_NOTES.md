# Upstream findings to report (Phase 1b work)

Running log of bugs/issues found in dependencies during Phase 1b, with status.
Owner: Renato (decides when/how to report).

## 1. GPU bug in `callmip_forcing.jl` on `kp/diff-multi-cols` (found 2026-08-18, FIXED locally)

**For:** Kevin Phan (ph-kev), branch `kp/diff-multi-cols` (companion to PR #1826).

**Symptom:** `callmip_run_demo.jl` with `CLIMACOMMS_DEVICE=CUDA` crashes in the snow
model's surface-temperature broadcast:

```
ERROR: ClimaCoreCUDAExt.ThisHost and ClimaCore.DataLayouts.ThisThreadPool do not
overlap, so they cannot be put in the same DataScope
```

(first hit inside `ClimaLand.Snow.solve_for_surface_temp_at_a_point` — but the defect
is upstream of the snow model.)

**Cause:** in `prescribed_forcing_callmip`, the per-column tower height is built as
`ClimaCore.Fields.array2field(heights, surface_space)` with `heights::Vector{FT}` — a
host `Array` wrapped into a Field on a GPU space. Every downstream broadcast that mixes
`atmos.h` with device fields then throws. Invisible on CPU.

**Fix (verified on 2× V100, CPU/GPU final states agree to ~11 digits):**

```julia
device_array = ClimaComms.array_type(ClimaComms.device(surface_space))(heights)
atmos_h = length(heights) == 1 ? heights[1] :
          ClimaCore.Fields.array2field(device_array, surface_space)
```

(+ `import ClimaComms` at the top.) Same pattern applies to ANY future per-column
Field built from a host vector — worth a helper or a doc note on PR #1826.

**Also useful for #1826 thread (kmdeck's question):** first timing points on
Tesla V100 (sm_70, CUDA runtime pinned 12.9 — CUDA 13 dropped Volta), 2-column
ensemble, integrated LandModel, Δt=450 s: GPU ≈ 125 SYPD vs CPU ≈ 780 SYPD —
launch-latency-bound at N=2, as expected. Full N-scaling (10-day windows, JIT
excluded, real 21-site forcing cycled to larger N):

| N | V100 SYPD | V100 column-SYPD | CPU 1-core SYPD | CPU column-SYPD |
|---|---|---|---|---|
| 2 | 145 | 291 | 1099 | 2199 |
| 10 | 150 | 1497 | 505 | 5049 |
| 21 | 151 | 3178 | 289 | 6070 |
| 30 | 147 | 4421 | 212 | 6368 |
| 100 | 151 | 15087 | 69 | 6862 |
| 300 | 146 | 43659 | 24 | 7183 |

V100 SYPD is FLAT to 300 columns (launch-latency-bound; column-throughput linear,
unsaturated at 43.7k col-SYPD). CPU saturates ~7k col-SYPD. Crossover ≈ 50
columns/process: below it a single CPU core beats a V100; above it the GPU wins and
keeps winning linearly — very promising for big-ensemble batching (sites×members) or
newer GPUs. (Caveat: sm_70 V100 on CUDA 12.9; A100+/CUDA 13 should look better.)

## 2. Stale `callmip_phase1` artifact in CliMA/ClimaArtifacts (found 2026-08-18, worked around)

The published artifact (tree `c8014d3…`, created 2026-06-30, PR #172) predates the
callmip-org/Phase1 site swap of 2026-07-17: it contains AU-How flux (now a validation
site) and is missing US-Var (calibration site). Fix: rerun
`callmip_phase1/create_artifact.jl` at current upstream (`d6d36ad` or later) and bump
`Artifacts.toml` in ClimaLand (we already carry the refreshed tree `9b4bd0d8…` locally;
see DATA_MANIFEST.md).

## 3. `callmip_phase1_forcing` should grow to 31 sites (2026-08-18)

The published forcing artifact has only the 21 calibration sites. Phase 1b also
requires met for the 10 ME-org validation sites (protocol Table 2). We built the
31-site tree `a597745c…` locally (see DATA_MANIFEST.md); the upstream ClimaArtifacts
refresh should include them (note: 3 validation sites are hourly, AU zones are +9.5,
provenance is OzFlux/LaThuile for some — README should say so).
