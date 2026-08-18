# CalLMIP Phase 1b — Data Manifest and Verification

Verified 2026-08-18 (Claude, session with Renato). Update this file whenever data
locations, hashes, or upstream contents change.

## Artifacts

| Artifact | git-tree-sha1 | Resolved path (Overrides.toml) | Status |
|---|---|---|---|
| `callmip_phase1_forcing` | `a0c004071c73d047095a01faedf017571e5237ef` | `/net/sampo/data1/ClimaArtifacts/artifacts/callmip_phase1_forcing` | ✅ hash verified 2026-08-18; 21 met files, correct site set |
| `callmip_phase1` (published) | `c8014d3bbd838fcaa01bd2a69523b3fa98f28c5c` | `/net/sampo/data1/ClimaArtifacts/artifacts/callmip_phase1` | ⚠️ STALE (see below); kept for Phase 1a reproduction only |
| `callmip_phase1` (refreshed, used by this branch) | `9b4bd0d80ad1a1ac54ab906ee999af5b952f0e7a` | `/net/sampo/data1/renatob/callmip_artifacts/callmip_phase1_20260818` | ✅ built 2026-08-18 from callmip-org/Phase1 @ `d6d36ad` (27 files, 3.5 MB) |

Overrides live in `~/.julia/artifacts/Overrides.toml` (backup:
`Overrides.toml.bak-20260818`). Both published tree hashes were recomputed with
`Pkg.GitTools.tree_hash` and match ClimaArtifacts `OutputArtifacts.toml` exactly.

### ⚠️ Stale-artifact finding (2026-08-18)

The published `callmip_phase1` artifact (created 2026-06-30, ClimaArtifacts PR #172) was
built **before** callmip-org swapped the Phase 1b site set on **2026-07-17**
(commits: "Delete Data/Phase1b/AU-How_…_Flux.nc" then "Add files via upload" = US-Var):

- artifact contains `AU-How_daily_aggregated_2003-2014_FLUXNET2015_Flux.nc`
  → AU-How is now a **validation** site (protocol Table 2) — must NOT be used for calibration;
- artifact is **missing** `US-Var_daily_aggregated_2001-2014_FLUXNET2015_Flux.nc`
  → US-Var IS a calibration site (protocol Table 1 / current GitHub Data/Phase1b).

Resolution: this branch bumps `Artifacts.toml` `[callmip_phase1]` to the refreshed hash
`9b4bd0d8…`. **Follow-up (open):** refresh the artifact upstream in CliMA/ClimaArtifacts
(rerun `callmip_phase1/create_artifact.jl`) so PR #1784/#the Phase 1b PR can use a
published artifact; also stage the refreshed copy to the shared
`/net/sampo/data1/ClimaArtifacts/artifacts/` (needs write access — dir owned by tapio group).

## Met forcing windows (from filenames; PLUMBER2 `*_Met.nc`, half-hourly, local time)

| Site | Met window | Site | Met window |
|---|---|---|---|
| CA-Qfo | 2004–2010 | RU-Fyo | 2003–2014 |
| CH-Dav | 1997–2014 | US-MMS | 1999–2014 |
| DE-Gri | 2004–2014 | US-NR1 | 1999–2014 |
| DE-Hai | 2000–2012 | US-SRG | 2009–2014 |
| DE-Tha | 1998–2014 | US-SRM | 2004–2014 |
| DK-Sor | 1997–2014 | US-Ton | 2001–2014 |
| FI-Hyy | 1996–2014 | US-Var | 2001–2014 |
| FR-Pue | 2000–2014 | US-Whs | 2008–2014 |
| IT-Lav | 2005–2014 | US-Wkg | 2005–2014 |
| IT-MBo | 2003–2012 | NL-Loo | 1997–2013 |
| IT-Noe | 2004–2014 | | |

- **Union:** 1996–2014 · **Intersection: 2009–2010 only** (US-SRG start × CA-Qfo end)
  → this is why the calibration uses a padded union axis (PLAN.md D3).
- Flux obs files carry the same year labels as met (per-file verification of actual
  first/last valid days: TODO in Phase 4, log results here).
- Met files verified to carry: Tair, Wind, Qair, Psurf, SWdown, LWdown, CO2air, Precip,
  VPD, RH, LAI(±alternative), reference_height, latitude, longitude, elevation.
  Time axis: `seconds since <start>-01-01`, **local standard time**, no UTC-offset
  attribute (checked US-SRG) → UTC offsets must come from our site table (Phase 2).
- Flux files carry: NEE, Qle, Qh + `<var>_uc` + missing-% per day. Daily, gap-days = NaN.

## Non-site forcing

- TRENDY CO2: `Data/Non-site-specific_forcing/CO2_1700_2024_TRENDYv2025.txt`
  (in both flux-artifact versions), for 1850-spin-up/transient CO2.
- Met `CO2air` variable available per site for the observed window.

## Validation-site met (10 sites) — NOT YET AVAILABLE LOCALLY

AU-ASM, AU-How, AU-Stp, AU-Tum, CH-Cha, US-FPe, US-GLE, US-Ha1, US-Me2, US-UMB.
Only distributed via ME-org ("Phase 1b Validation" experiment); requires Renato's ME-org
login (blocker for Phase 7; flagged in PLAN.md). Note: AU-How *flux* obs exist in the
stale artifact but protocol provides no validation-site flux — do not use.

## Other local copies (pre-artifact, kept for reference)

- `/net/sampo/data1/renatob/DK-Sor_1997-2014_FLUXNET2015_Met.nc` (= artifact copy)
- `/net/sampo/data1/renatob/DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc`
  (= Phase 1a-test file)
- `/net/sampo/data1/renatob/callmip_phase1a_dksor/` (Phase 1a working dir)
