# CalLMIP Phase 1b — Data Manifest and Verification

Verified 2026-08-18 (Claude, session with Renato). Update this file whenever data
locations, hashes, or upstream contents change.

## Artifacts

| Artifact | git-tree-sha1 | Resolved path (Overrides.toml) | Status |
|---|---|---|---|
| `callmip_phase1_forcing` | `a0c004071c73d047095a01faedf017571e5237ef` | `/net/sampo/data1/ClimaArtifacts/artifacts/callmip_phase1_forcing` (per Renato: use sampo; hash-identical copy exists at `/home/renatob/data/callmip_data/callmip_phase1_forcing`) | ✅ hash verified 2026-08-18; 21 met files, correct site set |
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
- Flux obs files carry the same year labels as met. All 21 refreshed flux files verified
  2026-08-18: variables `{NEE,Qle,Qh}_daily` + `{…}_uc_daily` + time/latitude/longitude;
  NEE units `gC/m2/d`; daily counts consistent with year spans (leap days included).
- Met files verified to carry: Tair, Wind, Qair, Psurf, SWdown, LWdown, CO2air, Precip,
  VPD, RH, LAI(±alternative), reference_height, latitude, longitude, elevation.
- **US-MMS is HOURLY** (140,256 = 5,844 d × 24); the other 20 sites are half-hourly.
  Handled by linear midpoint insertion at load time (`resample_to_30min` in
  `forcing_phase1b.jl`) — identical to what the TVI's linear interpolation of the
  hourly nodes would produce, so physically a no-op. Found 2026-08-18 when padding
  hit a missing 30-min key; the n-half-hours column below is the raw file count.
- **LAI product↔variable mapping varies by site** (`source` attribute): MODIS sits
  under `LAI` at FR-Pue, IT-Lav, IT-MBo, IT-Noe, RU-Fyo and all US sites; under
  `LAI_alternative` at the 8 northern/central European + Canadian sites. D7 (MODIS
  everywhere, Renato 2026-08-18) is therefore implemented by attribute
  (`modis_lai_var`), never by variable name.
- **US-SRG tower height is 3.25 m** (not 3.2 — earlier table entry came from a
  rounded print; caught by the `load_site` cross-check).
  Time axis: `seconds since <start>-01-01`, **local standard time**, no UTC-offset
  attribute (checked US-SRG) → UTC offsets must come from our site table (Phase 2).
- Flux files carry: NEE, Qle, Qh + `<var>_uc` + missing-% per day. Daily, gap-days = NaN.

### Flux valid-day coverage (% of file days with non-NaN obs; verified 2026-08-18)

Protocol filters (>25% gapfilled → NaN; no uncertainty → NaN) make NEE sparse at some
sites. Totals across 21 sites: **NEE 16,625 · Qle 60,485 · Qh 67,352 valid daily obs**.

| Site | NEE% | Qle% | Qh% | | Site | NEE% | Qle% | Qh% |
|---|---|---|---|---|---|---|---|---|
| CA-Qfo | 38.5 | 86.1 | 84.2 | | NL-Loo | 17.5 | 61.1 | 80.4 |
| CH-Dav | **1.4** | 43.3 | 56.0 | | RU-Fyo | 21.9 | 74.9 | 83.2 |
| DE-Gri | 26.8 | 76.0 | 87.3 | | US-MMS | 17.3 | 80.6 | 81.8 |
| DE-Hai | **1.5** | 22.6 | 64.5 | | US-NR1 | 24.0 | 87.3 | 85.7 |
| DE-Tha | 26.6 | 92.1 | 89.8 | | US-SRG | 31.5 | 89.0 | 96.9 |
| DK-Sor | 46.4 | 47.6 | 49.7 | | US-SRM | 15.9 | 87.9 | 96.6 |
| FI-Hyy | 18.7 | 25.4 | 26.3 | | US-Ton | **0.3** | 59.0 | 66.5 |
| FR-Pue | 18.4 | 40.2 | 26.7 | | US-Var | **1.5** | 79.6 | 89.3 |
| IT-Lav | **1.8** | 20.6 | 32.7 | | US-Whs | 19.9 | 91.0 | 97.0 |
| IT-MBo | 9.9 | 40.7 | 49.4 | | US-Wkg | 8.3 | 84.9 | 93.8 |
| IT-Noe | 7.2 | 34.0 | 37.9 | | | | | |

**Implications for calibration (Phase 4):** per-site-per-month valid-day masks are
mandatory (Phase 1a `MIN_VALID_DAYS` machinery); NEE contributes almost nothing at
US-Ton/CH-Dav/DE-Hai/IT-Lav/US-Var — misfit at those sites is energy-flux-dominated.
Same data limitation applies to every CalLMIP group.

## Calibration site table (extracted from met files 2026-08-18; UTC offsets VERIFIED)

UTC offsets verified two ways from the SWdown diurnal cycle in each met file:
(a) top-half centroid of the mean summer diurnal cycle (biased morning-ward at
convective sites: US-NR1 −7.79, US-MMS −5.56), and (b) **sunrise/sunset midpoint**
(unbiased; decisive). Method (b) matches the FLUXNET political offsets at all 21 sites
within ±0.15 h — the table below is safe to hard-code. Naive `round(lon/15)` would be
WRONG at 4 sites: FR-Pue, NL-Loo (+1 not 0), RU-Fyo (+3 not +2), US-MMS (−5 not −6).

| Site | lat | lon | UTC offset | z_ref (m) | elev (m) | met start | n half-hours |
|---|---|---|---|---|---|---|---|
| CA-Qfo | 49.6925 | −74.3421 | −5 | 24.0 | 382 | 2004-01-01 | 122736 |
| CH-Dav | 46.8153 | 9.8559 | +1 | 35.0 | 1639 | 1997-01-01 | 315552 |
| DE-Gri | 50.9495 | 13.5125 | +1 | 3.0 | 385 | 2004-01-01 | 192864 |
| DE-Hai | 51.0792 | 10.4530 | +1 | 43.5 | 430 | 2000-01-01 | 227952 |
| DE-Tha | 50.9636 | 13.5669 | +1 | 42.0 | 380 | 1998-01-01 | 298032 |
| DK-Sor | 55.4859 | 11.6446 | +1 | 57.0 | 40 | 1997-01-01 | 315552 |
| FI-Hyy | 61.8475 | 24.2950 | +2 | 23.0 | 181 | 1996-01-01 | 333120 |
| FR-Pue | 43.7414 | 3.5958 | +1 | 11.0 | 270 | 2000-01-01 | 262992 |
| IT-Lav | 45.9562 | 11.2813 | +1 | 33.0 | 1353 | 2005-01-01 | 175296 |
| IT-MBo | 46.0147 | 11.0458 | +1 | 2.5 | 1550 | 2003-01-01 | 175344 |
| IT-Noe | 40.6062 | 8.1512 | +1 | 3.0 | 25 | 2004-01-01 | 192864 |
| NL-Loo | 52.1666 | 5.7436 | +1 | 27.0 | 25 | 1997-01-01 | 298032 |
| RU-Fyo | 56.4615 | 32.9221 | +3 | 48.0 | 265 | 2003-01-01 | 210384 |
| US-MMS | 39.3232 | −86.4131 | −5 | 48.0 | 275 | 1999-01-01 | 140256 |
| US-NR1 | 40.0329 | −105.5464 | −7 | 26.0 | 3050 | 1999-01-01 | 280512 |
| US-SRG | 31.7894 | −110.8277 | −7 | 3.25 | 1291 | 2009-01-01 | 105168 |
| US-SRM | 31.8214 | −110.8660 | −7 | 6.4 | 1120 | 2004-01-01 | 192864 |
| US-Ton | 38.4316 | −120.9660 | −8 | 23.0 | 177 | 2001-01-01 | 245424 |
| US-Var | 38.4133 | −120.9507 | −8 | 3.0 | 129 | 2001-01-01 | 245424 |
| US-Whs | 31.7438 | −110.0522 | −7 | 4.0 | 1370 | 2008-01-01 | 122736 |
| US-Wkg | 31.7365 | −109.9419 | −7 | 6.4 | 1531 | 2005-01-01 | 175296 |

Convention: UTC = local − offset (matches `hour_offset_to_period` usage in
`callmip_forcing.jl` / `read_callmip_met`).

## Non-site forcing

- TRENDY CO2: `Data/Non-site-specific_forcing/CO2_1700_2024_TRENDYv2025.txt`
  (in both flux-artifact versions), for 1850-spin-up/transient CO2.
- Met `CO2air` variable available per site for the observed window.

## Local copy note (2026-08-18)

`/home/renatob/data/callmip_data/` holds byte-identical copies of the two PUBLISHED
artifacts (forcing tree hash re-verified = `a0c0040…`; forcing override now points
here). ⚠️ Its `callmip_phase1/Data/Phase1b` is the STALE flux set (AU-How present,
US-Var missing) — do not use; the refreshed flux artifact `9b4bd0d8…` remains the
one wired in.

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
