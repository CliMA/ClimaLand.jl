# 1850 spin-up of the LandModel with carbon pools (side branch, not for merging)

This branch records the code and the steps that produced the spun-up 1850 state used as
the initial condition of TRENDY-style transient runs (S0, S1, S2 over 1850–2023) and
packaged as a ClimaArtifacts initial condition. It is the provenance of that state; it is
not meant to be merged (AGENTS.md: spin-up scripts go with the artifact's recipe or on a
side branch).

## Code

- ClimaLand main at b7829b5fc, plus the change to `src/shared_utilities/drivers.jl` on this
  branch: two options of `prescribed_forcing_crujra`, `fill_missing` (missing CRU-JRA
  values are taken from the nearest valid nodes) and `zenith_scaled_shortwave` (the
  6-hourly shortwave is redistributed within each window with the cosine of the solar
  zenith angle, keeping its mean). The runs use both.
- The scripts in this folder, as run (only the default path of `CO2_CSV` in
  `run_trendy.jl` changed, to this folder). `spinup_tools.jl` is the later version used by
  the second offline solve (A5); the runs used an earlier one that differs only in the
  interpolation of `solve_soc`, which these runs do not apply (`APPLY_SOC=0`). The first
  solve (A3) used the earlier interpolation (linear in log SOC; within about 1%).
- Julia 1.13.3 with the `.buildkite` environment, Float64, one A100 per run
  (~17–18 minutes per simulated year).

## Configuration

- The global LandModel on the 1° longitude–latitude box (`global_box_domain`,
  `mask_threshold = 0.99`), 15 soil layers to 15 m, with prognostic LAI
  (`ZhouOptimalLAIModel`) and carbon pools (`PrognosticCarbonModel`), so litter from the
  pools feeds soil organic carbon (SOC).
- Forcing: CRU-JRA, with the 1901–1920 cycle recycled. Nominal year n is forced by
  forcing year 1901 + mod(n − 1850, 20).
- CO2 fixed at its 1850 annual mean, 286.8 ppm (`co2_annual_1700_2025.csv`, built by
  `co2_series.py` from the Law Dome spline, Mauna Loa and NOAA global means).
- Parameters: main's defaults, except the DAMM values of `respiration_override.toml`
  from C3 on: `soilCO2_reference_rate` 6.7e-7 kg C m⁻³ s⁻¹ and `michaelis_constant`
  0.70 kg C m⁻³, the values proposed as defaults in CliMA/ClimaLand.jl#1968.

## Steps

Each run sets the environment variables in its row and calls
`julia --project=.buildkite experiments/trendy_spinup_1850/run_trendy.jl` (the header of
`run_trendy.jl` lists the options). All runs set `FORCING=crujra`, `RECYCLE=1`,
`CYCLE_ORIGIN=1850` and `CO2_FIXED_YEAR=1850`. With `SPINUP=pools_soc`, a run accumulates
the mean pool fluxes and, per layer, the litter and the microbial respiration on a grid of
SOC values, and at its end writes `spun_<year>.hdf5` with the vegetation pools at their
equilibrium; `APPLY_SOC=0` keeps the SOC.

| step | years | from | settings | result |
|---|---|---|---|---|
| C1 | 1810–1829 | default IC (`rosetta_spunup_ic`), `REPAIR_IC=1` | default DAMM, `SPINUP=pools_soc`, `APPLY_SOC=0` | `spun_1829`: pools at equilibrium (wood 176 Pg C) |
| C2 | 1830–1849 | C1 `spun_1829` | default DAMM, `SPINUP=pools_soc`, `APPLY_SOC=0`, `SOC_ALTS` | respiration curves (`accum_1849`), `spun_1849` |
| A3 | offline | C2 `accum_1849`, `spun_1849` | `damm_choice.jl … 0 20.06 16.83` (targets from `targets.jl`) | V_ref ×20, kM ×16.8 match SoilGrids' top metre; `state_damm` (SOC 0–1 m 2,104, 0–2 m 2,869 Pg C) and the override TOML |
| C3 | 1850–1869 | A3 `state_damm` | override | NBP 3.09 Pg C/yr |
| C4 | 1870–1889 | C3 `state_1869` | override, `SPINUP=pools_soc`, `APPLY_SOC=0` | NBP 2.71 Pg C/yr; curves `accum_1889`, `spun_1889` |
| A5 | offline | C4 `accum_1889`, `spun_1889` | `damm_choice.jl … 0 1 1` (the override values) | `state_damm` (SOC 0–1 m 2,156 Pg C) |
| C5 | 1890–1909 | A5 `state_damm` | override, `SPINUP=pools_soc`, `APPLY_SOC=0` | **`spun_1909`, the 1850 state**; NBP 2.14 Pg C/yr |

`o2_soc_change.jl` computes the SOC drift by layer, band and column between two states.
`state_to_netcdf.jl` writes a state as the NetCDF file of the ClimaArtifacts initial
condition (every field of `Y`, read back and checked bit for bit), and
`soc_spinup_figure.jl` compares the SOC of a state (HDF5 or that NetCDF file) with
SoilGrids.

## The state at the end (C5, `spun_1909`)

Land north of 60°S, SoilGrids as the model reads it (on its layers and land points):

| Pg C | 0–1 m | SoilGrids | 0–2 m | SoilGrids | full column |
|---|---|---|---|---|---|
| global | 2,162 | 2,112 | 3,062 | 2,960 | 4,223 |
| north of 50°N | 597 | 724 | 665 | 1,008 | 695 |
| 30–50°N | 631 | 439 | 693 | 590 | 807 |
| 30°S–30°N | 819 | 815 | 1,576 | 1,136 | 2,588 |
| 60–30°S | 115 | 134 | 129 | 226 | 133 |

The state is not at a steady state everywhere. Over C5, NBP is 2.14 Pg C/yr with no
trend: SOC rises by 0.3 Pg C/yr in the top metre and by 1.76 Pg C/yr below 1 m, almost
all in tropical layers at 1–4 m, where root litter enters faster than respiration,
capped by its saturated rate and by O2, removes it. Column drift is below 10 g C m⁻² yr⁻¹
on 85% of the land (97% north of 50°N, 74% in the tropics). In transient runs, the S0
control removes the constant part. The vertical profile is shallower than SoilGrids'
(about 1,600 Pg C in the top 10 cm against about 500).
