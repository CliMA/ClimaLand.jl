# Transient runs from 1850 (TRENDY style)

`transient_from_1850.jl` runs the global `LandModel` with prognostic LAI
(`ZhouOptimalLAIModel`) and carbon pools (`PrognosticCarbonModel`) from a spun-up 1850
state through 2023, in the scenarios of the TRENDY protocol that apply to a model without
land-use change:

| scenario | CO2 | climate |
|---|---|---|
| S0 (control) | fixed at 1850 | 1901–1920 recycled |
| S1 | observed | 1901–1920 recycled |
| S2 | observed | 1901–1920 recycled before 1901, then observed |

The land sink is S_LAND = NBP(S2) − NBP(S0), with S0 removing the drift of the initial
state. S1 − S0 is the CO2 effect and S2 − S1 the climate effect.

## Inputs

- The ClimaArtifacts artifact `climaland_1850_spunup_ic`, in `ARTIFACT_DIR` until it is
  registered: `climaland_1850_spunup_ic.nc`, the full state at the end of an 1850 spin-up
  with CRU-JRA forcing (its README describes the spin-up and the drift left in the
  state), and `co2_annual_1700_2025.csv`. The state is read node by node, so the domain
  must be the 1° box domain with 15 layers. `IC` can instead be a checkpoint of an
  earlier run, to branch from it.
- Forcing: CRU-JRA (`FORCING=crujra`), with its missing values filled and its shortwave
  scaled with the zenith angle (`prescribed_forcing_crujra` options), and the drivers
  updated every step. `FORCING=era5` runs from 1979 on ERA5 instead.
- CO2: the annual global means of the artifact (Law Dome spline to 1958, Mauna Loa
  1959–1979, NOAA global means from 1980), placed at mid-year and interpolated linearly.

## Running

Each run writes, in `RUN_DIR`, a state per year (`states/state_<year>.hdf5`), monthly
diagnostics, and `annual_budget.csv`: per year, the global GPP, autotrophic and
heterotrophic respiration, litter input and NBP (Pg C), the carbon stocks (Pg C), and the
closure of the budget. A rerun in the same `RUN_DIR` resumes after the last saved year,
and a job stops before `WALLTIME_HOURS`, so runs can be chained.

```
SCENARIO=S2 ARTIFACT_DIR=<artifact folder> RUN_DIR=S2 WALLTIME_HOURS=6 \
    julia --project=.buildkite experiments/integrated/transient_from_1850/transient_from_1850.jl
```

S1 and S2 are identical before 1901, so S1 can start from S2's state at the end of 1900
(`SCENARIO=S1 START_YEAR=1901 IC=S2/states/state_1900.hdf5`), and an ERA5 variant from
S2's state at the end of 1978 (`FORCING=era5 START_YEAR=1979 IC=S2/states/state_1978.hdf5`);
both still read the CO2 series from `ARTIFACT_DIR`.
On one A100 a simulated year takes about 17 minutes, so a full run takes about 50 GPU
hours. CRU-JRA ends at 2023-12-31 18:00, so the year 2023 of a run on CRU-JRA has 364.75
days.
