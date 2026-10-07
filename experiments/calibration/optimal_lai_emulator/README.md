# Offline emulator of the optimal-LAI model

An offline version of `ZhouOptimalLAIModel`, used to screen its structure and calibrate
its parameters against MODIS LAI in minutes rather than with global runs. It calls the
model's own functions (`compute_A0_and_χ`, `compute_L_steady_target`, `optimal_chi`,
`climate_tree_share`, ...) on hourly low-resolution ERA5 for 2008 (recycled) at the land
points of its 8° grid, at equilibrium: the yearly totals are the annual means of their
rates, and LAI runs hourly over the year after a 60-day warm-up. It is scored against
MODIS LAI (Yuan et al. 2011, 2008) averaged over the 1° land cells around each point,
over natural vegetation (CLM5 crop fraction ≤ 0.5): 226 points, cos-latitude weighted.
The CLM5 tree share is a target of the calibration, never an input.

- `data.jl`: forcing, observations and the per-point climate (potential GPP, χ, PET).
- `model.jl`: climate features, the configuration (`Config`), `emulate` and scores.
- `snow.jl`: the degree-day snow store and the moist growing season.
- `calibrate.jl`: Nelder–Mead.
- `run.jl`: loads everything, scores the defaults of `toml/default_parameters.toml`,
  and holds the calibration (weak priors, global-bias penalty, checkerboard
  cross-validation).

```
CLM_SURFDATA=/path/to/surfdata_0.9x1.25_16pfts__CMIP6_simyr2000_c170616.nc \
    julia --project=.buildkite -i experiments/calibration/optimal_lai_emulator/run.jl
```

The surface data is the file that `artifacts/clm_tree_share/create_clm_tree_share.jl`
downloads. With the defaults, the emulator gives an annual-mean LAI RMSE of 0.73, a
monthly RMSE of 0.85 and a bias of −0.01; with the same structure, global long runs have
had RMSEs about 0.83 times those of the emulator.
