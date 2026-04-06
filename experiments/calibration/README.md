# Land Calibration Pipeline

The land calibration pipeline features calibration of the parameters of the full
land model or bucket model. The observations are seasonal averages of `lhf`,
`shf`, `lwu`, and `swu` spanning twenty years (1979-3 to 2024-8)  with a
latitude-weighted scalar covariance matrix.

# How do I run a calibration?

1. Pick (or create) a configuration file under
   `experiments/calibration/configs/`. `run_calibration.jl` loads
   `configs/energy_fluxes.jl` by default. To use a different config, set the
   `CALIBRATION_CONFIG` environment variable before running the script, e.g.
   `CALIBRATION_CONFIG=gpp.jl bash experiments/calibration/run_calibration.sh`.
   When `TEST_CALIBRATION` is set, the single-parameter `configs/test.jl` is
   loaded instead.
2. Load the `climacommon` version appropriate for your cluster (for example,
   `module load climacommon/2025_02_25` on Derecho). The loaded version is
   forwarded to the compute jobs automatically.
3. Run the calibration script on the cluster:
    - Execute `bash experiments/calibration/run_calibration.sh [OUTPUT_DIR]` in
      the terminal. `OUTPUT_DIR` is optional; on Derecho you typically want to
      pass a scratch path (e.g. `/glade/derecho/scratch/$USER/calibration_gpp`).
    - You can start this script from the login node. Slurm/PBS will handle job
      submission for the forward models.
    - You may want to use `tmux` to keep a persistent session on the cluster.

For example, launching er.jl on Derecho
Go into tmux:

```
tmux new -t calibration
```

(later on, attach via `tmux attach -t calibration`, and detach via `ctrl+b, release, press d`)
Start the calibration:

```
module load climacommon/2025_02_25
export HDF5_USE_FILE_LOCKING=FALSE
export CALIBRATION_CONFIG=er.jl
bash experiments/calibration/run_calibration.sh /glade/derecho/scratch/$USER/calibration_er
```

# Debugging

ClimaCalibrate tries to keep `climacommon` up to date. If errors result from
Julia version mismatches, use the `climacommon` versions below. These versions
are compatible with ClimaCalibrate v0.2.2.

- `ClimaGPUBackend`: `climacommon/2026_02_18`
- `DerechoBackend`: `climacommon/2025_02_25`
- `CaltechHPCBackend`: `climacommon/2024_10_09`
- `GCPBackend`: `climacommon` is not supported
