# Land Calibration Pipeline

The land calibration pipeline features calibration of the parameters of the full
land model or bucket model. The observations are seasonal averages of `lhf`,
`shf`, `lwu`, and `swu` spanning twenty years (1979-3 to 2024-8)  with a
latitude-weighted scalar covariance matrix.

# How do I run a calibration?

1. In `run_calibration.jl`, update the calibration configuration.
    - On Derecho, you may want to set `output_dir` to a scratch directory
      (e.g. "/glade/derecho/scratch").
2. Run the calibration script on the cluster.
    - Execute `bash experiments/calibration/run_calibration.sh` in the terminal.
    - You can start this script from the login node. Slurm/PBS will handle job
      submission for the forward models.
    - You may want to use `tmux` to keep a persistent session on the cluster.

# Debugging

ClimaCalibrate tries to keep `climacommon` up to date. If errors result from
Julia version mismatches, use the `climacommon` versions below. These versions
are compatible with ClimaCalibrate v0.2.2.

- `ClimaGPUBackend`: `climacommon/2026_02_18`
- `DerechoBackend`: `climacommon/2025_02_25`
- `CaltechHPCBackend`: `climacommon/2024_10_09`
- `GCPBackend`: `climacommon` is not supported
