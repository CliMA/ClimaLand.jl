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
2. Run the calibration script on the cluster:
    - Execute `bash experiments/calibration/run_calibration.sh [OUTPUT_DIR]` in
      the terminal. `OUTPUT_DIR` is optional; on Derecho you typically want to
      pass a scratch path (e.g. `/glade/derecho/scratch/$USER/calibration_gpp`).
    - You can start this script from the login node. Slurm/PBS will handle job
      submission for the forward models.
    - You may want to use `tmux` to keep a persistent session on the cluster.

Before invoking the script you still need to `module load climacommon/<version>`
to put Julia on the orchestrator's PATH. The version forward-model PBS jobs
load is pinned separately in `run_calibration.jl` via the `modules` keyword on
the backend constructor — edit that line directly to change it.

For example, launching er.jl on Derecho.

Note Derecho has multiple login nodes (`derecho1`–`derecho7`); tmux sessions are
per-node, so SSH to a specific one (e.g. `ssh derecho7.hpc.ucar.edu`) so you can
reattach later from the same node.

Go into tmux (`-s` sets the session name; `-t` is for target-session groups and
won't do what you want here):

```
tmux new -s calibration
```

(later on, attach via `tmux attach -t calibration`, and detach via `ctrl+b, release, press d`)
Start the calibration *inside* the tmux session:

```
module load climacommon/2025_02_25
export CALIBRATION_CONFIG=er.jl
bash experiments/calibration/run_calibration.sh /glade/derecho/scratch/$USER/calibration_er
```
