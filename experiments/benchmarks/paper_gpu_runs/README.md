# Generating Weak and Strong Scaling Data

This folder contains scripts that output data that can be used to produce
weak and strong scaling plots.

## Weak Scaling

To generate weak scaling data, run the `weak_scaling_gpu_runs.jl` script while using
`ClimaComms.CUDADevice`

For example:

```bash
export CLIMACOMMS_DEVICE="CUDA"
julia --project=.buildkite experiments/benchmarks/paper_gpu_runs/weak_scaling_gpu_runs.jl
```

This will print the scaling results for the Snow, Bucket, Richards, Canopy, Soil, Soil-Canopy,
and Soil-Canopy-Snow models to stdout. The results are grouped by model, with a tuple
on each new line of the form (num_columns, time_per_step). The results for each model are
wrapped in additional lines that allow users to copy and paste into a PGF plot.

Sample Output:
```bash
\addplot[
    color=white,
    mark=pentagon*,
    opacity = 0
    ]
    coordinates {
      (1536,0.008838768)
      (3456,0.009402341)
      (13824,0.011664405)
      (55296,0.02421938)
      (203136,0.063300975)
      (777600,0.21756136)
      (3110400,0.79437906)
    };
\end{loglogaxis}
```

## Strong Scaling

The process of generating the strong scaling plots varies depending on which HPC system is
used. A bash script, `scale.sh`, is provided as an example. The comments within the script
describe its usage. Results are printed to stdout as tuples in the form
(effective_resolution, number_of_gpus, time_per_step)

Sample Output:

```bash
(0.5 , 4 , 0.014412601)
(0.25 , 4 , 0.054200597)
(0.125 , 4 , 0.20701422)
```
