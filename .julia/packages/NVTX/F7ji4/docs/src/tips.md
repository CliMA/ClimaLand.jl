# Nsight Systems tips

## MPI

Nsight Systems can profile MPI operations, and can be used in two ways:

### Profiler outside launcher

```
nsys profile <launcher> <program>
```

In this approach, `nsys profile` calls the MPI launcher (e.g. `mpiexec`), which in turn calls the script:
```
nsys profile --trace=nvtx,mpi --mpi-impl=openmpi mpiexec -n 2 julia --project mpi.jl
```
This will generate one report for the whole run, but only works when all MPI processes are on the one machine.

### Launcher outside profiler

```
<launcher> nsys profile <program>
```

Alternatively the launcher can be used to start multiple copies of the profiler, which will generate a report for each MPI process. You will need to specify a unique name for each file, which can be done by `--output=` with using `%q` to interpolate an environment variable. Most launchers will set the `MPI_COMM_WORLD` rank of each process as an environment variable (e.g. `PMI_RANK` or `OMPI_COMM_WORLD_RANK`).

```
mpiexec -n 2 nsys profile --trace=nvtx,mpi --mpi-impl=openmpi --output=report.%q{OMPI_COMM_WORLD_RANK} julia --project mpi.jl
```

This file can be opened as a "multi-report view".

## Reducing overhead

The profiler itself has some overhead. Typically it will make use of another thread, but if there insufficient hardware resources available, it may result in occasional pauses in the program being profiled. 

When using a HPC cluster scheduler such as Slurm, this can be reduced by allocating an additional CPU core per MPI process for the profiler (e.g. via the `--cpus-per-task` option in Slurm). To ensure that these are scheduled correctly, it is best to "bind" the CPU cores per task. If launching using `srun`, then use the `--cpu-bind=cores`; if launching using Open MPI `mpiexec`, use `--map-by node:PE=$cpus_per_task --bind-to core`.

## Querying profiles

If you need to query the profile data programmatically, the easiest option is to first convert to make use of the SQLite format, either by the `--export` argument, or using the `nsys export` command.

The format of the resulting file is documented in the [SQLite Schema Reference](https://docs.nvidia.com/nsight-systems/UserGuide/index.html#exporter-sqlite-schema).