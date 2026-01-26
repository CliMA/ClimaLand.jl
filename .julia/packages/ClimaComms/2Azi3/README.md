# ClimaComms.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://CliMA.github.io/ClimaComms.jl/dev)
[![GitHub Actions status](https://github.com/CliMA/ClimaComms.jl/actions/workflows/OS-Tests.yml/badge.svg)](https://github.com/CliMA/ClimaComms.jl/actions/workflows/OS-Tests.yml)
[![Buildkite status](https://badge.buildkite.com/e3cbade62b514474b9f6abca474d58e760b9cb7a2545e46ad0.svg?branch=main)](https://buildkite.com/clima/climacomms-ci/builds?branch=main)

`ClimaComms.jl` is a small package to work with diverse computing devices and
environments (e.g., CPUs/GPUs, single-process/MPI). `ClimaComms.jl` provides
objects to represent such devices and environments and functions to interact
with them in a unified way. `ClimaComms.jl` is used extensively through the
`CliMA` ecosystem to control where and how simulations are run (e.g., on one
CPU, or on several GPUs using with MPI).

`ClimaComms.jl` supports the following `Device`s:
- `CPUSingleThreaded`
- `CPUMultiThreaded` (not actively used)
- `CUDADevice`
and `Context`es (i.e., environments for distributed computing):
- `SingletonCommsContext`
- `MPICommsContext`

Refer to the [documentation](https://CliMA.github.io/ClimaComms.jl/dev) for
more details.

