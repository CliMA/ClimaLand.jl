ClimaComms.jl Release Notes
========================

main
-------

v0.6.9
-------
- Added a device-agnostic API for querying available memory [PR 117](https://github.com/CliMA/ClimaComms.jl/pull/117).

v0.6.8
-------
- Extended `@threaded` to work with multiple iterators and lazy iterators (e.g., `enumerate`, `zip`, and `Iterators.partition`), and modified the `threaded` function to make it equivalent to `@threaded` [PR 115](https://github.com/CliMA/ClimaComms.jl/pull/115).

v0.6.7
-------
- Extended `@threaded` to work on GPU devices, with block sizes automatically determined by the CUDA occupancy API, and added the ability to control thread coarsening across all devices [PR 111](https://github.com/CliMA/ClimaComms.jl/pull/111).

v0.6.6
-------
- Replaced `MPIFileLogger` with `FileLogger` and added an `OnlyRootLogger` logger that silences non-root processes [PR 104](https://github.com/CliMA/ClimaComms.jl/pull/104).

v0.6.5
-------

- Support for adapt was added, so that users can convert between CPU and GPU
  devices, and contexts containing cpu and gpu devices [PR 103](https://github.com/CliMA/ClimaComms.jl/pull/103).

- New MPI logging tools were added, `MPIFileLogger` and `MPILogger` [PR 100](https://github.com/CliMA/ClimaComms.jl/pull/100).

v0.6.4
-------

- Add device-flexible `@assert` was added [PR 86](https://github.com/CliMA/ClimaComms.jl/pull/86).

v0.6.3
-------

- Bugfix: `cuda_sync` was missing in the extension and, as a result, `ClimaComms.@cuda_sync` was not actually synchronizing. We've also removed the abstract fallback, so that we instead method-error if we pass a cuda device when the cuda extension does not exist.

v0.6.2
-------

- We added a device-agnostic `allowscalar(f, ::AbstractDevice, args...; kwargs...)` to further assist in making CUDA an extension.

v0.6.1
-------

- Macros have been refactored to hopefully fix some code loading issues.

v0.6.0
-------

- ![][badge-💥breaking] `ClimaComms` does no longer try to guess the correct
  compute device: the default is now CPU. To control which device to use,
  use the `CLIMACOMMS_DEVICE` environment variable.
- ![][badge-💥breaking] `CUDA` and `MPI` are now extensions in `ClimaComms`. To
  use `CUDA`/`MPI`, `CUDA.jl`/`MPI.jl` have to be loaded. A convenience macro
  `ClimaComms.@import_required_backends` checks what device/context could be
  used and conditionally loads `CUDA.jl`/`MPI.jl`. It is recommended to change
  ```julia
  import ClimaComms
  ```
  to 
  ```julia
  import ClimaComms
  ClimaComms.@import_required_backends
  ```
  This has to be done before calling `ClimaComms.context()`.

[badge-💥breaking]: https://img.shields.io/badge/💥BREAKING-red.svg
