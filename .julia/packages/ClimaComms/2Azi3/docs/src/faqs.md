# Frequently Asked Questions

## How do I run my simulation on a GPU?

Set the environment variable `CLIMACOMMS_DEVICE` to `CUDA`. This can be
accomplished in your Julia script with (at the top)
```julia
ENV["CLIMACOMMS_DEVICE"] = "CUDA"
```
or calling
```julia
export CLIMACOMMS_DEVICE="CUDA"
```
in your shell (outside of Julia, no spaces).

## My simulation does not start and crashes with a `MPI` error. I don't want to run with `MPI`. What should I do?

`ClimaComms` tries to be smart and select the best configuration for your run.
Sometimes, it fails with an error message like the following.
```
cmd=init pmi_version=2 pmi_subversion=0
--------------------------------------------------------------------------
PMI2_Init failed to intialize.  Return code: 14
--------------------------------------------------------------------------
--------------------------------------------------------------------------
The application appears to have been direct launched using "srun",
but OMPI was not built with SLURM's PMI support and therefore cannot
execute. There are several options for building PMI support under
SLURM, depending upon the SLURM version you are using:

  version 16.05 or later: you can use SLURM's PMIx support. This
  requires that you configure and build SLURM --with-pmix.

  Versions earlier than 16.05: you must use either SLURM's PMI-1 or
  PMI-2 support. SLURM builds PMI-1 by default, or you can manually
  install PMI-2. You must then build Open MPI using --with-pmi pointing
  to the SLURM PMI library location.

Please configure as appropriate and try again.
--------------------------------------------------------------------------
*** An error occurred in MPI_Init_thread
*** on a NULL communicator
*** MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort,
***    and potentially your MPI job)
```

In this case, you can force `ClimaComms` to ignore `MPI`
with
```julia
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
```
at the top of your Julia script or by calling
```julia
export CLIMACOMMS_CONTEXT="SINGLETON"
```
in your shell (outside of Julia, no spaces).

## My code is saying something about `ClimaComms.@import_required_backends`, what does that mean?

When you are using the environment variables to control the execution of your
script, `ClimaComms` can detect that some important packages are not loaded. For
example, `ClimaComms` will emit an error if you set `CLIMACOMMS_DEVICE="CUDA"`
but do not import `CUDA.jl` in your code.

`ClimaComms` provides a macro, [`ClimaComms.@import_required_backends`](@ref),
that you can add at the top of your scripts to automatically load the required
packages when needed. Note, the packages have to be in your Julia environment,
so you might install packages like ` MPI.jl` and `CUDA.jl`.

## How can I see the MPI state and verify that MPI is set up correctly?

To inspect the current MPI state, use the `summary` function: `print(summary(context))`

The output varies depending on your communication context type:

- `SingletonCommsContext`: Displays basic context information and device type
- `MPICommsContext`: Shows detailed information including each node's rank.

When using GPU acceleration with `CUDADevice`, the summary additionally includes the device type and UUID.

To test that MPI and CUDA are set up correctly, see [this guide](https://github.com/CliMA/slurm-buildkite?tab=readme-ov-file#testing-cuda-and-mpi-modules).
