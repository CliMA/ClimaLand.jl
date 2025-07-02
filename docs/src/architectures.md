# Running on GPU or with MPI

ClimaLand.jl is designed to run on either CPU or GPU, and to be compatible with MPI.
This section will explain how to select the architecture you want to run your simulation on.

The CliMA ecosystem uses [ClimaComms.jl](https://github.com/CliMA/ClimaComms.jl)
to control which architecture to run on, and related internals are handled within that package.
This isolation makes specifying the architecture to run on relatively simple.

ClimaComms.jl introduces two concepts that we utilize to define the architecture
to run a simulation on: the `device` and the `context`.

The `device` specifies whether we want to run on CPU or GPU.
The `context` contains information about the device, and also defines whether we
want to run on an individual processor or distributedly with MPI.
Since the device is contained within the context's data structure, when we construct
these objects in ClimaLand we simply refer to them as the context.

Below, we have descriptions for how to specify the context when trying to run a
simulation in different setups.

For more detailed information, please see the ClimaComms.jl [documentation](https://clima.github.io/ClimaComms.jl/dev/).

## Running on CPU (without MPI)

```
context = ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())
```

## Running on GPU (without MPI)

```
context = ClimaComms.SingletonCommsContext(ClimaComms.CUDADevice())
```

## Running with MPI

On CPU:
```
context = ClimaComms.MPICommsContext(ClimaComms.CPUMultiThreaded())
```

On GPU:
```
context = ClimaComms.MPICommsContext(ClimaComms.CUDADevice())
```

## Default context and device

We can also set the context using the following default constructor. This is useful for runs where
we want to run on CPU and don't have MPI available, for example simple local runs.
```
context = ClimaComms.context()
```

This constructor will by default run on CPU without MPI, unless the relevant environment variables are set.
ClimaComms checks two environment variables to determine which device and context to use: `CLIMACOMMS_DEVICE` and
`CLIMACOMMS_CONTEXT`. Setting these environment variables provides an alternative to defining the device and context
explicitly, as was done above.

### CLIMACOMMS_DEVICE

To run on CPU (note that this is the default behavior, so this ENV variable doesn't need to be explicitly set):
```
ENV["CLIMACOMMS_DEVICE"] = "CPU"
```

To run on GPU:
```
ENV["CLIMACOMMS_DEVICE"] = "CUDA"
```

### CLIMACOMMS_CONTEXT

To run without MPI (note that this is the default behavior, so this ENV variable doesn't need to be explicitly set):
```
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
```

To run with MPI:
```
ENV["CLIMACOMMS_CONTEXT"] = "MPI"
```
