# Running on GPU or with MPI

ClimaLand.jl is designed to run on either CPU or GPU, and to be compatible with MPI.
This section will explain how to select the architecture you want to run your simulation on.

The CliMA ecosystem uses [ClimaComms.jl](https://github.com/CliMA/ClimaComms.jl)
to control which architecture to run on, and related internals are handled within that package.
This isolation makes specifying the architecture to use relatively simple.

ClimaComms.jl introduces two concepts that we utilize to define the architecture
to run a simulation on: the `device` and the `context`.

The `device` specifies whether we want to run on CPU or GPU.
The `context` contains information about the device, and also defines whether we
want to run on an individual processor or distributedly with MPI.
Since the device is contained within the context's data structure, when we construct
these objects in ClimaLand we simply refer to them as the context.

Below, we have descriptions for how to specify the context when trying to run a
simulation in different setups. There are two ways to do this: internally within a
script/REPL, or in the terminal via environment variables.

For more detailed information, please see the ClimaComms.jl [documentation](https://clima.github.io/ClimaComms.jl/stable/).

## Choosing the architecture within a script/REPL

To select the architecture to use within a script or the REPL, we will
explicitly set the device and context to use. The context is then
automatically used for the rest of the session.

### Running on CPU (without MPI)

```
context = ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())
```

### Running on GPU (without MPI)

```
context = ClimaComms.SingletonCommsContext(ClimaComms.CUDADevice())
```

### Running with MPI

On CPU:
```
context = ClimaComms.MPICommsContext(ClimaComms.CPUSingleThreaded())
```

On GPU:
```
context = ClimaComms.MPICommsContext(ClimaComms.CUDADevice())
```

### Default context and device

We can also set the context using the following default constructor. This is useful for runs where
we want to run on CPU and don't have MPI available, perhaps for example in simple local runs.
```
context = ClimaComms.context()
```

This constructor will by default run on CPU without MPI, unless the relevant environment variables are set.
See the section below for more details about controlling the architecture with environment variables.

## Choosing the architecture via environment variables

ClimaComms checks two environment variables to determine which device and context to use: `CLIMACOMMS_DEVICE` and
`CLIMACOMMS_CONTEXT`. Setting these environment variables provides an alternative to defining the device and context
explicitly, as was done above.

### Running on CPU (without MPI)

```
ENV["CLIMACOMMS_DEVICE"] = "CPU"
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
```

### Running on GPU (without MPI)

```
ENV["CLIMACOMMS_DEVICE"] = "CUDA"
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
```

### Running with MPI

On CPU:
```
ENV["CLIMACOMMS_DEVICE"] = "CPU"
ENV["CLIMACOMMS_CONTEXT"] = "MPI"
```

On GPU:
```
ENV["CLIMACOMMS_DEVICE"] = "CUDA"
ENV["CLIMACOMMS_CONTEXT"] = "MPI"
```
