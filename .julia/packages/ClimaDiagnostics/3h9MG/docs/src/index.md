# `ClimaDiagnostics.jl`

`ClimaDiagnostics.jl` is a module that adds diagnostics to your
[`CliMA`](https://github.com/CliMA) simulations. Diagnostics are variables that
are computed from your state when specific conditions are met (typically at set
intervals of time) throughout the run and typically saved to disk.

`ClimaDiagnostics.jl` provides the infrastructure to schedule, compute, reduce,
and output such variables.

If you are a user of a package that is already using `ClimaDiagnostics.jl`, you
can jump to the User guide page.

If you are a developer interested in adding support for `ClimaDiagnostics.jl` in
your package or learning about the internal design of this package, please read
the Developer guide page.

## Features

- Define diagnostics as function of the integrator state and the cache;
- Accumulate diagnostics over period of times with associative binary temporal
  reductions (eg, averages);
- Work with calendar dates (eg, monthly averages);
- Allow users to define arbitrary new diagnostics;
- Trigger diagnostics on arbitrary conditions;
- Save output to HDF5 or NetCDF files, or a dictionary in memory;
- Work with lazy expressions (such as the ones produced by
  [LazyBroadcast.jl](https://github.com/CliMA/LazyBroadcast.jl)).

