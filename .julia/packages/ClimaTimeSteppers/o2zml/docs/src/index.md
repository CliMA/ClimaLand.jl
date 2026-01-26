# ClimaTimeSteppers.jl

ClimaTimeSteppers.jl is a suite of ordinary differential equation (ODE) solvers for use as
time-stepping methods in a partial differential equation (PDE) solver, such as
[ClimateMachine.jl](https://github.com/CliMA/ClimateMachine.jl). They are specifically
written to support distributed and GPU computation, while minimising the memory footprint.

ClimaTimeSteppers.jl is built on top of
[DiffEqBase.jl](https://github.com/SciML/DiffEqBase.jl/), and aims to be compatible with
the [DifferentialEquations.jl ecosystem](https://diffeq.sciml.ai/latest/).

```@docs
ClimaTimeSteppers
```
