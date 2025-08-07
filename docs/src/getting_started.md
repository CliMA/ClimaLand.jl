# Getting Started

## Installation of Julia and ClimaLand

ClimaLand is provided as a Julia package, so it requires having Julia installed. Information about Julia packages is available on the [Julia website](https://julialang.org/packages/).

First, download and install Julia by following the instructions at [https://julialang.org/downloads/](https://julialang.org/downloads/).
Then, you can launch a [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) by running `julia --project` and install the
ClimaLand.jl package by running the following within the REPL:

```julia
using Pkg
Pkg.add(ClimaLand)
using ClimaLand
```

## Running a simple soil model

Now that we have Julia installed, let's set up and run a simple ClimaLand simulation!
You can follow along by copying the following code segments into your REPL.

First we import the necessary Julia packages:
```julia
import SciMLBase
import ClimaParams as CP
import ClimaTimeSteppers as CTS
using ClimaLand
using ClimaLand.Domains
using ClimaLand.Soil
import ClimaLand.Parameters as LP
```

Choose a floating point precision, and get the parameter set, which holds constants used across CliMA models.
```julia
FT = Float32
earth_param_set = LP.LandParameters(FT);
```

We will run this simulation on a column domain with 1 meter depth, at a random lat/lon location.
```julia
zmax = FT(0)
zmin = FT(-1.0)
longlat = FT.((34.1, -118.1))
domain = Domains.Column(; zlim = (zmin, zmax), nelements = 10, longlat);
```

For this simple case, set up prescribed analytical atmosphere and radiation forcing.
```julia
atmos, radiation = prescribed_analytic_forcing(FT);
```

Now, we can create the [`EnergyHydrology`](@ref ClimaLand.Soil.EnergyHydrology) model.
This constructor uses default parameters and parameterizations, but these can also be
overwritten, which we'll demonstrate in later tutorials.
```julia
soil = Soil.EnergyHydrology{FT}(
    domain,
    (; atmos, radiation),
    earth_param_set,
);
```

Now that the model is set up, we can initialize our variables.
```julia
Y, p, coords = initialize(soil);
```

Set some initial conditions for the prognostic variables.
```julia
Y.soil.ϑ_l .= FT(0.24);
Y.soil.θ_i .= FT(0.0);
T = FT(290.15);
ρc_s = Soil.volumetric_heat_capacity.(
    Y.soil.ϑ_l,
    Y.soil.θ_i,
    soil.parameters.ρc_ds,
    soil.parameters.earth_param_set,
);
Y.soil.ρe_int .=
    Soil.volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        T,
        soil.parameters.earth_param_set,
    );
```

We choose the initial and final simulation times and timestep in seconds.
```julia
t0 = Float64(0)
tf = Float64(60 * 60 * 72);
dt = Float64(1000.0);
```

Set the cache values corresponding to the initial conditions of the state Y.
```julia
set_initial_cache! = make_set_initial_cache(soil);
set_initial_cache!(p, Y, t0);
```

Choose a timestepper and set up the ODE problem.
```julia
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 1,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

exp_tendency! = make_exp_tendency(soil);
imp_tendency! = make_imp_tendency(soil);
jacobian! = ClimaLand.make_jacobian(soil);

jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);
prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLand.dss!,
    ),
    Y,
    (t0, tf),
    p,
);
```

Now we can solve the problem, i.e. run the simulation.
```julia
sol = SciMLBase.solve(prob, ode_algo; dt);
```

That's it! For more complex examples and pretty plots,
please see the tutorial section of our docs.
