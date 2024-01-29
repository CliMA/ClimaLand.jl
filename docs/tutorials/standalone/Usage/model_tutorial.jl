# The `AbstractModel` framework allows users to define land
# component models (e.g. for snow, soil, vegetation, carbon...)
# which can be run in standalone mode, or as part of a land
# surface model with many components. In order to achieve this
# flexibility, we require a standard interface, which is what
# `AbstractModel`s provide. The interface is designed to work
# with an external package for the time-stepping of
# ODEs, `ClimaTimesteppers.jl`,
# with `ClimaCore.jl`, for the spatial discretization of
# PDEs, and with `ClimaLand.jl`, for designing and
# running multi-component land surface models. For a developer of a new
# land model component, using `AbstractModel`s as shown
# below is the first step towards building a model which
# can be run in standalone or with other components in an
# integrated land surface model.


# This tutorial introduces some of the functionality of the
# `AbstractModel` interface functions and types. We demonstrate
# how to use a `Model <: AbstractModel` structure to define
# a set of equations, and explain a few core methods
# which must be defined for your `Model` type in order to
# run a simulation.

# # General setup
# We assume you are solving a system of the form
# of a set of PDEs or ODEs. Additional algebraic equations
# can be accomodated as well, but only in addition to
# variables advanced using differential equations.

# Spatially discretized PDEs reduce to a system of ODEs,
# so we can assume an ODE system in what follows without a
# loss of generality. When using `AbstractModel`s,
# you should use
# `ClimaCore` to discretize your PDE, as applicable.


# Your model defines a system of equations of the following form:

# ``
# \frac{d \vec{Y}}{d t} = \vec{f}(\vec{Y}, \vec{x}, t; \mbox{params} \ldots)
# ``

# The variables that are stepped forward
# via a differential equation are referred to as prognostic variables, and
# are stored in `` \vec{Y} ``. Generically, we will speak of the
# functions `` \vec{f} `` as tendencies; these can be
# functions of the prognostic state, of space `` \vec{x} ``, and of time
# `` t ``, as well as of other parameters. Note that quantities such as
# boundary conditions, source terms, etc, will appear within these
# tendency functions

# # The cache ("auxiliary variables")

# There are typically quantities, which depend on the state vector `` \vec{Y} ``,
# location, time, and other parameters, which are expensive to compute,
# needed multiple times in the tendency computation, or require "a lot" of  memory to
# store (e.g., most variables in global runs). Allocating memory "on-the-fly" is
# typically time-consuming.  In these cases, it is far better
# to compute a quantity once and store in a variable where memory has been pre-allocated.
# The location where memory is allocated is called the model cache.

# Denoting the cache as `` \vec{p} ``,
# your equations may be rewritten as:

# ``
# \frac{d \vec{Y}}{d t} = \vec{f}(\vec{Y}, \vec{p}, \vec{x}, t; \mbox{params} \ldots)
# ``

# ``
# \vec{p}(\vec{x}, t) = \vec{g}(\vec{Y}, \vec{x}, t; \mbox{params} \ldots)
# ``

# The variables `` \vec{p} `` at the current timestep are functions of the prognostic
# state, space,
# time, and parameters. These variables are referred to as auxiliary variables
# (or cache variables). **Their main purpose is for storing the value of a quantity
# in a pre-allocated spot in memory, to avoid computing something expensive many times
# per time-step, or to avoid allocating memory each timestep. From a mathematical
# point of view, they represent intermediate quantities computed in each tendency.**
# A model purely consisting of algebraic equations, with no prognostic variables,
# is not supported (`` \vec{Y} `` cannot be zero dimensional).

# In order to define this set of equations, in a manner which is consistent
# with the `AbstractModel` interface (used by `ClimaLand.jl`) and
# time-stepping algorithms (`OrdinaryDiffEq.jl` for the present),
# the following must be provided.

# # The Model

# All `ClimaLand` component models are concrete instances of `AbstractModel`s. The reason
# for grouping them in such a way is because they all have shared required
# functionality, as we will see, and can make use of common default behavior.

# The model structure holds all
# of the information needed to create the full right hand side function, including
# parameters (which can be functions of space and time), boundary conditions, and
# physical equations.

# The purpose of our `AbstractModel` interface is that it allows you to run land
# component models in
# standalone mode and in an LSM mode without a change in interface. 

# As a simple demonstration of use, we'll build a model now which solves
# Richards Equation assuming a prescribed flux at the surface,
# and zero flux at the bottom of the column.

# Note that some model equations are stiff and require a very small timestep
# if stepped explicitly in time. Some model equations are amenable to "imex"
# timestepping, where some tendency functions are stepped implicitly,
# and some are stepped explicitly. Tagging a tendency function as "explicit"
# or "implicit" hardcodes something about the timestepping, and as such,
# conflates the idea of the model (which defines the equations) and the independent
# idea of a simulation (which solves the equations). However, we decided we did not need
# to support the flexibility of solving any set of equations in any way, as we
# are focused on land surface modeling in particular.
# In this example, we will tag the tendency as
# an explicitly time-stepped tendency. A follow-on tutorial will explain
# how to define an implicit tendency and tendency Jacobian. 

# Let's first import some needed packages.
import ClimaTimeSteppers as CTS
using SciMLBase
using Plots
using ClimaCore
using ClimaLand

# Import the functions we are extending for our model:
import ClimaLand:
    name,
    make_exp_tendency,
    make_compute_exp_tendency,
    make_update_aux,
    make_update_boundary_fluxes,
    prognostic_vars,
    prognostic_types,
    prognostic_domain_names,
    auxiliary_vars,
    auxiliary_types,
    auxiliary_domain_names


# The model should contain everything you need to create the
# tendency function. In this case, that is some parameters,
# the surface flux boundary value, the floating point precision,
# and the domain of the model
# (single column, global run, etc..).
struct RichardsTutorialModel{FT, D} <: AbstractModel{FT}
    "van Genuchten model parameters"
    vGmodel::ClimaLand.Soil.vanGenuchten{FT}
    "Porosity [unitless]"
    ν::FT
    "Residual water fraction [unitless]"
    θ_r::FT
    "Saturated hydraulic conductiity [m/s]"
    Ksat::FT
    "Surface flux, used as boundary condition [m/s]"
    F_sfc::FT
    "Domain of the model"
    domain::D
end;

# For reasons that will be clear momentarily, let's also define the name of the model:
ClimaLand.name(model::RichardsTutorialModel) = :soil;

# # Explicit tendency

# Here is where we need to specify the equations of motion. The prognostic variables
# for Richards equation consist of the volumetric water content at each location
# in the domain, `θ`.
# The differential equations are:

# ``
# \frac{\partial ϑ_l}{\partial t} = - ∇⋅[-K(θ) ∇(ψ(θ)+z)],
# ``

# where K(θ) is the hydraulic conductivity, and ψ(θ) is the matric potential.
# We now create the function which makes the `compute_exp_tendency!` function:

function ClimaLand.make_compute_exp_tendency(model::RichardsTutorialModel)
    function compute_exp_tendency!(dY, Y, p, t)
        gradc2f = ClimaCore.Operators.GradientC2F()
        interpc2f = ClimaCore.Operators.InterpolateC2F()
        FT = FTfromY(Y)
        divf2c = ClimaCore.Operators.DivergenceF2C(
            top = ClimaCore.Operators.SetValue(
                ClimaCore.Geometry.WVector.(model.F_sfc),
            ),
            bottom = ClimaCore.Operators.SetValue(
                ClimaCore.Geometry.WVector.(FT(0)),
            ),
        )

        @. dY.soil.θ =
            -(divf2c(-interpc2f(p.soil.K) * gradc2f(p.soil.ψ + p.soil.z)))
    end
    return compute_exp_tendency!
end;

# A couple of notes: the vector `` \vec{dY} `` contains the evaluation of the tendency function
# for each variable in `` \vec{Y} ``. It is updated in place (so no extra allocations are needed). Note that
# Y is not a simple array. It is a `ClimaCore` `FieldVector`, which allow us to impose some
# organizational structure on the state while still behaving like arrays in some ways.
# We use the symbol returned by `name(model)` to create the naming hierarchy in `Y`, `dY`, `p`. This 
# is useful for multi-component models.

# The arguments of `compute_exp_tendency!` are generic for any time-stepping algorithm.
# The `compute_exp_tendency!`
# function is only created once. If there are time-varying forcing terms appearing, for example, the
# forcing functions must be stored in `model` and passed in that way.

# # The prognostic state vector `` \vec{Y} `` and cache `` \vec{p} ``

# We have given the state vector `` \vec{Y} `` a particular structure, and don't expect the user
# to build this themselves.
# In order to have the structure `Y` (and `p`) correctly created, the model developer needs to define
# the names of the prognostic and auxiliary variables, as well as their types (often a floating point
# scalar), and where in the domain they are defined. For example, the volumetric water content
# is a scalar (type FT), with name θ, and it is defined throughout the subsurface of the domain.

ClimaLand.prognostic_vars(::RichardsTutorialModel) = (:θ,);
ClimaLand.prognostic_types(::RichardsTutorialModel{FT}) where {FT} = (FT,);
ClimaLand.prognostic_domain_names(::RichardsTutorialModel) = (:subsurface,);


# The auxiliary variables for this model are the hydraulic conductivity, matric potential,
# boundary fluxes, and heights of each level in the domain. All of these are scalars, and
# some are defined throughout the soil volume, or subsurface, while some are defined only on
# a surface (at the top or bottom of the domain).

ClimaLand.auxiliary_vars(::RichardsTutorialModel) =
    (:K, :ψ, :top_flux, :bottom_flux, :z)
ClimaLand.auxiliary_types(::RichardsTutorialModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT);
ClimaLand.auxiliary_domain_names(::RichardsTutorialModel) =
    (:subsurface, :subsurface, :surface, :surface, :subsurface);

# # Updating the cache
# We next need to define how we update the auxiliary variables. These are split between two
# functions, `update_aux!`, and `update_boundary_fluxes!`. For standalone component models,
# these could be combined into a single function, and indeed they could also be part of the tendency
# function itself. 
function ClimaLand.make_update_aux(model::RichardsTutorialModel)
    function update_aux!(p, Y, t)
        p.soil.z .=
            ClimaCore.Fields.coordinate_field(model.domain.space.subsurface).z # technically this does not need to update each step
        @. p.soil.K = ClimaLand.Soil.hydraulic_conductivity(
            model.vGmodel,
            model.Ksat,
            ClimaLand.Soil.effective_saturation(model.ν, Y.soil.θ, model.θ_r),
        )
        @. p.soil.ψ = ClimaLand.Soil.matric_potential(
            model.vGmodel,
            ClimaLand.Soil.effective_saturation(model.ν, Y.soil.θ, model.θ_r),
        )
    end
    return update_aux!
end;

function ClimaLand.make_update_boundary_fluxes(model::RichardsTutorialModel)
    function update_boundary_fluxes!(p, Y, t)
        FT = ClimaLand.FTfromY(Y)
        p.soil.top_flux .= model.F_sfc
        p.soil.bottom_flux .= FT(0)
    end
    return update_boundary_fluxes!
end;

# The default tendency function in ClimaLand for any AbstractModel carries out the following:

# ```julia
# function make_exp_tendency(model::AbstractModel)
#     update_aux! = make_update_aux(model)
#     update_boundary_fluxes! = make_update_boundary_fluxes(model)
#     compute_exp_tendency! = make_compute_exp_tendency(model)
#     function exp_tendency!(dY,Y,p,t)
#         update_aux!(p,Y,t)
#         update_boundary_fluxes!(p,Y,t)
#         compute_exp_tendency!(dY,Y,p,t)
#     end
#     return exp_tendency!
# end;
# ```

# Therefore, each time we need the tendency, we first update auxiliary variables, then update boundary
# fluxes, and then compute the tendency itself.

# Why do we do this? It would be straightforward, and arguably a lot simpler,
# to update the cache `p` within `compute_exp_tendency!`itself.
# The reason why we introduce these other functions is because we want to be able to
# combine standalone "component" models, like this one, with others, to create
# land surface models. For example, if we would like to run a land surface model with
# the soil and the canopy, the canopy auxiliary variables (e.g. interception of water
# and snow, transmitted radiation) affect the boundary fluxes of the soil.
# In this case, we must update auxiliary variables for all components, before computing
# boundary conditions and tendency functions.
# Please see the (LSM tutorial) for further explanation.

# More complex cases might require the evaluation of external data. For this, use the
# `TimeVaryingInput` interface. You can wrap functions, 1D/2D data in `TimeVaryingInput` to
# obtain an object that know how to evaluate that data on the model time (e.g., by
# performing linear interpolation). Then, in your model, you can just call
# `evaluate!(destination, itp, time)` to evaluate the itp on the given `time` and write the
# result to `dest` (typically a `Field`). With this common interface, you do not have to
# worry about the detail of the underlying data.

# # Running a simulation
# Create a model instance.
FT = Float32
vGmodel = ClimaLand.Soil.vanGenuchten(; α = 2.3f0, n = 2.0f0)
Ksat = FT(4.0e-7)
ν = 0.5f0
θ_r = 0.0f0
F_sfc = FT(-3.0e-8)
domain = ClimaLand.Domains.Column(; zlim = (-1.0f0, 0.0f0), nelements = 10)
soil = RichardsTutorialModel{Float32, typeof(domain)}(
    vGmodel,
    ν,
    θ_r,
    Ksat,
    F_sfc,
    domain,
);

# Create the initial state structure, using the default method. This step
# creates the vector Y and cache p, but initializes them with zeros.
Y, p, cds = initialize(soil);

# Note that `Y` has the structure we planned on in our `compute_exp_tendency!` function, for `x`,
Y.soil

# The same is true for `p`:
p.soil

# Here we now update `Y` in place with initial conditions of our choosing.

Y.soil.θ = 0.25f0;

# Set initial cache variable values, and inspect values:

set_initial_cache! = make_set_initial_cache(soil);
set_initial_cache!(p, Y, 0.0);
@show p.soil.K

@show p.soil.ψ

@show p.soil.top_flux

# Next up is to create the `exp_tendency!` function:
exp_tendency! = make_exp_tendency(soil);

# # Running the simulation

# Set the initial and end times, timestep:

t0 = 0.0;
tf = 7 * 24 * 3600.0;
dt = 1800.0;

# Select the timestepping algorithm we want to use from CTS.jl.
timestepper = CTS.RK4()
ode_algo = CTS.ExplicitAlgorithm(timestepper)

# SciMLBase problem statement using CTS.jl internals:
prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(T_exp! = exp_tendency!),
    Y,
    (t0, tf),
    p,
);

# Solve command:
sol = SciMLBase.solve(prob, ode_algo; dt = dt, saveat = dt);

# The solution is stored in `sol.u[k].soil.θ`, where k ranges over the number
# of timesteps.
