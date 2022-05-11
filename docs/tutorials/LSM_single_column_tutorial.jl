# The `AbstractModel` 
# [tutorial](https://clima.github.io/ClimaLSM.jl/dev/generated/model_tutorial)
# describes how a user can run
# simulations of a physical system governed by differential equations.
# In this framework, the user must define a model type for their problem,
# which contains all of the information required to set up the system of
# equations. By extending the methods for `make_rhs(model)`,
# `prognostic_variables(model)`, etc, the information
# stored in the `model` is used to make the system of equations.
# Given initial conditions, these equations can then be stepped forward in
# time using the time-stepper of your choice (we are set up to use
# OrdinaryDiffEq.jl currently).

# The benefit of this framework is that it can be used for both individual
# components of an LSM (soil, snow, rivers, canopy biophysics, carbon...)
# **as well as the LSM itself**. Here we explain how a simple single column two
# component model can be set up using this software interface. Additionally,
# we demonstrate here the use of the auxiliary or cache variables,
# which were mentioned but not needed in the Henon-Heiles problem
# solved in the prior tutorial.

# We'll first demonstrate how to set up two components in
# standalone mode, before spending time explaining the LSM setup.
# In our example, we have a component which accounts for soil hydrology
# via the Richardson-Richards (RR) equation.  Our second component is a surface
# water model without lateral flow (standing water, as in a pond). For
# more details on these models, and how they were set up,
# please feel free to look at the source
# code [here](https://github.com/CliMA/ClimaLSM.jl/blob/main/src/Soil/Soil.jl)
# and
# [here](https://github.com/CliMA/ClimaLSM.jl/blob/main/src/SurfaceWater/Pond.jl).
# This tutorial focuses on using the `AbstractModel`s framework to set up
# the equations, rather than on running simulations.

# First, let's load the required modules:
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: LSMSingleColumnDomain, Column
using ClimaLSM.Soil
using ClimaLSM.Pond

FT = Float64;

# # The individual component models I - Soil Hydrology
# The RR equation for the volumetric water content of soil is given by

# ``
# \frac{\partial ϑ}{\partial t} = -∇ ⋅ (-K∇(ψ+z)) + S(x,y,z, t)
# ``

# In order to solve this, one must specify:
# - boundary conditions,
# - relevant parameters (closure models for `K` and `ψ`),
# - a domain and a spatial discretization scheme,
# - additional source terms `S`, if applicable,
# - a time-stepping algorithm,
# - initial conditions.

# We make the distinction between the spatially discretized
# equations (for which you need parameters, boundary conditions, source terms,
# and domain/
# discretization scheme information in order to write down and evaluate), and the
# simulation you want to run (for which you need the equations,
# initial conditions,
# a time span, and a time-stepping scheme in order to specify completely).

# Here, we'll focus on what you need to write the equations. In the design
# of all CliMA systems, everything you need to write the equations is stored
# in the model structure itself, so that we can call
# `make_ode_function(model)` and get back a function which computes the time
# derivative of the prognostic variables,
# which the ODE timestepper needs to advance the state forward
# in time.

# For the RR equation, we can create this as follows. First, we specify parameters:
ν = FT(0.495);
K_sat = FT(0.0443 / 3600 / 100); # m/s
S_s = FT(1e-3); #inverse meters
vg_n = FT(2.0);
vg_α = FT(2.6); # inverse meters
vg_m = FT(1) - FT(1) / vg_n;
θ_r = FT(0);
soil_ps = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r);

# Next, let's define the spatial domain and discretization:
zmax = FT(0);
zmin = FT(-1);
nelems = 20;
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems);

# And boundary conditions and source terms (none currently):
top_flux_bc = FT(0.0)
bot_flux_bc = FT(0.0)
sources = ()
boundary_fluxes = FluxBC{FT}(top_flux_bc, bot_flux_bc);

# With this information, we can make our model:
soil = Soil.RichardsModel{FT}(;
    parameters = soil_ps,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = sources,
);
# We also can create the soil prognostic and auxiliary
# `ClimaCore.Field.FieldVector`s using the default method for `initialize`,
Y_soil, p_soil, coords_soil = initialize(soil);
# and we can set up the ode function using the default as well,
soil_ode! = make_ode_function(soil);
# which computes, for the column domain,

# ``
# -\frac{∂ }{∂z} (-K\frac{∂(ψ+z)}{∂ z})
# ``


# for each value of ϑ on the mesh of our `soil_domain`.

# Note that the soil model does include both hydraulic `K` and pressure
# head `ψ` in the auxiliary vector, so the fields `p_soil.soil.K` and
# `p_soil.soil.ψ`
# are present. These are automatically updated first in each call
# to `soil_ode!`, as follows:
# ```julia
# function soil_ode!(dY, Y, p, t)
#          update_aux!(p,Y,t)
#          rhs!(dY, Y, p, t)
# end
# ```

# where `update_aux!` updates `K`, and `ψ`, in `p`, in place,
# and `rhs!` computes the divergence of the Darcy flux, using the updated
# `p`, and then updates `dY`
# in place with the computed values.
# For this reason,
# the `p` vector does not need to be set to
# some initial condition consistent
# with Y_soil.soil.ϑ(t=0) before starting a simulation, though initial
# conditions must be given for `Y`.

# Note also that we have defined methods `make_rhs` and `make_update_aux`, which only
# take the `model` as argument, and which return the functions `update_aux!`
# and `rhs!`, [here](https://github.com/CliMA/ClimaLSM.jl/blob/9e2b8df6d8d5b7f878a5b0f02e2bf80fa66aec33/src/Soil/Soil.jl#L292)
# and [here](https://github.com/CliMA/ClimaLSM.jl/blob/9e2b8df6d8d5b7f878a5b0f02e2bf80fa66aec33/src/Soil/Soil.jl#L400).

# Lastly, the coordinates returned by `initialize` contain
# the z-coordinates of the centers of the finite difference layers used
# for spatial discretization of the PDE.

# # The individual component models II - Surface Water
# The pond model has a single
# variable, the pond height η, which satisfies the ODE:

# ``
# \frac{∂ η}{∂ t} = -(P - I) = R,
# ``

# where P is the precipitation, I the infiltration into the soil, and R
# is the runoff. Note that P, I < 0 indicates flow in the -ẑ direction.

# To write down the pond equations, we need to specify
# - P
# - I
# which are akin to boundary conditions. In standalone mode, 
#  one would need to pass in **prescribed** functions of time and store
# them inside our pond model, since again, the pond model structure
# must contain everything needed to make the ode function:

precipitation(t::T) where {T} = t < T(20) ? -T(1e-5) : T(0.0) # m/s

infiltration(t::T) where {T} = -T(1e-6) #m/s
pond_model = Pond.PondModel{FT}(;
    runoff = PrescribedRunoff{FT}(precipitation, infiltration),
);

# Here, `PrescribedRunoff` is the structure holding the prescribed
# driving functions for `P` and `I`.

# Again we can initialize the state vector and auxiliary vectors:
Y_pond, p_pond, coords_pond = initialize(pond_model);

# We can make the ode function in the same way, for
# stepping the state forward in time:
pond_ode! = make_ode_function(pond_model);

# The `pond_ode!` function works in the same way as for the soil model:
# ```julia
# function pond_ode!(dY, Y, p, t)
#          update_aux!(p,Y,t)
#          rhs!(dY, Y, p, t)
# end
# ```
# but the `update_aux!` does not alter `p` at all in this case.
# The pond model does not have auxiliary variables, so
# `p_pond` is empty.



# The coordinates here are relatively meaningless - we
# are solving for the pond height at a point in space on the surface of the
# Earth, and by default
# this assigns a `Point` domain, with a coordinate of `z_sfc = 0`.
# In a simulation with horizontal
# resolution, the `coordinates` returned would be the `(x,y,z=z_sfc(x,y))`
# coordinates of the surface, which are more useful.


# # An LSM with pond and soil:
# The LSM model must contain everything needed to write
# down the joint system of equations

# ``
# \frac{\partial \eta}{\partial t} = -(P(t) - I(ϑ, η, P)) = R,
# ``

# ``
# \frac{\partial ϑ}{\partial t} = -∇ ⋅ (-K∇(ψ+z)) + S
# ``

# ``
# -K ∇(ψ+z)|_{z = zmax}  ⋅ ẑ = I(ϑ, η, P)
# ``

# ``
# -K ∇(ψ+z)|_{z = zmin}  ⋅ ẑ = 0.0.
# ``

# These two components
# interact via the infiltration term `I`. Infiltration is a boundary
# condition for the soil, and affects the source term for the surface
# water equation. Infiltration depends on precipitation, the soil
# moisture state, and the pond height.

# As in the standalone cases, defining our model
# requires specifying
# - parameters,
# - domains, discretizations
# - precipitation,
# - boundary conditions,
# - sources in the soil equation, if any. 

# First, let's make our LSM domain, which now contains 
# information about the subsurface domain and the surface domain. For
# a single column, this means specifying the boundaries of the soil domain
# and the number of elements.
lsm_domain = LSMSingleColumnDomain(; zlim = (zmin, zmax), nelements = nelems);
# The surface domain is again just a `Point` with
# `z_sfc = zmax`.
lsm_domain.surface
# The subsurface domain is a column from zmin to zmax:
lsm_domain.subsurface

# Let's now collect the needed arguments for the soil and pond
# models:

soil_args = (parameters = soil_ps, domain = lsm_domain.subsurface, sources = ());
surface_water_args = (domain = lsm_domain.surface,);
# Atmospheric drivers don't "belong" to either component alone:
land_args = (precip = precipitation,);
land = LandHydrology{FT}(;
    land_args = land_args,
    soil_model_type = Soil.RichardsModel{FT},
    soil_args = soil_args,
    surface_water_model_type = Pond.PondModel{FT},
    surface_water_args = surface_water_args,
);

# Here, `LandHydrology` is a type of `AbstractModel` which has a surface
# water model (Pond or otherwise) and a soil model (RR, or perhaps otherwise).
# Note that we pass in the type of the soil and surface water model -
# these could be more complex, e.g. a river model with lateral
# flow could be used in place of the `Pond`. We could also add in
# a snow component.

# Now, note that we did not specify the infiltration function, like we did
# in standalone pond mode, nor did we specify boundary conditions for the
# soil model. Yet, before we stressed that the model needs to have everything
# required to write down and evaluate the time derivative of the ODEs.
# So, how does this work?

# Here, the LSM model **constructor** is given the information needed to make
# both the soil model and the pond model. Then, it is like running
# the pond and soil model in standalone mode, in series, **except**
# we have defined methods internally for computing the boundary condition
# and pond source term correctly, based on `I`, instead of using
# prescribed values passed in.
# The LSM constructor creates the correct
# `boundary_fluxes` object for the soil model, and the correct
# `infiltration` object
# for the pond model under the hood.

# To advance the state of the joint system (ϑ, η)
# from time `t` to
# time `t+Δt`, we must compute the infiltration at `t`.
# This value is stored
# in `p.soil_infiltration`, and reflects a proper use of the auxiliary
# or cache state: storing a quantity which we would rather compute once
# and store, rather than compute twice, once in the soil ode function,
# and once in the pond ode function. This guarantees the same value is
# used for both equations. In pseudo code, we have:
# ```julia
# function make_update_aux(land)
#          soil_update_aux! = make_update_aux(land.soil)
#          surface_update_aux! = make_update_aux(land.surface_water)
#          interactions_update_aux! = make_update_aux(land, land.soil, land.surface_water)
#          function update_aux!(p,Y,t)
#                   surface_update_aux!(p,Y,t) # does nothing to `p`
#                   soil_update_aux!(p,Y,t) # updates p.soil.K and p.soil.ψ
#                   interactions_update_aux!(p,Y,t) # updates p.soil_infiltration
#          end
#          return update_aux!
# end
# ```

# and similarily for the `rhs!` functions:
# ```julia
# function make_rhs(land)
#          soil_rhs! = make_update_aux(land.soil)
#          surface_rhs! = make_update_aux(land.surface_water)
#          function rhs!(dY,Y,p,t)
#                   surface_rhs(dY,Y,p, t), # computes dY.surface.η
#                   soil_rhs!(dY,Y,p,t) # computes dY.soil.ϑ
#          end
#          return rhs!
# end
# ```

# The `ode_function!` for the
# land model is then again just
# ```julia
# function ode_function!(dY, Y, p, t)
#          update_aux!(p,Y,t)
#          rhs!(dY, Y, p, t)
# end
# ```

# In the above, we showed explicitly what occurs by hardcoding the
# `rhs!`, `update_aux!` with names for `soil` and `surface_water`.
# In reality, this is done by
# looping over the components of the land model, meaning that we can
# use the same code internally for land models with different components.

# A similar composition occurs for initializing the state itself:
# Calling `initialize(land)` does
# four things:
# - initialize(land.soil) 
# - initialize(land.surface_water) 
# - initializes interaction terms, like p.soil_infiltration
# - append these into Y, p, and coords:
Y, p, coords = initialize(land);
# We have volumetric liquid water fraction:
propertynames(Y.soil)
# and surface height of the pond:
propertynames(Y.surface_water)
# as well as auxiliary variables for the soil:
propertynames(p.soil)
# and nothing for surface water:
propertynames(p.surface_water)
# and the shared interaction term
propertynames(p)
# and finally, coordinates - useful for visualization of solutions:
coords.subsurface
# and the coordinates of the surface variables:
coords.surface

# And we can make the ode function as before:
land_ode! = make_ode_function(land);

# Next up would be to set initial conditions, choose a timestepping scheme, and
# run your simulation.

# # Advantages and disadvantages

# Some advantages to our interface design are as follows:
# - a developer only needs to learn a few concepts (`rhs!`, prognostic vs. aux variables, `update_aux!`, `initialize`, domains) to make a model which can be run in standalone or work with other components.
# - likewise, a user only needs to learn one interface to run all models, regardless of if they are standalone components or LSMs with multiple componnents.
# - the `ode_function! `is completely seperate from the timestepping scheme used, so any scheme can be used (with the exception of mixed implicit/explicit schemes, which we can't handle yet).
# - although we wrote it here in a ``hardwired`` fashion for surface water and soil, the `update_aux!`, `rhs!`, etc. many methods for LSM models generalize to any number and mix of components. One just needs to write a new model type (e.g. `BiophysicsModel <: AbstractModel` for a vegetation and carbon component model) and the appropriate `interaction` methods for that model.
# - the order in which the components are treated in the rhs or in update aux does not matter. What matters is that auxiliary/cache variables are updated first, and within this update, interactions are updated last. We assume that the rhs! function for a component only needs the entire `p` and `Y.component` in making this statement. Similarily, updating the aux variables of a single component does not require interaction variables. Yhis is also the same as saying they can be run in *standalone* mode. 
# - the code is also modular in terms of swapping out a simple component model for a more complex version.

# Possible disadvantages to our interface design:
# - Even in standalone model, variables are accessed in a nested way: Y.soil, p.soil, etc, which is excessive.
# - To accomodate the fact that some components involve PDEs, a developer for purely ODE based component does need to at least handle `ClimaCore.Field.FieldVector`s.
# - standalone models need to play by the rules of `AbstractModel`s, and LSMs need to play by the rules of `ClimaLSM.jl`.
