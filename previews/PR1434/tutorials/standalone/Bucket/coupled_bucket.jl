# # Setting up a Coupled Simulation

# For more information about the bucket model,
# please see [the bucket model tutorial](https://clima.github.io/ClimaLand.jl/stable/generated/bucket_tutorial/).

# This tutorial shows how to set up a simulation for a coupled simulation. More detail for coupled runs can be
# found in the ClimaCoupler.jl [documentation](https://clima.github.io/ClimaCoupler.jl/stable/). In
# preparation for understanding this tutorial, we recommend also reading the [intro to multi-component models tutorial](https://clima.github.io/ClimaLand.jl/stable/generated/LSM_single_column_tutorial/) as well as being familiar
# with multiple dispatch programming in Julia.

# # Background

# Recall that in order to drive the system in standalone mode,
# the user must provide
# prescribed functions of time for the water volume flux in precipitation,
#  for the net downward shortwave and longwave
# radiative energy fluxes,
# for the atmospheric temperature `T_a`,
# wind speed `u_a` (m/s), specific humidity `q_a`, and air density
# `ρ_a` (kg/m^3) at a reference height `h_a` (m).

# Turbulent surface fluxes are computed by the bucket model at each step of the
# simulation, using the land surface properties
# as well as the prescribed atmospheric properties, according to Monin-Obukhov theory.
# These fluxes, as well as the net radiation, are stored in the `auxiliary` state
# of the bucket model: `p.bucket.turbulent_fluxes.lhf`, `p.bucket.turbulent_fluxes.shf`,
# `p.bucket.turbulent_fluxes.vapor_flux`, `p.bucket.R_n`, where they are accessible
# when boundary conditions are required in the ODE functions (right hand side) of the
# prognostic equations. Similarily, the precipitation rates are provided from
# prescribed conditions and stored in `p.drivers.P_liq`, `p.drivers.P_snow`.

# In a coupled simulation, this changes. The coupler computes turbulent surface fluxes
# based on information (prognostic state, parameters) passed to it by both the atmosphere and land models.
# Net radiation
# is computed within the atmosphere model, using the prognostic land surface temperature and the land surface
# albedo, and passed back to the land model via the coupler. These details are important, but from
# the point of view of the land
# model, we only need to know that the coupler accesses land model variables to compute fluxes,
# and that the coupler passes these fluxes back to the land model.

# In our current setup, "passed back to the land model via the coupler" means that the coupler
# accesses the auxiliary state of the land model and modifies it, at each step in the simulation, so that
# it holds the current net radiation, precipitation, and turbulent surface fluxes (`p.bucket.turbulent_fluxes`,
# `p.bucket.R_n`, `p.drivers.P_liq`, `p.drivers.P_snow`).
# These quantities are then still available in the ODE functions
# of the prognostic equations for the bucket model, as in the standalone case.

# In order for the land model to be able to run both in standalone mode, and a coupled mode,
# within a single interface, we make use of multiple dispatch.

# # Turbulent Surface Fluxes and Radiation

# Let's review how turbulent surface fluxes and radiation are computed by the land model.
# The user first creates the prescribed atmosphere and prescribed radiation drivers. In
# pseudo code, this might look something like:

# ```julia
# prescribed_atmos = PrescribedAtmosphere{FT}(*driver data passed in here*)
# prescribed_radiation = PrescribedRadiativeFluxes{FT}(*driver data passed in here*)
# ```

# These are stored in the [BucketModel](https://clima.github.io/ClimaLand.jl/stable/APIs/Bucket/#ClimaLand.Bucket.BucketModel) object,
# along with [BucketParameters](https://clima.github.io/ClimaLand.jl/stable/APIs/Bucket/#ClimaLand.Bucket.BucketParameters).
# In order to compute turbulent surface fluxes, we call [turbulent_fluxes!](https://clima.github.io/ClimaLand.jl/stable/APIs/shared_utilities/#ClimaLand.turbulent_fluxes!),
# with arguments including `prescribed_atmos`. Since this argument is of the type `PrescribedAtmosphere`, the method of `turbulent_fluxes` which is executed is one which computes the turbulent surface fluxes
# using MOST. We have a similar function for [net_radiation](https://clima.github.io/ClimaLand.jl/stable/APIs/shared_utilities/#ClimaLand.net_radiation) and which computes the net radiation based on the prescribed downwelling radiative fluxes, stored in an argument
# `prescribed_radiation`, which is of type `PrescribedRadiation`.

# In the coupled case, we want different behavior. We have defined new _coupled_ types to
# use instead of the "prescribed" types:

# ```julia
# struct CoupledAtmosphere{FT} <: AbstractAtmosphericDrivers{FT} end
# struct CoupledRadiativeFluxes{FT} <: AbstractRadiativeDrivers{FT} end
# ```

# Then, we have defined a new method for `turbulent_fluxes!` and `net_radiation!` which dispatch for these types,
# and updates dest with the fluxes that the coupler has updated `p.bucket.turbulent_fluxes` and `p.bucket.R_n` with.
# In pseudo code:
#
# ```julia
# function ClimaLand.turbulent_fluxes!(
#    dest,
#    atmos::CoupledAtmosphere,
#    model::BucketModel,
#    p)
#    dest .= (
#         lhf = p.bucket.turbulent_fluxes.lhf,
#         shf = p.bucket.turbulent_fluxes.shf,
#         vapor_flux = p.bucket.turbulent_fluxes.vapor_flux,
#     )
# end
# ```

# similarily:

# ```julia
# function ClimaLand.net_radiation!(
#     dest,
#     radiation::CoupledRadiativeFluxes{FT},
#     model::BucketModel{FT},
#     p)
#     dest .= p.bucket.R_n
# end
# ```

# These methods update the values stored in the auxiliary state `p`. Importantly, these functions are
# called by the bucket model
# each time step **after** the coupler has already computed these values
# (or extracted them from another model) and modified `p`!

# # Surface air density
# Within the right hand side/ODE function calls for the bucket model, we need both the
# surface air density (for computing specific humidity at the surface). In standalone runs,
# we call the function [`surface_air_density`](https://clima.github.io/ClimaLand.jl/stable/APIs/shared_utilities/#ClimaLand.surface_air_density),
# When the `atmos` type is `PrescribedAtmosphere`, this function uses the atmospheric state and surface
# temperature to estimate the surface air density assuming an ideal gas and hydrostatic balance
# and by extrapolating from the air density at the lowest level of the atmosphere.

# In the coupled case, we need to extend these functions with a `CoupledAtmosphere` method:
# ```julia
# function ClimaLand.surface_air_density(
#     atmos::CoupledAtmosphere,
#     model::BucketModel,
#     p)
#     return p.bucket.ρ_sfc
# end
# ```

# Again, this functions is called in the ODE function of the bucket model *after* the coupler
# has updated the values of `p` with the correct values at that timestep.
