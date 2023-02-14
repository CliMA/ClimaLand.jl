# # Setting up a Coupled Simulation

# For more information about the bucket model,
# please see [the bucket model tutorial](https://clima.github.io/ClimaLSM.jl/dev/generated/bucket_tutorial/).

# This tutorial shows how to set up a simulation for a coupled simulation. More detail for coupled runs can be
# found in the ClimaCoupler.jl [documentation](https://clima.github.io/ClimaCoupler.jl/dev/). In
# preparation for understanding this tutorial, we recommend also reading the [intro to multi-component models tutorial](https://clima.github.io/ClimaLSM.jl/dev/generated/LSM_single_column_tutorial/) as well as being familiar
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
# of the bucket model: `p.bucket.turbulent_energy_flux`, `p.bucket.evaporation`, `p.bucket.R_n`, where they are accessible
# when boundary conditions are required in the ODE functions (right hand side) of the
# prognostic equations.

# In a coupled simulation, this changes. As we will see, the coupler computes turbulent surface fluxes
# based on information (prognostic state, parameters) passed to it by both the atmosphere and land models.
# Net radiation 
# is computed within the atmosphere model, using the prognostic land surface temperature and the land surface
# albedo, and passed back to the land model via the coupler. Similarily, precipitation is computed within the
# atmosphere model, and passed back to the land model via the coupler. These details are important, but from
# the point of view of the land
# model, we only need to know that the coupler accesses land model variables to compute fluxes,
# and that the coupler passes these fluxes back to the land model. 

# In our current setup, "passed back to the land model via the coupler" means that the coupler
# accesses the auxiliary state of the land model and modifies it, at each step in the simulation, so that
# it holds the current net radiation, precipitation, and turbulent surface fluxes (`p.bucket.turbulent_energy_flux`,
# `p.bucket.evaporation`, `p.bucket.R_n`, `p.bucket.P_liq`, `p.bucket.P_snow`). These quantities are then still available in the ODE functions
# of the prognostic equations for the bucket model, as in the standalone case.

# In order for the land model to be able to run both in standalone mode, and a coupled mode,
# within a single interface, we make use of multiple dispatch.

# To begin, let's import some necessary abstract and concrete types, as well as methods.
using ClimaLSM
using ClimaLSM.Drivers: AbstractAtmosphericDrivers, AbstractRadiativeDrivers
using ClimaLSM.Bucket: BucketModel, BucketModelParameters

import ClimaLSM.Bucket:
    surface_fluxes,
    surface_air_density,
    liquid_precipitation,
    snow_precipitation

FT = Float64;
# # Turbulent Surface Fluxes and Radiation

# Let's review how turbulent surface fluxes and radiation are computed by the land model.
# The user first creates the prescribed atmosphere and prescribed radiation drivers. In
# pseudo code, this might look something like:

# ` prescribed_atmos = PrescribedAtmosphere{FT}(*driver functions passed in here*)`
# ` prescribed_radiation = PrescribedRadiativeFluxes{FT}(*driver functions passed in here*) `

# These are stored in the [BucketModel](https://clima.github.io/ClimaLSM.jl/dev/APIs/Bucket/#ClimaLSM.Bucket.BucketModel) object,
# along with [BucketParameters](https://clima.github.io/ClimaLSM.jl/dev/APIs/Bucket/#ClimaLSM.Bucket.BucketParameters).
# In order to compute turbulent surface fluxes, we call [surface_fluxes](https://clima.github.io/ClimaLSM.jl/dev/APIs/Bucket/#ClimaLSM.Bucket.surface_fluxes),
# with arguments including `prescribed_atmosphere`. Since this argument is of the type `PrescribedAtmosphere`, the method of `surface_fluxes` which is executed is one which computes the turbulent surface fluxes
# using MOST. We have a similar function for [net_radiation](https://clima.github.io/ClimaLSM.jl/dev/APIs/Bucket/#ClimaLSM.Bucket.net_radiation) and which computes the net radiation based on the prescribed downwelling radiative fluxes, stored in an argument
# `prescribed_radiation`, which is of type `PrescribedRadiation`.

# In the coupled case, we want different behavior. Inside coupler source code, we define new ``coupled`` types to
# use instead of the "prescribed" types:

struct CoupledAtmosphere{FT} <: AbstractAtmosphericDrivers{FT} end
struct CoupledRadiativeFluxes{FT} <: AbstractRadiativeDrivers{FT} end

# Then, we define a new method for `surface_fluxes` and `net_radiation` which dispatch for these types:
function ClimaLSM.Bucket.surface_fluxes(
    atmos::CoupledAtmosphere{FT},
    p,
    t,
    _...,
) where {FT <: AbstractFloat}
    return (
        turbulent_energy_flux = p.bucket.turbulent_energy_flux,
        evaporation = p.bucket.evaporation,
    )
end


function ClimaLSM.Bucket.net_radiation(
    radiation::CoupledRadiativeFluxes{FT},
    p,
    t,
    _...,
) where {FT <: AbstractFloat}
    return p.bucket.R_n
end
# Essentially, these methods simply returns the values stored in the auxiliary state `p`. Importantly, these functions are
# called by the bucket model 
# each time step **after** the coupler has already computed these values
# (or extracted them from another model) and modifed `p`!

# # Precipitation and surface air density
# Within the right hand side/ODE function calls for the bucket model, we need both the liquid/snow precipitation
# and the surface air density (for computing specific humidity at the surface). In standalone runs,
# we call the functions [`surface_air_density`](https://clima.github.io/ClimaLSM.jl/dev/APIs/Bucket/#ClimaLSM.Bucket.surface_air_density),
#  [`liquid_precipitation`](https://clima.github.io/ClimaLSM.jl/dev/APIs/Bucket/#ClimaLSM.Bucket.liquid_precipitation), and
# [`snow_precipitation`](https://clima.github.io/ClimaLSM.jl/dev/APIs/Bucket/#ClimaLSM.Bucket.snow_precipitation).
# When the `atmos` type is `PrescribedAtmosphere`, these
# return the prescribed values for these quantities.

# In the coupled case, we need to extend these functions with a `CoupledAtmosphere` method:
function ClimaLSM.Bucket.surface_air_density(atmos::CoupledAtmosphere, p, _...)
    return p.bucket.ρ_sfc
end

function ClimaLSM.Bucket.liquid_precipitation(atmos::CoupledAtmosphere, p, _)
    return p.bucket.P_liq
end

function ClimaLSM.Bucket.snow_precipitation(atmos::CoupledAtmosphere, p, _)
    return p.bucket.P_snow
end

# Again, these functions are called in the ODE function of the bucket model *after* the coupler
# has updated the values of `p` with the correct values at that timestep.

# Two other notes: (1) the coupler code actually extends the auxiliary state `p` of the bucket model,
# adding in the fields `P_liq`, `P_snow`, and `ρ_sfc`. These fields do not exist in `p` if the standalone
# bucket model, and (2) the surface air density is computed assuming an ideal gas and hydrostatic balance
# and by extrapolating from the air density at the lowest level of the atmosphere, which is why it is
# handled by the coupler.
