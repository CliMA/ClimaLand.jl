# # Setting up a Coupled Simulation

# For more information about the bucket model,
# please see [the bucket model tutorial](@ref "Introduction to the Land Bucket Model").

# This tutorial verbally describes how to set up a simulation for a coupled simulation. More detail for coupled runs can be
# found in the ClimaCoupler.jl [documentation](https://clima.github.io/ClimaCoupler.jl/dev/).

# # Background

# Recall that in order to drive the system in standalone mode,
# the user must provide
# prescribed functions of time for the water volume flux in precipitation,
#  for the net downward shortwave and longwave
# radiative energy fluxes,
# for the atmospheric temperature `T_a`,
# wind speed `u_a` (m/s), specific humidity `q_a`, and air density
# `ρ_a` (kg/m^3) at a reference height `h_a` (m). The values
# of these quantities at any given time are stored in `p.drivers`.

# Turbulent surface fluxes are computed by the bucket model at each step of the
# simulation, using the land surface properties
# as well as the prescribed atmospheric properties, according to Monin-Obukhov theory.
# These fluxes, as well as the net radiation, are stored in the cache
# of the bucket model: `p.bucket.turbulent_fluxes.lhf`, `p.bucket.turbulent_fluxes.shf`,
# `p.bucket.turbulent_fluxes.vapor_flux`, `p.bucket.R_n`, where they are accessible
# when boundary conditions are required in the ODE functions (right hand side) of the
# prognostic equations. Similarily, the precipitation rates are provided from
# prescribed conditions and stored in `p.drivers.P_liq`, `p.drivers.P_snow`.

# # Coupled simulations

# In a coupled simulation, this changes in two key ways. First, the coupler now computes
# the turbulent surface fluxes and the net radiative flux for the bucket, but it still does so
# by calling the functions (`turbulent_fluxes!, net_radiation!`) defined in ClimaLand.
# This means that the coupler avoids needing to have any internal knowledge about how these
# fluxes are computed. In order for this to work, though, the fields in `p.drivers` need to be
# updated. The second key change then is that the coupler updates these fields directly with
# the atmopsheric state provided to it by the atmosphere model. 
# As in standalone runs, the cache is updated with the correct fluxes.
