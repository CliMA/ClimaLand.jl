# # Background
# When solving the partial differential equations describing the storage of water, ice, and energy in soil,
# boundary conditions must be set for the soil liquid water and energy.
# The tutorial on [boundary conditions for the soil model](@ref "Boundary conditions for the soil model") describes
# the various boundary condition options we currently support. Many of these are suitable for use when
# replicating laboratory experiments or other standalone soil settings - for example, a particular experiment
# may fix the temperature or water content at the top of a soil column. However, when we wish to model the
# soil interacting with the atmosphere, in standalone or integrated models, the particular boundary condition
# we use is one that computes the sensible, latent, and radiative fluxes in addition to evaporation and sublimation;
# this option is called the `AtmosDrivenFluxBC`. This boundary condition type includes the atmospheric and radiative
# forcings, the runoff parameterization, and finally, a Tuple which indicates which components of the land are present.
# The last argument is what we will describe here. For more information on how to supply the prescribed forcing
# data, please see the tutorial linked [here](@ref "Using atmospheric and radiative drivers").

# # Adjusting the boundary conditions for the soil model when run as part of an integrated land model

# The presence of other land components (a canopy, snow) affects the boundary conditions of the soil
# and changes them relative to what they would be if only bare soil was interacting with the atmosphere.
# For example, the canopy and snow absorb radiation that might otherwise be absorbed by the soil, the
# canopy intercepts precipitation, and the snow can melt and contribute to soil infiltration. Because of this,
# we need a way to indicate to the soil model that the boundary conditions must be computed in a different way.
# In all cases, the same `AtmosDrivenFluxBC` is used with the exception of the field `prognostic_land_components`,
# which is a Tuple of the symbols associated with each component in the land model. For example,
# if you are simulating the soil model in standalone mode, you would set

# `prognostic_land_components = (:soil,)`

# `top_soil_boundary_condition = ClimaLand.Soil.AtmosDrivenFluxBC(atmos_forcing, radiative_forcing, runoff_parameterization, prognostic_land_components)`

# while, if you are simulating both a prognostic canopy and soil model, you would set

# `prognostic_land_components = (:canopy, :soil)`

# `top_soil_boundary_condition = ClimaLand.Soil.AtmosDrivenFluxBC(atmos_forcing, radiative_forcing, runoff_parameterization, prognostic_land_components)`

# Note that the land components are *always* in alphabetical order, and the symbols associated with each
# land model component are *always* the same as the name of that component, i.e.,

# `ClimaLand.name(snow_model) = :snow`

# `ClimaLand.name(canopy_model) = :canopy`

# `ClimaLand.name(soilco2_model) = :soilco2`

# `ClimaLand.name(soil_model) = :soil`

# Currently, the only supported options are: `(:soil,), (:canopy, :soil, :soilco2), (:snow, :soil)`.

# # How it works
# When the soil model updates the boundary conditions, it calls the function
# `soil_boundary_fluxes!(bc, ClimaLand.TopBoundary(), soil_model, Δz_top, Y, p, t)`, where `bc`
# is the boundary condition being used at the top of the domain.
# This function has multiple methods depending on the boundary condition type `typeof(bc)`. When that type is `AtmosDrivenFluxBC`,
# the function then calls `soil_boundary_fluxes!(bc, Val(bc.prognostic_land_components), soil, Y, p, t)`.
# Again multiple dispatch is used, this time to compute the fluxes correctly according to the value of the `prognostic_land_components`,
# i.e., based on which components are present.
