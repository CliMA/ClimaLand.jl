# # Background
# When solving the system of equations describing the storage of water and energy in snow,
# "boundary conditions" must be set. These are only true boundary conditions
# if one considers the equations for water and energy in snow to be discretized PDEs,
# but they do represent the fluxes at upper boundary (the snow/atmosphere interface) and at
# the lower boundary (the snow/soil interface), so we refer to them as boundary conditions
# in order to use a consistent notation with the soil model.

# At present, the only type of boundary condition the snow model supports is
# `AtmosDrivenSnowBC`, which is an attempt at a compact way of indicating that with this
# choice, the snow model exchanges water and energy with the atmosphere via sensible, latent, and radiative heat fluxes, and by
# precipitation and sublimation/evaporation. With this choice, the snow model also has exchange fluxes with the soil (melt water,
# and a conductive ground heat flux). This boundary condition type includes the atmospheric and radiative
# forcings, and, importantly, a Tuple which indicates which components of the land are present.
# The last argument is what we will describe here. For more information on how to supply the prescribed forcing
# data, please see the tutorial linked [here](@ref "Using atmospheric and radiative drivers"). For an explanation of the
# same design implementation for the Soil model, please see [here](@ref "Boundary conditions for the soil model").

# # Adjusting the boundary conditions for the snow model when run as part of an integrated land model

# The presence of other land components (a canopy, the soil) affects the boundary conditions of the snow
# and changes them relative to what they would be if only bare snow was interacting with the atmosphere.
# For example, the canopy absorbs radiation that might otherwise be absorbed by the snow, the
# canopy intercepts precipitation, and the soil and snow interact thermally at the interface between them.
# Because of this, we need a way to indicate to the snow model that the boundary conditions must be computed in a different way
# depending on the components. We do this via the field in the boundary condition of type `AtmosDrivenSnowBC` via
# a field `prognostic_land_components`,
# which is a Tuple of the symbols associated with each component in the land model. For example,
# if you are simulating the snow model in standalone mode, you would set

# `prognostic_land_components = (:snow,)`

# `boundary_condition = ClimaLand.Snow.AtmosDrivenSnowBC(atmos_forcing, radiative_forcing, prognostic_land_components)`

# while, if you are simulating both a prognostic snow and soil model, you would set

# `prognostic_land_components = (:snow, :soil)`

# `boundary_condition = ClimaLand.Snow.AtmosDrivenSnowBC(atmos_forcing, radiative_forcing, prognostic_land_components)`

# Note that the land components are *always* in alphabetical order, and the symbols associated with each
# land model component are *always* the same as the name of that component, i.e.,

# `ClimaLand.name(snow_model) = :snow`

# `ClimaLand.name(canopy_model) = :canopy`

# `ClimaLand.name(soilco2_model) = :soilco2`

# `ClimaLand.name(soil_model) = :soil`

# Currently, the only supported options are: `(:snow,)`, and `(:snow, :soil)`.

# # How it works
# When the snow model computes its boundary fluxes, it calls the function
# `snow_boundary_fluxes!(bc, snow_model, Y, p, t)`, where `bc`
# is the boundary condition of type `AtmosDrivenSnowBC`.
# This function then calls `snow_boundary_fluxes!(bc, Val(bc.prognostic_land_components), snow, Y, p, t)`.
# Here multiple dispatch is used to compute the fluxes correctly according to the value of the `prognostic_land_components`,
# i.e., based on which components are present.

# # Some important notes
# When the snow model is run in standalone mode (`prognostic_land_components = (:snow,)`), the ground heat flux is approximate
# as zero. When the soil and snow model are run together, this flux is computed and affects both the snow and soil.
