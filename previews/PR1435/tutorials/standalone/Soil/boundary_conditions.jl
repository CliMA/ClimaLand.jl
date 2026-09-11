# # Boundary conditions for the soil model

# In general, you must supply two boundary conditions for each PDE being
# solved. These are passed to the model as a NamedTuple
# of the form `(; top = top_bc, bottom = bottom_bc)`, where both `top_bc`
# and `bottom_bc` are of type `ClimaLand.AbstractBC`.

# Flux boundary conditions are always passed as the (scalar)
# z-component of the flux `f`, i.e. F⃗ = f ẑ.

# # Boundary conditions for Richards equation
# 1. `FreeDrainage <: AbstractWaterBC`: this only can be used at the bottom of the domain.

# 2. `WaterFluxBC <: AbstractWaterBC`: this accepts a prescribed analytic function of the cache `p` and simulation time `t` which returns the flux value. e.g: `WaterFluxBC((p,t) -> 0.0)`.

# 3. `MoistureStateBC <: AbstractWaterBC`: this accepts a prescribed analytic function of the cache `p` and simulation time `t` which returns the value of `ϑ_l` at the boundary . e.g: `MoistureStateBC((p,t) -> 0.2)`.

# 4. `RichardsAtmosDrivenFluxBC <: AbstractWaterBC`: this requires a single argument of abstract type `AbstractTimeVaryingInput`. Under the hood, this specifies the precipitation as a function of space and time (using data read in from a file, or an analytic function) and applies this a flux BC, optionally accounting for surface/subsurface runoff.

# # Boundary conditions for the soil heat equation

# 1. `HeatFluxBC <: AbstractHeatBC`: this accepts a prescribed analytic function of the cache `p` and simulation time `t` which returns the flux value. e.g: `HeatFluxBC((p,t) -> 0.0)`.

# 2. `TemperatureStateBC <: AbstractHeatBC`: this accepts a prescribed analytic function of the cache `p` and simulation time `t` which returns the value of `T` at the boundary . e.g: `TemperatureStateBC((p,t) -> 273.15)`.

# # Boundary conditions for the soil heat + water equations (EnergyHydrology model)
# The full soil model requires boundary conditions for both Richards equation and the soil heat equation.

# 1. `WaterHeatBC <: AbstractEnergyHydrologyBC`: In many cases, the two boundary conditions can be treated independently. The `WaterHeatBC` requires a boundary condition of abstract type `AbstractWaterBC` and one of type `AbstractHeatBC`, for example, `top = WaterHeatBC(; water = MoistureBC(ϑ_l(p,t)), heat = TemperatureBC(T(p,t)))`.

# 2. `AtmosDrivenFluxBC <: AbstractEnergyHydrologyBC`: This is an example of a set of boundary conditions for the full soil model which cannot be decomposed into two independent boundary conditions for water and heat. In this case, we compute the turbulent surface fluxes with the atmosphere, obtaining a sensible heat, latent heat, and water vapor flux. We also take into account the net radiation at the surface and any precipitation or runoff. This is the BC type used in most land simulations.
