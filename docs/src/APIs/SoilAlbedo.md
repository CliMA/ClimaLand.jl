## Parameterizing Soil Albedo

# We currently support two parameterizations for the soil albedo. These
# options are defined using Julia types, using simplified code, as follows.
# We first introduce the abstract type `AbstractSoilAlbedoParameterization`,
# and then the two albedo parameterizations are concrete examples of that
# abstract type:

```julia
ConstantTwoBandSoilAlbedo <: AbstractSoilAlbedoParameterization
CLMTwoBandSoilAlbedo <: AbstractSoilAlbedoParameterization
```

# Each of these types define a method for `update_albedo!`, which is
# called in the soils tendency functions, and each stores the
# parameters required to compute the albedo for that parameterizations.
# For example, the parameterization based off of CLM's approach requires
# the wet and dry soil albedos in the PAR and NIR wavelength bands. For
# each band, the actual albedo is a linear combination of the wet and dry
# values, with a weight depending on the soil moisture. The method
# of update_albedo! for this parameterization computes this linear
# combination using a helper function `albedo_from_moisture` and ultimately
# sets

```julia
@. p.soil.PAR_albedo =
    albedo_from_moisture(S_sfc, PAR_albedo_dry, PAR_albedo_wet)
@. p.soil.NIR_albedo =
    albedo_from_moisture(S_sfc, NIR_albedo_dry, NIR_albedo_wet)
````
where `S_sfc` is the effective saturation at the surface.

### Creating a new albedo parameterization

# Suppose you want to define a new parameterization which models albedo
# as a linear combination of the albedos of quartz, organic matter, minerals
# and water, treating NIR and PAR albedos the same.

# First, create the type:

```julia
struct SoilAlbedoFromComposition{FT <: AbstractFloat} <: AbstractSoilAlbedoParameterization
    α_quartz::FT
    α_minerals::FT
    α_om::FT
    α_water::FT
end
```

# And then create the method. For now, don't worry about the other arguments
# and their types:

```julia
function update_albedo!(
    bc::AtmosDrivenFluxBC,
    albedo::SoilAlbedoFromComposition,
    p,
    soil_domain,
    model_parameters,
    )
    # unpack parameters of the albedo model
    (; α_quartz, α_minerals, α_om, α_water) = albedo
    # unpack composition parameters from the soil parameters
    (; ν_ss_om, ν_ss_quartz) = model_parameters
    # Extract the values at the top level
    ν_ss_om_sfc = ClimaLand.Domains.top_center_to_surface(ν_ss_om)
    ν_ss_quartz_sfc = ClimaLand.Domains.top_center_to_surface(ν_ss_quartz)
    # Extract the top layer's volumetric water content
    θ_l_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.θ_l)

    # Compute the linear combination
    @. p.soil.PAR_albedo = α_water * θ_l_sfc + (1-θ_l_sfc)* (α_quartz * ν_ss_quartz_sfc + α_om * ν_ss_om_sfc + α_minerals * (1-ν_ss_om - ν_ss_quartz))
    p.soil.NIR_albedo .= p.soil.PAR_albedo 
end
```