export SoilCarbonLitterInput

"""
    SoilCarbonLitterInput{FT} <: Soil.Biogeochemistry.AbstractCarbonSource{FT}

Couples soil organic carbon to the carbon pools of a `Canopy.PrognosticCarbonModel`:
`dSOC/dt = I_litter - Sm`, the litter shed by the pools, `p.soil_litter_input`
(kg C m^-3 s^-1, set by `update_soil_litter_input!`), minus the microbial respiration
`Sm` of `MicrobeProduction`. Without this source, SOC is held at its initial condition.
"""
struct SoilCarbonLitterInput{FT} <:
       Soil.Biogeochemistry.AbstractCarbonSource{FT} end

"""
    ClimaLand.source!(dY, src::SoilCarbonLitterInput, Y, p, params)

Adds the litter input minus the microbial respiration to the SOC tendency.
"""
NVTX.@annotate function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::SoilCarbonLitterInput,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    params,
)
    @. dY.soilco2.SOC += p.soil_litter_input - p.soilco2.Sm
    return nothing
end

"""
    soilco2_sources(biomass)

Sources of a soil CO2 model coupled to a canopy with the biomass model `biomass`:
microbial respiration, plus the litter input when the canopy carries carbon pools.
"""
soilco2_sources(::Canopy.AbstractBiomassModel{FT}) where {FT} =
    (Soil.Biogeochemistry.MicrobeProduction{FT}(),)
soilco2_sources(::Canopy.PrognosticCarbonModel{FT}) where {FT} =
    (Soil.Biogeochemistry.MicrobeProduction{FT}(), SoilCarbonLitterInput{FT}())

"""
    update_soil_litter_input!(p, Y, t, land)

Sets `p.soil_litter_input` (kg C m^-3 s^-1), the litter of the canopy carbon pools
distributed in the soil: leaf and stem litter on an exponential profile with e-folding
depth `soil_litter_depth`, root litter on the root distribution. Each profile is
normalized by its column integral, so the soil receives exactly the litter the pools
shed. Does nothing for a canopy without carbon pools.
"""
update_soil_litter_input!(p, Y, t, land) =
    update_soil_litter_input!(p, land, land.canopy.biomass)

update_soil_litter_input!(p, land, ::Canopy.AbstractBiomassModel) = nothing

function update_soil_litter_input!(
    p,
    land,
    biomass::Canopy.PrognosticCarbonModel,
)
    z = land.soil.domain.fields.z
    (; rooting_depth) = biomass
    (; soil_litter_depth) = biomass.parameters
    (; L_leaf, L_stem, L_root) = p.canopy.biomass.carbon
    surface_profile = @. lazy(Canopy.root_distribution(z, soil_litter_depth))
    root_profile = @. lazy(Canopy.root_distribution(z, rooting_depth))
    # Surface scratch fields, included in lsm_aux_vars
    surface_norm = p.scratch1
    root_norm = p.scratch2
    ClimaCore.Operators.column_integral_definite!(surface_norm, surface_profile)
    ClimaCore.Operators.column_integral_definite!(root_norm, root_profile)
    @. p.soil_litter_input =
        (L_leaf + L_stem) * surface_profile / surface_norm +
        L_root * root_profile / root_norm
    return nothing
end
