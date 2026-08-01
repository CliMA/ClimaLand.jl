export SoilCarbonLitterInput

"""
    update_soil_litter_input!(p, Y, t, land)

Updates `p.soil_litter_input` (kg C m^-3 s^-1) with the canopy litter reaching
the soil carbon pool, distributed through the column.

Leaf and stem litter enters near the surface and root litter follows the root
profile, both normalised by their own column integral so that the integral of
`soil_litter_input` returns the surface litter flux exactly, independent of the
vertical discretisation:

    I_litter(z) = (L_leaf + L_stem)·s(z)/∫s + L_root·r(z)/∫r

with `r` the root distribution and `s` the same exponential shape at the much
shallower `soil_litter_depth`. The specification places leaf and stem litter "in
the top layer"; a shallow exponential is the continuum form of that, and unlike
writing into a single cell it does not depend on the layer thickness or on
which end of the column is indexed first.

The litter fluxes themselves are already computed and cached by
`Canopy.update_carbon_fluxes!`, so this only redistributes them.
"""
update_soil_litter_input!(p, Y, t, land) =
    update_soil_litter_input!(p, Y, t, land, land.canopy.biomass)

# No carbon pools means no litter; the soil carbon pool then simply decays.
function update_soil_litter_input!(
    p,
    Y,
    t,
    land,
    biomass::Canopy.AbstractBiomassModel,
)
    p.soil_litter_input .= 0
    return nothing
end

function update_soil_litter_input!(
    p,
    Y,
    t,
    land,
    biomass::Canopy.PrognosticCarbonModel{FT},
) where {FT}
    z = land.soil.domain.fields.z
    rooting_depth = biomass.rooting_depth
    litter_depth = biomass.parameters.soil_litter_depth

    root_shape = @. lazy(Canopy.root_distribution(z, rooting_depth))
    surface_shape = @. lazy(Canopy.root_distribution(z, litter_depth))
    # Normalise each shape by its own column integral, so a truncated or coarse
    # column still receives the full litter flux rather than a fraction of it.
    root_norm = p.scratch1
    surface_norm = p.scratch2
    ClimaCore.Operators.column_integral_definite!(root_norm, root_shape)
    ClimaCore.Operators.column_integral_definite!(surface_norm, surface_shape)

    L_surface =
        @. lazy(p.canopy.biomass.carbon.L_leaf + p.canopy.biomass.carbon.L_stem)
    L_root = p.canopy.biomass.carbon.L_root
    @. p.soil_litter_input =
        L_surface * Canopy.root_distribution(z, litter_depth) /
        max(surface_norm, eps(FT)) +
        L_root * Canopy.root_distribution(z, rooting_depth) /
        max(root_norm, eps(FT))
    return nothing
end

"""
    SoilCarbonLitterInput{FT} <: Soil.Biogeochemistry.AbstractCarbonSource{FT}

Adds canopy litter to the soil organic carbon pool.

Together with the `Sm` debit in `MicrobeProduction` this closes the soil carbon
balance: `dSOC/dt = I_litter(z) - Sm`. Without the litter input SOC only decays,
which is why the two belong together.

The source reads `p.soil_litter_input`, filled by `update_soil_litter_input!`,
in the same way `RootExtraction` reads `p.root_extraction`. When the canopy
carries no carbon pools that field is zero, so including this source is
harmless.
"""
struct SoilCarbonLitterInput{FT} <:
       Soil.Biogeochemistry.AbstractCarbonSource{FT} end

"""
    ClimaLand.source!(dY, src::SoilCarbonLitterInput, Y, p, params)

Adds the litter input to the soil organic carbon tendency.
"""
NVTX.@annotate function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::SoilCarbonLitterInput,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    params,
)
    dY.soilco2.SOC .+= p.soil_litter_input
    return nothing
end
