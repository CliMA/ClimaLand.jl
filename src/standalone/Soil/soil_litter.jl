"""
    AbstractSoilLitter

The abstract type for the litter (dead leaves, needles and twigs) lying on the soil
surface, which sits between the soil skin and the top soil layer.
"""
abstract type AbstractSoilLitter end

"""
    NoLitter <: AbstractSoilLitter

A bare soil surface: the skin is connected to the top soil layer by the half-cell
conduction resistance of the soil only.
"""
struct NoLitter <: AbstractSoilLitter end

"""
    SoilLitter{FT, F} <: AbstractSoilLitter

A litter layer on the soil surface, which insulates the soil from the skin. Its heat
capacity is neglected (that of 3 cm of litter is about 3% of that of a 5 cm soil layer),
so it adds the thermal resistance `thickness / κ` in series with the half-cell
resistance of the top soil layer between the skin and the top layer (Napoly et al.,
2017, Eq. A10). The skin, which exchanges radiation and turbulent fluxes with the air,
is then the top of the litter; the litter does not change the soil evaporation.

Napoly, A., et al. (2017), The interactions between soil–biosphere–atmosphere (ISBA)
land surface model multi-energy balance (MEB) option in SURFEXv8 – Part 2: Introduction
of a litter formulation and model evaluation for local-scale forest sites, Geosci. Model
Dev., 10, 1621–1644, https://doi.org/10.5194/gmd-10-1621-2017.

$(DocStringExtensions.FIELDS)
"""
struct SoilLitter{
    FT <: AbstractFloat,
    F <: Union{FT, ClimaCore.Fields.Field},
} <: AbstractSoilLitter
    "The thickness of the litter layer (m); zero where there is no litter"
    thickness::F
    "The thermal conductivity of the litter (W/m/K)"
    κ::FT
end

SoilLitter{FT}(thickness::F, κ) where {FT <: AbstractFloat, F} =
    SoilLitter{FT, F}(thickness, FT(κ))

"""
    skin_conduction_resistance(litter::AbstractSoilLitter, Δz_top, κ_top)

Return the thermal resistance (m² K/W) between the soil skin and the center of the top
soil layer, at a distance `Δz_top` (m) below the soil surface, with the thermal
conductivity `κ_top` (W/m/K) of the top layer: the half-cell resistance of the soil,
plus that of the `litter` on the surface, if any.
"""
skin_conduction_resistance(::NoLitter, Δz_top, κ_top) = @. lazy(Δz_top / κ_top)
function skin_conduction_resistance(litter::SoilLitter, Δz_top, κ_top)
    (; thickness, κ) = litter
    return @. lazy(Δz_top / κ_top + thickness / κ)
end
