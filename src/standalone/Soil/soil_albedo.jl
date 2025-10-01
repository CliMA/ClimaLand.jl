export update_albedo!, CLMTwoBandSoilAlbedo, ConstantTwoBandSoilAlbedo

abstract type AbstractSoilAlbedoParameterization end

"""
    ConstantTwoBandSoilAlbedo{
        SF <: Union{AbstractFloat, ClimaCore.Fields.Field},
    } <: AbstractSoilAlbedoParameterization

A parameterization for soil albedo: the soil albedo is 
defined in two bands (PAR and NIR), can spatially vary or
be set to scalar, but is assumed not to vary in time (or
with water content).
"""
struct ConstantTwoBandSoilAlbedo{
    SF <: Union{AbstractFloat, ClimaCore.Fields.Field},
} <: AbstractSoilAlbedoParameterization
    "Soil PAR Albedo"
    PAR_albedo::SF
    "Soil NIR Albedo"
    NIR_albedo::SF
end

function ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo::FT = FT(0.2),
    NIR_albedo::FT = FT(0.4),
) where {FT}
    ConstantTwoBandSoilAlbedo{FT}(PAR_albedo, NIR_albedo)
end


"""
    update_albedo!(bc::AtmosDrivenFluxBC, albedo::ConstantTwoBandSoilAlbedo, p, soil_domain, model_parameters)

Updates PAR and NIR albedo using the temporally-constant parameters 
provided in `albedo`; these values may be spatially varying.

For the temporally-constant albedo model, there is no need to update 
the values each step, or have an allocated spot in the cache for them. This
can be optimized in the future.
"""
function update_albedo!(
    bc::AtmosDrivenFluxBC,
    albedo::ConstantTwoBandSoilAlbedo,
    p,
    soil_domain,
    model_parameters,
)
    @. p.soil.PAR_albedo = albedo.PAR_albedo
    @. p.soil.NIR_albedo = albedo.NIR_albedo
end

"""
     CLMTwoBandSoilAlbedo{
        FT <: AbstractFloat,
        SF <: Union{FT, ClimaCore.Fields.Field},
    } <: AbstractSoilAlbedoParameterization

A parameterization for soil albedo: the soil albedo is
defined in two bands (PAR and NIR), and can spatially vary or
be set to scalar. However, it varies temporally due to a
dependence on soil water content at the surface,
via the effective saturation S(θ\\_sfc):
α = α\\_wet*S + α\\_dry*(1-S)

We use a value for θ\\_sfc averaged over the depth `albedo_calc_top_thickness`.
If the model resolution is such that the first layer is thicker than this depth,
the value from the first layer is used.

CLM reference: Lawrence, P.J., and Chase, T.N. 2007. Representing a MODIS consistent land surface in the Community Land Model
(CLM 3.0). J. Geophys. Res. 112:G01023. DOI:10.1029/2006JG000168.
"""
struct CLMTwoBandSoilAlbedo{
    FT <: AbstractFloat,
    SF <: Union{FT, ClimaCore.Fields.Field},
} <: AbstractSoilAlbedoParameterization
    "Soil PAR Albedo dry"
    PAR_albedo_dry::SF
    "Soil NIR Albedo dry"
    NIR_albedo_dry::SF
    "Soil PAR Albedo wet"
    PAR_albedo_wet::SF
    "Soil NIR Albedo wet"
    NIR_albedo_wet::SF
    "Thickness of top of soil used in albedo calculations (m)"
    albedo_calc_top_thickness::FT
end

function CLMTwoBandSoilAlbedo{FT}(;
    PAR_albedo_dry,
    NIR_albedo_dry,
    PAR_albedo_wet,
    NIR_albedo_wet,
    albedo_calc_top_thickness = FT(0.02),
) where {FT}
    return CLMTwoBandSoilAlbedo{FT, typeof(PAR_albedo_dry)}(
        PAR_albedo_dry,
        NIR_albedo_dry,
        PAR_albedo_wet,
        NIR_albedo_wet,
        albedo_calc_top_thickness,
    )
end


"""
    albedo_from_moisture(S_sfc::FT, θ_int::FT, albedo_wet::FT)

Calculates pointwise albedo for any band as a function of soil surface moisture given
the dry and wet albedo values for that band using the CLM parameterization.

CLM reference: Lawrence, P.J., and Chase, T.N. 2007. Representing a MODIS consistent land surface in the Community Land Model
(CLM 3.0). J. Geophys. Res. 112:G01023. DOI:10.1029/2006JG000168.
"""
function albedo_from_moisture(
    S_sfc::FT,
    albedo_dry::FT,
    albedo_wet::FT,
) where {FT}
    return (1 - S_sfc) * albedo_dry + S_sfc * albedo_wet
end


"""
    update_albedo!(bc::AtmosDrivenFluxBC, albedo::CLMTwoBandSoilAlbedo, p, soil_domain, model_parameters)

Calculates and updates PAR and NIR albedo as a function of volumetric soil water content at
the top of the soil. If the soil layers are larger than the specified `albedo_calc_top_thickness`,
the water content of the top layer is used in the calclulation. For the PAR and NIR bands,

α\\_band = α\\_{band,dry} * (1 - S\\_e) +  α\\_{band,wet} * (S\\_e)

where S\\_e is the relative soil wetness above some depth, `albedo_calc_top_thickness`. This
is a modified version of Equation (1) of:

Braghiere, R. K., Wang, Y., Gagné-Landmann, A., Brodrick, P. G., Bloom, A. A., Norton,
A. J., et al. (2023). The importance of hyperspectral soil albedo information for improving
Earth system model projections. AGU Advances, 4, e2023AV000910. https://doi.org/10.1029/2023AV000910

where effective saturation is used in place of volumetric soil water content. The dry and wet
albedo values come from a global soil color map and soil color to albedo map from CLM.

CLM reference: Lawrence, P.J., and Chase, T.N. 2007. Representing a MODIS consistent land surface in the Community Land Model
(CLM 3.0). J. Geophys. Res. 112:G01023. DOI:10.1029/2006JG000168.
"""
function update_albedo!(
    bc::AtmosDrivenFluxBC,
    albedo::CLMTwoBandSoilAlbedo,
    p,
    soil_domain,
    model_parameters,
)
    (;
        PAR_albedo_dry,
        NIR_albedo_dry,
        PAR_albedo_wet,
        NIR_albedo_wet,
        albedo_calc_top_thickness,
    ) = albedo
    FT = eltype(soil_domain.fields.Δz_top)
    # checks if there is at least 1 layer centered within the top soil depth
    if soil_domain.fields.Δz_min < albedo_calc_top_thickness
        # We compute ∫H_θ_dz / N, where N =∫H_dz
        # is a normalization,
        # and then compute θ_sfc = ∫[(H θ)/N]dz.
        # N is integral of 1 from (surface-albedo_calc_top_thickness) to surface
        N = p.soil.sfc_scratch
        # zero all centers lower than boundary, set everything above to one
        @. p.soil.sub_sfc_scratch = ClimaLand.heaviside(
            albedo_calc_top_thickness + sqrt(eps(FT)),
            soil_domain.fields.z_sfc - soil_domain.fields.z,
        )
        ClimaCore.Operators.column_integral_definite!(N, p.soil.sub_sfc_scratch)
        # zeros all effective saturation at levels centered lower than boundary
        @. p.soil.sub_sfc_scratch =
            ClimaLand.heaviside(
                albedo_calc_top_thickness + sqrt(eps(FT)),
                soil_domain.fields.z_sfc - soil_domain.fields.z,
            ) * p.soil.θ_l / N
        # Now use scratch space to compute normed integral
        θ_sfc = p.soil.sfc_scratch
        ClimaCore.Operators.column_integral_definite!(
            θ_sfc,
            p.soil.sub_sfc_scratch,
        )
    else
        # in the case where no layer is centered above boundary, use the values of the top layer
        θ_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.θ_l)
    end
    ν_sfc = ClimaLand.Domains.top_center_to_surface(model_parameters.ν)
    θ_r_sfc = ClimaLand.Domains.top_center_to_surface(model_parameters.θ_r)
    S_sfc = @. lazy(effective_saturation(ν_sfc, θ_sfc, θ_r_sfc))
    @. p.soil.PAR_albedo =
        albedo_from_moisture(S_sfc, PAR_albedo_dry, PAR_albedo_wet)
    @. p.soil.NIR_albedo =
        albedo_from_moisture(S_sfc, NIR_albedo_dry, NIR_albedo_wet)
end

"""
    update_albedo!(bc::AbstractEnergyHydrologyBC, _...)

Does nothing for boundary conditions where albedo is not used.
"""
function update_albedo!(bc::AbstractEnergyHydrologyBC, _...) end


"""
    ClimaLand.surface_albedo(
        model::EnergyHydrology{FT},
        Y,
        p,
    ) where {FT}

Returns the surface albedo field of the
`EnergyHydrology` soil model.
"""
function ClimaLand.surface_albedo(model::EnergyHydrology{FT}, Y, p) where {FT}
    return @. lazy((p.soil.PAR_albedo + p.soil.NIR_albedo) / 2)
end
