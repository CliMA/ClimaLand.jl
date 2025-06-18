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

"""
    update_albedo!(bc::AtmosDrivenFluxBC, albedo::ConstantTwoBandSoilAlbedo, p, soil_domain, model_parameters)

Updates PAR and NIR albedo using the constant parameters provided in `albedo`.
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
defined in two bands (PAR and NIR), can spatially vary or
be set to scalar, and varies with water content at the surface
θ_sfc, according to CLM:

α = min(α_wet + Δ(θ_sfc), α_dry), where
Δ(θ_sfc) = max(0.11 - 0.40*θ_sfc,0).

We use a value for θ_sfc averaged over the depth `albedo_calc_top_thickness`.
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
    albedo_from_moisture(θ_sfc::FT, albedo_dry::FT, albedo_wet::FT)

Calculates pointwise albedo for any band as a function of soil surface moisture given
the dry and wet albedo values for that band.

CLM reference: Lawrence, P.J., and Chase, T.N. 2007. Representing a MODIS consistent land surface in the Community Land Model
(CLM 3.0). J. Geophys. Res. 112:G01023. DOI:10.1029/2006JG000168.
"""
function albedo_from_moisture(
    θ_sfc::FT,
    albedo_dry::FT,
    albedo_wet::FT,
) where {FT}
    Δ = max(FT(0.11) - FT(0.4) * θ_sfc, FT(0))
    return min(albedo_wet + Δ, albedo_dry)
end


"""
    update_albedo!(bc::AtmosDrivenFluxBC, albedo::CLMTwoBandSoilAlbedo, p, soil_domain, model_parameters)

Calculates and updates PAR and NIR albedo as a function of volumetric soil water content at
the top of the soil. If the soil layers are larger than the specified `albedo_calc_top_thickness`,
the water content of the top layer is used in the calclulation. .The dry and wet
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
        # and then computing ∫[(H θ)/N]dz.
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
        normed_∫H_θ_dz = p.soil.sfc_scratch
        ClimaCore.Operators.column_integral_definite!(
            normed_∫H_θ_dz,
            p.soil.sub_sfc_scratch,
        )
        θ_sfc = normed_∫H_θ_dz
    else
        # in the case where no layer is centered above boundary, use the values of the top layer
        θ_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.θ_l)
    end
    @. p.soil.PAR_albedo =
        albedo_from_moisture(θ_sfc, PAR_albedo_dry, PAR_albedo_wet)
    @. p.soil.NIR_albedo =
        albedo_from_moisture(θ_sfc, NIR_albedo_dry, NIR_albedo_wet)
end

"""
update_albedo!(bc::AbstractWaterBC, _, p, soil_domain, model_parameters)

    Does nothing for boundary conditions where albedo is not used.
"""
function update_albedo!(
    bc::BC,
    _,
    p,
    soil_domain,
    model_parameters,
) where {BC <: AbstractEnergyHydrologyBC} end

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
