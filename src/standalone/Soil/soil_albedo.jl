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
dependence on soil water content,
θ_sfc, according to CLM:

α = min(α_wet + Δ(θ_sfc), α_dry), where
Δ(θ_sfc) = max(θ_int - dαdθ*θ_sfc,0),
where θ_int = 0.11 and dαdθ = 0.4.

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
    "Slope of albedo vs moisture curve, dαdθ (unitless)"
    dαdθ::FT
    "Intercept parameter of slope of albedo vs moisture curve, θ_int (unitless)"
    θ_int::FT
    "Thickness of top of soil used in albedo calculations (m)"
    albedo_calc_top_thickness::FT
end

function CLMTwoBandSoilAlbedo{FT}(;
    PAR_albedo_dry,
    NIR_albedo_dry,
    PAR_albedo_wet,
    NIR_albedo_wet,
    dαdθ = FT(0.4),
    θ_int = FT(0.11),
    albedo_calc_top_thickness = FT(0.02),
) where {FT}
    return CLMTwoBandSoilAlbedo{FT, typeof(PAR_albedo_dry)}(
        PAR_albedo_dry,
        NIR_albedo_dry,
        PAR_albedo_wet,
        NIR_albedo_wet,
        dαdθ,
        θ_int,
        albedo_calc_top_thickness,
    )
end


"""
    albedo_from_moisture(θ_sfc::FT, θ_int::FT,dαdθ::FT,  albedo_dry::FT, albedo_wet::FT)

Calculates pointwise albedo for any band as a function of soil surface moisture given
the dry and wet albedo values for that band using the CLM parameterization.

CLM reference: Lawrence, P.J., and Chase, T.N. 2007. Representing a MODIS consistent land surface in the Community Land Model
(CLM 3.0). J. Geophys. Res. 112:G01023. DOI:10.1029/2006JG000168.
"""
function albedo_from_moisture(
    θ_sfc::FT,
    θ_int::FT,
    dαdθ::FT,
    albedo_dry::FT,
    albedo_wet::FT,
) where {FT}
    Δ = max(θ_int - dαdθ * θ_sfc, FT(0))
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
        dαdθ,
        θ_int,
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
    @. p.soil.PAR_albedo =
        albedo_from_moisture(θ_sfc, θ_int, dαdθ, PAR_albedo_dry, PAR_albedo_wet)
    @. p.soil.NIR_albedo =
        albedo_from_moisture(θ_sfc, θ_int, dαdθ, NIR_albedo_dry, NIR_albedo_wet)
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
