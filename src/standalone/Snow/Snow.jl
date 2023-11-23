module Snow

using DocStringExtensions
import ...Parameters as LSMP
export SnowParameters

"""
    SnowParameters{FT <: AbstractFloat, PSE}
A struct for storing parameters of the `SnowModel`.
$(DocStringExtensions.FIELDS)
"""
struct SnowParameters{FT <: AbstractFloat, PSE}
    "Density of snow (kg/m^3)"
    ρ_snow::FT
    "Roughness length over snow for momentum (m)"
    z_0m::FT
    "Roughness length over snow for scalars (m)"
    z_0b::FT
    "Albedo of snow (unitless)"
    α_snow::FT
    "Emissivity of snow (unitless)"
    ϵ_snow::FT
    "Volumetric holding capacity of water in snow (unitless)"
    θ_r::FT
    "Hydraulic conductivity of wet snow (m/s)"
    Ksat::FT
    "Critical threshold for snow cover fraction (m)"
    fS_c::FT
    "Thermal conductivity of ice (W/m/K)"
    κ_ice::FT
    "Timestep of the model"
    Δt::Any
    "Clima-wide parameters"
    earth_param_set::PSE
end

"""
   SnowParameters{FT}(Δt;
                      ρ_snow = FT(200),
                      z_0m = FT(0.0024),
                      z_0b = FT(0.00024),
                      α_snow = FT(0.8),
                      ϵ_snow = FT(0.99),
                      θ_r = FT(0.08),
                      Ksat = FT(1e-3),
                      fS_c = FT(0.2),
                      κ_ice = FT(2.21),
                      earth_param_set::PSE) where {FT, PSE}

An outer constructor for `SnowParameters` which supplies defaults for
all arguments but `earth_param_set`.
"""
function SnowParameters{FT}(
    Δt;
    ρ_snow = FT(200),
    z_0m = FT(0.0024),
    z_0b = FT(0.00024),
    α_snow = FT(0.8),
    ϵ_snow = FT(0.99),
    θ_r = FT(0.08),
    Ksat = FT(1e-3),
    fS_c = FT(0.2),
    κ_ice = FT(2.21),
    earth_param_set::PSE,
) where {FT <: AbstractFloat, PSE}
    return SnowParameters{FT, PSE}(
        ρ_snow,
        z_0m,
        z_0b,
        α_snow,
        ϵ_snow,
        θ_r,
        Ksat,
        fS_c,
        κ_ice,
        Δt,
        earth_param_set,
    )
end


include("./snow_parameterizations.jl")

end
