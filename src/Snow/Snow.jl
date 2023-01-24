module Snow

using UnPack
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
    z0_m::FT
    "Roughness length over snow for scalars (m)"
    z0_b::FT
    "Displacement height over snow (m)"
    d::FT
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
    Δt::FT
    "Clima-wide parameters"
    earth_param_set::PSE
end

"""
   SnowParameters{FT}(Δt::FT;
                      ρ_snow = FT(200),
                      z0_m = FT(0.0024),
                      z0_b = FT(0.00024),
                      d = FT(0),
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
    Δt::FT;
    ρ_snow = FT(200),
    z0_m = FT(0.0024),
    z0_b = FT(0.00024),
    d = FT(0),
    α_snow = FT(0.8),
    ϵ_snow = FT(0.99),
    θ_r = FT(0.08),
    Ksat = FT(1e-3),
    fS_c = FT(0.2),
    κ_ice = FT(2.21),
    earth_param_set::PSE,
) where {FT, PSE}
    return SnowParameters{FT, PSE}(
        ρ_snow,
        z0_m,
        z0_b,
        d,
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
