module Snow

using DocStringExtensions
import ...Parameters as LP
export SnowParameters

"""
    SnowParameters{FT <: AbstractFloat, PSE}
A struct for storing parameters of the `SnowModel`.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SnowParameters{FT <: AbstractFloat, PSE}
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

include("./snow_parameterizations.jl")

end
