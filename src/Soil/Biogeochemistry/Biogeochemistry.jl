module Biogeochemistry
using ClimaLSM
using DocStringExtensions
export DETECTModelParameters
"""
    DETECTParameters{FT <: AbstractFloat}
A struct for storing parameters of the `DETECTModel`.
"""
struct DETECTModelParameters{FT <: AbstractFloat}
    "Pressure at the surface of the soil (Pa)"
    P_sfc::FT
    "Root mass-base respiration rate at 10°C and mean environmental conditions (kg C m⁻³ s⁻¹)"
    Rb::FT
    "The effect of soil water content (θ) on root respiration (unitless)"
    α1r::FT
    "The effect of antecedent θ on root respiration (unitless)"
    α2r::FT
    "The interactive effect of θ and antecedent θ on root respiration (unitless)"
    α3r::FT
    "Value of Vₘₐₓ at 10°C and mean environmental conditions (kg C m⁻³ s⁻¹)"
    Vb::FT
    "The effect of soil water content (θ) on microbial respiration (unitless)"
    α1m::FT
    "The effect of antecedent θ on microbial respiration (unitless)"
    α2m::FT
    "The interactive effect of θ and antecedent θ on microbial respiration (unitless)"
    α3m::FT
    "Michaelis-Menten half-saturation constant (kg C m⁻³ s⁻¹)"
    Km::FT
    "Microbial carbon-use efficiency (kg C kg⁻¹ C⁻¹)"
    CUE::FT
    "Fraction of soil organic C that is soluble (-)"
    soluble_fraction::FT
    "Diffusivity of soil C substrate in liquid (unitless)"
    D_liq::FT
    "Temperature sensitivity parameter, somewhat analogous to an energy activation (Kelvin)"
    Estar::FT
    "Temperature sensitivity-related parameter (Kelvin)"
    T_ref::FT
    "The effect of antecedent soil temperature on root and microbial respiration (unitless)"
    α4::FT
    "Reference temperature (Kelvin)"
    T_ref_soil::FT
    "Absolute value of the slope of the line relating log(ψ) versus log(θ) (unitless)"
    α5::FT
    "Soil porosity (m³ m⁻³)"
    ν::FT
    "Air-filled porosity at soil water potential of -100 cm H₂O (~ 10 Pa)"
    θ_a100::FT
    "Diffusion coefficient for CO₂ in air at standard temperature and pressure (m² s⁻¹)"
    D_ref::FT
    "Standard pressure (Pa)"
    P_ref::FT
    "Parameter related to the pore size distribution of the soil (unitless)"
    b::FT
end

"""
    DETECTModelParameters(;
                           P_sfc::FT,
                           Rb::FT,
                           Vb::FT,
                           α1r::FT,
                           α2r::FT,
                           α3r::FT,
                           α1m::FT,
                           α2m::FT,
                           α3m::FT,
                           α4::FT,
                           α5::FT,
                           ν::FT,
                           b::FT, 
                           θ_a100::FT,
                           D_ref::FT,
                           P_ref::FT,
                           T_ref::FT,
                           T_ref_soil::FT
                           Km::FT,
                           CUE::FT,
                           Estar::FT,
                           D_liq::FT,
                           soluble_fraction::FT
                           ) where {FT}

An outer constructor for creating the parameter struct of the `DETECTModel`,
    based on keyword arguments.
"""
function DETECTModelParameters(;
    P_sfc::FT,
    Rb::FT,
    Vb::FT,
    α1r::FT,
    α2r::FT,
    α3r::FT,
    α1m::FT,
    α2m::FT,
    α3m::FT,
    α4::FT,
    α5::FT,
    ν::FT,
    b::FT,
    θ_a100::FT,
    D_ref::FT,
    P_ref::FT,
    T_ref::FT,
    T_ref_soil::FT,
    Km::FT,
    CUE::FT,
    Estar::FT,
    D_liq::FT,
    soluble_fraction::FT,
) where {FT}
    return DETECTModelParameters{FT}(
        P_sfc,
        Rb,
        α1r,
        α2r,
        α3r,
        Vb,
        α1m,
        α2m,
        α3m,
        Km,
        CUE,
        soluble_fraction,
        D_liq,
        Estar,
        T_ref,
        α4,
        T_ref_soil,
        α5,
        ν,
        θ_a100,
        D_ref,
        P_ref,
        b,
    )
end

include("./co2_parameterizations.jl")

end # module
