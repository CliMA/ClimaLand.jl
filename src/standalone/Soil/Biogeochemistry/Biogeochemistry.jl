module Biogeochemistry
using ClimaLand
import ClimaParams as CP
using LinearAlgebra
using DocStringExtensions
using ClimaCore
using LazyBroadcast: lazy
import ...Parameters as LP
using ClimaCore.MatrixFields
import ClimaCore.MatrixFields: @name, ⋅
using NVTX
import ClimaCore: Fields, Operators, Geometry, Spaces

import ClimaLand.Domains: AbstractDomain
import ClimaLand:
    AbstractExpModel,
    make_update_aux,
    make_update_implicit_aux,
    make_compute_exp_tendency,
    make_compute_imp_tendency,
    make_update_boundary_fluxes,
    make_update_implicit_boundary_fluxes,
    make_compute_jacobian,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    auxiliary_types,
    prognostic_domain_names,
    auxiliary_domain_names,
    TopBoundary,
    BottomBoundary,
    AbstractBC,
    boundary_flux!,
    AbstractSource,
    source!,
    boundary_vars,
    boundary_var_types,
    boundary_var_domain_names
export SoilCO2ModelParameters,
    SoilCO2Model,
    PrescribedMet,
    MicrobeProduction,
    SoilCO2FluxBC,
    SoilO2FluxBC,
    AtmosCO2StateBC,
    AtmosO2StateBC,
    SoilCO2StateBC,
    AbstractSoilDriver,
    SoilDrivers

"""
    SoilCO2ModelParameters{FT <: AbstractFloat, PSE}

A struct for storing parameters of the `SoilCO2Model`.

All of these parameters are currently treated as global constants.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SoilCO2ModelParameters{FT <: AbstractFloat, PSE}
    "Diffusion coefficient for CO₂ in air at standard temperature and pressure (m² s⁻¹)"
    D_ref::FT
    "Diffusion coefficient for O₂ in air at standard temperature and pressure (m² s⁻¹)"
    D_ref_o2::FT
    "Diffusivity of soil C substrate in liquid (unitless)"
    D_liq::FT
    # DAMM — centered Arrhenius: Vmax = V_ref_sx * exp(-Ea_sx/R * (1/T - 1/T_ref_sx))
    "Maximum respiration rate at T_ref_sx (kg C m-3 s-1)"
    V_ref_sx::FT
    "Reference temperature for the centered Arrhenius form (K)"
    T_ref_sx::FT
    "Activation energy (J mol-1)"
    Ea_sx::FT
    "Michaelis constant (kg C m-3)"
    kM_sx::FT
    "Michaelis constant for O2 (m3 m-3)"
    kM_o2::FT
    "Volumetric fraction of O₂ in atmospheric air (reference value for boundary condition), dimensionless"
    O2_f_atm::FT
    "Diffusion coefficient of oxygen in air, dimensionless"
    D_oa::FT
    "Fraction of soil carbon that is considered soluble, dimensionless"
    p_sx::FT
    "Molar mass of carbon (kg/mol)"
    M_C::FT
    "Molar mass of oxygen (kg/mol)"
    M_O2::FT
    # Henry's law parameters for air-water partitioning (Sander, 2015)
    "Henry's law constant for CO2 at 298K (mol/(m³·Pa))"
    K_H_co2_298::FT
    "Temperature coefficient for CO2 Henry's law (K)"
    dln_K_H_co2_dT::FT
    "Henry's law constant for O2 at 298K (mol/(m³·Pa))"
    K_H_o2_298::FT
    "Temperature coefficient for O2 Henry's law (K)"
    dln_K_H_o2_dT::FT
    "Reference temperature for Henry's law (K), Sander (2015)"
    T_ref_henry::FT
    "Temperature exponent for free-air gas diffusivity correction (dimensionless), Ryan et al. (2018)"
    T_exp_diffusivity::FT
    "Physical constants used Clima-wide"
    earth_param_set::PSE
end

## For interfacing with ClimaParams


"""
    SoilCO2ModelParameters(toml_dict::CP.ParamDict)

SoilCO2ModelParameters provides a constructor using the TOML dict.
Keywords arguments can be used to directly override any parameters.
"""
function SoilCO2ModelParameters(
    toml_dict::CP.ParamDict;
    D_ref = toml_dict["CO2_diffusion_coefficient"],
    D_ref_o2 = toml_dict["O2_diffusion_coefficient_air"],
)
    name_map = (;
        :soil_C_substrate_diffusivity => :D_liq,
        :soilCO2_reference_rate => :V_ref_sx,
        :soilCO2_reference_temperature => :T_ref_sx,
        :soilCO2_activation_energy => :Ea_sx,
        :michaelis_constant => :kM_sx,
        :O2_michaelis_constant => :kM_o2,
        :O2_volume_fraction => :O2_f_atm,
        :oxygen_diffusion_coefficient => :D_oa,
        :soluble_soil_carbon_fraction => :p_sx,
        :molar_mass_carbon => :M_C,
        :molar_mass_oxygen => :M_O2,
        # Henry's law constants from Sander (2015), Atmos. Chem. Phys., 15, 4399-4981
        :CO2_henry_k298 => :K_H_co2_298,
        :CO2_henry_Tcoeff => :dln_K_H_co2_dT,
        :O2_henry_k298 => :K_H_o2_298,
        :O2_henry_Tcoeff => :dln_K_H_o2_dT,
        :henry_reference_temperature => :T_ref_henry,
        :gas_diffusivity_temperature_exponent => :T_exp_diffusivity,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    earth_param_set = LP.LandParameters(toml_dict)
    return SoilCO2ModelParameters{FT, typeof(earth_param_set)}(;
        earth_param_set,
        D_ref,
        D_ref_o2,
        parameters...,
    )
end

"""
    AbstractSoilBiogeochemistryModel{FT} <: ClimaLand.AbstractImExModel{FT}

An abstract model type for soil biogeochemistry models.
"""
abstract type AbstractSoilBiogeochemistryModel{FT} <:
              ClimaLand.AbstractImExModel{FT} end

"""
    SoilCO2Model

A model for simulating the production and transport of CO₂ and O₂ in the soil with dynamic
source and diffusion terms.
$(DocStringExtensions.FIELDS)
"""
struct SoilCO2Model{FT, PS, D, BC, S, DT} <:
       AbstractSoilBiogeochemistryModel{FT}
    "the parameter set"
    parameters::PS
    "the soil domain, using ClimaCore.Domains"
    domain::D
    "the boundary conditions, of type NamedTuple"
    boundary_conditions::BC
    "A tuple of sources, each of type AbstractSource"
    sources::S
    "Drivers"
    drivers::DT
end


"""
SoilCO2Model{FT}(
        domain::ClimaLand.AbstractDomain,
        drivers::DT,
        toml_dict::CP.ParamDict;
        parameters::SoilCO2ModelParameters{FT} = SoilCO2ModelParameters(toml_dict),
        sources::Tuple = (MicrobeProduction{FT}(),),
    ) where {FT, DT}

A constructor for `SoilCO2Model`.
Defaults are provided for the parameters and sources.
These can be overridden by providing the appropriate keyword arguments.

Boundary conditions are set automatically:
- Top: AtmosCO2StateBC(earth_param_set, M_C) for CO2, AtmosO2StateBC(earth_param_set, M_O2, O2_f_atm) for O2
- Bottom: Zero flux for both CO2 and O2
"""
function SoilCO2Model{FT}(
    domain::ClimaLand.AbstractDomain,
    drivers::DT,
    toml_dict::CP.ParamDict;
    parameters::SoilCO2ModelParameters{FT} = SoilCO2ModelParameters(toml_dict),
    sources::Tuple = (MicrobeProduction{FT}(),),
) where {FT, DT}
    boundary_conditions = (
        top = (
            co2 = AtmosCO2StateBC(parameters.earth_param_set, parameters.M_C),
            o2 = AtmosO2StateBC(
                parameters.earth_param_set,
                parameters.M_O2,
                parameters.O2_f_atm,
            ),
        ),
        bottom = (
            co2 = SoilCO2FluxBC((p, t) -> 0.0),
            o2 = SoilO2FluxBC((p, t) -> 0.0),
        ),
    )
    args = (parameters, domain, boundary_conditions, sources, drivers)
    SoilCO2Model{FT, typeof.(args)...}(args...)
end

ClimaLand.name(model::SoilCO2Model) = :soilco2
ClimaLand.prognostic_vars(::SoilCO2Model) = (:CO2, :O2, :SOC)
ClimaLand.prognostic_types(::SoilCO2Model{FT}) where {FT} = (FT, FT, FT)
ClimaLand.prognostic_domain_names(::SoilCO2Model) =
    (:subsurface, :subsurface, :subsurface)

ClimaLand.auxiliary_vars(model::SoilCO2Model) = (
    :D,
    :D_o2,
    :Sm,
    :T,
    :θ_eff,        # Effective porosity for CO2 (air + dissolved in water)
    :θ_eff_o2,     # Effective porosity for O2 (air + dissolved in water)
    # CO2 boundary vars (top_bc, bottom_bc, top_bc_wvec, bottom_bc_wvec)
    ClimaLand.boundary_vars(
        model.boundary_conditions.top.co2,
        ClimaLand.TopBoundary(),
    )...,
    ClimaLand.boundary_vars(
        model.boundary_conditions.bottom.co2,
        ClimaLand.BottomBoundary(),
    )...,
    # O2 boundary vars (top_bc_o2, bottom_bc_o2, top_bc_o2_wvec, bottom_bc_o2_wvec)
    ClimaLand.boundary_vars(
        model.boundary_conditions.top.o2,
        ClimaLand.TopBoundary(),
    )...,
    ClimaLand.boundary_vars(
        model.boundary_conditions.bottom.o2,
        ClimaLand.BottomBoundary(),
    )...,
    :bidiag_matrix_scratch,
    :full_bidiag_matrix_scratch,
    :topBC_scratch, # useful for most cases, so we always include it
)


ClimaLand.auxiliary_types(model::SoilCO2Model{FT}) where {FT} = (
    FT,
    FT,
    FT,
    FT,
    FT,  # θ_eff
    FT,  # θ_eff_o2
    # CO2 boundary var types
    ClimaLand.boundary_var_types(
        model,
        model.boundary_conditions.top.co2,
        ClimaLand.TopBoundary(),
    )...,
    ClimaLand.boundary_var_types(
        model,
        model.boundary_conditions.bottom.co2,
        ClimaLand.BottomBoundary(),
    )...,
    ClimaLand.boundary_var_types(
        model,
        model.boundary_conditions.top.o2,
        ClimaLand.TopBoundary(),
    )...,
    ClimaLand.boundary_var_types(
        model,
        model.boundary_conditions.bottom.o2,
        ClimaLand.BottomBoundary(),
    )...,
    MatrixFields.BidiagonalMatrixRow{Geometry.Covariant3Vector{FT}},
    MatrixFields.BidiagonalMatrixRow{Geometry.Covariant3Vector{FT}},
    ClimaCore.Geometry.Covariant3Vector{FT},
)
ClimaLand.auxiliary_domain_names(model::SoilCO2Model) = (
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,  # θ_eff
    :subsurface,  # θ_eff_o2
    # CO2 boundary var domain names
    ClimaLand.boundary_var_domain_names(
        model.boundary_conditions.top.co2,
        ClimaLand.TopBoundary(),
    )...,
    ClimaLand.boundary_var_domain_names(
        model.boundary_conditions.bottom.co2,
        ClimaLand.BottomBoundary(),
    )...,
    ClimaLand.boundary_var_domain_names(
        model.boundary_conditions.top.o2,
        ClimaLand.TopBoundary(),
    )...,
    ClimaLand.boundary_var_domain_names(
        model.boundary_conditions.bottom.o2,
        ClimaLand.BottomBoundary(),
    )...,
    :subsurface_face,
    :subsurface_face,
    :subsurface_face,
)

function ClimaLand.make_update_implicit_boundary_fluxes(model::SoilCO2Model)
    function update_boundary_fluxes!(p, Y, t)
        Δz_top = model.domain.fields.Δz_top
        # Update CO2 boundary fluxes
        boundary_flux!(
            p.soilco2.top_bc,
            model.boundary_conditions.top.co2,
            TopBoundary(),
            Δz_top,
            Y,
            p,
            t,
        )
        # Update O2 boundary fluxes
        boundary_flux!(
            p.soilco2.top_bc_o2,
            model.boundary_conditions.top.o2,
            TopBoundary(),
            Δz_top,
            Y,
            p,
            t,
        )
    end
    return update_boundary_fluxes!
end

"""
    make_compute_exp_tendency(model::SoilCO2Model)

An extension of the function `make_compute_exp_tendency`, for the soilco2 equation.
This function creates and returns a function which computes the entire
right hand side of the PDE for `CO2`, `O2`, and `SOC`, and updates `dY.soilco2.CO2`,
`dY.soilco2.O2`, and `dY.soilco2.SOC` in place with those values.
These quantities will be stepped explicitly.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_compute_exp_tendency(model::SoilCO2Model{FT}) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        dY.soilco2.CO2 .= FT(0)
        dY.soilco2.O2 .= FT(0)
        dY.soilco2.SOC .= FT(0)

        for src in model.sources
            source!(dY, src, Y, p, model.parameters) # explicit sources
        end
    end
    return compute_exp_tendency!
end

function ClimaLand.make_compute_imp_tendency(model::SoilCO2Model{FT}) where {FT}
    function compute_imp_tendency!(dY, Y, p, t)
        top_flux_bc = p.soilco2.top_bc
        bottom_flux_bc = p.soilco2.bottom_bc
        @. p.soilco2.top_bc_wvec = Geometry.WVector(top_flux_bc)
        @. p.soilco2.bottom_bc_wvec = Geometry.WVector(bottom_flux_bc)

        top_flux_bc_o2 = p.soilco2.top_bc_o2
        bottom_flux_bc_o2 = p.soilco2.bottom_bc_o2
        @. p.soilco2.top_bc_o2_wvec = Geometry.WVector(top_flux_bc_o2)
        @. p.soilco2.bottom_bc_o2_wvec = Geometry.WVector(bottom_flux_bc_o2)

        interpc2f = Operators.InterpolateC2F()
        gradc2f_C = Operators.GradientC2F()
        gradc2f_O2 = Operators.GradientC2F()
        divf2c_C = Operators.DivergenceF2C(
            top = Operators.SetValue(p.soilco2.top_bc_wvec),
            bottom = Operators.SetValue(p.soilco2.bottom_bc_wvec),
        )
        divf2c_O2 = Operators.DivergenceF2C(
            top = Operators.SetValue(p.soilco2.top_bc_o2_wvec),
            bottom = Operators.SetValue(p.soilco2.bottom_bc_o2_wvec),
        )

        # CO₂ diffusion using air-equivalent concentration for stable saturated soil behavior
        # D from Ryan et al. already accounts for air-filled porosity via tortuosity
        @. dY.soilco2.CO2 =
            -divf2c_C(
                -interpc2f(p.soilco2.D) *
                gradc2f_C(max(Y.soilco2.CO2, 0) / p.soilco2.θ_eff),
            )

        # O₂ diffusion: compute ∇·[D_O2 * ∇ρ_O2] on mass concentration,
        # then convert to O2 tendency using θ_eff_o2 for stability in saturated soils
        @. dY.soilco2.O2 =
            -divf2c_O2(
                -interpc2f(p.soilco2.D_o2) *
                gradc2f_O2(max(Y.soilco2.O2, 0) / p.soilco2.θ_eff_o2),
            )

        # SOC has no implicit piece
        @. dY.soilco2.SOC = 0.0
    end
    return compute_imp_tendency!
end

function ClimaLand.make_update_implicit_aux(model::SoilCO2Model)
    update_imp_aux!(p, Y, t) = nothing
    return update_imp_aux!
end

function ClimaLand.make_update_boundary_fluxes(model::SoilCO2Model)
    function update_bf!(p, Y, t)
        Δz_top = model.domain.fields.Δz_top
        Δz_bottom = model.domain.fields.Δz_bottom

        # Update CO2 boundary fluxes
        boundary_flux!(
            p.soilco2.top_bc,
            model.boundary_conditions.top.co2,
            TopBoundary(),
            Δz_top,
            Y,
            p,
            t,
        )
        boundary_flux!(
            p.soilco2.bottom_bc,
            model.boundary_conditions.bottom.co2,
            BottomBoundary(),
            Δz_bottom,
            Y,
            p,
            t,
        )

        # Update O2 boundary fluxes
        boundary_flux!(
            p.soilco2.top_bc_o2,
            model.boundary_conditions.top.o2,
            TopBoundary(),
            Δz_top,
            Y,
            p,
            t,
        )
        boundary_flux!(
            p.soilco2.bottom_bc_o2,
            model.boundary_conditions.bottom.o2,
            BottomBoundary(),
            Δz_bottom,
            Y,
            p,
            t,
        )
    end
    return update_bf!
end


"""
    AbstractCarbonSource{FT} <: ClimaLand.AbstractSource{FT}

An abstract type for soil CO2 sources. The model is currently
heterotrophic-only, with a single concrete source `MicrobeProduction`
(microbial decomposition of soil organic matter).
"""
abstract type AbstractCarbonSource{FT} <: ClimaLand.AbstractSource{FT} end

"""
    MicrobeProduction{FT} <: AbstractCarbonSource{FT}

Struct for the microbe production of CO2, appearing as a source
term in the differential equation.
"""
struct MicrobeProduction{FT} <: AbstractCarbonSource{FT} end

"""
    ClimaLand.source!(dY::ClimaCore.Fields.FieldVector,
                          src::MicrobeProduction,
                          Y::ClimaCore.Fields.FieldVector,
                          p::NamedTuple,
                          params)

A method which extends the ClimaLand source! function for the
case of microbe production of CO2 in soil and consumption of O2
SOC is held constant (initialized from data, no tendency).

Physics:
- CO2 production from microbial respiration (kg C m⁻³ s⁻¹)
- O2 consumption with correct stoichiometry: C + O₂ → CO₂
  For every 12 kg C respired, 32 kg O₂ is consumed (ratio = 32/12 = 8/3)
"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::MicrobeProduction,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    params,
)
    dY.soilco2.CO2 .+= p.soilco2.Sm

    FT = eltype(p.soilco2.Sm)
    # Stoichiometry of aerobic respiration C + O₂ → CO₂: 1 mol O₂ consumed per
    # mol C respired, i.e. (M_O2 / M_C) kg O₂ per kg C respired. Sm is in
    # kg C m⁻³ s⁻¹ and dY.soilco2.O2 is in kg O₂ m⁻³ s⁻¹, so this mass-ratio
    # factor (= 32/12 = 8/3) is required for the O₂ sink to balance the C source.
    stoich_O2_per_C = FT(params.M_O2 / params.M_C)
    T_soil = p.soilco2.T  # soil temperature (K)
    P_sfc = p.drivers.P   # atmospheric pressure (Pa)

    # Extra Michaelis-Menten attenuation drives O₂ consumption to ~O2_f² near
    # zero (Sm already carries an MM_o2 ~O2_f factor), keeping O₂ from being
    # pulled below zero by the explicit source step.
    O2_f_lim = FT(1e-4)
    @. dY.soilco2.O2 -=
        stoich_O2_per_C *
        p.soilco2.Sm *
        (
            o2_fraction_from_concentration(
                max(Y.soilco2.O2, 0) / p.soilco2.θ_eff_o2,
                T_soil,
                P_sfc,
                params,
            ) / (
                o2_fraction_from_concentration(
                    max(Y.soilco2.O2, 0) / p.soilco2.θ_eff_o2,
                    T_soil,
                    P_sfc,
                    params,
                ) + O2_f_lim
            )
        )
end

"""
    AbstractSoilDriver

An abstract type for drivers of soil CO2 production and diffusion.
These are soil temperature, soil moisture,
root carbon, soil organic matter and microbe carbon, and atmospheric pressure.
Soil temperature and moisture, as well as soc, vary in space (horizontally and vertically) and time.
Atmospheric pressure vary in time (defined at the surface only, not with depth).
"""
abstract type AbstractSoilDriver end

"""
    SoilDrivers

A container which passes in the soil drivers to the biogeochemistry
model. These drivers are either of type Prescribed (for standalone mode)
or Prognostic (for running with a prognostic model for soil temp and moisture).

$(DocStringExtensions.FIELDS)
"""
struct SoilDrivers{
    FT,
    MET <: AbstractSoilDriver,
    ATM <: AbstractAtmosphericDrivers{FT},
}
    "Soil temperature and moisture drivers - Prescribed or Prognostic"
    met::MET
    "Prescribed or coupled atmospheric variables"
    atmos::ATM
end

"""
    PrescribedMet <: AbstractSoilDriver

A container which holds the prescribed functions for soil temperature
and moisture.

This is meant for use when running the biogeochemistry model in standalone mode,
without a prognostic soil model.

$(DocStringExtensions.FIELDS)
"""
struct PrescribedMet{
    FT,
    F1 <: Function,
    F2 <: Function,
    F <: Union{AbstractFloat, ClimaCore.Fields.Field},
} <: AbstractSoilDriver
    "The temperature of the soil, of the form f(z::FT,t) where FT <: AbstractFloat"
    temperature::F1
    "Soil moisture, of the form f(z::FT,t) FT <: AbstractFloat"
    volumetric_liquid_fraction::F2
    "Soil porosity (m³ m⁻³)"
    ν::F
    "Air-filled porosity at soil water potential of -100 cm H₂O (~ 10 Pa)"
    θ_a100::F
    "Absolute value of the slope of the line relating log(ψ) versus log(S) (unitless)"
    b::F
end

function PrescribedMet{FT}(
    temperature::Function,
    volumetric_liquid_fraction::Function,
    ν::F,
    θ_r::F,
    hcm::CF,
) where {FT <: AbstractFloat, F, CF}
    if F <: AbstractFloat
        θ_a100 =
            ClimaLand.Soil.inverse_matric_potential(hcm, -FT(1)) * (ν - θ_r) +
            θ_r
        b = ClimaLand.Soil.approximate_ψ_S_slope(hcm)
    else
        θ_a100 = @. ClimaLand.Soil.inverse_matric_potential(hcm, -FT(1)) *
           (ν - θ_r) + θ_r
        b = @. ClimaLand.Soil.approximate_ψ_S_slope(hcm)
    end

    return PrescribedMet{
        FT,
        typeof(temperature),
        typeof(volumetric_liquid_fraction),
        F,
    }(
        temperature,
        volumetric_liquid_fraction,
        ν,
        θ_a100,
        b,
    )
end

"""
    soil_temperature(driver::PrescribedMet, p, Y, t, z)

Returns the soil temperature at location (z) and time (t) for the prescribed
soil case.
"""
function soil_temperature(driver::PrescribedMet, p, Y, t, z)
    return driver.temperature.(z, t)
end

"""
    soil_moisture(driver::PrescribedMet, p, Y, t, z)

Returns the soil moisture at location (z) and time (t) for the prescribed
soil case.
"""
function soil_moisture(driver::PrescribedMet, p, Y, t, z)
    return driver.volumetric_liquid_fraction.(z, t)
end

"""
    soil_ice(driver::PrescribedMet, p, Y, t, z)

Returns zero ice content for prescribed soil case (standalone mode has no ice).
"""
function soil_ice(driver::PrescribedMet, p, Y, t, z)
    FT = eltype(z)
    return FT(0)
end

"""
    make_update_aux(model::SoilCO2Model)

An extension of the function `make_update_aux`, for the soilco2 equation.
This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.
This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_update_aux(model::SoilCO2Model)
    function update_aux!(p, Y, t)
        FT = eltype(Y.soilco2.CO2)

        params = model.parameters
        z = model.domain.fields.z
        T_soil = soil_temperature(model.drivers.met, p, Y, t, z)
        θ_l = soil_moisture(model.drivers.met, p, Y, t, z)
        θ_i = soil_ice(model.drivers.met, p, Y, t, z)
        Csom = Y.soilco2.SOC
        P_sfc = p.drivers.P
        ν = model.drivers.met.ν
        θ_a100 = model.drivers.met.θ_a100
        b = model.drivers.met.b

        @. p.soilco2.T = T_soil
        @. p.soilco2.D =
            co2_diffusivity(T_soil, θ_l + θ_i, P_sfc, θ_a100, b, ν, params)
        @. p.soilco2.D_o2 =
            o2_diffusivity(T_soil, θ_l + θ_i, P_sfc, θ_a100, b, ν, params)

        # Compute Henry's law factors (temperature-dependent)
        R = LP.gas_constant(params.earth_param_set)
        K_H_co2_298 = params.K_H_co2_298
        dln_K_H_co2_dT = params.dln_K_H_co2_dT
        K_H_o2_298 = params.K_H_o2_298
        dln_K_H_o2_dT = params.dln_K_H_o2_dT
        T_ref_henry = params.T_ref_henry

        # Compute effective porosities (θ_l only, not θ_i - gas dissolves in liquid water)
        @. p.soilco2.θ_eff = effective_porosity(
            volumetric_air_content(θ_l + θ_i, ν),
            θ_l,
            beta_gas(
                henry_constant(
                    K_H_co2_298,
                    dln_K_H_co2_dT,
                    T_soil,
                    T_ref_henry,
                ),
                R,
                T_soil,
            ),
        )

        @. p.soilco2.θ_eff_o2 = effective_porosity(
            volumetric_air_content(θ_l + θ_i, ν),
            θ_l,
            beta_gas(
                henry_constant(K_H_o2_298, dln_K_H_o2_dT, T_soil, T_ref_henry),
                R,
                T_soil,
            ),
        )
        (; D_oa) = params
        @. p.soilco2.Sm = microbe_source(
            T_soil,
            θ_l,
            max(Csom, 0),
            o2_availability(
                o2_fraction_from_concentration(
                    max(Y.soilco2.O2, 0) / p.soilco2.θ_eff_o2,
                    T_soil,
                    P_sfc,
                    params,
                ),
                volumetric_air_content(θ_l + θ_i, ν),
                D_oa,
            ),
            params,
        ) # in case Csom < 0 due to numerical issues
    end
    return update_aux!
end

# Boundary condition types/methods for Soil CO2
"""
    SoilCO2FluxBC <: ClimaLand.AbstractBC

A container holding the CO2 flux boundary condition,
which is a function `f(p,t)`, where `p` is the auxiliary state
vector.
"""
struct SoilCO2FluxBC{F <: Function} <: ClimaLand.AbstractBC
    bc::F
end

"""
    ClimaLand.boundary_flux!(bc_field,
        bc::SoilCO2FluxBC,
        boundary::ClimaLand.AbstractBoundary,
        Δz::ClimaCore.Fields.Field,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    )

A method of ClimaLand.boundary_flux which updates the soilco2
flux (kg C m⁻² s⁻¹) in the case of a prescribed flux BC at either the top
or bottom of the domain.
"""
function ClimaLand.boundary_flux!(
    bc_field,
    bc::SoilCO2FluxBC,
    boundary::ClimaLand.AbstractBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    FT = eltype(Δz)
    bc_field .= FT.(bc.bc(p, t))
end

"""
    SoilCO2StateBC <: ClimaLand.AbstractBC

A container holding the CO2 state boundary condition (kg C m⁻³ air-equivalent),
which is a function `f(p,t)`, where `p` is the auxiliary state
vector.
"""
struct SoilCO2StateBC{F <: Function} <: ClimaLand.AbstractBC
    bc::F
end

"""
    ClimaLand.boundary_flux!(bc_field,
    bc::SoilCO2StateBC,
    boundary::ClimaLand.TopBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
    )

A method of ClimaLand.boundary_flux which returns the soilco2
flux in the case of a prescribed state BC at
top of the domain.

The prescribed state should be the air-equivalent CO2 concentration
(kg C/m³ air), consistent with the diffused variable.
"""
function ClimaLand.boundary_flux!(
    bc_field,
    bc::SoilCO2StateBC,
    boundary::ClimaLand.TopBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    FT = eltype(Δz)

    # We need to project center values onto the face space
    D_c = ClimaLand.Domains.top_center_to_surface(p.soilco2.D)
    C_c = ClimaLand.Domains.top_center_to_surface(Y.soilco2.CO2)
    θ_sfc = ClimaLand.Domains.top_center_to_surface(p.soilco2.θ_eff)
    C_bc = FT.(bc.bc(p, t))
    @. bc_field = ClimaLand.diffusive_flux(D_c, C_bc, C_c / θ_sfc, Δz)
    levels =
        ClimaCore.Spaces.nlevels(ClimaCore.Spaces.face_space(axes(p.soilco2.D)))
    local_geometry_faceN = ClimaCore.Fields.level(
        Fields.local_geometry_field(
            ClimaCore.Spaces.face_space(axes(p.soilco2.D)),
        ),
        levels - ClimaCore.Utilities.half,
    )
    @. p.soilco2.dfluxBCdY =
        ClimaLand.Soil.covariant3_unit_vector(local_geometry_faceN) *
        (D_c / θ_sfc / Δz)
end

"""
    ClimaLand.boundary_flux!(bc_field,
        bc::SoilCO2StateBC,
        boundary::ClimaLand.BottomBoundary,
        Δz::ClimaCore.Fields.Field,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    )

A method of ClimaLand.boundary_flux which returns the soilco2
flux in the case of a prescribed state BC at
bottom of the domain.

The prescribed state should be the air-equivalent CO2 concentration
(kg C/m³ air), consistent with the diffused variable.
"""
function ClimaLand.boundary_flux!(
    bc_field,
    bc::SoilCO2StateBC,
    ::ClimaLand.BottomBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    FT = eltype(Δz)
    D_c = ClimaLand.Domains.bottom_center_to_surface(p.soilco2.D)
    C_c = ClimaLand.Domains.bottom_center_to_surface(Y.soilco2.CO2)
    θ_sfc = ClimaLand.Domains.bottom_center_to_surface(p.soilco2.θ_eff)
    C_bc = FT.(bc.bc(p, t))
    @. bc_field = ClimaLand.diffusive_flux(D_c, C_c / θ_sfc, C_bc, Δz)
end

"""
    AtmosCO2StateBC <: ClimaLand.AbstractBC

Set the CO2 concentration to the atmospheric one.
Stores physical constants needed for the boundary flux calculation.

$(DocStringExtensions.FIELDS)
"""
struct AtmosCO2StateBC{FT <: AbstractFloat} <: ClimaLand.AbstractBC
    "Universal gas constant (J/(mol·K))"
    R::FT
    "Molar mass of carbon (kg/mol)"
    M_C::FT
end

"""
    AtmosCO2StateBC(earth_param_set, M_C::FT) where {FT}

Constructor for AtmosCO2StateBC that gets parameters from LandParameters and model parameters.
"""
function AtmosCO2StateBC(earth_param_set, M_C::FT) where {FT}
    R = LP.gas_constant(earth_param_set)
    return AtmosCO2StateBC{FT}(R, M_C)
end

"""
    ClimaLand.boundary_flux!(bc_field,
    bc::AtmosCO2StateBC,
    boundary::ClimaLand.TopBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
    )

A method of ClimaLand.boundary_flux which returns the soilco2 flux in the case when the
atmospheric CO2 is used at top of the domain. Uses air-equivalent CO2 concentration
for stable behavior in saturated soils.
"""
function ClimaLand.boundary_flux!(
    bc_field,
    bc::AtmosCO2StateBC,
    boundary::ClimaLand.TopBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    D_c = ClimaLand.Domains.top_center_to_surface(p.soilco2.D)
    C_c = ClimaLand.Domains.top_center_to_surface(Y.soilco2.CO2)
    θ_sfc = ClimaLand.Domains.top_center_to_surface(p.soilco2.θ_eff)
    T_sfc = p.drivers.T
    P_sfc = p.drivers.P
    R = bc.R
    M_C = bc.M_C
    FT = eltype(T_sfc)
    @. bc_field = ClimaLand.diffusive_flux(
        D_c,
        p.drivers.c_co2 * P_sfc * M_C / (R * T_sfc),
        max(C_c / θ_sfc, 0),
        Δz,
    )
    levels =
        ClimaCore.Spaces.nlevels(ClimaCore.Spaces.face_space(axes(p.soilco2.D)))
    local_geometry_faceN = ClimaCore.Fields.level(
        Fields.local_geometry_field(
            ClimaCore.Spaces.face_space(axes(p.soilco2.D)),
        ),
        levels - ClimaCore.Utilities.half,
    )
    @. p.soilco2.dfluxBCdY =
        ClimaLand.Soil.covariant3_unit_vector(local_geometry_faceN) *
        (D_c / θ_sfc / Δz)
end

"""
    SoilO2FluxBC <: ClimaLand.AbstractBC

A container holding the O2 flux boundary condition,
which is a function `f(p,t)`, where `p` is the auxiliary state
vector.
"""
struct SoilO2FluxBC{F <: Function} <: ClimaLand.AbstractBC
    bc::F
end

"""
    ClimaLand.boundary_flux!(bc_field,
        bc::SoilO2FluxBC,
        boundary::ClimaLand.AbstractBoundary,
        Δz::ClimaCore.Fields.Field,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    )

A method of ClimaLand.boundary_flux which updates the soil o2
flux (kg C m⁻² s⁻¹) in the case of a prescribed flux BC at either the top
or bottom of the domain.
"""
function ClimaLand.boundary_flux!(
    bc_field,
    bc::SoilO2FluxBC,
    boundary::ClimaLand.AbstractBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    FT = eltype(Δz)
    bc_field .= FT.(bc.bc(p, t))
end

"""
    AtmosO2StateBC{FT} <: ClimaLand.AbstractBC

Set the O2 mass concentration to the atmospheric one.
Stores physical constants needed for the boundary flux calculation.

$(DocStringExtensions.FIELDS)
"""
struct AtmosO2StateBC{FT <: AbstractFloat} <: ClimaLand.AbstractBC
    "Universal gas constant (J/(mol·K))"
    R::FT
    "Molar mass of oxygen (kg/mol)"
    M_O2::FT
    "Atmospheric O2 volumetric fraction (dimensionless)"
    O2_f_atm::FT
end

"""
    AtmosO2StateBC(earth_param_set, M_O2::FT, O2_f_atm::FT) where {FT}

Constructor for AtmosO2StateBC that gets parameters from LandParameters and model parameters.
"""
function AtmosO2StateBC(earth_param_set, M_O2::FT, O2_f_atm::FT) where {FT}
    R = LP.gas_constant(earth_param_set)
    return AtmosO2StateBC{FT}(R, M_O2, O2_f_atm)
end

# CO2 BC names; the default applies for SoilCO2FluxBC
boundary_vars(
    ::Union{AtmosCO2StateBC, SoilCO2StateBC},
    ::ClimaLand.TopBoundary,
) = (:top_bc, :top_bc_wvec, :dfluxBCdY)
boundary_var_domain_names(
    ::Union{AtmosCO2StateBC, SoilCO2StateBC},
    ::ClimaLand.TopBoundary,
) = (:surface, :surface, :surface)
function boundary_var_types(
    model::SoilCO2Model{FT},
    ::Union{AtmosCO2StateBC, SoilCO2StateBC},
    ::ClimaLand.TopBoundary,
) where {FT}
    (
        FT,
        ClimaCore.Geometry.WVector{FT},
        ClimaCore.Geometry.Covariant3Vector{FT},
    )
end
# O2 BC names. The default will *not* work because we would have duplicate "top_bc" etc names. so O2 BC types always need new methods
boundary_vars(::AtmosO2StateBC, ::ClimaLand.TopBoundary) =
    (:top_bc_o2, :top_bc_o2_wvec, :dfluxBCdY_o2)
boundary_var_domain_names(::AtmosO2StateBC, ::ClimaLand.TopBoundary) =
    (:surface, :surface, :surface)
function boundary_var_types(
    model::SoilCO2Model{FT},
    ::AtmosO2StateBC,
    ::ClimaLand.TopBoundary,
) where {FT}
    (
        FT,
        ClimaCore.Geometry.WVector{FT},
        ClimaCore.Geometry.Covariant3Vector{FT},
    )
end
boundary_vars(
    ::Union{SoilO2FluxBC, AtmosO2StateBC},
    ::ClimaLand.BottomBoundary,
) = (:bottom_bc_o2, :bottom_bc_o2_wvec)
boundary_vars(::Union{SoilO2FluxBC, AtmosO2StateBC}, ::ClimaLand.TopBoundary) =
    (:top_bc_o2, :top_bc_o2_wvec)

"""
    ClimaLand.boundary_flux!(bc_field,
    bc::AtmosO2StateBC,
    boundary::ClimaLand.TopBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
    )

A method of ClimaLand.boundary_flux which returns the O2 diffusive flux in the case when the
atmospheric O2 (0.21 volumetric fraction) is used at top of the domain.
"""
function ClimaLand.boundary_flux!(
    bc_field,
    bc::AtmosO2StateBC,
    boundary::ClimaLand.TopBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
)
    FT = eltype(Δz)
    D_o2 = ClimaLand.Domains.top_center_to_surface(p.soilco2.D_o2)
    O2_c = ClimaLand.Domains.top_center_to_surface(Y.soilco2.O2)
    θ_sfc = ClimaLand.Domains.top_center_to_surface(p.soilco2.θ_eff_o2)
    T_sfc = p.drivers.T
    P_sfc = p.drivers.P
    O2_f_atm = bc.O2_f_atm
    @. bc_field = ClimaLand.diffusive_flux(
        D_o2,
        O2_f_atm * P_sfc * bc.M_O2 / (bc.R * T_sfc),
        max(O2_c / θ_sfc, 0),
        Δz,
    )
    levels =
        ClimaCore.Spaces.nlevels(ClimaCore.Spaces.face_space(axes(p.soilco2.D)))
    local_geometry_faceN = ClimaCore.Fields.level(
        Fields.local_geometry_field(
            ClimaCore.Spaces.face_space(axes(p.soilco2.D)),
        ),
        levels - ClimaCore.Utilities.half,
    )
    @. p.soilco2.dfluxBCdY_o2 =
        ClimaLand.Soil.covariant3_unit_vector(local_geometry_faceN) *
        (D_o2 / θ_sfc / Δz)
end

function ClimaLand.get_drivers(model::SoilCO2Model)
    return (model.drivers.atmos,)
end

Base.broadcastable(ps::SoilCO2ModelParameters) = tuple(ps)

function ClimaLand.make_compute_jacobian(model::SoilCO2Model{FT}) where {FT}
    function compute_jacobian!(
        jacobian::MatrixFields.FieldMatrixWithSolver,
        Y,
        p,
        dtγ,
        t,
    )
        (; matrix) = jacobian

        # Create divergence operator
        divf2c_op = Operators.DivergenceF2C()
        divf2c_matrix = MatrixFields.operator_matrix(divf2c_op)
        # Create gradient operator, and set gradient at boundaries to 0
        gradc2f_op = Operators.GradientC2F(
            top = Operators.SetGradient(Geometry.WVector(FT(0))),
            bottom = Operators.SetGradient(Geometry.WVector(FT(0))),
        )
        gradc2f_matrix = MatrixFields.operator_matrix(gradc2f_op)
        # Create interpolation operator
        interpc2f_op = Operators.InterpolateC2F(
            bottom = Operators.Extrapolate(),
            top = Operators.Extrapolate(),
        )

        # The derivative of the residual with respect to the prognostic variable
        ∂CO2res∂CO2 = matrix[@name(soilco2.CO2), @name(soilco2.CO2)]
        @. p.soilco2.bidiag_matrix_scratch =
            gradc2f_matrix() ⋅
            MatrixFields.DiagonalMatrixRow(1 / p.soilco2.θ_eff)
        @. p.soilco2.full_bidiag_matrix_scratch =
            MatrixFields.DiagonalMatrixRow(interpc2f_op(p.soilco2.D)) ⋅
            p.soilco2.bidiag_matrix_scratch
        if haskey(p.soilco2, :dfluxBCdY)
            topBC_op = Operators.SetBoundaryOperator(
                top = Operators.SetValue(p.soilco2.dfluxBCdY),
                bottom = Operators.SetValue(
                    Geometry.Covariant3Vector(zero(FT)),
                ),
            )
            @. p.soilco2.topBC_scratch = topBC_op(
                Geometry.Covariant3Vector(zero(interpc2f_op(p.soilco2.D))),
            )
            @. p.soilco2.full_bidiag_matrix_scratch +=
                -MatrixFields.LowerDiagonalMatrixRow(p.soilco2.topBC_scratch)
        end
        @. ∂CO2res∂CO2 =
            float(dtγ) *
            (divf2c_matrix() ⋅ p.soilco2.full_bidiag_matrix_scratch) - (I,)

        ∂O2res∂O2 = matrix[@name(soilco2.O2), @name(soilco2.O2)]
        @. p.soilco2.bidiag_matrix_scratch =
            gradc2f_matrix() ⋅
            MatrixFields.DiagonalMatrixRow(1 / p.soilco2.θ_eff_o2)
        @. p.soilco2.full_bidiag_matrix_scratch =
            MatrixFields.DiagonalMatrixRow(interpc2f_op(p.soilco2.D_o2)) ⋅
            p.soilco2.bidiag_matrix_scratch
        if haskey(p.soilco2, :dfluxBCdY_o2)
            topBC_op_o2 = Operators.SetBoundaryOperator(
                top = Operators.SetValue(p.soilco2.dfluxBCdY_o2),
                bottom = Operators.SetValue(
                    Geometry.Covariant3Vector(zero(FT)),
                ),
            )
            @. p.soilco2.topBC_scratch = topBC_op_o2(
                Geometry.Covariant3Vector(zero(interpc2f_op(p.soilco2.D_o2))),
            )
            @. p.soilco2.full_bidiag_matrix_scratch +=
                -MatrixFields.LowerDiagonalMatrixRow(p.soilco2.topBC_scratch)
        end
        @. ∂O2res∂O2 =
            float(dtγ) *
            (divf2c_matrix() ⋅ p.soilco2.full_bidiag_matrix_scratch) - (I,)
    end
    return compute_jacobian!
end

include("./co2_parameterizations.jl")

end # module
