module Biogeochemistry
using ClimaLand
import ClimaParams as CP
using DocStringExtensions
using ClimaCore
import ...Parameters as LP
import ClimaCore: Fields, Operators, Geometry, Spaces

import ClimaLand.Domains: AbstractDomain
import ClimaLand:
    AbstractExpModel,
    make_update_aux,
    make_compute_exp_tendency,
    make_update_boundary_fluxes,
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
    source!
export SoilCO2ModelParameters,
    SoilCO2Model,
    PrescribedMet,
    MicrobeProduction,
    SoilCO2FluxBC,
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
    "Diffusivity of soil C substrate in liquid (unitless)"
    D_liq::FT
    # DAMM
    "Pre-exponential factor (kg C m-3 s-1)"
    α_sx::FT
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
)
    name_map = (;
        :soil_C_substrate_diffusivity => :D_liq,
        :soilCO2_pre_exponential_factor => :α_sx,
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
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    earth_param_set = LP.LandParameters(toml_dict)
    return SoilCO2ModelParameters{FT, typeof(earth_param_set)}(;
        earth_param_set,
        D_ref,
        parameters...,
    )
end

"""
    AbstractSoilBiogeochemistryModel{FT} <: ClimaLand.AbstractExpModel{FT}

An abstract model type for soil biogeochemistry models.
"""
abstract type AbstractSoilBiogeochemistryModel{FT} <:
              ClimaLand.AbstractExpModel{FT} end

"""
    SoilCO2Model

A model for simulating the production and transport of CO₂ in the soil with dynamic
source and diffusion terms.

ClimaLand v1: SoilCO2 is still under testing; in particular, in global runs,
an instability appears in some columns, and the prognostic equation does not
enforce the positivity of CO2.

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
            co2 = AtmosCO2StateBC(
                parameters.earth_param_set,
                parameters.M_C,
            ),
            o2 = AtmosO2StateBC(
                parameters.earth_param_set,
                parameters.M_O2,
                parameters.O2_f_atm,
            ),
        ),
        bottom = (
            co2 = SoilCO2FluxBC((p, t) -> 0.0),
            o2 = SoilCO2FluxBC((p, t) -> 0.0),
        ),
    )
    args = (parameters, domain, boundary_conditions, sources, drivers)
    SoilCO2Model{FT, typeof.(args)...}(args...)
end

ClimaLand.name(model::SoilCO2Model) = :soilco2
ClimaLand.prognostic_vars(::SoilCO2Model) = (:CO2, :O2_f, :SOC)
ClimaLand.prognostic_types(::SoilCO2Model{FT}) where {FT} = (FT, FT, FT)
ClimaLand.prognostic_domain_names(::SoilCO2Model) =
    (:subsurface, :subsurface, :subsurface)

ClimaLand.auxiliary_vars(model::SoilCO2Model) = (
    :D,
    :D_o2,
    :O2,
    :O2_avail,
    :Sm,
    :θ_a,
    :T,
    :θ_eff,        # Effective porosity for CO2 (air + dissolved in water)
    :θ_eff_o2,     # Effective porosity for O2 (air + dissolved in water)
    :CO2_air_eq,   # Air-equivalent CO2 concentration for diffusion
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
    :top_bc_o2,
    :bottom_bc_o2,
    :top_bc_o2_wvec,
    :bottom_bc_o2_wvec,
)


ClimaLand.auxiliary_types(model::SoilCO2Model{FT}) where {FT} = (
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,  # θ_eff
    FT,  # θ_eff_o2
    FT,  # CO2_air_eq
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
    # O2 boundary var types
    FT,
    FT,
    Geometry.WVector{FT},
    Geometry.WVector{FT},
)
ClimaLand.auxiliary_domain_names(model::SoilCO2Model) = (
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,  # θ_eff
    :subsurface,  # θ_eff_o2
    :subsurface,  # CO2_air_eq
    # CO2 boundary var domain names
    ClimaLand.boundary_var_domain_names(
        model.boundary_conditions.top.co2,
        ClimaLand.TopBoundary(),
    )...,
    ClimaLand.boundary_var_domain_names(
        model.boundary_conditions.bottom.co2,
        ClimaLand.BottomBoundary(),
    )...,
    # O2 boundary var domain names
    :surface,
    :surface,
    :surface,
    :surface,
)

function make_update_boundary_fluxes(model::SoilCO2Model)
    function update_boundary_fluxes!(p, Y, t)
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
    return update_boundary_fluxes!
end

"""
    make_compute_exp_tendency(model::SoilCO2Model)

An extension of the function `make_compute_exp_tendency`, for the soilco2 equation.
This function creates and returns a function which computes the entire
right hand side of the PDE for `CO2`, `O2_f`, and `SOC`, and updates `dY.soilco2.CO2`,
`dY.soilco2.O2_f`, and `dY.soilco2.SOC` in place with those values.
These quantities will be stepped explicitly.

For O2_f (volumetric fraction), we convert to O2 mass concentration using ideal gas law,
apply diffusion, then convert back to O2_f tendency.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_compute_exp_tendency(model::SoilCO2Model)
    function compute_exp_tendency!(dY, Y, p, t)
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
        # Multiply D by θ_a because diffusion only occurs through air-filled pores
        @. dY.soilco2.CO2 =
            -divf2c_C(-interpc2f(p.soilco2.D * p.soilco2.θ_a) * gradc2f_C(p.soilco2.CO2_air_eq))

        FT = eltype(Y.soilco2.CO2)
        T_soil = p.soilco2.T
        P_sfc = p.drivers.P
        R = FT(LP.gas_constant(model.parameters.earth_param_set))
        M_O2 = FT(model.parameters.M_O2)

        # O₂ diffusion: compute ∇·[D_O2 * θ_a * ∇ρ_O2] on mass concentration,
        # then convert to O2_f tendency using θ_eff_o2 for stability in saturated soils
        @. dY.soilco2.O2_f =
            -divf2c_O2(
                -interpc2f(p.soilco2.D_o2 * p.soilco2.θ_a) *
                gradc2f_O2(p.soilco2.O2),
            ) *
            R *
            T_soil / max(p.soilco2.θ_eff_o2 * P_sfc * M_O2, eps(FT))

        # SOC has no diffusion, only source/sink from microbial activity
        @. dY.soilco2.SOC = 0.0

        for src in model.sources
            source!(dY, src, Y, p, model.parameters)
        end

    end
    return compute_exp_tendency!
end

"""
    AbstractCarbonSource{FT} <: ClimaLand.AbstractSource{FT}

An abstract type for soil CO2 sources. There are two sources:
roots and microbes, in struct RootProduction and MicrobeProduction.
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
case of microbe production of CO2 in soil, consumption of O2_f (volumetric O2 fraction),
and consumption of SOC.

Physics:
- CO2 production from microbial respiration (kg C m⁻³ s⁻¹)
- O2 consumption with correct stoichiometry: C + O₂ → CO₂
  For every 12 kg C respired, 32 kg O₂ is consumed (ratio = 32/12 = 8/3)
- SOC consumption equals CO2 production to conserve carbon mass
"""
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::MicrobeProduction,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    params,
)
    dY.soilco2.CO2 .+= p.soilco2.Sm

    M_C = eltype(p.soilco2.Sm)(params.M_C)
    R = eltype(p.soilco2.Sm)(LP.gas_constant(params.earth_param_set))
    T_soil = p.soilco2.T  # soil temperature (K)
    P_sfc = p.drivers.P   # atmospheric pressure (Pa)

    # Use θ_eff_o2 for stability in saturated soils
    θ_eff_o2 = p.soilco2.θ_eff_o2

    @. dY.soilco2.O2_f -=
        (R * T_soil) / max(M_C * θ_eff_o2 * P_sfc, eps(eltype(p.soilco2.Sm))) *
        p.soilco2.Sm

    @. dY.soilco2.SOC -= p.soilco2.Sm
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
    return FT(0) .* z  # Return field of zeros matching z
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
        Csom = max.(Y.soilco2.SOC, FT(0))  # Clamp SOC to non-negative
        P_sfc = p.drivers.P
        ν = model.drivers.met.ν
        θ_a100 = model.drivers.met.θ_a100
        b = model.drivers.met.b

        @. p.soilco2.T = T_soil

        # Compute θ_a using total water (liquid + ice)
        θ_w = θ_l .+ θ_i
        @. p.soilco2.θ_a = volumetric_air_content(θ_w, ν)

        @. p.soilco2.D =
            co2_diffusivity.(T_soil, θ_w, P_sfc, θ_a100, b, ν, params)
        @. p.soilco2.D_o2 =
            co2_diffusivity.(T_soil, θ_w, P_sfc, θ_a100, b, ν, params)

        # Compute Henry's law factors (temperature-dependent)
        R = FT(LP.gas_constant(params.earth_param_set))
        K_H_co2_298 = params.K_H_co2_298
        dln_K_H_co2_dT = params.dln_K_H_co2_dT
        K_H_o2_298 = params.K_H_o2_298
        dln_K_H_o2_dT = params.dln_K_H_o2_dT

        # Compute effective porosities (θ_l only, not θ_i - gas dissolves in liquid water)
        @. p.soilco2.θ_eff = effective_porosity(
            p.soilco2.θ_a,
            θ_l,
            beta_gas(henry_constant(K_H_co2_298, dln_K_H_co2_dT, T_soil), R, T_soil),
        )
        @. p.soilco2.θ_eff_o2 = effective_porosity(
            p.soilco2.θ_a,
            θ_l,
            beta_gas(henry_constant(K_H_o2_298, dln_K_H_o2_dT, T_soil), R, T_soil),
        )

        # Compute air-equivalent CO2 concentration for diffusion
        @. p.soilco2.CO2_air_eq = Y.soilco2.CO2 / max(p.soilco2.θ_eff, eps(FT))

        @. p.soilco2.O2 =
            o2_concentration(Y.soilco2.O2_f, T_soil, P_sfc, params)

        (; D_oa) = params
        @. p.soilco2.O2_avail =
            o2_availability(Y.soilco2.O2_f, p.soilco2.θ_a, D_oa)

        @. p.soilco2.Sm =
            microbe_source.(T_soil, θ_l, Csom, p.soilco2.O2_avail, params)
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
flux (kg C m⁻² s⁻¹) in the case of a prescribed flux BC at either the top
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

A container holding the CO2 state boundary condition (kg C m⁻³ air-equivalent),
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
    θ_a_c = ClimaLand.Domains.top_center_to_surface(p.soilco2.θ_a)
    C_c = ClimaLand.Domains.top_center_to_surface(p.soilco2.CO2_air_eq)
    C_bc = FT.(bc.bc(p, t))
    # Multiply D by θ_a because diffusion only occurs through air-filled pores
    @. bc_field = ClimaLand.diffusive_flux(D_c * θ_a_c, C_bc, C_c, Δz)
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
    θ_a_c = ClimaLand.Domains.bottom_center_to_surface(p.soilco2.θ_a)
    C_c = ClimaLand.Domains.bottom_center_to_surface(p.soilco2.CO2_air_eq)
    C_bc = FT.(bc.bc(p, t))
    # Multiply D by θ_a because diffusion only occurs through air-filled pores
    @. bc_field = ClimaLand.diffusive_flux(D_c * θ_a_c, C_c, C_bc, Δz)
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
    R = FT(LP.gas_constant(earth_param_set))
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
    θ_a_c = ClimaLand.Domains.top_center_to_surface(p.soilco2.θ_a)
    D_c = ClimaLand.Domains.top_center_to_surface(p.soilco2.D)
    C_c = ClimaLand.Domains.top_center_to_surface(p.soilco2.CO2_air_eq)
    T_soil_top = ClimaLand.Domains.top_center_to_surface(p.soilco2.T)
    P_sfc = p.drivers.P
    R = bc.R
    M_C = bc.M_C
    C_bc = @. p.drivers.c_co2 * P_sfc * M_C / (R * T_soil_top)
    # Multiply D by θ_a because diffusion only occurs through air-filled pores
    @. bc_field = ClimaLand.diffusive_flux(D_c * θ_a_c, C_bc, C_c, Δz)
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
    R = FT(LP.gas_constant(earth_param_set))
    return AtmosO2StateBC{FT}(R, M_O2, O2_f_atm)
end

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
    θ_a_top = ClimaLand.Domains.top_center_to_surface(p.soilco2.θ_a)
    D_o2 = ClimaLand.Domains.top_center_to_surface(p.soilco2.D_o2)
    O2_c = ClimaLand.Domains.top_center_to_surface(p.soilco2.O2)  # Current O2 mass concentration in air at top (kg O2/m³ air)

    T_soil_top = ClimaLand.Domains.top_center_to_surface(p.soilco2.T)
    P_sfc = p.drivers.P

    R = bc.R
    M_O2 = bc.M_O2
    O2_f_atm = bc.O2_f_atm

    @. bc_field = ClimaLand.diffusive_flux(
        D_o2 * θ_a_top,
        O2_f_atm * P_sfc * M_O2 / (R * T_soil_top),
        O2_c,
        Δz,
    )
end

function ClimaLand.get_drivers(model::SoilCO2Model)
    return (model.drivers.atmos,)
end

Base.broadcastable(ps::SoilCO2ModelParameters) = tuple(ps)

include("./co2_parameterizations.jl")

end # module
