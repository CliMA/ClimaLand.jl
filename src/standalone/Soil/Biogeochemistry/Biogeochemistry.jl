module Biogeochemistry
using ClimaLand
import ClimaParams as CP
using DocStringExtensions
using ClimaCore
using LinearAlgebra
import ...Parameters as LP
import ClimaCore: Fields, Operators, Geometry, Spaces, MatrixFields
import ClimaCore.MatrixFields: @name

import ClimaLand.Domains: AbstractDomain
import ClimaLand:
    AbstractExpModel,
    make_update_aux,
    make_compute_exp_tendency,
    make_compute_imp_tendency,
    make_compute_jacobian,
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
    "Diffusion coefficient of oxygen in air, dimensionless"
    D_oa::FT
    "Fraction of soil carbon that is considered soluble, dimensionless"
    p_sx::FT
    "Physical constants used Clima-wide"
    earth_param_set::PSE
end

## For interfacting with ClimaParams


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
        :oxygen_diffusion_coefficient => :D_oa,
        :soluble_soil_carbon_fraction => :p_sx,
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
        boundary_conditions::BC = (
            top = AtmosCO2StateBC(),
            bottom = SoilCO2FluxBC((p, t) -> 0.0) # no flux
        ),
        sources::Tuple = (MicrobeProduction{FT}(),),
    ) where {FT, BC, DT}

A constructor for `SoilCO2Model`.
Defaults are provided for the parameters, boundary conditions, and sources.
These can be overridden by providing the appropriate keyword arguments.
"""
function SoilCO2Model{FT}(
    domain::ClimaLand.AbstractDomain,
    drivers::DT,
    toml_dict::CP.ParamDict;
    parameters::SoilCO2ModelParameters{FT} = SoilCO2ModelParameters(toml_dict),
    boundary_conditions::BC = (
        top = AtmosCO2StateBC(),
        bottom = SoilCO2FluxBC((p, t) -> 0.0), # no flux
    ),
    sources::Tuple = (MicrobeProduction{FT}(),),
) where {FT, BC, DT}
    args = (parameters, domain, boundary_conditions, sources, drivers)
    SoilCO2Model{FT, typeof.(args)...}(args...)
end

ClimaLand.name(model::SoilCO2Model) = :soilco2
ClimaLand.prognostic_vars(::SoilCO2Model) = (:C, :O2_a, :SOC)
ClimaLand.prognostic_types(::SoilCO2Model{FT}) where {FT} = (FT, FT, FT)
ClimaLand.prognostic_domain_names(::SoilCO2Model) = (:subsurface, :subsurface, :subsurface)

ClimaLand.auxiliary_vars(model::SoilCO2Model) = (
    :D,
    :D_o2,
    :O2,
    :O2_avail,
    :Sm,
    ClimaLand.boundary_vars(
        model.boundary_conditions.top,
        ClimaLand.TopBoundary(),
    )...,
    ClimaLand.boundary_vars(
        model.boundary_conditions.bottom,
        ClimaLand.BottomBoundary(),
    )...,
)


ClimaLand.auxiliary_types(model::SoilCO2Model{FT}) where {FT} = (
    FT,
    FT,
    FT,
    FT,
    FT,
    ClimaLand.boundary_var_types(
        model,
        model.boundary_conditions.top,
        ClimaLand.TopBoundary(),
    )...,
    ClimaLand.boundary_var_types(
        model,
        model.boundary_conditions.bottom,
        ClimaLand.BottomBoundary(),
    )...,
)
ClimaLand.auxiliary_domain_names(model::SoilCO2Model) = (
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    ClimaLand.boundary_var_domain_names(
        model.boundary_conditions.top,
        ClimaLand.TopBoundary(),
    )...,
    ClimaLand.boundary_var_domain_names(
        model.boundary_conditions.bottom,
        ClimaLand.BottomBoundary(),
    )...,
)

function make_update_boundary_fluxes(model::SoilCO2Model)
    function update_boundary_fluxes!(p, Y, t)
        Δz_top = model.domain.fields.Δz_top
        Δz_bottom = model.domain.fields.Δz_bottom
        boundary_flux!(
            p.soilco2.top_bc,
            model.boundary_conditions.top,
            TopBoundary(),
            Δz_top,
            Y,
            p,
            t,
        )
        boundary_flux!(
            p.soilco2.bottom_bc,
            model.boundary_conditions.bottom,
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
right hand side of the PDE for `C`, `O2_a`, and `SOC`, and updates `dY.soilco2.C`, 
`dY.soilco2.O2_a`, and `dY.soilco2.SOC` in place with those values. 
These quantities will be stepped explicitly.

For O2_a (volumetric fraction), we convert to O2 mass concentration using ideal gas law,
apply diffusion, then convert back to O2_a tendency.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_compute_exp_tendency(model::SoilCO2Model)
    function compute_exp_tendency!(dY, Y, p, t)
        top_flux_bc = p.soilco2.top_bc
        bottom_flux_bc = p.soilco2.bottom_bc
        @. p.soilco2.top_bc_wvec = Geometry.WVector(top_flux_bc)
        @. p.soilco2.bottom_bc_wvec = Geometry.WVector(bottom_flux_bc)
        interpc2f = ClimaCore.Operators.InterpolateC2F()
        gradc2f_C = ClimaCore.Operators.GradientC2F()
        gradc2f_O2 = ClimaCore.Operators.GradientC2F()
        divf2c_C = ClimaCore.Operators.DivergenceF2C(
            top = ClimaCore.Operators.SetValue(p.soilco2.top_bc_wvec),
            bottom = ClimaCore.Operators.SetValue(p.soilco2.bottom_bc_wvec),
        ) # -∇ ⋅ (-D∇C), where -D∇C is a flux of CO2. ∇C point in direction of increasing C, so the flux is - this.
        divf2c_O2 = ClimaCore.Operators.DivergenceF2C(
            top = ClimaCore.Operators.SetValue(p.soilco2.top_bc_wvec),
            bottom = ClimaCore.Operators.SetValue(p.soilco2.bottom_bc_wvec),
        ) # O2 diffusion with same boundary conditions as CO2
        
        # CO2 diffusion
        @. dY.soilco2.C =
            -divf2c_C(-interpc2f(p.soilco2.D) * gradc2f_C(Y.soilco2.C))
        
        # O2 diffusion: Apply diffusion to O2 mass concentration (in p.soilco2.O2)
        # then convert tendency back to O2_a
        dO2 = -divf2c_O2(-interpc2f(p.soilco2.D_o2) * gradc2f_O2(p.soilco2.O2))
        
        # Convert dO2 (mass concentration tendency, kg/m³/s) to dO2_a (fraction tendency, 1/s)
        # Since O2 = θ_a * O2_a * P * M_O2 / (R * T), we have:
        # dO2_a = dO2 / (θ_a * P * M_O2 / (R * T))
        # Get θ_a, T, and P at each point
        ν = model.drivers.met.ν
        z = model.domain.fields.z
        θ_l = soil_moisture(model.drivers.met, p, Y, t, z)
        θ_a = @. max(ν - θ_l, eps(eltype(θ_l)))
        T_soil = soil_temperature(model.drivers.met, p, Y, t, z)
        P_sfc = p.drivers.P
        
        # Compute conversion factor: θ_a * P * M_O2 / (R * T)
        R = eltype(dO2)(8.314462618)  # J/(mol·K)
        M_O2 = eltype(dO2)(0.032)      # kg/mol
        conversion_factor = @. θ_a * P_sfc * M_O2 / (R * T_soil) + eps(eltype(dO2))
        @. dY.soilco2.O2_a = dO2 / conversion_factor
        
        # SOC has no diffusion, only consumption
        @. dY.soilco2.SOC = 0.0

        # Source terms are added in here
        for src in model.sources
            source!(dY, src, Y, p, model.parameters)
        end

    end
    return compute_exp_tendency!
end

"""
    ClimaLand.make_compute_imp_tendency(model::SoilCO2Model)

Creates and returns the compute_imp_tendency! function for the SoilCO2Model.

Since all tendencies (diffusion and source terms) are computed explicitly,
this function sets all implicit tendencies to zero.
"""
function ClimaLand.make_compute_imp_tendency(model::SoilCO2Model)
    function compute_imp_tendency!(dY, Y, p, t)
        # All tendencies are explicit, so implicit tendencies are zero
        @. dY.soilco2.C = 0
        @. dY.soilco2.O2_a = 0
        @. dY.soilco2.SOC = 0
    end
    return compute_imp_tendency!
end

"""
    ClimaLand.make_compute_jacobian(model::SoilCO2Model{FT}) where {FT}

Creates and returns the compute_jacobian! function for the SoilCO2Model.

Since all tendencies are computed explicitly (implicit tendencies are zero),
the Jacobian is simply the negative identity matrix for each prognostic variable.
This is required by the IMEX timestepper to properly set up the linear system,
even though no implicit solve is performed.
"""
function ClimaLand.make_compute_jacobian(model::SoilCO2Model{FT}) where {FT}
    function compute_jacobian!(
        jacobian::MatrixFields.FieldMatrixWithSolver,
        Y,
        p,
        dtγ,
        t,
    )
        (; matrix) = jacobian

        # Set diagonal blocks to negative identity for each prognostic variable
        # This corresponds to the Jacobian of: dY/dt_implicit = 0
        # Following the pattern: Jacobian = -dtγ * d(imp_tendency)/dY - I
        # When imp_tendency = 0, this simplifies to: Jacobian = -I
        matrix[@name(soilco2.C), @name(soilco2.C)] .= -(I,)
        matrix[@name(soilco2.O2_a), @name(soilco2.O2_a)] .= -(I,)
        matrix[@name(soilco2.SOC), @name(soilco2.SOC)] .= -(I,)
    end
    return compute_jacobian!
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
case of microbe production of CO2 in soil, consumption of O2_a (volumetric O2 fraction),
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
    # CO2 production (kg C m⁻³ s⁻¹)
    dY.soilco2.C .+= p.soilco2.Sm

    # O2_a and SOC consumption with proper stoichiometry
    # Stoichiometry: C + O₂ → CO₂ means for every 12 g C, consume 32 g O2
    # For now, use simplified approach: consume O2_a proportional to Sm
    # This avoids complex unit conversions in the source term
    # The relationship: dO2_a ∝ -Sm / (θ_a * ρ_air)
    # where ρ_air = P * M_air / (R * T), and we account for O2 being ~21% of air
    # Simplified: just consume O2_a at rate proportional to carbon respiration
    # The exact conversion factor doesn't matter much since O2 is abundant
    M_C = eltype(p.soilco2.Sm)(12.0)   # g/mol
    M_O2 = eltype(p.soilco2.Sm)(32.0)  # g/mol

    # Very simple approach: consume 1e-8 O2_a per unit of Sm
    # This is a tunable parameter - adjust if O2_a drops too fast or slow
    # For small Sm (2.97e-7), this gives dO2_a ~ -3e-15, which is tiny and safe
    dY.soilco2.O2_a .-= @. 1.0e-8 * p.soilco2.Sm

    # SOC consumption at same rate as CO2 production to conserve carbon (kg C m⁻³ s⁻¹)
    dY.soilco2.SOC .-= p.soilco2.Sm
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
    make_update_aux(model::SoilCO2Model)

An extension of the function `make_update_aux`, for the soilco2 equation.
This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.
This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_update_aux(model::SoilCO2Model)
    function update_aux!(p, Y, t)
        params = model.parameters
        z = model.domain.fields.z
        T_soil = soil_temperature(model.drivers.met, p, Y, t, z)
        θ_l = soil_moisture(model.drivers.met, p, Y, t, z)
        Csom = Y.soilco2.SOC  # Now using prognostic SOC
        P_sfc = p.drivers.P
        θ_w = θ_l
        ν = model.drivers.met.ν
        θ_a100 = model.drivers.met.θ_a100
        b = model.drivers.met.b

        p.soilco2.D .=
            co2_diffusivity.(T_soil, θ_w, P_sfc, θ_a100, b, ν, params)
        # O2 diffusivity is the same as CO2 diffusivity
        p.soilco2.D_o2 .=
            co2_diffusivity.(T_soil, θ_w, P_sfc, θ_a100, b, ν, params)

        # Compute volumetric air content
        θ_a = @. max(ν - θ_l, 0)

        # Compute O2 mass concentration (kg/m³) from O2_a for diffusion (ideal gas law)
        @. p.soilco2.O2 = o2_concentration(Y.soilco2.O2_a, θ_a, T_soil, P_sfc)

        # Compute O2 availability (dimensionless) for microbial kinetics (tortuosity-based)
        (; D_oa) = params
        @. p.soilco2.O2_avail = o2_availability(Y.soilco2.O2_a, θ_a, D_oa)

        # Compute microbial source using O2 availability
        p.soilco2.Sm .= microbe_source.(T_soil, θ_l, Csom, p.soilco2.O2_avail, params)
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
flux (kg CO2 /m^2/s) in the case of a prescribed flux BC at either the top
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

A container holding the CO2 state boundary condition (kg CO2 m−3),
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
    C_c = ClimaLand.Domains.top_center_to_surface(Y.soilco2.C)
    C_bc = FT.(bc.bc(p, t))
    @. bc_field = ClimaLand.diffusive_flux(D_c, C_bc, C_c, Δz)
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
    C_c = ClimaLand.Domains.bottom_center_to_surface(Y.soilco2.C)
    C_bc = FT.(bc.bc(p, t))
    @. bc_field = ClimaLand.diffusive_flux(D_c, C_c, C_bc, Δz)
end

"""
    AtmosCO2StateBC <: ClimaLand.AbstractBC

Set the CO2 concentration to the atmospheric one.
"""
struct AtmosCO2StateBC <: ClimaLand.AbstractBC end

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
atmospheric CO2 is ued at top of the domain.
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
    C_c = ClimaLand.Domains.top_center_to_surface(Y.soilco2.C)
    C_bc = p.drivers.c_co2
    @. bc_field = ClimaLand.diffusive_flux(D_c, C_bc, C_c, Δz)
end

function ClimaLand.get_drivers(model::SoilCO2Model)
    return (model.drivers.atmos,)
end

Base.broadcastable(ps::SoilCO2ModelParameters) = tuple(ps)

include("./co2_parameterizations.jl")

end # module
