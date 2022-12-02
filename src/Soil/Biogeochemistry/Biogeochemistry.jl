module Biogeochemistry
using ClimaLSM
using UnPack
using DocStringExtensions
using ClimaCore
import ..Parameters as LSMP
import ClimaCore: Fields, Operators, Geometry, Spaces

import ClimaLSM.Domains: AbstractDomain
import ClimaLSM:
    AbstractModel,
    make_update_aux,
    make_rhs,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    auxiliary_types,
    TopBoundary,
    BottomBoundary,
    boundary_flux,
    AbstractBC,
    AbstractSource,
    source!
export soilco2Parameters, SoilCO2Model, PrescribedSoil, RootProduction, MicrobeProduction, SoilCO2FluxBC, SoilCO2StateBC

"""
    soilco2Parameters{FT <: AbstractFloat}

A struct for storing parameters of the `SoilCO2Model`.
$(DocStringExtensions.FIELDS)
"""
struct SoilCO2ModelParameters{FT <: AbstractFloat}
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
    SoilCO2ModelParameters(;
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

An outer constructor for creating the parameter struct of the `SoilCO2Model`,
    based on keyword arguments.
"""
function SoilCO2ModelParameters(;
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
    return SoilCO2ModelParameters{FT}(
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

"""
    AbstractSoilBiogeochemistryModel{FT} <: ClimaLSM.AbstractModel{FT}

An abstract model type for soil biogeochemistry models.
"""
abstract type AbstractSoilBiogeochemistryModel{FT} <: ClimaLSM.AbstractModel{FT} end

"""
    SoilCO2Model

A model for simulating the production and transport of CO₂ in the soil with dynamic
source and diffusion terms.
$(DocStringExtensions.FIELDS)
"""
struct SoilCO2Model{FT, PS, D, BC, S, DT} <: AbstractSoilBiogeochemistryModel{FT}
    "the parameter set"
    parameters::PS 
    "the soil domain, using ClimaCore.Domains"
    domain::D
    "the boundary conditions, of type AbstractBoundaryConditions"
    boundary_conditions::BC 
    "A tuple of sources, each of type AbstractSource"
    sources::S
    " Drivers"
    driver::DT
end


"""
SoilCO2Model{FT}(;
        parameters::soilco2Parameters{FT},
        domain::ClimaLSM.AbstractDomain,
        boundary_conditions::NamedTuple,
        sources::Tuple,
    ) where {FT}

A constructor for `SoilCO2Model`.
"""
function SoilCO2Model{FT}(;
    parameters::soilco2Parameters{FT},
    domain::ClimaLSM.AbstractDomain,
    boundary_conditions::BC,
    sources::Tuple,
    drivers::DT
) where {FT,BC, DT}
    args = (parameters, domain, boundary_conditions, sources, drivers)
    SoilCO2Model{FT, typeof.(args)...}(args...)
end

ClimaLSM.name(model::SoilCO2Model) = :soilco2;
ClimaLSM.domain(model::SoilCO2Model) = :subsurface

# 3. Prognostic and Auxiliary variables 

ClimaLSM.prognostic_vars(::SoilCO2Model) = (:C,) # pCO2 in soil, [ppm] # stored in Y.soilco2.C
ClimaLSM.prognostic_types(::SoilCO2Model{FT}) where {FT} = (FT,)
ClimaLSM.auxiliary_vars(::SoilCO2Model) = (:D, :Sₘ, :Sᵣ) # Diffusivity, Source (microbe + root) # stored in p.soilco2.D
ClimaLSM.auxiliary_types(::SoilCO2Model{FT}) where {FT} = (FT, FT, FT)


# 4. RHS

"""
    make_rhs(model::SoilCO2Model)

An extension of the function `make_rhs`, for the soilco2 equation.
This function creates and returns a function which computes the entire
right hand side of the PDE for `C`, and updates `dY.soil.C` in place
with that value.
This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_rhs(model::SoilCO2Model)
    function rhs!(dY, Y, p, t)
        z = ClimaCore.Fields.coordinate_field(model.domain.space).z
        Δz_top, Δz_bottom = get_Δz(z)

        top_flux_bc = boundary_flux(
            model.boundary_conditions.CO2.top,
            TopBoundary(),
            Δz_top,
            Y,
            p,
            t,
        )
        bot_flux_bc = boundary_flux(
            model.boundary_conditions.CO2.bottom,
            BottomBoundary(),
            Δz_bottom,
            Y,
            p,
            t,
        )

        interpc2f = ClimaCore.Operators.InterpolateC2F()
        gradc2f_C = ClimaCore.Operators.GradientC2F() # set a BC on state
        divf2c_C = ClimaCore.Operators.DivergenceF2C(
            top = ClimaCore.Operators.SetValue(ClimaCore.Geometry.WVector.(top_flux_bc)),
            bottom = ClimaCore.Operators.SetValue(ClimaCore.Geometry.WVector.(bot_flux_bc)),
        ) # -∇ ⋅ (-D∇C), where -D∇C is a flux of C02. ∇C point in direction of increasing C, so the flux is - this.
        @. dY.soilco2.C =
	                 -divf2c_C(-interpc2f(p.soilco2.D)*gradc2f_C(Y.soilco2.C))

        # Source terms are added in here
        for src in model.sources
            source!(dY, src, Y, p, model.parameters)
        end

    end
    return rhs!
end

"""
    AbstractCarbonSource

An abstract type for soil CO2 sources. There are two sources:
roots and microbes, in struct RootProduction and MicrobeProduction.
"""
abstract type AbstractCarbonSource <: ClimaLSM.AbstractSource end

"""
    RootProduction <: AbstractCarbonSource

Struct for the root production of CO2, which appears as a source
term in the differential equation.
"""
struct RootProduction <: AbstractCarbonSource end

"""
    MicrobeProduction <: AbstractCarbonSource

Struct for the microbe production of CO2, appearing as a source
term in the differential equation.
"""
struct MicrobeProduction <: AbstractCarbonSource end

"""
    ClimaLSM.source!(dY, src::RootProduction, Y, p, params)

A method which extends the ClimaLSM source! function for the 
case of root production of CO2 in soil.
"""
function ClimaLSM.source!(dY, src::RootProduction, Y, p, params)
    dY.soilco2.C .+= p.soilco2.Sᵣ
end

"""
    ClimaLSM.source!(dY, src::MicrobeProduction, Y, p, params)
   
A method which extends the ClimaLSM source! function for the
case of microbe production of CO2 in soil.
"""
function ClimaLSM.source!(dY, src::MicrobeProduction, Y, p, params)
    dY.soilco2.C .+= p.soilco2.Sₘ
end

# 5. Auxiliary variables

"""
    AbstractSoilDriver

An abstract type for drivers of soil CO2 production.
These are soil temperature, soil moisture,
root carbon, soil organic matter and microbe carbon.
All varying in space (horizontally and vertically) and time.
"""
abstract type AbstractSoilDriver end

"""
    PrescribedSoil <: AbstractSoilDriver

A container which holds the prescribed functions for soil temperature
and moisture. 

This is meant for use when running the biogeochemistry model in standalone mode,
without a prognostic soil model.

$(DocStringExtensions.FIELDS)
"""
struct PrescribedSoil <: AbstractSoilDriver
    "The temperature of the soil, of the form f(t::FT,z::FT) where {FT <: AbstractFloat}"
    temperature::Function 
    "Soil moisture, of the form f(t::FT,z::FT) where {FT <: AbstractFloat}"
    volumetric_liquid_fraction::Function
    "Antecedant soil temperature, of the form f(t::FT,z::FT) where {FT <: AbstractFloat}"
    antecedent_temperature::Function
    "Antecedant soil moisture acting on microbe, of the form f(t::FT,z::FT) where {FT <: AbstractFloat}"
    antecedent_volumetric_liquid_fraction_m::Function
    "Antecedant soil moisture acting on roots, of the form f(t::FT,z::FT) where {FT <: AbstractFloat}"
    antecedent_volumetric_liquid_fraction_r::Function
    "Carbon content of root in soil, of the form f(t::FT,z::FT) where {FT <: AbstractFloat}"
    root_carbon::Function
    "Carbon content of soil organic matter, of the form f(t::FT,z::FT) where {FT <: AbstractFloat}"
    soil_organic_carbon::Function
    "Carbon content of microbes in soil, of the form f(t::FT,z::FT) where {FT <: AbstractFloat}"
    microbe_carbon::Function
end

"""
    antecedent_soil_temperature(driver::PrescribedSoil, p, Y, t, z)

Returns the antecedent soil temperature at location (z) and time (t)
for the prescribed soil case.
"""
function antecedent_soil_temperature(driver::PrescribedSoil, p, Y, t, z)
    return driver.antecedent_temperature.(z, t)
end

"""
    antecedent_soil_moisture_microbes(driver::PrescribedSoil, p, Y, t, z)

Returns the antecedent soil moisture (for microbes) at location (z) and time (t) for the prescribed
soil case.
"""
function antecedent_soil_moisture_microbes(driver::PrescribedSoil, p, Y, t, z)
    return driver.antecedent_volumetric_liquid_fraction_m.(z, t)
end

"""
    antecedent_soil_moisture_roots(driver::PrescribedSoil, p, Y, t, z)

Returns the antecedent soil moisture (for roots) at location (z) and time (t) for the prescribed
soil case.
"""
function antecedent_soil_moisture_roots(driver::PrescribedSoil, p, Y, t, z)
    return driver.antecedent_volumetric_liquid_fraction_r.(z, t)
end

"""
    soil_temperature(driver::PrescribedSoil, p, Y, t, z)

Returns the soil temperature at location (z) and time (t) for the prescribed
soil case.
"""
function soil_temperature(driver::PrescribedSoil, p, Y, t, z)
    return driver.temperature.(z, t)
end

"""
    soil_moisture(driver::PrescribedSoil, p, Y, t, z)

Returns the soil moisture at location (z) and time (t) for the prescribed
soil case.
"""
function soil_moisture(driver::PrescribedSoil, p, Y, t, z)
    return driver.volumetric_liquid_fraction.(z, t)
end

"""
    soil_root_carbon(driver::PrescribedSoil, p, Y, t, z)

Returns the soil root carbon at location (z) and time (t) for the prescribed
soil case.
"""
function soil_root_carbon(driver::PrescribedSoil, p, Y, t, z)
    return driver.root_carbon.(z, t) # for now, z only, but in time, z and t
end

"""
    soil_som_C(driver::PrescribedSoil, p, Y, t, z)

Returns the carbon soil organic matter (SOM) at location (z) and time (t) for the prescribed
soil case.
"""
function soil_SOM_C(driver::PrescribedSoil, p, Y, t, z)
    return driver.soil_organic_carbon.(z, t)
end

"""
    soil_microbe_carbon(driver::PrescribedSoil, p, Y, t, z)

Returns the soil microbe carbon at location (z) and time (t) for the prescribed
soil case.
"""
function soil_microbe_carbon(driver::PrescribedSoil, p, Y, t, z)
    return driver.microbe_carbon.(z, t)
end

"""
    make_update_aux(model::SoilCO2Model)

An extension of the function `make_update_aux`, for the soilco2 equation. 
This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.
This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_update_aux(model::SoilCO2Model)
    function update_aux!(p, Y, t) 
    params = model.parameters 
    z = ClimaCore.Fields.coordinate_field(model.domain.space).z
    T_soil = soil_temperature(model.driver, p, Y, t, z) 
	θ_l = soil_moisture(model.driver, p, Y, t, z)
    θ_ant_roots = antecedent_soil_moisture_roots(model.driver, p, Y, t, z)
    θ_ant_microbe = antecedent_soil_moisture_microbes(model.driver, p, Y, t, z)
    T_ant_soil = antecedent_soil_temperature(model.driver, p, Y, t, z)

	Cr = soil_root_carbon(model.driver, p, Y, t, z)
	Csom = soil_SOM_C(model.driver, p, Y, t, z)
	Cmic = soil_microbe_carbon(model.driver, p, Y, t, z)

	@. p.soilco2.D = co2_diffusivity(
                                        T_soil,
                                        θ_w,
                                        params,
                                        )

	@. p.soilco2.Sᵣ = root_source(
                                     T_soil,
                                     T_ant_soil,
                                     θ_l,
                                     θ_ant_roots,
                                     Cr,
                                     params,
                                     ) 

	@. p.soilco2.Sₘ = microbe_source(
                                        T_soil,
                                        T_ant_soil,
                                        θ_l,
                                        θ_ant_microbe,
                                        Csom,
                                        Cmic,
                                        params,
                                        )
    end
    return update_aux!
end

"""
    AbstractSoilCO2BC{FT} <: ClimaLSM. AbstractBC{FT}

An abstract type for co2-specific types of boundary conditions.
"""
abstract type AbstractSoilCO2BC{FT} <: ClimaLSM.AbstractBC{FT} end

"""
    SoilCO2FluxBC{FT} <: AbstractSoilCO2BC{FT}

A container holding the CO2 flux boundary condition,
which is a function `f(p,t)`, where `p` is the auxiliary state
vector.
"""
struct SoilCO2FluxBC{FT} <: AbstractSoilCO2BC{FT}
    bc::Function
end

"""
    ClimaLSM.boundary_flux(
        bc::SoilCO2FluxBC{FT},
        boundary::ClimaLSM.AbstractBoundary,
        Δz::ClimaCore.Fields.Field,
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::ClimaCore.Fields.Field where {FT}

A method of ClimaLSM.boundary_flux which returns the soilco2
flux (kg CO2 /m^2/s) in the case of a prescribed flux BC at either the top
or bottom of the domain.
"""
function ClimaLSM.boundary_flux(
    bc::SoilCO2FluxBC{FT},
    boundary::ClimaLSM.AbstractBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::ClimaCore.Fields.Field where {FT}
    return bc.bc(p, t) .+ ClimaCore.Fields.zeros(axes(Δz))
end

"""
    SoilCO2StateBC{FT} <: AbstractSoilCO2BC{FT}

A container holding the CO2 state boundary condition (kg CO2 m−3),
which is a function `f(p,t)`, where `p` is the auxiliary state
vector. 
"""
struct SoilCO2StateBC{FT} <: AbstractSoilCO2BC{FT}
    bc::Function
end

"""
    ClimaLSM.boundary_flux(
    bc::SoilCO2StateBC{FT},
    boundary::ClimaLSM.TopBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
    )::ClimaCore.Fields.Field where {FT}

A method of ClimaLSM.boundary_flux which returns the soilco2
flux in the case of a prescribed state BC at 
 top of the domain.
"""
function ClimaLSM.boundary_flux(
    bc::SoilCO2StateBC{FT},
    boundary::ClimaLSM.TopBoundary,
    Δz::ClimaCore.Fields.Field,
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::ClimaCore.Fields.Field where {FT}
    p_len = Spaces.nlevels(axes(p.soilco2.D))
    D_c = Fields.level(p.soilco2.D, p_len)
    C_c = Fields.level(Y.soilco2.C, p_len)
    C_bc = bc.bc(p, t)
    return ClimaLSM.diffusive_flux(D_c, C_bc, C_c, Δz)
end

"""
    ClimaLSM.boundary_flux(
        bc::SoilCO2StateBC{FT},
        boundary::ClimaLSM.BottomBoundary,
        Δz::ClimaCore.Fields.Field,
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::ClimaCore.Fields.Field where {FT}

A method of ClimaLSM.boundary_flux which returns the soilco2
flux in the case of a prescribed state BC at 
bottom of the domain.
"""
function ClimaLSM.boundary_flux(
    bc::SoilCO2StateBC,
    ::ClimaLSM.BottomBoundary,
    Δz,
    Y,
    p,
    t,
)::ClimaCore.Fields.Field
    D_c = Fields.level(p.soilco2.D, 1)
    C_c = Fields.level(Y.soilco2.C, 1)
    C_bc = bc.bc(p, t)
    return ClimaLSM.diffusive_flux(D_c, C_c, C_bc, Δz)
end

include("./co2_parameterizations.jl") # is this needed?

end # module
