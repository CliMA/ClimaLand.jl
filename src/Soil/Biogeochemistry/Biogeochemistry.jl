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
    AbstractBC
export DETECTParameters, DETECTModel, PrescribedSoil, RootProduction, MicrobeProduction, SoilCO2FluxBC, SoilCO2StateBC

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




#################################

"""
    DETECTModel

A model for simulating the production and transport of CO₂ in the soil with dynamic
source and diffusion terms.
"""
struct DETECTModel{FT, PS, D, BC, S, DT} <: ClimaLSM.AbstractModel{FT}
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
DETECTModel{FT}(;
        parameters::DETECTParameters{FT},
        domain::ClimaLSM.AbstractDomain,
        boundary_conditions::NamedTuple,
        sources::Tuple,
    ) where {FT}

A constructor for `DETECTModel`.
"""
function DETECTModel{FT}(;
    parameters::DETECTParameters{FT},
    domain::ClimaLSM.AbstractDomain,
    boundary_conditions::BC,
    sources::Tuple,
    drivers::DT
) where {FT,BC, DT}
    args = (parameters, domain, boundary_conditions, sources, drivers)
    DETECTModel{FT, typeof.(args)...}(args...)
end

ClimaLSM.name(model::DETECTModel) = :DETECT;
ClimaLSM.domain(model::DETECTModel) = :subsurface

# 3. Prognostic and Auxiliary variables 

ClimaLSM.prognostic_vars(::DETECTModel) = (:C,) # pCO2 in soil, [ppm] # stored in Y.DETECT.C
ClimaLSM.prognostic_types(::DETECTModel{FT}) where {FT} = (FT,)
ClimaLSM.auxiliary_vars(::DETECTModel) = (:D, :Sₘ, :Sᵣ) # Diffusivity, Source (microbe + root) # stored in p.DETECT.D
ClimaLSM.auxiliary_types(::DETECTModel{FT}) where {FT} = (FT, FT, FT)


# 4. RHS

"""
    make_rhs(model::DETECTModel)

An extension of the function `make_rhs`, for the DETECT equation.
This function creates and returns a function which computes the entire
right hand side of the PDE for `C`, and updates `dY.soil.C` in place
with that value.
This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_rhs(model::DETECTModel)
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
        @. dY.DETECT.C =
	                 -divf2c_C(-interpc2f(p.DETECT.D)*gradc2f_C(Y.DETECT.C))

        # Source terms are added in here
        for src in model.sources
            source!(dY, src, Y, p, model.parameters)
        end

    end
    return rhs!
end

abstract type AbstractCarbonSource end
struct RootProduction end
struct MicrobeProduction end

function source!(dY, src::RootProduction, Y, p, params)
    dY.DETECT.C .+= p.DETECT.Sᵣ
end

function source!(dY, src::MicrobeProduction, Y, p, params)
    dY.DETECT.C .+= p.DETECT.Sₘ
end

# 5. Auxiliary variables

abstract type AbstractSoilDriver end

struct PrescribedSoil <: AbstractSoilDriver
    # functions of z and t
    temperature::Function # (t,z) -> exp(-z)*sin(t), or e.g. a spline fit to data
    volumetric_liquid_fraction::Function
    antecedent_temperature::Function
    antecedent_volumetric_liquid_fraction_m::Function
    antecedent_volumetric_liquid_fraction_r::Function

    # functions of z only (also t in the future)
    rootC::Function
    SOMC::Function
    microbeC::Function
end

function antecedent_soil_temperature(driver::PrescribedSoil, p, Y, t, z)
    return driver.antecedent_temperature.(z, t)
end

function antecedent_soil_moisture_microbes(driver::PrescribedSoil, p, Y, t, z)
    return driver.antecedent_volumetric_liquid_fraction_m.(z, t)
end

function antecedent_soil_moisture_roots(driver::PrescribedSoil, p, Y, t, z)
    return driver.antecedent_volumetric_liquid_fraction_r.(z, t)
end

function soil_temperature(driver::PrescribedSoil, p, Y, t, z)
    return driver.temperature.(z, t)
end

function soil_moisture(driver::PrescribedSoil, p, Y, t, z)
    return driver.volumetric_liquid_fraction.(z, t)
end

function soil_rootC(driver::PrescribedSoil, p, Y, t, z)
    return driver.rootC.(z, t) # for now, z only, but in time, z and t
end

function soil_SOMC(driver::PrescribedSoil, p, Y, t, z)
    return driver.SOMC.(z, t)
end

function soil_microbeC(driver::PrescribedSoil, p, Y, t, z)
    return driver.microbeC.(z, t)
end

"""
    make_update_aux(model::DETECTModel)

An extension of the function `make_update_aux`, for the DETECT equation. 
This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.
This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_update_aux(model::DETECTModel)
    function update_aux!(p, Y, t) 
        @unpack P_sfc, 
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
        b
        = model.parameters
        params = model.parameters # do we need to unpack params above?
        z = ClimaCore.Fields.coordinate_field(model.domain.space).z
        T_soil = soil_temperature(model.driver, p, Y, t, z) 
	θ_l = soil_moisture(model.driver, p, Y, t, z)
        θ_ant_roots = antecedent_soil_moisture_roots(model.driver, p, Y, t, z)
        θ_ant_microbe = antecedent_soil_moisture_microbes(model.driver, p, Y, t, z)
        T_ant_soil = antecedent_soil_temperature(model.driver, p, Y, t, z)

	Cr = soil_rootC(model.driver, p, Y, t, z)
	Csom = soil_SOMC(model.driver, p, Y, t, z)
	Cmic = soil_microbeC(model.driver, p, Y, t, z)

	@. p.DETECT.D = co2_diffusivity(
                                        T_soil,
                                        θ_w,
                                        params,
                                        )

	@. p.DETECT.Sᵣ = root_source(
                                     T_soil,
                                     T_ant_soil,
                                     θ_l,
                                     θ_ant_roots,
                                     Cr,
                                     params,
                                     ) 

	@. p.DETECT.Sₘ = microbe_source(
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


struct SoilCO2FluxBC{FT} <: AbstractSoilCO2BC{FT}
    bc::Function
end

"""
ClimaLSM.boundary_flux

test
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

struct SoilCO2StateBC{FT} <: AbstractSoilCO2BC{FT}
    bc::Function
end

function ClimaLSM.boundary_flux(
    bc::SoilCO2StateBC,
    ::ClimaLSM.TopBoundary,
    Δz,
    Y,
    p,
    t,
)::ClimaCore.Fields.Field
    p_len = Spaces.nlevels(axes(p.DETECT.D))
    D_c = Fields.level(p.DETECT.D, p_len)
    C_c = Fields.level(Y.DETECT.C, p_len)
    C_bc = bc.bc(p, t)
    return ClimaLSM.diffusive_flux(D_c, C_bc, C_c, Δz)
end

function ClimaLSM.boundary_flux(
    bc::SoilCO2StateBC,
    ::ClimaLSM.BottomBoundary,
    Δz,
    Y,
    p,
    t,
)::ClimaCore.Fields.Field
    D_c = Fields.level(p.DETECT.D, 1)
    C_c = Fields.level(Y.DETECT.C, 1)
    C_bc = bc.bc(p, t)
    return ClimaLSM.diffusive_flux(D_c, C_c, C_bc, Δz)
end

include("./co2_parameterizations.jl") # is this needed?

end # module
