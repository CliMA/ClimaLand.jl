# This script aim to run DETECT as a standalone (not inside ClimaLSM.jl)
# I use ClimaLSM.jl tutorial and Katherine ppt 

# Required deps, from tutorial
using OrdinaryDiffEq: ODEProblem, solve, RK4
using SciMLBase
using Plots
using ClimaCore
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains
import ClimaLSM: name, make_rhs, prognostic_vars, prognostic_types
import ClimaLSM.Domains: coordinates, Column
using UnPack

include("./DETECTModel_auxiliary.jl")

# 1. Parameters  

"""
    DETECTParameters{FT <: AbstractFloat}
A struct for storing parameters of the `DETECTModel`.
"""
struct DETECTParameters{FT <: AbstractFloat}
  # parameters (drivers?) that vary with depth and time
  #  "Soil moisture (m³ m⁻³) - vary in time and depth"
  #  θ::FT # driver?
  #  "Soil temperature (Kelvin) - vary in time and depth"
  #  Tₛ::FT # driver?
  #  "Antecedent soil moisture for root - vary in time and depth"
  #  θₐᵣ::FT # ?
  #  "Antecedent soil moisture for microbe - vary in time and depth"
  #  θₐₘ::FT
  #  "Antecedent soil temperature - vary in time and depth"
  #  Tₛₐ::FT

  # parameters that vary with depth
    "Pressure (kPa)"
    Pz::FT
    "Root biomass carbon at depth z (!! should be a function of z) (mg C cm⁻²)"
    Cᵣ::FT
    "Soil organic C content at depth z (!! should be a function of z) (mg C cm⁻²)"
    Csom::FT
    "Microbial C pool at depth z (!! should be a function of z) (mg C cm⁻²)"
    Cmic::FT
#Cmic(z) = exp(-z/a)
#function microbial_carbon(params, z)
#    return exp(-z/params.a)
#end

  # root submodel parameters
    #"Total root biomass C in a 1 m deep by 1 cm² soil column (mg C cm⁻²)" --> this is Cᵣ, see section S1 paper
    #Rᶜ::FT
    "Root mass-base respiration rate at 10°C and mean environmental conditions (mg C cm⁻³ h⁻¹)"
    Rᵦ::FT
    "The effect of soil water content (θ) on root respiration (unitless)"
    α₁ᵣ::FT
    "The effect of antecedent θ on root respiration (unitless)"
    α₂ᵣ::FT
    "The interactive effect of θ and antecedent θ on root respiration (unitless)"
    α₃ᵣ::FT

  # microbial submodel parameters
    "Total soil organic C in a 1 m deep by 1 cm² soil column (mg C cm⁻²)"
    Sᶜ::FT
    "Total microbial biomass C in a 1 m deep by 1 cm² column of soil (mg C cm⁻²)"
    Mᶜ::FT
    "Value of Vₘₐₓ at 10°C and mean environmental conditions (mg C cm⁻³ h⁻¹)"
    Vᵦ::FT
    "The effect of soil water content (θ) on microbial respiration (unitless)"
    α₁ₘ::FT
    "The effect of antecedent θ on microbial respiration (unitless)"
    α₂ₘ::FT
    "The interactive effect of θ and antecedent θ on microbial respiration (unitless)"
    α₃ₘ::FT
    "Michaelis-Menten half-saturation constant (mg C cm⁻³ h⁻¹)"
    Kₘ::FT
    "Microbial carbon-use efficiency (mg C mg⁻¹ C⁻¹)"
    CUE::FT
    "Fraction of soil organic C that is soluble (-)"
    pf::FT
    "Diffusivity of soil C substrate in liquid (unitless)"
    Dₗᵢ::FT

  # shared parameters between root/microbial submodels
    "Temperature sensitivity parameter, somewhat analogous to an energy activation (Kelvin)"
    E₀ₛ::FT 
    "Temperature sensitivity-related parameter (Kelvin)"
    T₀::FT
    "The effect of antecedent soil temperature on root and microbial respiration (unitless)"
    α₄::FT
    "Reference temperature (Kelvin)"
    Tᵣₑ::FT

  # soil CO₂ diffusivity submodel parameters
    "Absolute value of the slope of the line relating log(ψ) versus log(θ) (unitless)"
    α₅::FT
    "Soil bulk density (g cm⁻³)"
    BD::FT
    "Air-filled porosity at soil water potential of -100 cm H₂O (~ 10 kPa) (%)"
    ϕ₁₀₀::FT
    "Particle density"
    PD::FT
    "Diffusion coefficient for CO₂ in air at standard temperature and pressure (m² s⁻¹)"
    Dstp::FT
    "Standard pressure (kPa)"
    P₀::FT
    "Parameter related to the pore size distribution of the soil (unitless)"
    b::FT
end

function DETECTParameters(;
    Pz::FT,
    Cᵣ::FT,
    Csom::FT,
    Cmic::FT,
    Rᵦ::FT,
    α₁ᵣ::FT,
    α₂ᵣ::FT,
    α₃ᵣ::FT,
    Sᶜ::FT,
    Mᶜ::FT,
    Vᵦ::FT,
    α₁ₘ::FT,
    α₂ₘ::FT,
    α₃ₘ::FT,
    Kₘ::FT,
    CUE::FT,
    pf::FT,
    Dₗᵢ::FT,
    E₀ₛ::FT, 
    T₀::FT,
    α₄::FT,
    Tᵣₑ::FT,
    α₅::FT,
    BD::FT,
    ϕ₁₀₀::FT,
    PD::FT,
    Dstp::FT,
    P₀::FT,
    b::FT,
) where {FT}
    return DETECTParameters{FT}(Pz, Cᵣ, Csom, Cmic, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, pf, Dₗᵢ, E₀ₛ, T₀, α₄, α₅, BD, ϕ₁₀₀, PD, Dstp, P₀, b)
end


# 2. Model 

"""
    DETECTModel

A model for simulating the production and transport of CO₂ in the soil with dynamic
source and diffusion terms.
"""
struct DETECTModel{FT, PS, D, BC, S, DT} <: AbstractModel{FT}
    "the parameter set"
    parameters::PS # in constructor you could enforce that it is ::DETECTParameters{FT}
    "the soil domain, using ClimaCore.Domains"
    domain::D
    "the boundary conditions, of type AbstractSoilBoundaryConditions"
    boundary_conditions::BC # maybe also have an FT
    "A tuple of sources, each of type AbstractSoilSource"
    sources::S
    " Drivers"
    driver::DT
end

ClimaLSM.name(model::DETECTModel) = :DETECT;


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
        
	top_flux_bc, bot_flux_bc =
            boundary_fluxes(model.boundary_conditions, p, t)

        interpc2f = ClimaCore.Operators.InterpolateC2F()
        gradc2f_C = ClimaCore.Operators.GradientC2F()
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
    temperature::Function # (t,z) -> exp(-z)*sin(t), or e.g. a spline fit to data
    volumetric_liquid_fraction::Function
    antecedent_temperature::Function
    antecedent_volumetric_liquid_fraction_m::Function
    antecedent_volumetric_liquid_fraction_r::Function

end

function antecedent_soil_temperature(driver::PrescribedSoil, p, Y, t, z)
    return driver.antecedent_temperature(t, z)
end

function antecedent_soil_moisture_microbes(driver::PrescribedSoil, p, Y, t, z)
    return driver.antecedent_volumetric_liquid_fraction_m(t, z)
end

function antecedent_soil_moisture_roots(driver::PrescribedSoil, p, Y, t, z)
    return driver.antecedent_volumetric_liquid_fraction_r(t, z)
end

function soil_temperature(driver::PrescribedSoil, p, Y, t, z)
    return driver.temperature(t, z)
end

function soil_moisture(driver::PrescribedSoil, p, Y, t, z)
    return driver.volumetric_liquid_fraction(t, z)
end


# function for θₐᵣ(θ, t, z)
# function for θₐₘ(θ, t, z)
# function for Tₛₐ(Tₛ, t, z)
# function for P(z)
# function for Cᵣ(z)
# function for Csom(z)
# function for Cmic(z)


#=
function soil_moisture(driver::PrognosticSoil, p, Y, t, z)
    return Y.soil.ϑ_l
end

function soil_temperature(driver::PrognosticSoil, p, Y, t, z)
    return p.soil.T
end
=#


"""
    make_update_aux(model::DETECTModel)

An extension of the function `make_update_aux`, for the DETECT equation. 
This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_update_aux(model::DETECTModel)
    function update_aux!(p, Y, t) 
        @unpack Pz, Cᵣ, Csom, Cmic, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, pf, Dₗᵢ, E₀ₛ, T₀, α₄, α₅, BD, ϕ₁₀₀, PD, Dstp, P₀, b = model.parameters
        z = ClimaCore.Fields.coordinate_field(model.domain.space).z
        Tₛ = soil_temperature(model.driver, p, Y, t, z)
	    θ = soil_moisture(model.driver, p, Y, t, z)
        θₐᵣ = antecedent_soil_moisture_roots(model.driver, p, Y, t, z)
        θₐₘ = antecedent_soil_moisture_microbes(model.driver, p, Y, t, z)
        Tₛₐ = antecedent_soil_temperature(model.driver, p, Y, t, z)

	@. p.DETECT.D = CO2_diffusivity(
					ϕ₁₀₀,
					b,
					diffusion_coefficient(Dstp, Tₛ, T₀, P₀, Pz
                ),
					air_filled_soil_porosity(BD, PD, θ),
					) 

	@. p.DETECT.Sᵣ = root_source(
				     Rᵦ,
				     Cᵣ,
				     root_θ_adj(θ, θₐᵣ, α₁ᵣ, α₂ᵣ, α₃ᵣ), # fᵣ
				     temp_adj( # g
					      Tᵣₑ,
					      T₀,
					      Tₛ,
				              energy_act(E₀ₛ, α₄, Tₛₐ), # E₀
					      ),
    				     )
	@. p.DETECT.Sₘ = microbe_source(
					Kₘ,
					Cmic,
					CUE,
					fVmax( # Vmax
					     Vᵦ,
					     microbe_θ_adj(θ, θₐₘ, α₁ₘ, α₂ₘ, α₃ₘ), # fₘ
					     temp_adj( # g
					              Tᵣₑ,
					              T₀,
					              Tₛ,
				                      energy_act(E₀ₛ, α₄, Tₛₐ), # E₀
					              ),
					       ),
					fCsol(Csom, Dₗᵢ, θ, pf), # Csol
					)

    end
    return update_aux!
end


struct FluxBC
    top::Function
    bottom::Function
end

function boundary_fluxes(f::FluxBC, p, t)
    return f.top(t), f.bottom(t)
end



