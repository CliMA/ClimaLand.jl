# code to run DETECTModel
using UnPack
using DiffEqCallbacks
using OrdinaryDiffEq: ODEProblem, solve, Euler
#import CLIMAParameters as CP
using ClimaCore
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.DETECT: DETECTModel, DETECTParameters, RootProduction, MicrobeProduction, PrescribedSoil, SoilCO2FluxBC, SoilCO2StateBC

#import ClimaLSM
#import ClimaLSM.Parameters as LSMP
#include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float32
include("./DETECTModel_params.jl")
nelems = 50 # number of layers in the vertical
zmin = FT(-1) # 0 to 1 m depth
zmax = FT(0.0)
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
top_bc = SoilCO2StateBC{FT}((p, t) -> FT(3000.0))
bot_bc = SoilCO2StateBC{FT}((p, t) -> FT(100.0))
sources = (RootProduction(),MicrobeProduction()) 
boundary_conditions = (; CO2 = (top = top_bc, bottom = bot_bc))
params = DETECTParameters{FT}(Pz, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, pf, Dₗᵢ, E₀ₛ, T₀, α₄, Tᵣₑ,	α₅, BD, ϕ₁₀₀, PD, Dstp, P₀, b)


θ(t, z) = FT(0.3) # exp(-z)*sin(t)*100 # for now

# plt = lineplot([cos, sin], -π/2, 2π)
# Tₛ(t) = 10*sin(.0001t) + 290
# plt = lineplot(Tₛ, 0, 86400) # 86400 seconds in a day...

Tₛ(t, z) = FT(303.0) #exp(-z)*sin(t) #
θₐᵣ(t, z) = FT(0.3) #exp(-z)*sin(t) #
θₐₘ(t, z) = FT(0.3) #exp(-z)*sin(t) #
Tₛₐ(t, z) = FT(303.0) #exp(-z)*sin(t)*100 + 273
Cᵣ(t, z) = FT(15.0) # FT(15.0*(z+1)/1)
Csom(t, z) = FT(25.0) # FT(25.0*(z+1)/1)
Cmic(t, z) = FT(1.0) # FT(1.0*(z+1)/1)

soil_drivers = PrescribedSoil(Tₛ, θ, Tₛₐ, θₐᵣ, θₐₘ, Cᵣ, Csom, Cmic)
model = DETECTModel{FT}(; parameters = params, domain = soil_domain, sources = sources, boundary_conditions= boundary_conditions, drivers = soil_drivers)

Y, p, coords = initialize(model)

# Initial conditions
function init_DETECT!(Y, z, params)
	function CO2_profile(
		z::FT,
		params::DETECTParameters{FT},
	) where {FT}
		@unpack Pz, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, pf, Dₗᵢ, E₀ₛ, T₀, α₄, Tᵣₑ, α₅, BD, ϕ₁₀₀, PD, Dstp, P₀, b = params
		C = FT(2900.0)*(z+FT(1))/FT(1) + 100 
		return FT(C)
	end
	Y.DETECT.C .= CO2_profile.(z, Ref(params))
end

init_DETECT!(Y, coords.z, model.parameters)

DETECT_ode! = make_ode_function(model)

t0 = FT(0)
tf = FT(100) # why not
dt = FT(1)
saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), saved_values) #?
prob = ODEProblem(DETECT_ode!, Y, (t0, tf), p)
sol = solve(prob, Euler(); dt = dt, callback = cb) # RK4() bugs because t is Float64? (forced by DifferentialEquations.jl?)

# t = sol.t
# state = [parent(sol.u[k].DETECT.C) for k in 1:length(sol.t)]
# source_r = [parent(saved_values.saveval[k].DETECT.Sᵣ) for k in 1:length(sol.t)]
# diffusivity = [parent(saved_values.saveval[k].DETECT.D) for k in 1:length(sol.t)]
# one idea for estimating flux
# Δz_srf = -parent(z)[end]
# estimated_surface_flux = - diffusivity .* (c_atm.(t) .- state) ./ Δz_srf # top of soil to atmosphere carbon at surface of soil

