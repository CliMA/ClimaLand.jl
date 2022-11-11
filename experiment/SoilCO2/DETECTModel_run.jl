# code to run DETECTModel
using UnPack
using DiffEqCallbacks
using OrdinaryDiffEq: ODEProblem, solve, RK4
#import CLIMAParameters as CP
using ClimaCore
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.DETECT: DETECTModel, DETECTParameters, RootProduction, MicrobeProduction, PrescribedSoil, FluxBC

#import ClimaLSM
#import ClimaLSM.Parameters as LSMP
#include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float32
include("./DETECTModel_params.jl")
nelems = 50 # number of layers in the vertical
zmin = FT(-5)
zmax = FT(0.0)
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
top_bc = (t) -> FT(0.0)
bot_bc = (t) -> FT(0.0) 
sources = (RootProduction(),MicrobeProduction()) 
bc = FluxBC(top_bc, bot_bc)
params = DETECTParameters{FT}(Pz, Cᵣ, Csom, Cmic, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, pf, Dₗᵢ, E₀ₛ, T₀, α₄, Tᵣₑ,	α₅, BD, ϕ₁₀₀, PD, Dstp, P₀, b)

θ(t, z) = FT(0.3) #exp(-z)*sin(t)*100 # for now
Tₛ(t, z) = FT(303) #exp(-z)*sin(t) #
θₐᵣ(t,z) = FT(0.3) #exp(-z)*sin(t) #
θₐₘ(t,z) = FT(0.3) #exp(-z)*sin(t) #
Tₛₐ(t,z) = FT(303) #exp(-z)*sin(t)*100 + 273

soil_drivers = PrescribedSoil(Tₛ, θ, Tₛₐ, θₐᵣ, θₐₘ)
model = DETECTModel{FT}(; parameters = params, domain = soil_domain, sources = sources, boundary_conditions= bc, drivers = soil_drivers)

Y, p, coords = initialize(model)

# Initial conditions
function init_DETECT!(Y, z, params)
	function CO2_profile(
		z::FT,
		params::DETECTParameters{FT},
	) where {FT}
		@unpack Pz, Cᵣ, Csom, Cmic, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, pf, Dₗᵢ, E₀ₛ, T₀, α₄, Tᵣₑ, α₅, BD, ϕ₁₀₀, PD, Dstp, P₀, b = params
		C = FT(2900.0)*(z+FT(5))/FT(5) + 100 
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
sol = solve(prob, RK4(); dt = dt, callback = cb) # do we want Euler or another algorithm for this?

t = sol.t
state = [parent(sol.u[k].DETECT.C) for k in 1:length(sol.t)]
source_r = [parent(saved_values.saveval[k].DETECT.Sᵣ) for k in 1:length(sol.t)]
diffusivity = [parent(saved_values.saveval[k].DETECT.D) for k in 1:length(sol.t)]
# one idea for estimating flux
Δz_srf = -parent(z)[end]
estimated_surface_flux = - diffusivity .* (c_atm.(t) .- state) ./ Δz_srf # top of soil to atmosphere carbon at surface of soil

