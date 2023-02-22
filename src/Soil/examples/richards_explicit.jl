using UnPack
using Statistics
using DiffEqCallbacks
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil


# include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float64

saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)

# parameters for sand
# ν = FT(0.495)
# K_sat = FT(0.0443 / 3600 / 100) # m/s
# S_s = FT(1e-3) #inverse meters
# vg_n = FT(2.0)
# vg_α = FT(2.6) # inverse meters
# vg_m = FT(1) - FT(1) / vg_n
# θ_r = FT(0)

# parameters for clay from Bonan 2019 table 8.3 van Genuchten
ν = FT(0.38)
K_sat = FT(0.20 / 3600 / 100) # m/s
vg_n = FT(1.09)
vg_α = FT(0.8) # inverse meters
θ_r = FT(0.068)
vg_m = FT(1) - FT(1) / vg_n
S_s = FT(1e-3) #inverse meters

zmax = FT(0)
zmin = FT(-10)
nelems = 50

soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
# saturated top boundary condition
top_bc = FluxBC((p, t) -> eltype(t)(0.0))
bot_bc = FluxBC((p, t) -> eltype(t)(0.0))
sources = ()
boundary_conds = (; water = (top = top_bc, bottom = bot_bc))
params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)

soil = Soil.RichardsModel{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_conds,
    sources = sources,
)

Y, p, coords = initialize(soil)

# specify ICs
function init_soil!(Ysoil, z, params)
    function hydrostatic_profile(
        z::FT,
        params::RichardsParameters{FT},
    ) where {FT}
        @unpack ν, vg_α, vg_n, vg_m, θ_r = params
        #unsaturated zone only, assumes water table starts at z_∇
        z_∇ = FT(-10)# matches zmin
        S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
        ϑ_l = S * (ν - θ_r) + θ_r
        return FT(ϑ_l)
    end
    Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
end

init_soil!(Y, coords.z, soil.parameters)

soil_ode! = make_ode_function(soil)

t0 = FT(0)
tf = FT(60)
dt = FT(1)
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), saved_values)
prob = ODEProblem(soil_ode!, Y, (t0, tf), p)

using BenchmarkTools
@btime solve(prob, Euler(); dt = dt, callback = cb)
sol = solve(prob, Euler(); dt = dt, callback = cb)
