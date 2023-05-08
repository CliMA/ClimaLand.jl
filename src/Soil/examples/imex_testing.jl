"""
Usage:
    Run this file using julia with the following args after the file name:
1. "imp" or "exp" - use the implicit or explicit solver
2. "true_picard" or "mod_picard" - specify the method for the solver to use
3. "iters_n" - use NewtonsMethod with `max_iters = n`
4. "dt_n" - use a timestep of `n` seconds

Ex: to run the implicit solver using modified Picard with 1 iteration at each
timestep, and a timestep of 1 second, run:
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_1 dt_1

Note that when using the explicit solver, the second and third
arguments are irrelevant.
"""

using ClimaCore
using ClimaCore: Operators, Spaces, Fields, Geometry
using UnPack
using LinearAlgebra
using DocStringExtensions
using NVTX
using Colors
using DiffEqBase
using BenchmarkTools
using JLD2
using Plots

import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Soil
using ClimaLSM.Domains: HybridBox, SphericalShell, Column


FT = Float64

to_scalar_coefs(vector_coefs) = map(vector_coef -> vector_coef.u₃, vector_coefs)

is_imex_CTS_algo(::CTS.IMEXAlgorithm) = true
is_imex_CTS_algo(::DiffEqBase.AbstractODEAlgorithm) = false

is_implicit(::ODE.OrdinaryDiffEqImplicitAlgorithm) = true
is_implicit(::ODE.OrdinaryDiffEqAdaptiveImplicitAlgorithm) = true
is_implicit(ode_algo) = is_imex_CTS_algo(ode_algo)

is_rosenbrock(::ODE.Rosenbrock23) = true
is_rosenbrock(::ODE.Rosenbrock32) = true
is_rosenbrock(::DiffEqBase.AbstractODEAlgorithm) = false
use_transform(ode_algo) =
    !(is_imex_CTS_algo(ode_algo) || is_rosenbrock(ode_algo))


if isinteractive()
    is_true_picard = false
else
    is_true_picard = (ARGS[2] == "true_picard")
end

# parameters for clay from Bonan 2019 supplemental program 8.2
ν = FT(0.495)
K_sat = FT(0.0443 / 3600 / 100) # m/s
vg_n = FT(1.43)
vg_α = FT(0.026 * 100) # inverse meters
vg_m = FT(1) - FT(1) / vg_n
θ_r = FT(0.124)
S_s = FT(1e-3) #inverse meters

zmax = FT(0)
zmin = FT(-1.5)
nelems = 150

# soil_domain = HybridBox(;
#     zlim = (-10.0, 0.0),
#     xlim = (0.0, 100.0),
#     ylim = (0.0, 100.0),
#     nelements = (10, 10, 10),
#     npolynomial = 1,
#     periodic = (true, true),
# )
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems);

top_bc = Soil.MoistureStateBC((p, t) -> eltype(t)(ν - 1e-3))
bot_bc = Soil.FreeDrainage()
sources = ()
boundary_fluxes = (; top = (water = top_bc,), bottom = (water = bot_bc,))
params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)

soil = Soil.RichardsModel{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = sources,
)

Y, p, coords = initialize(soil)
Y.soil.ϑ_l = FT(0.24)

if isinteractive()
    dt = FT(64)
else
    dt = parse(FT, replace(ARGS[4][4:end], "p" => "."))
end

t_start = FT(0)
t_end = FT(1e6) #FT(256)

if isinteractive()
    max_iters = 1
else
    max_iters = parse(Int64, ARGS[3][7:end])
end

# if isinteractive()
#     convergence_cond = CTS.MaximumRelativeError(1e-10)
# else
#     convergence_cond = CTS.MaximumRelativeError(parse(FT, ARGS[5][13:end]))
# end

# TODO when conv_checker was used before, the solver thought it had converged
#  when it hadn't. For now, don't use a conv checker, but maybe try again later
# conv_checker = CTS.ConvergenceChecker(component_condition = convergence_cond)
conv_checker = nothing
ode_algo = CTS.IMEXAlgorithm(
    CTS.ARS343(),
    CTS.NewtonsMethod(
        max_iters = max_iters,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)

if isinteractive()
    is_imp = true
else
    is_imp = ARGS[1] == "imp"
end

if is_imp # mixed implicit/explicit case
    transform = use_transform(ode_algo)

    W =
        is_true_picard ? IdentityW(Ref(zero(FT)), transform) :
        TridiagonalW(Y, transform)

    implicit_tendency! = make_implicit_tendency(soil)
    explicit_tendency! = make_explicit_tendency(soil)

    jac_kwargs = if use_transform(ode_algo)
        (; jac_prototype = W, Wfact_t = Wfact!)
    else
        (; jac_prototype = W, Wfact = Wfact!)
    end

    problem = ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = explicit_tendency!,
            T_imp! = ODEFunction(implicit_tendency!; jac_kwargs...),
            dss! = ClimaLSM.Soil.dss!,
        ),
        Y,
        (t_start, t_end),
        p,
    )
else # explicit-only case
    ode_function! = make_ode_function(soil)

    problem = ODE.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = ode_function!,
            T_imp! = nothing,
            dss! = ClimaLSM.Soil.dss!,
        ),
        Y,
        (t_start, t_end),
        p,
    )
end

integrator = init(
    problem,
    ode_algo;
    dt = dt,
    adaptive = false,
    progress = true,
    saveat = t_start:1:t_end,
)

sol = ODE.solve!(integrator)

sol_vals = parent(sol.u[end].soil.ϑ_l)
zs = parent(Fields.coordinate_field(axes(sol.u[end].soil.ϑ_l)))

isdir("imex-res") ? nothing : mkpath("imex-res")
filename = "imex-res/rre_sol"
for a in ARGS
    a = string(a)
    replace(a, "." => "p")
    global filename *= "-" * string(a)
end

plot(
    parent(sol.u[1].soil.ϑ_l)[51:150],
    zs[51:150],
    xlims = (0.2, 0.5),
    xticks = 0.2:0.05:0.5,
    yticks = -1.0:0.2:0.0,
    title = "implicit sol'n with dt=$dt at t=0",
    xlabel = "vol. water content (m^3 m^-3)",
    ylabel = "depth (m)",
)
savefig(filename * "-start")

t_end = Int(t_end)
plot(
    parent(sol.u[end].soil.ϑ_l)[51:150],
    zs[51:150],
    xlims = (0.2, 0.5),
    xticks = 0.2:0.05:0.5,
    yticks = -1.0:0.2:0.0,
    title = "implicit sol'n with dt=$dt at t=$t_end s",
    xlabel = "vol. water content (m^3 m^-3)",
    ylabel = "depth (m)",
)
savefig(filename * "-$t_end")


b = @benchmark ODE.solve!(integrator)

filename *= ".jld2"

jldsave(filename; ϑ_l = sol_vals, z = zs, benchmark = b, dt = dt)
@show filename

@assert all(x -> x > θ_r, sol_vals)
@assert all(x -> x < ν * 1.03, sol_vals)
