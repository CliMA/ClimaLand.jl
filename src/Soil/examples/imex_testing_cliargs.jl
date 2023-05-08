"""
Usage:
    Run this file using julia with the following args after the file name:
1. "imp" or "exp" - use the implicit or explicit solver
2. "true_picard" or "mod_picard" - specify the Jacobian approx for the solver to use
3. "iters_n" - use NewtonsMethod with `max_iters = n`
4. "dt_n" - use a timestep of `n` seconds
5. "top_flux_n" - use a flux BC of n at the top of the soil (neg. is into soil)

Ex: to run the implicit solver using modified Picard with 1 iteration at each
timestep, and a timestep of 1 second, run:
julia --project src/Soil/examples/imex_testing_cliargs.jl imp mod_picard iters_1 dt_1000 top_flux_-1e-7

Note that when using the explicit solver, the second and third
arguments are irrelevant.
"""

using ClimaCore
using ClimaCore: Operators, Spaces, Fields, Geometry
using LinearAlgebra
using DocStringExtensions
using NVTX
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

include("./TridiagonalJacobian.jl")
using .TridiagonalJacobian:
    TridiagonalW, make_Wfact, make_implicit_tendency, explicit_tendency!

FT = Float64

is_imex_CTS_algo(::CTS.IMEXAlgorithm) = true
is_imex_CTS_algo(::DiffEqBase.AbstractODEAlgorithm) = false

is_rosenbrock(::ODE.Rosenbrock23) = true
is_rosenbrock(::ODE.Rosenbrock32) = true
is_rosenbrock(::DiffEqBase.AbstractODEAlgorithm) = false
use_transform(ode_algo) =
    !(is_imex_CTS_algo(ode_algo) || is_rosenbrock(ode_algo))


# parse command line arguments
if isinteractive()
    explicit = false
    is_true_picard = false
    max_iters = 1
    dt = FT(10)
    flux_in = FT(-1e-7)
else
    explicit = ARGS[1] == "exp"
    is_true_picard = (ARGS[2] == "true_picard")
    max_iters = parse(Int64, ARGS[3][7:end])
    dt = parse(FT, replace(ARGS[4][4:end], "p" => "."))
    flux_in = parse(FT, ARGS[5][10:end])
end

# set up timestepper
stepper = CTS.ARS111()
norm_condition = CTS.MaximumError(Float64(1e-6))
conv_checker = CTS.ConvergenceChecker(; norm_condition)

ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = max_iters,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        # convergence_checker = conv_checker,
        # verbose = CTS.Verbose()
    ),
)

t_start = FT(0)
t_end = FT(1e6)

# van Genuchten parameters for clay (from Bonan 2019 supplemental program 8.2)
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

soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems);

# specify boundary conditions (currently constant fluxes)
top_bc = Soil.FluxBC((p, t) -> eltype(t)(flux_in))
flux_out = FT(0)
bot_bc = Soil.FluxBC((p, t) -> eltype(t)(flux_out))

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
@. Y.soil.ϑ_l = FT(0.24)


if !explicit # mixed implicit/explicit case
    transform = use_transform(ode_algo)

    W =
        is_true_picard ? IdentityW(Ref(zero(FT)), transform) :
        TridiagonalW(Y, transform)

    implicit_tendency! = make_implicit_tendency(soil)
    Wfact! = make_Wfact(soil)
    Wfact!(W, Y, p, dt, t_start)

    jac_kwargs = if use_transform(ode_algo)
        (; jac_prototype = W, Wfact_t = Wfact!)
    else
        (; jac_prototype = W, Wfact = Wfact!)
    end

    implicit_problem = ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = explicit_tendency!,
            T_imp! = ODEFunction(implicit_tendency!; jac_kwargs...),
            dss! = ClimaLSM.Soil.dss!,
        ),
        Y,
        (t_start, t_end),
        p,
    )

    integrator = init(
        implicit_problem,
        ode_algo;
        dt = dt,
        adaptive = false,
        progress = true,
        saveat = t_start:1:t_end,
    )
else # explicit-only case
    implicit_tendency! = make_implicit_tendency(soil)
    explicit_problem = ODE.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = implicit_tendency!,
            T_imp! = nothing,
            dss! = ClimaLSM.Soil.dss!,
        ),
        Y,
        (t_start, t_end),
        p,
    )
    integrator = init(
        explicit_problem,
        ode_algo;
        dt = dt,
        adaptive = false,
        progress = true,
        saveat = t_start:1000:t_end,
    )
end

ODE.solve!(integrator)

# calculate water mass balance over entire simulation
mass_end = sum(integrator.sol.u[end].soil.ϑ_l)
mass_start = sum(integrator.sol.u[1].soil.ϑ_l)
t_sim = integrator.sol.t[end] - integrator.sol.t[1]
# flux changes water content every timestep (assumes constant flux_in, flux_out)
mass_change_exp = -(flux_in - flux_out) * t_sim
mass_change_actual = mass_end - mass_start
relerr = abs(mass_change_actual - mass_change_exp) / mass_change_exp * 100 # %


sol_vals = parent(integrator.sol.u[end].soil.ϑ_l)
zs = parent(Fields.coordinate_field(axes(integrator.sol.u[end].soil.ϑ_l)))

isdir("imex_output") ? nothing : mkpath("imex_output")
filename = "imex_output/"
for a in ARGS
    a = string(a)
    replace(a, "." => "p")
    global filename *= string(a) * "-"
end
filename = chop(filename)

if explicit
    plot_name = "explicit"
else
    plot_name = "implicit"
end

if flux_in > FT(-1e-7)
    x_lims = (0.2, 0.7)
    xticks = 0.2:0.1:0.7
else
    x_lims = (0.2, 0.5)
    xticks = 0.2:0.05:0.5
end

tstart_vals = parent(integrator.sol.u[1].soil.ϑ_l)
plot(
    tstart_vals[51:150],
    zs[51:150],
    xlims = x_lims,
    xticks = xticks,
    yticks = -1.0:0.2:0.0,
    title = plot_name * " sol'n with dt=$dt at t=0, top flux=$flux_in",
    xlabel = "vol. water content (m^3 m^-3)",
    ylabel = "depth (m)",
)
savefig(filename * "-start")

t_length = length(integrator.sol.t)
t1_pos = Int64(round(t_length / 3))
t1 = Int64(integrator.sol.t[t1_pos])
t1_vals = parent(integrator.sol.u[t1_pos].soil.ϑ_l)
plot(
    t1_vals[51:150],
    zs[51:150],
    xlims = x_lims,
    xticks = xticks,
    yticks = -1.0:0.2:0.0,
    title = plot_name * " sol'n with dt=$dt at t=$t1 s, top flux=$flux_in",
    xlabel = "vol. water content (m^3 m^-3)",
    ylabel = "depth (m)",
)
savefig(filename * "-$t1")


t2_pos = Int64(round(2 * t_length / 3))
t2 = Int64(integrator.sol.t[t2_pos])
t2_vals = parent(integrator.sol.u[t2_pos].soil.ϑ_l)
plot(
    t2_vals[51:150],
    zs[51:150],
    xlims = x_lims,
    xticks = xticks,
    yticks = -1.0:0.2:0.0,
    title = plot_name * " sol'n with dt=$dt at t=$t2 s, top flux=$flux_in",
    xlabel = "vol. water content (m^3 m^-3)",
    ylabel = "depth (m)",
)
savefig(filename * "-$t2")

t_end = Int64(t_end)
tend_vals = parent(integrator.sol.u[end].soil.ϑ_l)
plot(
    tend_vals[51:150],
    zs[51:150],
    xlims = x_lims,
    xticks = xticks,
    yticks = -1.0:0.2:0.0,
    title = plot_name * " sol'n with dt=$dt at t=$t_end s, top flux=$flux_in",
    xlabel = "vol. water content (m^3 m^-3)",
    ylabel = "depth (m)",
)
savefig(filename * "-$t_end")


b = @benchmark ODE.solve!(integrator)

filename *= ".jld2"

jldsave(
    filename;
    tstart_vals = tstart_vals,
    t1_vals = t1_vals,
    t2_vals = t2_vals,
    tend_vals = tend_vals,
    z = zs,
    benchmark = b,
    dt = dt,
    mass_err = relerr,
)

@show filename
