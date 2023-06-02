using ClimaCore
using DiffEqBase
import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS
using StaticArrays
using Plots
using DelimitedFiles
using Statistics
using ArtifactWrappers

using ClimaLSM
using ClimaLSM.Soil
using ClimaLSM.Domains: Column
import ClimaLSM
dir = pkgdir(ClimaLSM)
include(joinpath(dir,"src/Soil/examples/TridiagonalJacobian.jl"))
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
stepper = CTS.ARS111()
norm_condition = CTS.MaximumError(Float64(1e-8))
conv_checker = CTS.ConvergenceChecker(; norm_condition)

ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 50,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
       # verbose = CTS.Verbose()
    ),
)


t_start = FT(0)
t_end = FT(1e6)
dt = FT(1e4) #1000 = 1e3, 10000 = 1e4
explicit = false

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

top_bc = Soil.MoistureStateBC((p, t) -> eltype(t)(ν - 1e-3))
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

if !explicit
    transform = use_transform(ode_algo)

    W = TridiagonalW(Y, transform)

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
else
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

sol = ODE.solve!(integrator)
N = length(sol.t)
ϑ_l = parent(sol.u[N].soil.ϑ_l)[:]
z = parent(coords.z)[:]
bonan_clay_dataset = ArtifactWrapper(
    @__DIR__,
    "richards_clay",
    ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/nk89znth59gcsdb65lnywnzjnuno3h6k.txt",
        filename = "clay_bonan_sp801_22323.txt",
    ),],
)
datapath = get_data_folder(bonan_clay_dataset)
data = joinpath(datapath, "clay_bonan_sp801_22323.txt")
ds_bonan = readdlm(data)
bonan_moisture = reverse(ds_bonan[:, 1])
bonan_z = reverse(ds_bonan[:, 2]) ./ 100.0
plot(ϑ_l, parent(z), label = "Clima")
plot!(bonan_moisture, bonan_z, label = "Bonan's Matlab code")
# dt = 1e-3 passes, dt = 1e-4 does not, but the solution doesnt have NaNs
@show(sqrt.(mean((bonan_moisture .- ϑ_l) .^ 2.0)) < FT(1e-3))
