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
        verbose = CTS.Verbose()
    ),
)


t_start = FT(0)
t_end = FT(1e6)
dt = FT(1e4) #1000 = 1e3, 10000 = 1e4
explicit =  false
truth = readdlm("./src/Soil/examples/explicit_solution.csv", ',')


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
ϑ_l_bc = ν - 1e-3
top_bc = Soil.MoistureStateBC((p, t) -> eltype(t)(ϑ_l_bc))
bot_bc = Soil.FluxBC((p, t) -> eltype(t)(0))
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
ϑ_l_top = [
    parent(integrator.sol.u[k].soil.ϑ_l)[end] for
    k in 1:length(integrator.sol.t)
]
ψ_top = ClimaLSM.Soil.pressure_head.(vg_α, vg_n, vg_m, θ_r, ϑ_l_top, ν, S_s)
K_top = @. ClimaLSM.Soil.hydraulic_conductivity(
    K_sat,
    vg_m,
    effective_saturation(ν, ϑ_l_top, θ_r),
)
ϑ_l_bot = [
    parent(integrator.sol.u[k].soil.ϑ_l)[1] for k in 1:length(integrator.sol.t)
]
K_bot = @. ClimaLSM.Soil.hydraulic_conductivity(
    K_sat,
    vg_m,
    effective_saturation(ν, ϑ_l_bot, θ_r),
)

ψ_bc = ClimaLSM.Soil.pressure_head(vg_α, vg_n, vg_m, θ_r, ϑ_l_bc, ν, S_s)
dz_top = 0.005
QN = @. -K_top / dz_top * (ψ_bc - ψ_top) - K_top;
Q0 = 0#-K_bot
# calculate water mass balance over entire simulation
implicit_sol = integrator.sol.u[end].soil.ϑ_l
mass_end = sum(implicit_sol)
mass_start = sum(integrator.sol.u[1].soil.ϑ_l)
# flux changes water content every timestep (assumes constant flux_in, flux_out)
#mass_change_exp = sum(@. (-(QN[2:end] - Q0[2:end]) * dt))
mass_change_exp = sum(@. (-(QN[2:end]) * dt))
mass_change_actual = mass_end - mass_start
relerr = abs(mass_change_actual - mass_change_exp) / mass_change_exp
rmse = sqrt(mean((truth .- parent(implicit_sol)).^2))
@info dt
@info rmse
@info relerr
#=
BC = ν - ϵ
dts = [10, 25, 100, 250, 1000, 2500, 1e4]
rmses = [4.0358176789400375e-6, 5.645720063701383e-6, 4.292680404253761e-5, 0.00010665890769993408, 0.00043309844085658716, 0.0011132606742645518,0.005153354873715869]
masserrs = [3.234895886551343e-6, 1.1431852057916104e-6, 2.6003855544650854e-7, 4.256542768948846e-7, 5.475376618370726e-8, 3.6814091263920326e-8, 3.989774071947699e-9]


p = Plots.plot()
p_twin = twinx(p)
Plots.plot!(p,dts, rmses, xaxis = :log10, yaxis=:log10, label = "RMSE", color = "red", linewidth = 3, xlabel = "Δt (s)", ylabel = "RMSE ϑ", xlim = [10,1e4])
Plots.plot!(p_twin, dts, masserrs,label = "Vol. error", color = "purple", linewidth = 3, ylabel = "|∑ϑ-∑ϑ(0)|/∫ΔFdt", xlabel = "")
Plots.savefig("./src/Soil/examples/state_bc.png")
=#




#=
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
=#
