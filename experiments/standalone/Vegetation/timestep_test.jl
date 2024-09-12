import SciMLBase
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
using CairoMakie
using Statistics
using Dates
using Insolation
using DelimitedFiles

# Load CliMA Packages and ClimaLand Modules:

using ClimaCore
import ClimaParams as CP
import ClimaTimeSteppers as CTS
using StaticArrays
using ClimaLand
using ClimaLand.Domains: Point
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
const FT = Float64;
earth_param_set = LP.LandParameters(FT);
f_root_to_shoot = FT(3.5)
plant_ν = FT(2.46e-4) # kg/m^2
n_stem = Int64(1)
n_leaf = Int64(1)
h_leaf = FT(9.5)
h_stem = FT(9)
compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2]
compartment_surfaces = [FT(0), h_stem, h_stem + h_leaf]
land_domain = Point(; z_sfc = FT(0.0))
include(
    joinpath(pkgdir(ClimaLand), "experiments/integrated/fluxnet/data_tools.jl"),
);
time_offset = 7
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree
atmos_h = FT(32)
site_ID = "US-MOz"
data_link = "https://caltech.box.com/shared/static/7r0ci9pacsnwyo0o9c25mhhcjhsu6d72.csv"

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments/integrated/fluxnet/met_drivers_FLUXNET.jl",
    ),
);

z0_m = FT(2)
z0_b = FT(0.2)

shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
);
ψ_soil0 = FT(0.0)

soil_driver = PrescribedSoil(
    FT;
    root_depths = SVector{10, FT}(-(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0),
    ψ = t -> ψ_soil0,
    α_PAR = FT(0.2),
    α_NIR = FT(0.4),
    T = t -> 298.0,
    ϵ = FT(0.99),
);

rt_params = TwoStreamParameters(
    FT;
    G_Function = ConstantGFunction(FT(0.5)),
    α_PAR_leaf = FT(0.1),
    α_NIR_leaf = FT(0.45),
    τ_PAR_leaf = FT(0.05),
    τ_NIR_leaf = FT(0.25),
    Ω = FT(0.69),
    λ_γ_PAR = FT(5e-7),
    λ_γ_NIR = FT(1.65e-6),
)

rt_model = TwoStreamModel{FT}(rt_params);

cond_params = MedlynConductanceParameters(FT; g1 = FT(141.0))

stomatal_model = MedlynConductanceModel{FT}(cond_params);


photo_params = FarquharParameters(FT, Canopy.C3(); Vcmax25 = FT(5e-5))

photosynthesis_model = FarquharModel{FT}(photo_params);

AR_params = AutotrophicRespirationParameters(FT)
AR_model = AutotrophicRespirationModel{FT}(AR_params);

f_root_to_shoot = FT(3.5)
SAI = FT(1.0)
RAI = FT(3f_root_to_shoot)
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)
rooting_depth = FT(1.0)
function root_distribution(z::T; rooting_depth = rooting_depth) where {T}
    return T(1.0 / rooting_depth) * exp(z / T(rooting_depth)) # 1/m
end

K_sat_plant = FT(1.8e-6)
ψ63 = FT(-4 / 0.0098)
Weibull_param = FT(4)
a = FT(0.05 * 0.0098)

conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)

retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a);

ν = FT(0.7)
S_s = FT(1e-2 * 0.0098)

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = ν,
    S_s = S_s,
    root_distribution = root_distribution,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
);


plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_surfaces = compartment_surfaces,
    compartment_midpoints = compartment_midpoints,
);
ac_canopy = FT(1e3)
energy_model = ClimaLand.Canopy.BigLeafEnergyModel{FT}(
    BigLeafEnergyParameters{FT}(ac_canopy),
)

canopy = ClimaLand.Canopy.CanopyModel{FT}(;
    parameters = shared_params,
    domain = land_domain,
    autotrophic_respiration = AR_model,
    radiative_transfer = rt_model,
    photosynthesis = photosynthesis_model,
    conductance = stomatal_model,
    energy = energy_model,
    hydraulics = plant_hydraulics,
    soil_driver = soil_driver,
    atmos = atmos,
    radiation = radiation,
);


Y, p, coords = ClimaLand.initialize(canopy)
exp_tendency! = make_exp_tendency(canopy)
imp_tendency! = make_imp_tendency(canopy)
jacobian! = make_jacobian(canopy);
jac_kwargs =
    (; jac_prototype = ClimaLand.ImplicitEquationJacobian(Y), Wfact = jacobian!);


ψ_leaf_0 = FT(-2e5 / 9800)
ψ_stem_0 = FT(-1e5 / 9800)

S_l_ini =
    inverse_water_retention_curve.(
        retention_model,
        [ψ_stem_0, ψ_leaf_0],
        ν,
        S_s,
    )

Y.canopy.hydraulics.ϑ_l.:1 .= augmented_liquid_fraction.(ν, S_l_ini[1])
Y.canopy.hydraulics.ϑ_l.:2 .= augmented_liquid_fraction.(ν, S_l_ini[2])


seconds_per_day = 3600 * 24.0
t0 = 150seconds_per_day
N_days = 0.5
tf = t0 + N_days * seconds_per_day
evaluate!(Y.canopy.energy.T, atmos.T, t0)
set_initial_cache! = make_set_initial_cache(canopy)
set_initial_cache!(p, Y, t0);

saveat = Array(t0:(3 * 3600):tf)
updateat = Array(t0:(3600 * 3):tf)
drivers = ClimaLand.get_drivers(canopy)
updatefunc = ClimaLand.make_update_drivers(drivers)
cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

timestepper = CTS.ARS111();
err = (FT == Float64) ? 1e-8 : 1e-4
convergence_cond = CTS.MaximumError(err)
conv_checker = CTS.ConvergenceChecker(norm_condition = convergence_cond)
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        # convergence_checker = conv_checker,
    ),
);

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLand.dss!,
    ),
    Y,
    (t0, tf),
    p,
);

# ref_dt = 6.0
# ref_sol =
#     SciMLBase.solve(prob, ode_algo; dt = ref_dt, callback = cb, saveat = saveat);
# ref_T = [parent(ref_sol.u[k].canopy.energy.T)[1] for k in 1:length(ref_sol.t)]

# Read in solution from saved delimited file (experiment run explicitly with dt = 0.1s)
savedir = joinpath(pkgdir(ClimaLand), "experiments/standalone/Vegetation");
ref_file = joinpath(savedir, "ref_T_dt0p1_0p5days.txt")
ref_T = vec(readdlm(ref_file, ','))

mean_err = []
p95_err = []
p99_err = []
dts = [3600.0, 5400.0, 7200.0, 9000.0, 10800.0, 12600.0, 13500.0, 14040.0]#[12.0, 24.0, 48.0, 100.0, 225.0, 450.0, 900.0, 1800.0, 3600.0, 7200.0, 14400.0]
sol_ends = []
T_states = []
times = []
for dt in dts
    @info dt
    saveat = Array(t0:(3 * 3600):tf)
    evaluate!(Y.canopy.energy.T, atmos.T, t0)
    updateat = Array(t0:(3600 * 0.5):tf)
    updatefunc = ClimaLand.make_update_drivers(drivers)
    cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

    @time sol =
        SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb, saveat = saveat)
    T = [parent(sol.u[k].canopy.energy.T)[1] for k in 1:length(sol.t)]
    ΔT = abs.(T .- ref_T)
    push!(mean_err, mean(ΔT))
    push!(p95_err, percentile(ΔT, 95))
    push!(p99_err, percentile(ΔT, 99))
    push!(sol_ends, T[end])
    push!(T_states, T)
    push!(times, sol.t)
end

# Create plot with statistics
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel = "Timestep (minutes)",
    ylabel = "Temperature (K)",
    xscale = log10,
    yscale = log10,
)
dts = dts ./ 60
lines!(ax, dts, FT.(mean_err), label = "Mean Error")
lines!(ax, dts, FT.(p95_err), label = "95th% Error")
lines!(ax, dts, FT.(p99_err), label = "99th% Error")
axislegend(ax)
save(joinpath(savedir, "errors.png"), fig)

# Create convergence plot
errors = abs.(sol_ends .- ref_T[end])
fig2 = Figure()
ax2 = Axis(
    fig2[1, 1],
    xlabel = "log(dt)",
    ylabel = "log(|T[end] - T_ref[end]|)",
    xscale = log10,
    yscale = log10,
)
scatter!(ax2, dts, FT.(errors))
lines!(ax2, dts, dts)
save(joinpath(savedir, "convergence.png"), fig2)

# Create states plot
fig3 = Figure()
ax3 = Axis(fig3[1, 1], xlabel = "time (min)", ylabel = "T (K)")
times = times ./ 60.0
for i in 1:length(times)
    lines!(ax3, times[i], T_states[i], label = "dt $(dts[i]) min")
end
axislegend(ax3, position = :lt)
save(joinpath(savedir, "states.png"), fig3)
