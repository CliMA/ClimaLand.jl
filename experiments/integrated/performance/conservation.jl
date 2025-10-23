import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaComms
ClimaComms.@import_required_backends
using CairoMakie
using Statistics
using Dates
import ClimaUtilities.TimeManager: ITime
import ClimaUtilities.OutputPathGenerator: generate_output_path

using ClimaLand
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams


climaland_dir = pkgdir(ClimaLand)
FT = Float64
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
start_date = DateTime("2008-03-01")
stop_date = start_date + Year(2)
Δt = 450.0
domain = ClimaLand.Domains.Column(;
    dz_tuple = FT.((3, 0.05)),
    nelements = 15,
    longlat = FT.((-62.0, 3.0)),
    zlim = FT.((-15, 0)),
)

surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
surface_space = domain.space.surface
# Forcing data - high resolution
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict,
    FT;
    max_wind_speed = 25.0,
    use_lowres_forcing = true,
)
forcing = (; atmos, radiation)

# Read in LAI from MODIS data
LAI =
    ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
land = LandModel{FT}(forcing, LAI, toml_dict, domain, Δt)

timestepper = CTS.ARS111()
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 3,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
)
exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land)
jacobian! = make_jacobian(land)
set_initial_cache! = make_set_initial_cache(land)
Y, p, cds = initialize(land)
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!)
set_ic! = ClimaLand.Simulations.make_set_initial_state_from_file(
    ClimaLand.Artifacts.saturated_land_ic_path(),
    land,
)
t0 = ITime(0, Second(1), start_date)
tf = ITime(Second(stop_date - start_date).value, Second(1), start_date)
set_ic!(Y, p, t0, land)
set_initial_cache!(p, Y, t0)
saveat = ITime(Δt)
n_saves = length(t0:saveat:tf)
sv = (;
    t = Array{Float64}(undef, n_saves),
    saveval = Array{NamedTuple}(undef, n_saves),
)
saving_affect! = ClimaLand.SavingAffect(sv, 0)
saving_initialize = (_, _, _, x) -> saving_affect!(x)
saving_cb = ClimaLand.IntervalBasedCallback(
    saveat,
    t0,
    ITime(Int(Δt), Second(1), start_date),
    saving_affect!;
    callback_start = t0,
    initialize = saving_initialize,
)

drivers = ClimaLand.get_drivers(land)
updatefunc = ClimaLand.make_update_drivers(drivers)
updateat = ITime(3600 * 3)
driver_cb = ClimaLand.DriverUpdateCallback(
    updatefunc,
    updateat,
    t0;
    dt = ITime(Int(Δt), Second(1), start_date),
)
prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLand.dss!,
    ),
    Y,
    (t0, tf),
    p,
)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt = ITime(Int(Δt), Second(1), start_date),
    callback = cb,
    adaptive = false,
    saveat = collect(t0:saveat:tf),
);

daily = float.(sol.t[2:end]) ./ 3600 ./ 24
savedir =
    generate_output_path("experiments/integrated/performance/conservation")
# Energy balance
total_energy_array =
    [parent(sv.saveval[k].total_energy)[1] for k in 1:length(sol.t)]
mass_change_actual = total_energy_array[3:end] .- total_energy_array[2];
mass_change_expected =
    [parent(sol.u[k].∫F_e_dt)[1] for k in 2:(length(sol.t) - 1)]
fig = Figure(size = (1500, 400))
ax = Axis(
    fig[1, 1],
    xlabel = "Days",
    ylabel = "Fractional Error",
    yscale = log10,
)


lines!(
    ax,
    daily[2:end],
    eps(FT) .+
    abs.((mass_change_actual - mass_change_expected) ./ mass_change_expected,),
    label = "Energy Balance",
)
total_energy_array =
    [parent(sv.saveval[k].total_energy)[1] for k in 1:2112:length(sol.t)]
mass_change_actual = total_energy_array[2:end] .- total_energy_array[1];
mass_change_expected =
    [parent(sol.u[k].∫F_e_dt)[1] for k in 2:2112:length(sol.t)]
lines!(
    ax,
    daily[1:2112:end][2:end],
    eps(FT) .+
    abs.(
        (mass_change_actual - mass_change_expected[2:end]) ./
        mass_change_expected[2:end],
    ),
    label = "Energy Balance Sparse Save",
)

axislegend(position = :lb)
CairoMakie.save(joinpath(savedir, "energy_conservation.png"), fig)

# Water balance
total_water_array =
    [parent(sv.saveval[k].total_water)[1] for k in 1:length(sol.t)]
mass_change_actual = total_water_array[3:end] .- total_water_array[2];
mass_change_expected =
    [parent(sol.u[k].∫F_vol_liq_water_dt)[1] for k in 2:(length(sol.t) - 1)]
fig2 = Figure(size = (1500, 400))
ax2 = Axis(
    fig2[1, 1],
    xlabel = "Days",
    ylabel = "Fractional Error",
    yscale = log10,
)


lines!(
    ax2,
    daily[2:end],
    eps(FT) .+
    abs.((mass_change_actual - mass_change_expected) ./ mass_change_expected,),
    label = "Water Balance",
)
total_water_array =
    [parent(sv.saveval[k].total_water)[1] for k in 1:2112:length(sol.t)]
mass_change_actual = total_water_array[2:end] .- total_water_array[1];
mass_change_expected =
    [parent(sol.u[k].∫F_vol_liq_water_dt)[1] for k in 2:2112:length(sol.t)]
lines!(
    ax2,
    daily[1:2112:end][2:end],
    eps(FT) .+
    abs.(
        (mass_change_actual - mass_change_expected[2:end]) ./
        mass_change_expected[2:end],
    ),
    label = "Water Balance Sparse Save",
)

axislegend(position = :lb)
CairoMakie.save(joinpath(savedir, "water_conservation.png"), fig2)
