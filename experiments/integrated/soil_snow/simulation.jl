import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
using ClimaCore
using Plots
using Statistics
using Dates
using NCDatasets

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Snow
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams

climaland_dir = pkgdir(ClimaLand)

FT = Float64
earth_param_set = LP.LandParameters(FT)
# Utility functions for reading in and filling fluxnet data
include(joinpath(climaland_dir, "experiments/integrated/fluxnet/data_tools.jl"))

# Column dimensions - separation of layers at the top and bottom of the column:
nelements = 10
zmin = FT(-2)
zmax = FT(0)
dz_bottom = FT(1.0)
dz_top = FT(0.04)

land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = (dz_bottom, dz_top),
)

# TIME STEPS:
t0 = FT(1800)
N_days_spinup = 0
N_days = N_days_spinup + 100
dt = FT(60)
tf = t0 + FT(3600 * 24 * N_days)

# Read in the parameters for the site
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/soil_snow/parameters_fluxnet.jl",
    ),
)
# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(climaland_dir, "experiments/integrated/soil_snow/met_drivers_fluxnet.jl"),
)


# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain
soil_ps = Soil.EnergyHydrologyParameters(
    FT;
    ν = soil_ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    K_sat = soil_K_sat,
    S_s = soil_S_s,
    θ_r,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
)

soil_args = (domain = soil_domain, parameters = soil_ps)
soil_model_type = Soil.EnergyHydrology{FT}

snow_domain = ClimaLand.Domains.obtain_surface_domain(land_domain);
ρ = 300.0
α = 0.8
snow_parameters = SnowParameters{FT}(
    dt;
    α_snow = α,
    ρ_snow = ρ,
    earth_param_set = earth_param_set,
);
snow_args = (; domain = snow_domain, parameters = snow_parameters);
snow_model_type = Snow.SnowModel
land_input = (atmos = atmos, radiation = radiation, runoff = ClimaLand.Soil.SiteLevelSurfaceRunoff())
land = ClimaLand.LandHydrologyModel{FT}(;
    land_args = land_input,
    soil_model_type = soil_model_type,
    soil_args = soil_args,
    snow_args = snow_args,
    snow_model_type = snow_model_type,
)
Y, p, cds = initialize(land)
exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land)
tendency_jacobian! = ClimaLand.make_tendency_jacobian(land);
jac_kwargs =
    (; jac_prototype = ImplicitEquationJacobian(Y), Wfact = tendency_jacobian!);
set_initial_cache! = make_set_initial_cache(land)


#Initial conditions
Y.soil.ϑ_l =
    drivers.SWC.status != absent ?
    drivers.SWC.values[1 + Int(round(t0 / DATA_DT))] : soil_ν / 2
Y.soil.θ_i = 0
T_0 =
    drivers.TS.status != absent ?
    drivers.TS.values[1 + Int(round(t0 / DATA_DT))] :
    drivers.TA.values[1 + Int(round(t0 / DATA_DT))]
ρc_s = volumetric_heat_capacity.(Y.soil.ϑ_l, Y.soil.θ_i, soil_ps.ρc_ds, earth_param_set)
Y.soil.ρe_int = volumetric_internal_energy.(Y.soil.θ_i, ρc_s, T_0, earth_param_set)

Y.snow.S .= 0.0
Y.snow.U .= 0.0
set_initial_cache!(p, Y, t0)

saveat = Array(t0:(3 * 3600):tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)

updateat = deepcopy(saveat)
updatefunc = ClimaLand.make_update_drivers(atmos, radiation)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
# TIME STEPPING
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 3,
        update_j = CTS.UpdateEvery(CTS.NewTimeStep),
    ),
);

# Problem definition and callbacks
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
sol = SciMLBase.solve(prob, ode_algo;    callback = cb, dt = dt, saveat = t0:3*3600:tf);

# Plotting
daily = sol.t ./ 3600 ./ 24
savedir = joinpath(climaland_dir, "experiments/integrated/soil_snow")


# Water content

plt1 = Plots.plot(size = (1500, 800))
Plots.plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end - 1] for k in 1:1:length(sol.t)],
    label = "6cm",
    xlim = [minimum(daily), maximum(daily)],
    xlabel = "Days",
    ylabel = "SWC [m/m]",
    color = "blue",
    margin = 10Plots.mm,
)
Plots.plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.θ_i)[end - 1] for k in 1:1:length(sol.t)],
    label = "6cm, ice",
)

Plots.plot!(plt1, seconds ./ 3600 ./ 24, drivers.SWC.values[:], label = "Data, ?cm")

plt2 = Plots.plot(
    seconds ./ 3600 ./ 24,
    P_liq .* (1e3 * 24 * 3600),
    label = "Liquid",
    ylabel = "Precipitation [mm/day]",
    xlim = [minimum(daily), maximum(daily)],
    margin = 10Plots.mm,
    ylim = [-1300, 0],
    size = (1500, 400),
)
Plots.plot!(
    plt2,
    seconds ./ 3600 ./ 24,
    P_snow .* (1e3 * 24 * 3600),
    label = "Snow",
)

plt3 = Plots.plot(size = (1500, 800))
Plots.plot!(
    plt3,
    daily,
    [parent(sol.u[k].snow.S)[1] for k in 1:1:length(sol.t)],
    label = "Model, SWE",
    margin = 10Plots.mm,
)
Plots.plot(plt3, plt2, plt1, layout = grid(3, 1))


Plots.savefig(joinpath(savedir, "soil_water_content.png"))

# Temp
plt3 = Plots.plot(size = (1500, 800))
Plots.plot!(
    plt3,
    sv.t ./ 24 ./ 3600,
    [parent(sv.saveval[k].soil.T)[end - 1] for k in 1:1:length(sv.t)],
    label = "Model, 6cm",
    xlim = [minimum(daily), maximum(daily)],
    xlabel = "Days",
    ylabel = "T [m/m]",
    color = "blue",
    margin = 10Plots.mm,
)

Plots.plot!(
    plt3,
    sv.t ./ 24 ./ 3600,
    [parent(sv.saveval[k].snow.T)[1] for k in 1:1:length(sv.t)],
    label = "Tsnow",
)
Plots.plot!(plt3, seconds ./ 3600 ./ 24, drivers.TS.values[:], label = "Data, ?cm")
Plots.plot!(plt3, seconds ./ 3600 ./ 24, drivers.TA.values[:], label = "Air")


# LHF
lhf = [
    parent(sv.saveval[k].soil.turbulent_fluxes.lhf)[1] for k in 1:length(sol.t)
]
Plots.plot(daily, lhf, label = "Model", xlim = extrema(daily))
Plots.plot!(seconds ./ 3600 ./ 24, LE, label = "Data")

# SHF
shf =
    [parent(sv.saveval[k].soil.sfc_conditions.shf)[1] for k in 1:length(sol.t)]
Plots.plot(daily, shf, label = "Model", xlim = extrema(daily))
Plots.plot!(seconds ./ 3600 ./ 24, H, label = "Data")

soil_Rn = [parent(sv.saveval[k].soil.R_n)[1] for k in 1:length(sol.t)]
soil_Rn_data = @. (SW_IN - SW_OUT + LW_IN - LW_OUT)
Plots.plot(daily, -soil_Rn, label = "Model", xlim = extrema(daily))
Plots.plot!(seconds ./ 3600 ./ 24, soil_Rn_data, label = "Data")
