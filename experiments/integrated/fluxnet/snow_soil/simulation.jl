## Some site parameters have been taken from
## Ma, S., Baldocchi, D. D., Xu, L., Hehn, T. (2007)
## Inter-Annual Variability In Carbon Dioxide Exchange Of An
## Oak/Grass Savanna And Open Grassland In California, Agricultural
## And Forest Meteorology, 147(3-4), 157-171. https://doi.org/10.1016/j.agrformet.2007.07.008 
## CLM 5.0 Tech Note: https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
# Bonan, G. Climate change and terrestrial ecosystem modeling. Cambridge University Press, 2019.

import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using CairoMakie
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
site_ID = "US-MOz"
time_offset = 7 # offset from UTC in hrs
lat = FT(38.7441) # degree
long = FT(-92.2000) # degree

earth_param_set = LP.LandParameters(FT)
# Utility functions for reading in and filling fluxnet data
include(joinpath(climaland_dir, "experiments/integrated/fluxnet/data_tools.jl"))

# Column dimensions - separation of layers at the top and bottom of the column:
nelements = 10
zmin = FT(-2)
zmax = FT(0)
dz_bottom = FT(1.0)
dz_top = FT(0.04)

domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = (dz_bottom, dz_top),
)

# TIME STEPS:
t0 = FT(1800)
N_days_spinup = 0
N_days = N_days_spinup + 360
dt = FT(900)
tf = t0 + FT(3600 * 24 * N_days)


data_link = "https://caltech.box.com/shared/static/7r0ci9pacsnwyo0o9c25mhhcjhsu6d72.csv"
# Height of sensor on flux tower
atmos_h = FT(32)

# Required to use the same met_drivers_FLUXNET script, even though we currently have no canopy
h_leaf = FT(0)
f_root_to_shoot = FT(0)
plant_ν = FT(0)

# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/met_drivers_FLUXNET.jl",
    ),
)
forcing = (; atmos, radiation)
prognostic_land_components = (:snow, :soil)
α_soil = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = FT(0.25),
    NIR_albedo = FT(0.25),
)
runoff = ClimaLand.Soil.SurfaceRunoff()
retention_parameters = (;
    ν = FT(0.5),
    K_sat = FT(0.45 / 3600 / 100),
    hydrology_cm = vanGenuchten{FT}(; α = FT(2), n = FT(2)),
    θ_r = FT(0.09),
)
composition_parameters =
    (; ν_ss_quartz = FT(0.2), ν_ss_om = FT(0.0), ν_ss_gravel = FT(0.4))
S_s = FT(1e-3)
z_0m = FT(0.01)
z_0b = FT(0.001)
emissivity = FT(0.98)
soil_model = Soil.EnergyHydrology(
    FT,
    domain,
    forcing,
    earth_param_set;
    prognostic_land_components,
    albedo = α_soil,
    runoff,
    retention_parameters,
    S_s,
    composition_parameters,
    z_0m,
    z_0b,
    emissivity,
)
α_snow = Snow.ConstantAlbedoModel(0.8)
density = Snow.MinimumDensityModel(300.0)
snow_model = Snow.SnowModel(
    FT,
    ClimaLand.Domains.obtain_surface_domain(domain),
    forcing,
    earth_param_set,
    dt;
    prognostic_land_components,
    α_snow,
    density,
)

land = ClimaLand.SoilSnowModel(; snow = snow_model, soil = soil_model)
Y, p, cds = initialize(land)
exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land)
jacobian! = ClimaLand.make_jacobian(land);
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);
set_initial_cache! = make_set_initial_cache(land)


#Initial conditions
# Find data index corresponds to t0
t0_idx = 1 + Int(round(t0 / DATA_DT))
Y.soil.ϑ_l =
    drivers.SWC.status != absent ? drivers.SWC.values[t0_idx] : soil_ν / 2
Y.soil.θ_i = 0
T_0 =
    drivers.TS.status != absent ? drivers.TS.values[t0_idx] :
    drivers.TA.values[t0_idx]
ρc_s =
    volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        soil_model.parameters.ρc_ds,
        earth_param_set,
    )
Y.soil.ρe_int =
    volumetric_internal_energy.(Y.soil.θ_i, ρc_s, T_0, earth_param_set)

Y.snow.S .= 0.0
Y.snow.S_l .= 0.0
Y.snow.U .= 0.0
set_initial_cache!(p, Y, t0)

saveat = Array(t0:dt:tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)

updateat = deepcopy(saveat)
model_drivers = ClimaLand.get_drivers(land)
updatefunc = ClimaLand.make_update_drivers(model_drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
# TIME STEPPING
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 3,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
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
sol = SciMLBase.solve(
    prob,
    ode_algo;
    callback = cb,
    dt = dt,
    saveat = collect(t0:dt:tf),
);

# Plotting
daily = sol.t ./ 3600 ./ 24
savedir = joinpath(climaland_dir, "experiments/integrated/fluxnet/snow_soil")


# Water content

fig = Figure(size = (1600, 1200), fontsize = 26)
ax1 = Axis(fig[2, 2], ylabel = "SWC", xlabel = "Days")
lines!(
    ax1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end - 2] for k in 1:1:length(sol.t)],
    label = "10cm",
)
lines!(
    ax1,
    daily,
    [
        parent(sol.u[k].soil.θ_i .+ sol.u[k].soil.ϑ_l)[end - 2] for
        k in 1:1:length(sol.t)
    ],
    label = "10cm, liq+ice",
)

lines!(
    ax1,
    seconds ./ 3600 ./ 24,
    drivers.SWC.values[:],
    label = "Data, Unknown Depth",
)
axislegend(ax1, position = :rt)

ax2 = Axis(fig[1, 2], ylabel = "Precipitation (mm/day)")
ylims!(ax2, [-1300, 0])
hidexdecorations!(ax2, ticks = false)
lines!(ax2, seconds ./ 3600 ./ 24, P_liq .* (1e3 * 24 * 3600), label = "Liquid")
lines!(ax2, seconds ./ 3600 ./ 24, P_snow .* (1e3 * 24 * 3600), label = "Snow")
axislegend(ax2, position = :rb)
ax3 = Axis(fig[2, 1], ylabel = "SWE (m)", xlabel = "Days")
lines!(ax3, daily, [parent(sol.u[k].snow.S)[1] for k in 1:1:length(sol.t)])

# Temp
ax4 = Axis(fig[1, 1], ylabel = "T (K)")
hidexdecorations!(ax4, ticks = false)
lines!(ax4, seconds ./ 3600 ./ 24, drivers.TA.values[:], label = "Data, Air")
lines!(
    ax4,
    sv.t ./ 24 ./ 3600,
    [parent(sv.saveval[k].soil.T)[end - 2] for k in 1:1:length(sv.t)],
    label = "Model 10 cm",
)

lines!(
    ax4,
    sv.t ./ 24 ./ 3600,
    [parent(sv.saveval[k].snow.T)[1] for k in 1:1:length(sv.t)],
    label = "Snow",
)
lines!(
    ax4,
    seconds ./ 3600 ./ 24,
    drivers.TS.values[:],
    label = "Data, Unknown depth",
)
axislegend(ax4, position = :rt)
CairoMakie.save(joinpath(savedir, "results.png"), fig)

# Assess conservation
_ρ_i = FT(LP.ρ_cloud_ice(earth_param_set))
_ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
fig = Figure(size = (1600, 1200), fontsize = 26)
ax1 = Axis(fig[2, 1], ylabel = "ΔEnergy (J/A)", xlabel = "Days")
function compute_atmos_energy_fluxes(p)
    e_flux_falling_snow =
        Snow.energy_flux_falling_snow(atmos, p, land.snow.parameters)
    e_flux_falling_rain =
        Snow.energy_flux_falling_rain(atmos, p, land.snow.parameters)

    return @. (1 - p.snow.snow_cover_fraction) * (
                  p.soil.turbulent_fluxes.lhf +
                  p.soil.turbulent_fluxes.shf +
                  p.soil.R_n +
                  e_flux_falling_rain
              ) +
              p.snow.snow_cover_fraction * (
                  p.snow.turbulent_fluxes.lhf +
                  p.snow.turbulent_fluxes.shf +
                  p.snow.R_n +
                  e_flux_falling_rain
              ) +
              e_flux_falling_snow
end

function compute_atmos_water_vol_fluxes(p)
    return @. p.drivers.P_snow +
              p.drivers.P_liq +
              (1 - p.snow.snow_cover_fraction) * (
                  p.soil.turbulent_fluxes.vapor_flux_liq +
                  p.soil.turbulent_fluxes.vapor_flux_ice
              ) +
              p.snow.snow_cover_fraction * p.snow.turbulent_fluxes.vapor_flux
end

function compute_energy_of_runoff(p)
    liquid_influx = @. p.snow.water_runoff * p.snow.snow_cover_fraction +
       (1 - p.snow.snow_cover_fraction) * p.drivers.P_liq
    e_flux_falling_rain =
        Soil.volumetric_internal_energy_liq.(p.drivers.T, earth_param_set) .*
        p.drivers.P_liq
    influx_energy = @. p.snow.energy_runoff * p.snow.snow_cover_fraction +
       (1 - p.snow.snow_cover_fraction) * e_flux_falling_rain
    runoff_fraction = @. 1 - ClimaLand.Soil.compute_infiltration_fraction(
        p.soil.infiltration,
        liquid_influx,
    )
    return runoff_fraction .* influx_energy
end

function compute_runoff(p)
    liquid_influx = @. p.snow.water_runoff * p.snow.snow_cover_fraction +
       (1 - p.snow.snow_cover_fraction) * p.drivers.P_liq
    runoff_fraction = @. 1 - ClimaLand.Soil.compute_infiltration_fraction(
        p.soil.infiltration,
        liquid_influx,
    )
    return runoff_fraction .* liquid_influx
end

ΔE_expected =
    cumsum(
        -1 .* [
            parent(
                compute_atmos_energy_fluxes(sv.saveval[k]) .-
                compute_energy_of_runoff(sv.saveval[k]) .-
                sv.saveval[k].soil.bottom_bc.heat,
            )[end] for k in 1:1:(length(sv.t) - 1)
        ],
    ) * (sv.t[2] - sv.t[1])
E_measured = [
    sum(sol.u[k].soil.ρe_int) + parent(sol.u[k].snow.U)[end] for
    k in 1:1:length(sv.t)
]
ΔW_expected =
    cumsum(
        -1 .* [
            parent(
                compute_atmos_water_vol_fluxes(sv.saveval[k]) .-
                compute_runoff(sv.saveval[k]) .-
                sv.saveval[k].soil.bottom_bc.water,
            )[end] for k in 1:1:(length(sv.t) - 1)
        ],
    ) * (sv.t[2] - sv.t[1])
W_measured = [
    sum(sol.u[k].soil.ϑ_l) +
    sum(sol.u[k].soil.θ_i) * _ρ_i / _ρ_l +
    parent(sol.u[k].snow.S)[end] for k in 1:1:length(sv.t)
]
lines!(
    ax1,
    daily[2:end],
    E_measured[2:end] .- E_measured[1],
    label = "Simulated",
)
lines!(ax1, daily[2:end], ΔE_expected, label = "Expected")
axislegend(ax1, position = :rt)

# Temp
ax4 = Axis(fig[1, 1], ylabel = "ΔWater (m)")
hidexdecorations!(ax4, ticks = false)
lines!(
    ax4,
    daily[2:end],
    W_measured[2:end] .- W_measured[1],
    label = "Simulated",
)

lines!(ax4, daily[2:end], ΔW_expected, label = "Expected")
axislegend(ax4, position = :rt)

ax3 = Axis(fig[2, 2], ylabel = "ΔE/E", xlabel = "Days", yscale = log10)
lines!(
    ax3,
    daily[2:end],
    abs.(E_measured[2:end] .- E_measured[1] .- ΔE_expected) ./ mean(E_measured),
)

ax2 = Axis(fig[1, 2], ylabel = "ΔW/W", yscale = log10)
hidexdecorations!(ax2, ticks = false)
lines!(
    ax2,
    daily[2:end],
    abs.(W_measured[2:end] .- W_measured[1] .- ΔW_expected) ./ mean(W_measured),
)

CairoMakie.save(joinpath(savedir, "results_conservation.png"), fig)
