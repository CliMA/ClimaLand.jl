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
using Dates

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Snow
import ClimaLand
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
import ClimaParams

using DelimitedFiles
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using CairoMakie, StatsBase

const FT = Float64
earth_param_set = LP.LandParameters(FT)
climaland_dir = pkgdir(ClimaLand)

site_ID = "US-MOz"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)

# Get the default values for this site's domain, location, and parameters
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))
(;
    soil_ν,
    soil_K_sat,
    soil_S_s,
    soil_vg_n,
    soil_vg_α,
    θ_r,
    ν_ss_quartz,
    ν_ss_om,
    ν_ss_gravel,
    z_0m_soil,
    z_0b_soil,
    soil_ϵ,
    soil_α_PAR,
    soil_α_NIR,
    Ω,
    χl,
    G_Function,
    α_PAR_leaf,
    λ_γ_PAR,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
    ϵ_canopy,
    ac_canopy,
    g1,
    Drel,
    g0,
    Vcmax25,
    SAI,
    f_root_to_shoot,
    K_sat_plant,
    ψ63,
    Weibull_param,
    a,
    conductivity_model,
    retention_model,
    plant_ν,
    plant_S_s,
    rooting_depth,
    n_stem,
    n_leaf,
    h_leaf,
    h_stem,
    h_canopy,
    z0_m,
    z0_b,
) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

compartment_midpoints =
    n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
compartment_surfaces = n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]

# Construct the ClimaLand domain to run the simulation on
domain = Column(; zlim = (zmin, zmax), nelements = nelements, dz_tuple)

# Set up the timestepping information for the simulation
start_date = DateTime(2010) + Hour(time_offset)
N_days = 360
end_date = start_date + Day(N_days)
dt = FT(900)

# Height of sensor on flux tower
atmos_h = FT(32)
forcing = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    earth_param_set,
    FT,
)

# Construct the soil model
prognostic_land_components = (:snow, :soil)
α_soil = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
)
runoff = ClimaLand.Soil.SurfaceRunoff()
retention_parameters = (;
    ν = soil_ν,
    K_sat = soil_K_sat,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    θ_r = θ_r,
)
composition_parameters = (; ν_ss_quartz, ν_ss_om, ν_ss_gravel)
soil_model = Soil.EnergyHydrology{FT}(
    domain,
    forcing,
    earth_param_set;
    prognostic_land_components,
    albedo = α_soil,
    runoff,
    retention_parameters,
    S_s = soil_S_s,
    composition_parameters,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
)

# Construct the snow model
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

# Construct the land model
land = ClimaLand.SoilSnowModel{FT}(; snow = snow_model, soil = soil_model)

# Initial conditions
set_ic!(Y, p, _, model) = FluxnetSimulations.set_fluxnet_ic!(
    Y,
    site_ID,
    start_date,
    time_offset,
    model,
)

saveat = Array(start_date:Second(dt):end_date)
sv = (;
    t = Array{DateTime}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)

updateat = deepcopy(saveat)

simulation = LandSimulation(
    start_date,
    end_date,
    dt,
    land;
    user_callbacks = (saving_cb,),
    set_ic! = set_ic!,
    updateat,
    solver_kwargs = (; saveat = saveat),
)
sol = solve!(simulation)

# Plotting
daily = FT.(sol.t) ./ 3600 ./ 24
savedir =
    joinpath(climaland_dir, "experiments/integrated/fluxnet/ozark_soilsnow")
mkpath(savedir)
comparison_data = FluxnetSimulations.get_comparison_data(site_ID, time_offset)

# Water content
seconds =
    [Second(dt).value for dt in (comparison_data.UTC_datetime .- start_date)]
fig = Figure(size = (1600, 1200), fontsize = 26)
ax1 = Axis(fig[2, 2], ylabel = "SWC", xlabel = "Days")
lines!(
    ax1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end - 2] for k in 1:1:length(FT.(sol.t))],
    label = "10cm",
)
lines!(
    ax1,
    daily,
    [
        parent(sol.u[k].soil.θ_i .+ sol.u[k].soil.ϑ_l)[end - 2] for
        k in 1:1:length(FT.(sol.t))
    ],
    label = "10cm, liq+ice",
)

lines!(
    ax1,
    seconds ./ 3600 ./ 24,
    comparison_data.swc,
    label = "Data, Unknown Depth",
)
axislegend(ax1, position = :rt)

ax2 = Axis(fig[1, 2], ylabel = "Precipitation (mm/day)")
ylims!(ax2, [-1300, 0])
hidexdecorations!(ax2, ticks = false)
lines!(
    ax2,
    seconds ./ 3600 ./ 24,
    comparison_data.precip .* (1e3 * 24 * 3600),
    label = "Total precip",
)
axislegend(ax2, position = :rb)
ax3 = Axis(fig[2, 1], ylabel = "SWE (m)", xlabel = "Days")
lines!(ax3, daily, [parent(sol.u[k].snow.S)[1] for k in 1:1:length(FT.(sol.t))])

# Temp
sv_times = Dates.value.(Second.(sv.t .- sv.t[1]))
ax4 = Axis(fig[1, 1], ylabel = "T (K)")
hidexdecorations!(ax4, ticks = false)
lines!(
    ax4,
    sv_times ./ 24 ./ 3600,
    [parent(sv.saveval[k].soil.T)[end - 2] for k in 1:1:length(sv_times)],
    label = "Model 10 cm",
)

lines!(
    ax4,
    sv_times ./ 24 ./ 3600,
    [parent(sv.saveval[k].snow.T)[1] for k in 1:1:length(sv_times)],
    label = "Snow",
)
lines!(
    ax4,
    seconds ./ 3600 ./ 24,
    comparison_data.tsoil,
    label = "Data, Unknown depth",
)
axislegend(ax4, position = :rt)
CairoMakie.save(joinpath(savedir, "results.png"), fig)

# Assess conservation
atmos = forcing.atmos
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
            )[end] for k in 1:1:(length(sv_times) - 1)
        ],
    ) * (sv_times[2] - sv_times[1])
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
            )[end] for k in 1:1:(length(sv_times) - 1)
        ],
    ) * (sv_times[2] - sv_times[1])
W_measured = [
    sum(sol.u[k].soil.ϑ_l) +
    sum(sol.u[k].soil.θ_i) * _ρ_i / _ρ_l +
    parent(sol.u[k].snow.S)[end] for k in 1:1:length(sv_times)
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
