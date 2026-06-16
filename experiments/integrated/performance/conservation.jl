import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaComms
ClimaComms.@import_required_backends
using CairoMakie
using Statistics
using Dates
using Printf
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
stop_date = start_date + Year(3)
Δt_seconds = 450.0
domain = ClimaLand.Domains.Column(;
    dz_tuple = FT.((3, 0.05)),
    nelements = 15,
    longlat = FT.((-92.0, 38.7441)),
    zlim = FT.((-15, 0)),
);

surface_domain = ClimaLand.Domains.obtain_surface_domain(domain);
surface_space = domain.space.surface;
# Forcing data - high resolution
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict,
    FT;
    max_wind_speed = 25.0,
    use_lowres_forcing = true,
);
forcing = (; atmos, radiation);

LAI = TimeVaryingInput((t) -> 1.0)
land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    domain,
    Δt_seconds;
    conservation = true,
);
saveat = Second(Int(Δt_seconds))
saving_cb = ClimaLand.NonInterpSavingCallback(start_date, stop_date, saveat);
sv = saving_cb.affect!.saved_values;
savedir =
    generate_output_path("experiments/integrated/performance/conservation")
diagnostics = ClimaLand.default_diagnostics(
    land,
    start_date,
    savedir;
    conservation = true,
    conservation_period = Month(1),
)
simulation = ClimaLand.Simulations.LandSimulation(
    start_date,
    stop_date,
    Δt_seconds,
    land;
    user_callbacks = (saving_cb,),
    updateat = Second(3600 * 3),
    solver_kwargs = (; saveat),
    diagnostics,
);
@time sol = ClimaLand.Simulations.solve!(simulation);

# Make plots using every timestep
yearly = float.(sol.t[2:end]) ./ 3600 ./ 24 ./ 364
# Energy balance
total_energy_array =
    [parent(sv.saveval[k].total_energy)[1] for k in 2:length(sol.t)];
mass_change_actual = total_energy_array[2:end] .- total_energy_array[1];
mass_change_expected =
    [parent(sol.u[k].∫F_vol_e_dt)[1] for k in 2:(length(sol.t) - 1)];
mean_energy = mean(total_energy_array)
fig = Figure(size = (1500, 400))
ax = Axis(
    fig[1, 1],
    xlabel = "Years",
    ylabel = "Fractional Error",
    yscale = log10,
)

lines!(
    ax,
    yearly[2:end],
    eps(FT) .+
    abs.((mass_change_actual - mass_change_expected) ./ mean_energy,),
    label = "Energy Balance",
)

CairoMakie.save(joinpath(savedir, "energy_conservation.png"), fig)

# Water balance
total_water_array =
    [parent(sv.saveval[k].total_water)[1] for k in 2:length(sol.t)]
mass_change_actual = total_water_array[2:end] .- total_water_array[1];
mass_change_expected =
    [parent(sol.u[k].∫F_vol_liq_water_dt)[1] for k in 2:(length(sol.t) - 1)];
mean_water = mean(total_water_array)
fig2 = Figure(size = (1500, 400))
ax2 = Axis(
    fig2[1, 1],
    xlabel = "Years",
    ylabel = "Fractional Error",
    yscale = log10,
)

lines!(
    ax2,
    yearly[2:end],
    eps(FT) .+ abs.((mass_change_actual - mass_change_expected) ./ mean_water,),
    label = "Water Balance",
)
CairoMakie.save(joinpath(savedir, "water_conservation.png"), fig2)

# Make plots using the diagnostics
times = diagnostics[1].output_writer.dict["wvpac_1M_inst"].keys
water_volume_per_area = diagnostics[1].output_writer.dict["wvpa_1M_inst"]
water_volume_per_area_change =
    diagnostics[1].output_writer.dict["wvpac_1M_inst"]
water_volume_0 = parent(water_volume_per_area[times[1]])[1]
water_volume_end = parent(water_volume_per_area[times[end]])[1]
mean_water_volume = (water_volume_0 + water_volume_end) / 2

energy_per_area = diagnostics[1].output_writer.dict["epa_1M_inst"]
energy_per_area_change = diagnostics[1].output_writer.dict["epac_1M_inst"]
energy_0 = parent(energy_per_area[times[1]])[1]
energy_end = parent(energy_per_area[times[end]])[1]
mean_energy = (energy_0 + energy_end) / 2

N = length(times)
energy_error = zeros(N - 1)
water_volume_error = zeros(N - 1)
for (i, t) in enumerate(times[1:(N - 1)])
    # error = nanmean[(X(t) - X(0) - Expected Change in X)]
    energy_error[i] =
        abs(
            eps(FT) + parent(energy_per_area[t])[1] - energy_0 -
            parent(energy_per_area_change[times[i + 1]])[1],
        ) ./ mean_energy
    water_volume_error[i] =
        abs(
            eps(FT) + parent(water_volume_per_area[t])[1] - water_volume_0 -
            parent(water_volume_per_area_change[times[i + 1]])[1],
        ) ./ mean_water_volume
end
titles = [" mean energy per area", " mean water volume per area"]
quantity_names = ["energy", "water"]
errors = [energy_error, water_volume_error]
typical_value =
    [@sprintf("%1.2le", mean_energy), @sprintf("%1.2le", mean_water_volume)]
units = ["J/m²", "m³/m²"]
for i in 1:2
    fig_cycle = CairoMakie.Figure(size = (600, 400))
    ax_d = Axis(
        fig_cycle[1, 1],
        xlabel = "Years",
        ylabel = " Mean Fractional Conservation Error",
        title = "$(titles[i]), typical value = $(typical_value[i]) $(units[i])",
        yscale = log10,
    )
    CairoMakie.lines!(
        ax_d,
        float.(times[1:(N - 1)]) ./ 24 ./ 3600 ./ 365,
        errors[i],
    )
    CairoMakie.save(
        joinpath(savedir, "$(quantity_names[i])_diagnostics.png"),
        fig_cycle,
    )
end
