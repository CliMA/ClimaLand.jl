import SciMLBase
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
using ClimaLand.Simulations
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams

climaland_dir = pkgdir(ClimaLand)
FT = Float64
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
start_date = DateTime("2008-03-01")
stop_date = DateTime("2010-03-01")
Δt = 450.0
domain = ClimaLand.Domains.Column(;
    dz_tuple = FT.((3, 0.05)),
    nelements = 15,
    longlat = FT.((-62.0, 3.0)),
    zlim = FT.((-15, 0)),
)

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
outdir =
    generate_output_path("experiments/integrated/performance/conservation_simulation")
diagnostics = ClimaLand.default_diagnostics(
    land,
    start_date,
    outdir;
    conservation = true,
    conservation_period = Day(10),
)
simulation = LandSimulation(start_date, stop_date, Δt, land; outdir, diagnostics);
solve!(simulation)

times = diagnostics[1].output_writer.dict["wvpac_10d_inst"].keys
water_volume_per_area = diagnostics[1].output_writer.dict["wvpa_10d_inst"]
water_volume_per_area_change = diagnostics[1].output_writer.dict["wvpac_10d_inst"]
water_volume_0 = parent(water_volume_per_area[times[1]])[1]
water_volume_end = parent(water_volume_per_area[times[end]])[1]
mean_water_volume =(water_volume_0 +water_volume_end) / 2

energy_per_area = diagnostics[1].output_writer.dict["epa_10d_inst"]
energy_per_area_change = diagnostics[1].output_writer.dict["epac_10d_inst"]
energy_0 = parent(energy_per_area[times[1]])[1]
energy_end = parent(energy_per_area[times[end]])[1]
mean_energy =(energy_0 +energy_end) / 2

N = length(times)
energy_error = zeros(N)
water_volume_error = zeros(N)
for (i, t) in enumerate(times[1:end-1])
    # error = nanmean[(X(t) - X(0) - Expected Change in X)]
    energy_error[i] = abs(parent(energy_per_area[t])[1] - energy_0 - parent(energy_per_area_change[times[i+1]])[1])
    water_volume_error[i] = abs(parent(water_volume_per_area[t])[1] - water_volume_0 - parent(water_volume_per_area_change[times[i+1]])[1])
end
titles =
    [" mean energy per area", " mean water volume per area"]
quantity_names = ["energy", "water"]
errors = [energy_error, water_volume_error]
typical_value =
    [@sprintf("%1.2le", mean_energy), @sprintf("%1.2le", mean_water_volume)]
units = ["J/m²", "m³/m²"]
for i in 1:2
    fig_cycle = CairoMakie.Figure(size = (600, 400))
    ax = Axis(
        fig_cycle[1, 1],
        xlabel = "Years",
        ylabel = " Mean Conservation Error [$(units[i])]",
        title = "$(titles[i]), typical value = $(typical_value[i]) $(units[i])",
    )
    CairoMakie.lines!(ax, float.(times) ./ 24 ./ 3600 ./ 365, errors[i])
    CairoMakie.save(
        joinpath(outdir, "$(quantity_names[i]).png"),
        fig_cycle,
    )
end
