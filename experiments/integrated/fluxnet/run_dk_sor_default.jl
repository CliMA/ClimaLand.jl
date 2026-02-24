"""
Run DK-Sor site with default parameters for 2004, plot daily NEE, LH, SH.
Reuses setup from calibrate_dk_sor.jl
"""

# ── 1. Environment & Imports ──────────────────────────────────────────────────
import ClimaLand
import SciMLBase
using ClimaCore
import ClimaComms
import ClimaParams as CP
using Dates
using Insolation

using ClimaLand
using ClimaLand: LandModel
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Snow
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.TimeManager: date

using CairoMakie
CairoMakie.activate!()
using Statistics
using Logging
using NCDatasets

# ── 2. Configuration ─────────────────────────────────────────────────────────
const FT = Float64
site_ID = "DK-Sor"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)
climaland_dir = pkgdir(ClimaLand)

met_nc_path = joinpath(climaland_dir, "DK_Sor", "DK-Sor_1997-2014_FLUXNET2015_Met.nc")
flux_nc_path = joinpath(climaland_dir, "DK_Sor", "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")

spinup_days = 60
dt = Float64(450)

sim_start_year = 2004
sim_end_year = 2006

(; time_offset, lat, long) = FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

start_date = DateTime(sim_start_year, 1, 1) - Day(spinup_days)
stop_date = DateTime(sim_end_year + 1, 1, 1)
spinup_date = start_date + Day(spinup_days)
println("Simulation: $start_date to $stop_date")
println("Spinup until: $spinup_date")

# ── 3. Domain, forcing, LAI ──────────────────────────────────────────────────
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))

land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
surface_space = land_domain.space.surface
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

toml_dict = LP.create_toml_dict(FT)
#toml_dict.data["soil_momentum_roughness_length"]["value"] = FT(1.5)
#toml_dict.data["soil_scalar_roughness_length"]["value"] = FT(1.5)
#toml_dict.data["snow_momentum_roughness_length"]["value"] = FT(1.5)
#toml_dict.data["snow_scalar_roughness_length"]["value"] = FT(1.5)

(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_netcdf(
    met_nc_path, lat, long, time_offset, atmos_h, start_date, toml_dict, FT,
)

met_ds = NCDataset(met_nc_path, "r")
lai_data = Float64.(coalesce.(met_ds["LAI"][1, 1, :], NaN))
lai_times = met_ds["time"][:]
close(met_ds)

lai_seconds = [Float64(Second(t - Dates.Hour(time_offset) - start_date).value) for t in lai_times]
valid_lai = .!isnan.(lai_data)
LAI = TimeVaryingInput(lai_seconds[valid_lai], lai_data[valid_lai])

# ── 4. Load observations ─────────────────────────────────────────────────────
flux_ds = NCDataset(flux_nc_path, "r")
flux_times_dt = flux_ds["time"][:]
nee_raw = Float64.(coalesce.(flux_ds["NEE_daily"][:], NaN))
qle_raw = Float64.(coalesce.(flux_ds["Qle_daily"][:], NaN))
qh_raw = Float64.(coalesce.(flux_ds["Qh_daily"][:], NaN))
close(flux_ds)

flux_dates = Date.(flux_times_dt)
spinup_date_d = Date(spinup_date)
stop_date_d = Date(stop_date)
valid_mask = (flux_dates .>= spinup_date_d) .& (flux_dates .<= stop_date_d) .&
             .!isnan.(nee_raw) .& .!isnan.(qle_raw) .& .!isnan.(qh_raw) .&
             (abs.(nee_raw) .< 1e10) .& (abs.(qle_raw) .< 1e10) .& (abs.(qh_raw) .< 1e10)

obs_dates = flux_dates[valid_mask]
nee_obs = nee_raw[valid_mask]
qle_obs = qle_raw[valid_mask]
qh_obs = qh_raw[valid_mask]
n_obs = length(obs_dates)
println("Obs days: $n_obs")

# ── 5. Build & run model with default parameters ─────────────────────────────
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
forcing_nt = (; atmos, radiation, ground = ClimaLand.PrognosticGroundConditions{FT}())
 biomass = ClimaLand.Canopy.PrescribedBiomassModel{FT}(
        land_domain,
        LAI,
        toml_dict;
        height = FT(25))
canopy = Canopy.CanopyModel{FT}(
    canopy_domain, forcing_nt, LAI, toml_dict;
    prognostic_land_components,
    photosynthesis = Canopy.PModel{FT}(canopy_domain, toml_dict),
    conductance = Canopy.PModelConductance{FT}(toml_dict),
    soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(land_domain, toml_dict),
    biomass
);

land = LandModel{FT}(
    (; atmos, radiation), LAI, toml_dict, land_domain, dt;
    prognostic_land_components, canopy,
);

function custom_set_ic!(Y, p, t, model)
    earth_param_set = ClimaLand.get_earth_param_set(model.soil)
    evaluate!(p.drivers.T, atmos.T, t)
    (; θ_r, ν, ρc_ds) = model.soil.parameters
    @. Y.soil.ϑ_l = θ_r + (ν - θ_r) .* FT(0.95)
    Y.soil.θ_i .= FT(0.0)
    ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l, Y.soil.θ_i, ρc_ds, earth_param_set,
    )
    Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
        Y.soil.θ_i, ρc_s, p.drivers.T, earth_param_set,
    )
    Y.snow.S .= FT(0)
    Y.snow.S_l .= FT(0)
    Y.snow.U .= FT(0)
    if model.canopy.energy isa ClimaLand.Canopy.BigLeafEnergyModel
        Y.canopy.energy.T .= p.drivers.T
    end
    n_stem = model.canopy.hydraulics.n_stem
    n_leaf = model.canopy.hydraulics.n_leaf
    for i in 1:(n_stem + n_leaf)
        Y.canopy.hydraulics.ϑ_l.:($i) .= model.canopy.hydraulics.parameters.ν
    end
    if !isnothing(model.soilco2)
        Y.soilco2.CO2 .= FT(0.000412)
        Y.soilco2.O2_f .= FT(0.21)
        # Prescribed SOC profile: 15 kgC/m³ at surface, 0.5 kgC/m³ at 1m depth
        # Exponential decay: SOC(z) = SOC_bot + (SOC_top - SOC_bot) * exp(z/τ)
        # where z is negative (depth below surface), τ chosen so SOC(-1) = 0.5
        SOC_top = FT(15.0)
        SOC_bot = FT(0.5)
        τ_soc = FT(1.0 / log(SOC_top / SOC_bot))  # ~0.294 m
        z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
        @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
    end
end

output_writer = ClimaDiagnostics.Writers.DictWriter();
diags = ClimaLand.default_diagnostics(
    land, start_date;
    output_writer, output_vars = :long, reduction_period = :daily,
);

simulation = LandSimulation(
    start_date, stop_date, dt, land;
    set_ic! = custom_set_ic!, updateat = Second(dt), diagnostics = diags,
);

println("Running simulation...")
io = open("logfile.txt", "w")
logger = ConsoleLogger(io)
with_logger(logger) do
    solve!(simulation);
end
close(io)
println("Simulation complete.")

# ── 6. Extract diagnostics ────────────────────────────────────────────────────
function extract_daily_diag(simulation, diag_name, target_dates)
    writer = nothing
    for d in simulation.diagnostics
        if haskey(d.output_writer.dict, diag_name)
            writer = d.output_writer
            break
        end
    end
    if isnothing(writer)
        available = String[]
        for d in simulation.diagnostics
            append!(available, collect(keys(d.output_writer.dict)))
        end
        error("Diagnostic '$diag_name' not found. Available: $(unique(available))")
    end

    (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, diag_name)
    model_dates_dt = times isa Vector{DateTime} ? times : date.(times)
    model_dates = Date.(model_dates_dt)

    println("  $diag_name: $(length(model_dates)) model dates, $(first(model_dates)) to $(last(model_dates))")

    model_dict = Dict{Date, Float64}()
    for (d, v) in zip(model_dates, data)
        model_dict[d] = Float64(v)
    end

    result = Float64[]
    for d in target_dates
        push!(result, get(model_dict, d, NaN))
    end
    n_match = count(!isnan, result)
    println("  Matched $n_match / $(length(target_dates)) obs dates")
    return result
end

gpp_model = extract_daily_diag(simulation, "gpp_1d_average", obs_dates)
er_model = extract_daily_diag(simulation, "er_1d_average", obs_dates)
lhf_model = extract_daily_diag(simulation, "lhf_1d_average", obs_dates)
shf_model = extract_daily_diag(simulation, "shf_1d_average", obs_dates)
Tair = extract_daily_diag(simulation, "tair_1d_average", obs_dates)
Tsoil = extract_daily_diag(simulation, "tsoil_1d_average", obs_dates)
Tcanopy = extract_daily_diag(simulation, "ct_1d_average", obs_dates)
scf = extract_daily_diag(simulation, "snowc_1d_average", obs_dates)

# NEE = ER - GPP (both in mol CO₂/m²/s), convert to gC/m²/d
nee_model_gC = (er_model .- gpp_model) .* 12.0 .* 86400.0

# ── 7. Extract depth-resolved soilco2 diagnostics ────────────────────────────
# Helper to get a daily timeseries for a depth-resolved diagnostic at a given layer
function extract_daily_diag_layer(simulation, diag_name, layer)
    writer = nothing
    for d in simulation.diagnostics
        if haskey(d.output_writer.dict, diag_name)
            writer = d.output_writer
            break
        end
    end
    isnothing(writer) && error("Diagnostic '$diag_name' not found")
    (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, diag_name; layer)
    model_dates_dt = times isa Vector{DateTime} ? times : date.(times)
    return Date.(model_dates_dt), Float64.(data)
end

# Get hr (heterotrophic respiration = top_bc flux, surface only)
hr_dates, hr_vals = extract_daily_diag_layer(simulation, "hr_1d_average", nothing)
# Convert mol CO₂/m²/s → gC/m²/d
hr_gC = hr_vals .* 12.0 .* 86400.0
println("  hr: $(length(hr_dates)) dates, range $(round(minimum(hr_gC), sigdigits=4)) to $(round(maximum(hr_gC), sigdigits=4)) gC/m²/d")

# Get depth-resolved fields: get number of layers from first snapshot
sco2_writer = let w = nothing
    for d in simulation.diagnostics
        if haskey(d.output_writer.dict, "sco2_1d_average")
            w = d.output_writer
            break
        end
    end
    w
end
first_field = first(values(sco2_writer["sco2_1d_average"]))
n_layers = size(parent(first_field), 1)
println("  Soil has $n_layers layers")

# Extract CO2 concentration, microbe source, and SOC at select layers
# layer 1 = bottom, n_layers = top (surface)
layer_ids = [n_layers, round(Int, n_layers * 3 ÷ 4), n_layers ÷ 2, 1]
layer_labels = String[]

sco2_by_layer = Dict{Int, Vector{Float64}}()
scms_by_layer = Dict{Int, Vector{Float64}}()
soc_by_layer = Dict{Int, Vector{Float64}}()

# Get the depth coordinates from the domain
z_centers = parent(land_domain.fields.z)
println("  Layer depths (center): ", [round(z_centers[i], digits=3) for i in layer_ids])

for lid in layer_ids
    d, v = extract_daily_diag_layer(simulation, "sco2_1d_average", lid)
    sco2_by_layer[lid] = v
    d, v = extract_daily_diag_layer(simulation, "scms_1d_average", lid)
    scms_by_layer[lid] = v
    d, v = extract_daily_diag_layer(simulation, "soc_1d_average", lid)
    soc_by_layer[lid] = v
    push!(layer_labels, "z=$(round(z_centers[lid], digits=2))m")
    println("  Layer $lid (z=$(round(z_centers[lid], digits=2))m): " *
            "CO2=$(round(minimum(sco2_by_layer[lid]), sigdigits=4))..$(round(maximum(sco2_by_layer[lid]), sigdigits=4)), " *
            "Sm=$(round(minimum(scms_by_layer[lid]), sigdigits=4))..$(round(maximum(scms_by_layer[lid]), sigdigits=4)), " *
            "SOC=$(round(minimum(soc_by_layer[lid]), sigdigits=4))..$(round(maximum(soc_by_layer[lid]), sigdigits=4))")
end

# Use hr_dates for the x-axis of soilco2 plots
spinup_date_d2 = Date(spinup_date)
post_spinup = hr_dates .>= spinup_date_d2
hr_dates_ps = hr_dates[post_spinup]
hr_gC_ps = hr_gC[post_spinup]
t_days_co2 = [Dates.value(d - spinup_date_d2) for d in hr_dates_ps]

# ── 8. Plot ───────────────────────────────────────────────────────────────────
savedir = joinpath(@__DIR__, "calibrate_dk_sor_output")
mkpath(savedir)

# Figure 1: fluxes vs obs
fig = Figure(size = (1200, 900))

ax1 = Axis(fig[1, 1]; ylabel = "NEE (gC/m²/d)", title = "DK-Sor 2004 — Default Parameters")
lines!(ax1, obs_dates, nee_obs; color = :black, linewidth = 1.5, label = "Obs")
lines!(ax1, obs_dates, nee_model_gC; color = :blue, linewidth = 1.5, label = "Model")
axislegend(ax1; position = :rt, framevisible = false)

ax2 = Axis(fig[2, 1]; ylabel = "Latent Heat (W/m²)")
lines!(ax2, obs_dates, qle_obs; color = :black, linewidth = 1.5, label = "Obs")
lines!(ax2, obs_dates, lhf_model; color = :blue, linewidth = 1.5, label = "Model")
axislegend(ax2; position = :rt, framevisible = false)

ax3 = Axis(fig[3, 1]; xlabel = "Day of year", ylabel = "Sensible Heat (W/m²)")
lines!(ax3, obs_dates, qh_obs; color = :black, linewidth = 1.5, label = "Obs")
lines!(ax3, obs_dates, shf_model; color = :blue, linewidth = 1.5, label = "Model")
axislegend(ax3; position = :rt, framevisible = false)

outpath = joinpath(savedir, "dk_sor_default_2004.png")
CairoMakie.save(outpath, fig)
println("Saved: $outpath")

# Figure 2: Soil CO2 diagnostics
colors = [:red, :orange, :blue, :purple]
fig2 = Figure(size = (1400, 1200))

ax_hr = Axis(fig2[1, 1]; ylabel = "HR (gC/m²/d)",
    title = "DK-Sor 2004 — Soil CO₂ Diagnostics (post-spinup)")
lines!(ax_hr, t_days_co2, hr_gC_ps; color = :brown, linewidth = 1.5)

ax_co2 = Axis(fig2[2, 1]; ylabel = "CO₂ conc (kgC/m³ air-eq)")
for (i, lid) in enumerate(layer_ids)
    vals = sco2_by_layer[lid][post_spinup]
    lines!(ax_co2, t_days_co2, vals; color = colors[i], linewidth = 1.2,
        label = layer_labels[i])
end
axislegend(ax_co2; position = :rt, framevisible = false)

ax_sm = Axis(fig2[3, 1]; ylabel = "Microbe source (kgC/m³/s)")
for (i, lid) in enumerate(layer_ids)
    vals = scms_by_layer[lid][post_spinup]
    lines!(ax_sm, t_days_co2, vals; color = colors[i], linewidth = 1.2,
        label = layer_labels[i])
end
axislegend(ax_sm; position = :rt, framevisible = false)

ax_soc = Axis(fig2[4, 1]; ylabel = "SOC (kgC/m³)", xlabel = "Days since Jan 1, 2004")
for (i, lid) in enumerate(layer_ids)
    vals = soc_by_layer[lid][post_spinup]
    lines!(ax_soc, t_days_co2, vals; color = colors[i], linewidth = 1.2,
        label = layer_labels[i])
end
axislegend(ax_soc; position = :rt, framevisible = false)

outpath2 = joinpath(savedir, "dk_sor_soilco2_2004.png")
CairoMakie.save(outpath2, fig2)
println("Saved: $outpath2")



# Figure 3: Flux diags
fig = Figure(size = (600, 900))

ax1 = Axis(fig[1, 1]; ylabel = "Snow cover", title = "DK-Sor 2004 — Default Parameters")
lines!(ax1, 1:n_obs, scf; color = :black, linewidth = 1.5, label = "Model")
axislegend(ax1; position = :rt, framevisible = false)

ax2 = Axis(fig[2, 1]; ylabel = "ΔT")
lines!(ax2, 1:n_obs, -1 .*(Tair .- Tsoil); color = :black, linewidth = 1.5, label = "Tsoil-Tair")
lines!(ax2, 1:n_obs, -1 .*(Tair .- Tcanopy); color = :blue, linewidth = 1.5, label = "Tcanopy-Tair")
axislegend(ax2; position = :rt, framevisible = false)

outpath = joinpath(savedir, "dk_sor_temps.png")
CairoMakie.save(outpath, fig)
println("Saved: $outpath")
