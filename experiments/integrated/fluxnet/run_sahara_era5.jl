"""
Run a single column in the central Sahara (lat=23.5°N, lon=5°E) with ERA5 forcing
and the same model settings as DK-Sor (canopy + snow + soil + soilco2).

4-month simulation (Jan–Apr 2008) using low-res ERA5.

Usage:
    julia --project=.buildkite experiments/integrated/fluxnet/run_sahara_era5.jl
"""

# ── 1. Imports ───────────────────────────────────────────────────────────────
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
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeManager: date

using CairoMakie
CairoMakie.activate!()
using Statistics

# ── 2. Configuration ────────────────────────────────────────────────────────
const FT = Float64
lat = FT(23.5)   # central Sahara
long = FT(5.0)

spinup_days = 30
dt = Float64(900)  # 15 minutes

start_date = DateTime(2008, 1, 1) - Day(spinup_days)
stop_date = DateTime(2009, 1, 1)
spinup_date = start_date + Day(spinup_days)
println("Simulation: $start_date to $stop_date")
println("Spinup until: $spinup_date")

# ── 3. Domain ────────────────────────────────────────────────────────────────
zmin = FT(-15.0)
zmax = FT(0.0)
nelements = 20
dz_tuple = FT.((3.0, 0.05))

land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
surface_space = land_domain.space.surface
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

# ── 4. Forcing (ERA5) ───────────────────────────────────────────────────────
toml_dict = LP.create_toml_dict(FT)
(; atmos, radiation) = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict,
    FT;
    use_lowres_forcing = true,
)

# LAI from MODIS (standard ClimaLand approach)
LAI = ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)

# ── 5. Build model (same settings as DK-Sor) ────────────────────────────────
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
forcing_nt = (; atmos, radiation, ground = ClimaLand.PrognosticGroundConditions{FT}())

# Canopy parameters (same as DK-Sor)
χl = FT(0.25)
α_PAR_leaf = FT(0.1)
α_NIR_leaf = FT(0.45)
τ_PAR_leaf = FT(0.05)
τ_NIR_leaf = FT(0.25)
Ω = FT(1)
rooting_depth = FT(0.3)

biomass = ClimaLand.Canopy.PrescribedBiomassModel{FT}(
    land_domain,
    LAI,
    toml_dict;
    rooting_depth,
    height = FT(1),
    SAI = FT(0.0),
    RAI = FT(0.0),
)

radiation_parameters = (;
    Ω,
    G_Function = CLMGFunction(χl),
    α_PAR_leaf,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
)
radiative_transfer = Canopy.TwoStreamModel{FT}(
    canopy_domain,
    toml_dict;
    radiation_parameters,
)

# Soil parameters (same as DK-Sor)
ν_ss_cf = FT(0.12)
ν_ss_sand = FT(0.47)
ν_ss_om = FT(0.03)
K_sat = FT(1e-5)
ν = FT(0.45)
θ_r = FT(0.07)
vg_n = FT(1.6)
vg_α = FT(1.6)
hydrology_cm = ClimaLand.Soil.vanGenuchten{FT}(; α = vg_α, n = vg_n)

retention_parameters = (; ν, hydrology_cm, θ_r, K_sat)
composition_parameters = (;
    ν_ss_om,
    ν_ss_quartz = ν_ss_sand,
    ν_ss_gravel = ν_ss_cf,
)

hydraulics = Canopy.PlantHydraulicsModel{FT}(
    canopy_domain, toml_dict;
    n_stem = 0,
    h_stem = FT(0),
    h_leaf = FT(1),
)

canopy = Canopy.CanopyModel{FT}(
    canopy_domain, forcing_nt, LAI, toml_dict;
    prognostic_land_components,
    photosynthesis = Canopy.PModel{FT}(canopy_domain, toml_dict),
    conductance = Canopy.PModelConductance{FT}(toml_dict),
    soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(
        land_domain, toml_dict; soil_params = (; ν, θ_r),
    ),
    biomass,
    radiative_transfer,
    hydraulics,
)

forcing = (; atmos, radiation)
S_s = FT(1e-3)
soil = ClimaLand.Soil.EnergyHydrology{FT}(
    land_domain, forcing, toml_dict;
    prognostic_land_components,
    retention_parameters,
    composition_parameters,
    S_s,
    additional_sources = (ClimaLand.RootExtraction{FT}(),),
)

land = LandModel{FT}(
    forcing, LAI, toml_dict, land_domain, dt;
    prognostic_land_components, canopy, soil,
)

# ── 6. Initial conditions (same as DK-Sor) ──────────────────────────────────
function custom_set_ic!(Y, p, t, model)
    earth_param_set = ClimaLand.get_earth_param_set(model.soil)
    ClimaUtilities.TimeVaryingInputs.evaluate!(p.drivers.T, atmos.T, t)

    (; θ_r, ν, ρc_ds) = model.soil.parameters
    @. Y.soil.ϑ_l = θ_r + (ν - θ_r) * FT(0.95)
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
        SOC_top = FT(15.0)
        SOC_bot = FT(0.5)
        τ_soc = FT(1.0 / log(SOC_top / SOC_bot))
        z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
        @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
    end
end

# ── 7. Diagnostics & Simulation ─────────────────────────────────────────────
output_writer = ClimaDiagnostics.Writers.DictWriter()
diags = ClimaLand.default_diagnostics(
    land, start_date;
    output_writer, output_vars = :long, reduction_period = :daily,
)

simulation = LandSimulation(
    start_date, stop_date, dt, land;
    set_ic! = custom_set_ic!, updateat = Second(dt), diagnostics = diags,
)

println("Running simulation...")
@time solve!(simulation)
println("Simulation complete.")

# ── 8. Debug: soilco2 state at end ──────────────────────────────────────────
let Y = simulation._integrator.u, p = simulation._integrator.p
    println("\n=== SoilCO2 debug at end of simulation ===")
    println("  Y.soilco2.CO2:  ", extrema(parent(Y.soilco2.CO2)))
    println("  Y.soilco2.O2_f: ", extrema(parent(Y.soilco2.O2_f)))
    println("  Y.soilco2.SOC:  ", extrema(parent(Y.soilco2.SOC)))
    println("  p.soilco2.Sm:   ", extrema(parent(p.soilco2.Sm)))
    println("  p.soilco2.T:    ", extrema(parent(p.soilco2.T)))
    println("  p.soilco2.θ_a:  ", extrema(parent(p.soilco2.θ_a)))
    println("  p.soilco2.D:    ", extrema(parent(p.soilco2.D)))
    println("  p.soilco2.O2_avail: ", extrema(parent(p.soilco2.O2_avail)))
    println("  p.soil.T:       ", extrema(parent(p.soil.T)))
    println("  p.soil.θ_l:     ", extrema(parent(p.soil.θ_l)))
    println("  Y.soil.θ_i:     ", extrema(parent(Y.soil.θ_i)))
    println("  Y.soil.ϑ_l:     ", extrema(parent(Y.soil.ϑ_l)))
end

# ── 9. Extract diagnostics ──────────────────────────────────────────────────
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
    println("  Matched $n_match / $(length(target_dates)) dates")
    return result
end

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

# Simulation date range (post-spinup)
spinup_date_d = Date(spinup_date)
stop_date_d = Date(stop_date)

# Get all model dates from a surface diagnostic
hr_dates, hr_vals_raw = extract_daily_diag_layer(simulation, "hr_1d_average", nothing)
hr_gC = hr_vals_raw .* 12.0 .* 86400.0
post_spinup = hr_dates .>= spinup_date_d
hr_dates_ps = hr_dates[post_spinup]
hr_gC_ps = hr_gC[post_spinup]

sim_dates = hr_dates_ps
t_days = [Dates.value(d - spinup_date_d) for d in hr_dates_ps]

# Surface fluxes
gpp_model = extract_daily_diag(simulation, "gpp_1d_average", sim_dates)
er_model = extract_daily_diag(simulation, "er_1d_average", sim_dates)
lhf_model = extract_daily_diag(simulation, "lhf_1d_average", sim_dates)
shf_model = extract_daily_diag(simulation, "shf_1d_average", sim_dates)
nee_model_gC = (er_model .- gpp_model) .* 12.0 .* 86400.0

println("  hr: $(length(hr_dates_ps)) dates, range $(round(minimum(hr_gC_ps), sigdigits=4)) to $(round(maximum(hr_gC_ps), sigdigits=4)) gC/m²/d")

# Depth-resolved soil CO2
n_layers = let w = nothing
    for d in simulation.diagnostics
        if haskey(d.output_writer.dict, "sco2_1d_average")
            w = d.output_writer; break
        end
    end
    first_field = first(values(w["sco2_1d_average"]))
    size(parent(first_field), 1)
end
println("  Soil has $n_layers layers")

layer_ids = [n_layers, round(Int, n_layers * 3 ÷ 4), n_layers ÷ 2, 1]
layer_labels = String[]
z_centers = parent(land_domain.fields.z)
println("  Layer depths (center): ", [round(z_centers[i], digits=3) for i in layer_ids])

sco2_ppm_by_layer = Dict{Int, Vector{Float64}}()
scms_by_layer = Dict{Int, Vector{Float64}}()
soc_by_layer = Dict{Int, Vector{Float64}}()

for lid in layer_ids
    d, v = extract_daily_diag_layer(simulation, "sco2_ppm_1d_average", lid)
    sco2_ppm_by_layer[lid] = v
    d, v = extract_daily_diag_layer(simulation, "scms_1d_average", lid)
    scms_by_layer[lid] = v
    d, v = extract_daily_diag_layer(simulation, "soc_1d_average", lid)
    soc_by_layer[lid] = v
    push!(layer_labels, "z=$(round(z_centers[lid], digits=2))m")
    println("  Layer $lid (z=$(round(z_centers[lid], digits=2))m): " *
            "CO2_ppm=$(round(minimum(sco2_ppm_by_layer[lid]), sigdigits=4))..$(round(maximum(sco2_ppm_by_layer[lid]), sigdigits=4)), " *
            "Sm=$(round(minimum(scms_by_layer[lid]), sigdigits=4))..$(round(maximum(scms_by_layer[lid]), sigdigits=4)) kgC/m³/s, " *
            "SOC=$(round(minimum(soc_by_layer[lid]), sigdigits=4))..$(round(maximum(soc_by_layer[lid]), sigdigits=4)) kgC/m³")
end

# O2 fraction
has_o2 = false
so2_by_layer = Dict{Int, Vector{Float64}}()
try
    for lid in layer_ids
        d, v = extract_daily_diag_layer(simulation, "so2_1d_average", lid)
        so2_by_layer[lid] = v
    end
    global has_o2 = true
catch e
    println("  O2 diagnostic not available: $e")
end

# ── 10. Plot ─────────────────────────────────────────────────────────────────
savedir = joinpath(@__DIR__, "sahara_era5_output")
mkpath(savedir)

# Figure 1: Surface fluxes
fig = Figure(size = (1200, 900))
ax1 = Axis(fig[1, 1]; ylabel = "NEE (gC/m²/d)",
    title = "Sahara ($(lat)°N, $(long)°E) — Jan–Apr 2008")
lines!(ax1, sim_dates, nee_model_gC; color = :blue, linewidth = 1.5, label = "Model NEE")
axislegend(ax1; position = :rt, framevisible = false)

ax2 = Axis(fig[2, 1]; ylabel = "Latent Heat (W/m²)")
lines!(ax2, sim_dates, lhf_model; color = :blue, linewidth = 1.5)

ax3 = Axis(fig[3, 1]; xlabel = "Date", ylabel = "Sensible Heat (W/m²)")
lines!(ax3, sim_dates, shf_model; color = :blue, linewidth = 1.5)

outpath = joinpath(savedir, "sahara_fluxes_2008.png")
CairoMakie.save(outpath, fig)
println("Saved: $outpath")

# Figure 2: Soil CO2 diagnostics (same layout as DK-Sor)
colors = [:red, :orange, :blue, :purple]
n_panels = 3 + (has_o2 ? 1 : 0) + 1  # HR + CO2 + Sm + (O2) + SOC
fig2 = Figure(size = (1400, 300 * n_panels))

ax_hr = Axis(fig2[1, 1]; ylabel = "HR (gC/m²/d)",
    title = "Sahara ($(lat)°N, $(long)°E) — Soil CO₂ Diagnostics (post-spinup)")
lines!(ax_hr, t_days, hr_gC_ps; color = :brown, linewidth = 1.5)

ax_co2 = Axis(fig2[2, 1]; ylabel = "CO₂ (ppm)")
for (i, lid) in enumerate(layer_ids)
    vals = sco2_ppm_by_layer[lid][post_spinup]
    lines!(ax_co2, t_days, vals; color = colors[i], linewidth = 1.2,
        label = layer_labels[i])
end
axislegend(ax_co2; position = :rt, framevisible = false)

ax_sm = Axis(fig2[3, 1]; ylabel = "Microbe source (mgC/m³/s)")
for (i, lid) in enumerate(layer_ids)
    vals = scms_by_layer[lid][post_spinup] .* 1e6
    lines!(ax_sm, t_days, vals; color = colors[i], linewidth = 1.2,
        label = layer_labels[i])
end
axislegend(ax_sm; position = :rt, framevisible = false)

panel_idx = 4
if has_o2
    ax_o2 = Axis(fig2[panel_idx, 1]; ylabel = "O₂ fraction (mol/mol)")
    for (i, lid) in enumerate(layer_ids)
        vals = so2_by_layer[lid][post_spinup]
        lines!(ax_o2, t_days, vals; color = colors[i], linewidth = 1.2,
            label = layer_labels[i])
    end
    axislegend(ax_o2; position = :rt, framevisible = false)
    panel_idx += 1
end

ax_soc = Axis(fig2[panel_idx, 1]; ylabel = "SOC (kgC/m³)",
    xlabel = "Days since Jan 1, 2008")
for (i, lid) in enumerate(layer_ids)
    vals = soc_by_layer[lid][post_spinup]
    lines!(ax_soc, t_days, vals; color = colors[i], linewidth = 1.2,
        label = layer_labels[i])
end
axislegend(ax_soc; position = :rt, framevisible = false)

outpath2 = joinpath(savedir, "sahara_soilco2_2008.png")
CairoMakie.save(outpath2, fig2)
println("Saved: $outpath2")

# Figure 3: Flux component breakdown
cshf = extract_daily_diag(simulation, "cshf_1d_average", sim_dates)
soilshf = extract_daily_diag(simulation, "soilshf_1d_average", sim_dates)
clhf = extract_daily_diag(simulation, "clhf_1d_average", sim_dates)
soillhf = extract_daily_diag(simulation, "soillhf_1d_average", sim_dates)
gpp_all = extract_daily_diag(simulation, "gpp_1d_average", sim_dates)
gs_model = extract_daily_diag(simulation, "gs_1d_average", sim_dates)
scf = extract_daily_diag(simulation, "snowc_1d_average", sim_dates)

fig3 = Figure(size = (600, 1500))

ax1 = Axis(fig3[1, 1]; ylabel = "Snow cover",
    title = "Sahara ($(lat)°N, $(long)°E) — Flux Diagnostics")
lines!(ax1, sim_dates, scf; color = :black, linewidth = 1.5)

ax2 = Axis(fig3[2, 1]; ylabel = "SHF components (W/m²)")
lines!(ax2, sim_dates, soilshf .* (1 .- scf); color = :green, linewidth = 1.5, label = "Soil")
lines!(ax2, sim_dates, cshf; color = :red, linewidth = 1.5, label = "Canopy")
axislegend(ax2; position = :rt, framevisible = false)

ax3 = Axis(fig3[3, 1]; ylabel = "LHF components (W/m²)")
lines!(ax3, sim_dates, soillhf .* (1 .- scf); color = :green, linewidth = 1.5, label = "Soil")
lines!(ax3, sim_dates, clhf; color = :red, linewidth = 1.5, label = "Canopy")
axislegend(ax3; position = :rt, framevisible = false)

ax4 = Axis(fig3[4, 1]; ylabel = "GPP (mol CO₂/m²/s)")
lines!(ax4, sim_dates, gpp_all; color = :green, linewidth = 1.5)

ax5 = Axis(fig3[5, 1]; ylabel = "gs (mol/m²/s)", xlabel = "Date")
lines!(ax5, sim_dates, clamp.(gs_model, 0, 0.1); color = :green, linewidth = 1.5)

outpath3 = joinpath(savedir, "sahara_fluxes_debug.png")
CairoMakie.save(outpath3, fig3)
println("Saved: $outpath3")

# Figure 4: Water budget
precip_dates, precip_vals = extract_daily_diag_layer(simulation, "precip_1d_average", nothing)
precip_mm = precip_vals .* 86400.0
et_dates, et_vals = extract_daily_diag_layer(simulation, "et_1d_average", nothing)
et_mm = et_vals .* 86400.0

swc_by_layer = Dict{Int, Vector{Float64}}()
si_by_layer = Dict{Int, Vector{Float64}}()
for lid in layer_ids
    d, v = extract_daily_diag_layer(simulation, "swc_1d_average", lid)
    swc_by_layer[lid] = v
    d, v = extract_daily_diag_layer(simulation, "si_1d_average", lid)
    si_by_layer[lid] = v
    println("  Layer $lid (z=$(round(z_centers[lid], digits=2))m): " *
            "θ_l=$(round(minimum(swc_by_layer[lid]), sigdigits=4))..$(round(maximum(swc_by_layer[lid]), sigdigits=4)), " *
            "θ_i=$(round(minimum(si_by_layer[lid]), sigdigits=4))..$(round(maximum(si_by_layer[lid]), sigdigits=4))")
end

fig4 = Figure(size = (1400, 1000))

ax_p = Axis(fig4[1, 1]; ylabel = "Precip (mm/d)",
    title = "Sahara ($(lat)°N, $(long)°E) — Water Budget (post-spinup)")
precip_ps = precip_mm[precip_dates .>= spinup_date_d]
et_ps = et_mm[et_dates .>= spinup_date_d]
t_days_w = [Dates.value(d - spinup_date_d) for d in precip_dates[precip_dates .>= spinup_date_d]]
lines!(ax_p, t_days_w, precip_ps; color = :blue, linewidth = 1.5, label = "Precip")
lines!(ax_p, t_days_w, et_ps; color = :red, linewidth = 1.5, label = "ET")
axislegend(ax_p; position = :rt, framevisible = false)

ax_swc = Axis(fig4[2, 1]; ylabel = "Soil Water Content (m³/m³)")
for (i, lid) in enumerate(layer_ids)
    vals = swc_by_layer[lid][precip_dates .>= spinup_date_d]
    lines!(ax_swc, t_days_w, vals; color = colors[i], linewidth = 1.2,
        label = layer_labels[i])
end
axislegend(ax_swc; position = :rt, framevisible = false)

ax_si = Axis(fig4[3, 1]; ylabel = "Soil Ice (m³/m³)",
    xlabel = "Days since Jan 1, 2008")
for (i, lid) in enumerate(layer_ids)
    vals = si_by_layer[lid][precip_dates .>= spinup_date_d]
    lines!(ax_si, t_days_w, vals; color = colors[i], linewidth = 1.2,
        label = layer_labels[i])
end
axislegend(ax_si; position = :rt, framevisible = false)

outpath4 = joinpath(savedir, "sahara_water_2008.png")
CairoMakie.save(outpath4, fig4)
println("Saved: $outpath4")

println("\nDone.")
