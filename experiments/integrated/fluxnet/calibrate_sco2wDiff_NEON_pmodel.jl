"""
Calibrate 4 DAMM soil CO₂ biogeochemistry parameters (α_sx, Ea_sx, kM_sx, kM_o2)
against NEON soil CO₂ concentration observations at the CPER site using
Unscented Kalman Inversion (UKI).

Depth: 502 only (~6 cm, matching model layer ~5.7 cm)
Metric: Daily mean time series (~345 values after 20-day spinup)
SOC: Fixed at 2.0 kg C/m³ (override default of 5.0)

This version uses the full LandModel (soil + canopy + snow + soilco2) with
global-style parameters: the Column domain with longlat enables spatial
parameter lookup from data files (soil properties, albedo, etc.) matching
what a global simulation would use for this pixel.

Usage:
    julia --project=.buildkite experiments/integrated/fluxnet/calibrate_sco2_NEON_pmodel.jl
"""

# ── 1. Environment & Imports ──────────────────────────────────────────────────
#import Pkg
#Pkg.activate(joinpath(@__DIR__, "..", "..", "..", ".buildkite"))

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
import ClimaUtilities.TimeManager: date
import ClimaLand.LandSimVis as LandSimVis

import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using CairoMakie
CairoMakie.activate!()
using Statistics
using Logging
import Random
using CSV
using DataFrames
using TOML

# ── 2. Configuration ─────────────────────────────────────────────────────────
const FT = Float64
site_ID = "NEON-cper"
savefolder = "20260211_initial_optCO2diffv3_pmodel"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)
climaland_dir = pkgdir(ClimaLand)

spinup_days = 365
SOC_init = FT(2.0)
dt = Float64(450)  # 7.5 minutes

# Dates - get location from FLUXNET metadata, but use ERA5 forcing
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset)
start_date = start_date - Day(spinup_days)  # include spinup period
spinup_date = start_date + Day(spinup_days)

# UKI settings
rng_seed = 1234
rng = Random.MersenneTwister(rng_seed)
N_iterations = 11

# ── 3. One-time Setup (domain, forcing, LAI) ─────────────────────────────────
# Use global domain settings: 15m depth, 15 vertical elements, stretched grid
zmin = FT(-15.0)
zmax = FT(0.0)
nelements = 15
dz_tuple = FT.((3.0, 0.05))

# Build domain with longlat - this enables spatial parameter lookup
land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
surface_space = land_domain.space.surface
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val)) #CHECK THAT

# Build ERA5 forcing (global-style)
toml_dict_base = LP.create_toml_dict(FT)
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict_base,
    FT,
)
#change some canopy parameters
toml_dict_base.data["canopy_d_coeff"]["value"] = FT(0.67)
toml_dict_base.data["canopy_z_0b_coeff"]["value"] = FT(0.013)
toml_dict_base.data["canopy_z_0m_coeff"]["value"] = FT(0.13)

# LAI from MODIS (global-style)
LAI = ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)

# ── 4. Load & Process NEON Observations ──────────────────────────────────────
csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(site_ID)
obs_df = CSV.read(csv_path, DataFrame)

# Extract depth-502 columns for plots 001–005
co2_cols_502 = [
    Symbol("soilCO2concentrationMean_001_502"),
    Symbol("soilCO2concentrationMean_002_502"),
    Symbol("soilCO2concentrationMean_003_502"),
    Symbol("soilCO2concentrationMean_004_502"),
    Symbol("soilCO2concentrationMean_005_502"),
]

# Row-wise mean across plots (skip missing/NaN)
function rowmean_skipinvalid(row, cols)
    vals = Float64[]
    for c in cols
        v = row[c]
        if !ismissing(v) && !isnan(Float64(v))
            push!(vals, Float64(v))
        end
    end
    return isempty(vals) ? NaN : mean(vals)
end

obs_df[!, :sco2_mean_502] =
    [rowmean_skipinvalid(row, co2_cols_502) for row in eachrow(obs_df)]

# Parse timestamps → DateTime
obs_df[!, :datetime] =
    DateTime.(string.(Int.(obs_df.timestamp_fmt)), dateformat"yyyymmddHHMM")

# Compute inter-sensor variance per row for noise estimate
function rowvar_skipinvalid(row, cols)
    vals = Float64[]
    for c in cols
        v = row[c]
        if !ismissing(v) && !isnan(Float64(v))
            push!(vals, Float64(v))
        end
    end
    return length(vals) >= 2 ? var(vals) : NaN
end

obs_df[!, :sco2_var_502] =
    [rowvar_skipinvalid(row, co2_cols_502) for row in eachrow(obs_df)]

# Add date column
obs_df[!, :date] = Date.(obs_df.datetime)

# Group by date, compute daily means (require ≥24 valid half-hours)
daily_df = combine(
    groupby(obs_df, :date),
    :sco2_mean_502 =>
        (x -> begin
            valid = filter(!isnan, x)
            length(valid) >= 24 ? mean(valid) : NaN
        end) => :daily_mean,
    :sco2_var_502 =>
        (x -> begin
            valid = filter(!isnan, x)
            isempty(valid) ? NaN : mean(valid)
        end) => :daily_var,
)

# Trim to after spinup
daily_df = filter(row -> row.date >= Date(spinup_date), daily_df)
daily_df = filter(row -> !isnan(row.daily_mean), daily_df)
sort!(daily_df, :date)

y_obs = Float64.(daily_df.daily_mean)
obs_dates = daily_df.date
n_obs = length(y_obs)
println("Observation vector length: $n_obs daily values")

# Noise covariance from inter-sensor variance
mean_sensor_var = mean(filter(!isnan, daily_df.daily_var))
noise_cov = mean_sensor_var * EKP.I
println("Noise variance (inter-sensor): $mean_sensor_var ppm²")

# ── 5. Determine Target Model Layer ──────────────────────────────────────────
# The z-coordinates of cell centers in the subsurface domain.
# Layer 1 = bottom, layer N = surface.
z_field = ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z
z_vals = parent(z_field)[:, 1]  # extract as vector
target_depth = -0.06  # 6 cm below surface
target_layer = argmin(abs.(z_vals .- target_depth))
println("Target layer index: $target_layer (z = $(z_vals[target_layer]) m, target = $target_depth m)")

# ── 6. Model Function ────────────────────────────────────────────────────────
function run_model(α_sx_val, Ea_sx_val, kM_sx_val, kM_o2_val, D_co2_val)
    # Write temporary TOML override for calibrated DAMM parameters
    toml_content = """
    [soilCO2_pre_exponential_factor]
    value = $α_sx_val
    type = "float"

    [soilCO2_activation_energy]
    value = $Ea_sx_val
    type = "float"

    [michaelis_constant]
    value = $kM_sx_val
    type = "float"

    [O2_michaelis_constant]
    value = $kM_o2_val
    type = "float"

    [CO2_diffusion_coefficient]
    value = $D_co2_val
    type = "float"
    """
    tmp_toml_path = tempname() * ".toml"
    open(tmp_toml_path, "w") do io
        write(io, toml_content)
    end

    toml_dict = LP.create_toml_dict(FT; override_files = [tmp_toml_path])
    
     # Prognostic components (full land model with snow and soilco2)
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

    # Forcing NamedTuple
    forcing = (; atmos, radiation)
    
    # Construct the P model manually since it is not a default
    photosynthesis = PModel{FT}(land_domain, toml_dict_base)
    conductance = PModelConductance{FT}(toml_dict_base)
    # Use the soil moisture stress function based on soil moisture only
    soil_moisture_stress =
        ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(land_domain, toml_dict_base)

    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)

    # Canopy model with custom components
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        canopy_domain,
        canopy_forcing,
        LAI,
        toml_dict_base;
        prognostic_land_components,
        photosynthesis,
        conductance,
        soil_moisture_stress,
    )

    # ── 7. Snow Model ────────────────────────────────────────────────────────────
    # Snow model with zenith angle dependent albedo
    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict_base)
    snow = Snow.SnowModel(
        FT,
        canopy_domain,
        forcing,
        toml_dict_base,
        dt;
        prognostic_land_components,
        α_snow,
        #scf,
    )
   

    # Full LandModel with global-style spatial parameter lookup
    # The LandModel constructor automatically sets up:
    # - Soil: Van Genuchten, composition, albedo, runoff from spatial data
    # - Canopy: default component models with spatial parameters
    # - Snow: snow model with default parameters
    # - SoilCO2: soil CO2 model (included via prognostic_land_components)
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        land_domain,
        dt;
        prognostic_land_components,
        snow,
        canopy,
    )

    # Custom set_ic! that overrides SOC to SOC_init (2.0 instead of default 5.0)
    base_set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        time_offset,
        land,
    )
    function custom_set_ic!(Y, p, t, model)
        base_set_ic!(Y, p, t, model)
        Y.soilco2.SOC .= SOC_init
    end

    # Diagnostics: sco2_ppm at halfhourly resolution
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    output_vars = ["swc", "tsoil", "si", "sco2", "soc", "so2","sco2_ppm"]
    diags = ClimaLand.default_diagnostics(
        land,
        start_date;
        output_writer = output_writer,
        output_vars,
        reduction_period = :halfhourly,
    )
    

    # Build and run simulation
    simulation = LandSimulation(
        start_date,
        stop_date,
        dt,
        land;
        set_ic! = custom_set_ic!,
        updateat = Second(dt),
        diagnostics = diags,
    )
    solve!(simulation)

    # Clean up temp file
    rm(tmp_toml_path; force = true)

    return simulation, diags
end

# ── 7. Forward Model G ───────────────────────────────────────────────────────
function G(α_sx_val, Ea_sx_val, kM_sx_val, kM_o2_val, D_co2_val)
    simulation, _ = run_model(α_sx_val, Ea_sx_val, kM_sx_val, kM_o2_val, D_co2_val)

    # Extract sco2_ppm at target layer
    (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
        simulation.diagnostics[1].output_writer,
        "sco2_ppm_30m_average";
        layer = target_layer,
    )

    # Convert times to DateTime
    model_dates = times isa Vector{DateTime} ? times : date.(times)

    # Build a DataFrame for daily averaging
    model_df = DataFrame(datetime = model_dates, sco2_ppm = Float64.(data))
    model_df[!, :date] = Date.(model_df.datetime)

    # Only keep dates after spinup
    model_df = filter(row -> row.date >= Date(spinup_date), model_df)

    # Daily means
    model_daily = combine(
        groupby(model_df, :date),
        :sco2_ppm => mean => :daily_mean,
    )
    sort!(model_daily, :date)

    # Match to observation dates
    model_dict = Dict(row.date => row.daily_mean for row in eachrow(model_daily))
    result = Float64[]
    for d in obs_dates
        if haskey(model_dict, d)
            push!(result, model_dict[d])
        else
            push!(result, NaN)
        end
    end

    return result
end

# ── 8. UKI Calibration ───────────────────────────────────────────────────────
println("\n=== Setting up UKI calibration ===")

# 4 constrained Gaussian priors
priors = [
    PD.constrained_gaussian("α_sx", 23835.0, 12000.0, 1000.0, 200000.0),
    PD.constrained_gaussian("Ea_sx", 61000.0, 10000.0, 40000.0, 80000.0),
    PD.constrained_gaussian("kM_sx", 0.005, 0.003, 1e-5, 0.1),
    PD.constrained_gaussian("kM_o2", 0.004, 0.002, 1e-5, 0.1),
    PD.constrained_gaussian("D_co2", 3e-5, 2e-5, 0.1e-5, 6e-5),
]
prior = PD.combine_distributions(priors)

# Create EKP object with Unscented process (ensemble size = 2n+1 = 11 for 5 params)
ensemble_kalman_process = EKP.EnsembleKalmanProcess(
    y_obs,
    noise_cov,
    EKP.Unscented(prior; α_reg = 1.0, update_freq = 1);
    rng,
)

# UKI loop
n_ens = EKP.get_N_ens(ensemble_kalman_process)
Logging.with_logger(SimpleLogger(devnull, Logging.Error)) do
    for i in 1:N_iterations
        println("\n--- UKI Iteration $i / $N_iterations ---")
        params_i = EKP.get_ϕ_final(prior, ensemble_kalman_process)
        println("  Parameter means: ", round.(mean(params_i, dims = 2)[:, 1], sigdigits = 4))

        G_ens = hcat(
            [G(params_i[:, j]...) for j in 1:n_ens]...,
        )
        EKP.update_ensemble!(ensemble_kalman_process, G_ens)
        println("  Ensemble updated.")
    end
end

# ── 9. Results ────────────────────────────────────────────────────────────────
final_params = EKP.get_ϕ_mean_final(prior, ensemble_kalman_process)
println("\n=== Calibration Complete ===")
println("Final parameter means:")
println("  α_sx  = ", round(final_params[1], sigdigits = 5))
println("  Ea_sx = ", round(final_params[2], sigdigits = 5))
println("  kM_sx = ", round(final_params[3], sigdigits = 4))
println("  kM_o2 = ", round(final_params[4], sigdigits = 4))
println("  D_co2 = ", round(final_params[5], sigdigits = 4))

# ── 10. Visualization ────────────────────────────────────────────────────────
savedir = joinpath(@__DIR__, "calibrate_sco2_$(site_ID)_output")
savedir = "/Users/evametz/Documents/PostDoc/Projekte/CliMA/Siteruns/Calibrations_$(site_ID)/$(savefolder)/"
mkpath(savedir)

# Write final parameters to file
params_file = joinpath(savedir, "final_parameters.txt")
open(params_file, "w") do io
    println(io, "=== Calibration Complete ===")
    println(io, "Final parameter means:")
    println(io, "  α_sx  = ", round(final_params[1], sigdigits = 5))
    println(io, "  Ea_sx = ", round(final_params[2], sigdigits = 5))
    println(io, "  kM_sx = ", round(final_params[3], sigdigits = 4))
    println(io, "  kM_o2 = ", round(final_params[4], sigdigits = 4))
    println(io, "  D_co2 = ", round(final_params[5], sigdigits = 4))
end
println("Saved final parameters to: $params_file")

# Parameter convergence + error plot
dim_size = sum(length.(EKP.batch(prior)))
fig1 = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500))
for i in 1:dim_size
    EKP.Visualize.plot_ϕ_over_iters(
        fig1[1, i],
        ensemble_kalman_process,
        prior,
        i,
    )
end
EKP.Visualize.plot_error_over_iters(
    fig1[1, dim_size + 1],
    ensemble_kalman_process,
)
CairoMakie.save(joinpath(savedir, "sco2_param_convergence.png"), fig1)
println("Saved: sco2_param_convergence.png")

# First vs last ensemble G vs observations
fig2 = CairoMakie.Figure(size = (1200, 500))
first_G_ensemble = EKP.get_g(ensemble_kalman_process, 1)
last_iter = EKP.get_N_iterations(ensemble_kalman_process)
last_G_ensemble = EKP.get_g(ensemble_kalman_process, last_iter)
n_ens = EKP.get_N_ens(ensemble_kalman_process)

ax = Axis(
    fig2[1, 1];
    title = "Soil CO₂ (ppm) at 6 cm: first vs last UKI iteration (n=$n_ens, iters 1 vs $last_iter)",
    xlabel = "Day index (after spinup)",
    ylabel = "Soil CO₂ (ppm)",
)

for g_col in eachcol(first_G_ensemble)
    lines!(ax, 1:length(g_col), g_col; color = (:red, 0.5), linewidth = 1)
end
for g_col in eachcol(last_G_ensemble)
    lines!(ax, 1:length(g_col), g_col; color = (:blue, 0.5), linewidth = 1)
end
lines!(
    ax,
    1:n_obs,
    y_obs;
    color = (:black, 0.8),
    linewidth = 2.5,
)

axislegend(
    ax,
    [
        LineElement(color = :red, linewidth = 2),
        LineElement(color = :blue, linewidth = 2),
        LineElement(color = :black, linewidth = 3),
    ],
    ["First ensemble", "Last ensemble", "NEON obs (502)"];
    position = :rt,
    framevisible = false,
)

CairoMakie.resize_to_layout!(fig2)
CairoMakie.save(joinpath(savedir, "sco2_G_first_and_last.png"), fig2)
println("Saved: sco2_G_first_and_last.png")
println("\nDone.")



# plot the best one

# Find the ensemble member with lowest RMSE vs observations
  rmse_per_member = [sqrt(mean((last_G_ensemble[:, j] .- y_obs).^2)) for j in 1:size(last_G_ensemble, 2)]
  best_idx = argmin(rmse_per_member)
  println("Best member: $best_idx, RMSE: $(round(rmse_per_member[best_idx], sigdigits=4))")
  println("Worst member: $(argmax(rmse_per_member)), RMSE: $(round(maximum(rmse_per_member), sigdigits=4))")

  # Get the corresponding parameters
  final_params_all = EKP.get_ϕ_final(prior, ensemble_kalman_process)
  best_params = final_params_all[:, best_idx]
  println("Best params: α_sx=$(round(best_params[1],sigdigits=5)), Ea_sx=$(round(best_params[2],sigdigits=5)),
  kM_sx=$(round(best_params[3],sigdigits=4)), kM_o2=$(round(best_params[4],sigdigits=4)), D_co2=$(round(best_params[5],sigdigits=4))")

  # Write best parameters to file
  best_params_file = joinpath(savedir, "best_parameters.txt")
  open(best_params_file, "w") do io
      println(io, "=== Best Ensemble Member ===")
      println(io, "Best member: $best_idx")
      println(io, "RMSE: $(round(rmse_per_member[best_idx], sigdigits=4)) ppm")
      println(io, "Worst member: $(argmax(rmse_per_member))")
      println(io, "Worst RMSE: $(round(maximum(rmse_per_member), sigdigits=4)) ppm")
      println(io, "\nBest parameters:")
      println(io, "  α_sx  = $(round(best_params[1], sigdigits=5))")
      println(io, "  Ea_sx = $(round(best_params[2], sigdigits=5))")
      println(io, "  kM_sx = $(round(best_params[3], sigdigits=4))")
      println(io, "  kM_o2 = $(round(best_params[4], sigdigits=4))")
      println(io, "  D_co2 = $(round(best_params[5], sigdigits=4))")
  end
  println("Saved best parameters to: $best_params_file")

  # Plot best member vs obs
  fig = Figure(size=(1000, 500))
  ax = Axis(fig[1,1]; xlabel="Day index", ylabel="Soil CO₂ (ppm)", title="Best ensemble member vs NEON obs (502)")
  lines!(ax, 1:n_obs, y_obs; color=:black, linewidth=2, label="NEON obs")
  lines!(ax, 1:n_obs, last_G_ensemble[:, best_idx]; color=:blue, linewidth=1.5, label="Best member ($(best_idx))")
  axislegend(ax; position=:rt)
  save(joinpath(savedir, "sco2_best_vs_obs.png"), fig)
  println("Saved: sco2_best_vs_obs.png")

# ── 11. Run final model with best parameters and create timeseries plots ─────
println("\n=== Generating timeseries plots with best parameters ===")

# Run model with best parameters to get diagnostics
simulation_best, _ = run_model(best_params...)

# Plot SWC daily mean at target layer (similar to sco2_best_vs_obs)
(times_swc, swc_data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
    simulation_best.diagnostics[1].output_writer,
    "swc_30m_average";
    layer = target_layer,
)
swc_dates = times_swc isa Vector{DateTime} ? times_swc : date.(times_swc)
swc_df = DataFrame(datetime = swc_dates, swc = Float64.(swc_data))
swc_df[!, :date] = Date.(swc_df.datetime)
swc_df = filter(row -> row.date >= Date(spinup_date), swc_df)
swc_daily = combine(groupby(swc_df, :date), :swc => mean => :daily_mean)
sort!(swc_daily, :date)

fig_swc = Figure(size = (1000, 500))
ax_swc = Axis(
    fig_swc[1, 1];
    xlabel = "Day index (after spinup)",
    ylabel = "SWC",
    title = "Soil water content at 6 cm (daily mean)",
)
lines!(ax_swc, 1:nrow(swc_daily), swc_daily.daily_mean; color = :blue, linewidth = 1.5, label = "SWC")
axislegend(ax_swc; position = :rt)
CairoMakie.save(joinpath(savedir, "swc_daily_best.png"), fig_swc)
println("Saved: swc_daily_best.png")

# Plot soil temperature daily mean at target layer
(times_tsoil, tsoil_data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
    simulation_best.diagnostics[1].output_writer,
    "tsoil_30m_average";
    layer = target_layer,
)
tsoil_dates = times_tsoil isa Vector{DateTime} ? times_tsoil : date.(times_tsoil)
tsoil_df = DataFrame(datetime = tsoil_dates, tsoil = Float64.(tsoil_data))
tsoil_df[!, :date] = Date.(tsoil_df.datetime)
tsoil_df = filter(row -> row.date >= Date(spinup_date), tsoil_df)
tsoil_daily = combine(groupby(tsoil_df, :date), :tsoil => mean => :daily_mean)
sort!(tsoil_daily, :date)

fig_tsoil = Figure(size = (1000, 500))
ax_tsoil = Axis(
    fig_tsoil[1, 1];
    xlabel = "Day index (after spinup)",
    ylabel = "Tsoil",
    title = "Soil temperature at 6 cm (daily mean)",
)
lines!(
    ax_tsoil,
    1:nrow(tsoil_daily),
    tsoil_daily.daily_mean;
    color = :red,
    linewidth = 1.5,
    label = "Tsoil",
)
axislegend(ax_tsoil; position = :rt)
CairoMakie.save(joinpath(savedir, "tsoil_daily_best.png"), fig_tsoil)
println("Saved: tsoil_daily_best.png")

# Plot O2 daily mean at target layer
(times_so2, so2_data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
    simulation_best.diagnostics[1].output_writer,
    "so2_30m_average";
    layer = target_layer,
)
so2_dates = times_so2 isa Vector{DateTime} ? times_so2 : date.(times_so2)
so2_df = DataFrame(datetime = so2_dates, so2 = Float64.(so2_data))
so2_df[!, :date] = Date.(so2_df.datetime)
so2_df = filter(row -> row.date >= Date(spinup_date), so2_df)
so2_daily = combine(groupby(so2_df, :date), :so2 => mean => :daily_mean)
sort!(so2_daily, :date)

fig_so2 = Figure(size = (1000, 500))
ax_so2 = Axis(
    fig_so2[1, 1];
    xlabel = "Day index (after spinup)",
    ylabel = "O2",
    title = "Soil O₂ at 6 cm (daily mean)",
)
lines!(
    ax_so2,
    1:nrow(so2_daily),
    so2_daily.daily_mean;
    color = :green,
    linewidth = 1.5,
    label = "O2",
)
axislegend(ax_so2; position = :rt)
CairoMakie.save(joinpath(savedir, "so2_daily_best.png"), fig_so2)
println("Saved: so2_daily_best.png")


