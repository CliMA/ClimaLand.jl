"""
Calibrate 4 DAMM soil CO₂ biogeochemistry parameters (α_sx, Ea_sx, kM_sx, kM_o2)
against NEON soil CO₂ concentration observations at the CPER site using
Ensemble Kalman Inversion (EKI).

Depth: 502 only (~6 cm, matching model layer ~5.7 cm)
Metric: Daily mean time series (~345 values after 20-day spinup)
SOC: Fixed at 2.0 kg C/m³ (override default of 5.0)

Usage:
    julia --project=.buildkite experiments/integrated/fluxnet/calibrate_sco2_cper.jl
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
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeManager: date

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
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)
climaland_dir = pkgdir(ClimaLand)

spinup_days = 20
SOC_init = FT(2.0)
dt = Float64(450)  # 7.5 minutes

# Dates
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset)
spinup_date = start_date + Day(spinup_days)

# EKI settings
rng_seed = 1234
rng = Random.MersenneTwister(rng_seed)
ensemble_size = 10
N_iterations = 5

# ── 3. One-time Setup (domain, forcing, LAI, site params) ────────────────────
(; dz_tuple, nelements, zmin, zmax) =
    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
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
) = FluxnetSimulations.get_parameters(FT, Val(site_ID_val))

# PFT: C4 grass
pft_pcts = [
    0.0, # NET_Temp
    0.0, # NET_Bor
    0.0, # NDT_Bor
    0.0, # BET_Trop
    0.0, # BET_Temp
    0.0, # BDT_Trop
    0.0, # BDT_Temp
    0.0, # BDT_Bor
    0.0, # BES_Temp
    0.0, # BDS_Temp
    0.0, # BDT_Bor
    0.0, # C3G_A
    0.0, # C3G_NA
    1.0, # C4G
]
(
    Ω,
    α_PAR_leaf,
    α_NIR_leaf,
    τ_PAR_leaf,
    τ_NIR_leaf,
    ϵ_canopy,
    χl,
    ac_canopy,
    g1,
    Vcmax25,
    f_root_to_shoot,
    K_sat_plant,
    ψ63,
    plant_ν,
    rooting_depth,
) = FT.(params_from_pfts(pft_pcts))

# Build domain
land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

# Compartments for plant hydraulics
compartment_midpoints =
    n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
compartment_surfaces = n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]

# Build forcing (once)
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

# LAI (once)
surface_space = land_domain.space.surface
LAI = ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long)
RAI = maxLAI * f_root_to_shoot

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
function run_model(α_sx_val, Ea_sx_val, kM_sx_val, kM_o2_val)
    # Write temporary TOML override
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
    """
    tmp_toml_path = tempname() * ".toml"
    open(tmp_toml_path, "w") do io
        write(io, toml_content)
    end

    toml_dict = LP.create_toml_dict(FT; override_files = [tmp_toml_path])

    # Prognostic components
    prognostic_land_components = (:canopy, :soil, :soilco2)

    # Soil model
    soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
        PAR_albedo = soil_α_PAR,
        NIR_albedo = soil_α_NIR,
    )
    runoff = ClimaLand.Soil.Runoff.SurfaceRunoff()
    retention_parameters = (;
        ν = soil_ν,
        θ_r,
        K_sat = soil_K_sat,
        hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    )
    composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)
    soil_forcing = (; atmos, radiation)
    soil = Soil.EnergyHydrology{FT}(
        land_domain,
        soil_forcing,
        toml_dict;
        prognostic_land_components,
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
        albedo = soil_albedo,
        runoff,
        retention_parameters,
        composition_parameters,
        S_s = soil_S_s,
        z_0m = z_0m_soil,
        z_0b = z_0b_soil,
        emissivity = soil_ϵ,
    )

    # Soil CO2 model
    co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
    drivers = Soil.Biogeochemistry.SoilDrivers(co2_prognostic_soil, atmos)
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(land_domain, drivers, toml_dict)

    # Canopy model
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
        ϵ_canopy,
    )
    conductance = Canopy.MedlynConductanceModel{FT}(canopy_domain, toml_dict; g1)
    photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25)
    photosynthesis =
        FarquharModel{FT}(canopy_domain, toml_dict; photosynthesis_parameters)
    hydraulics = Canopy.PlantHydraulicsModel{FT}(
        canopy_domain,
        toml_dict;
        n_stem,
        n_leaf,
        h_stem,
        h_leaf,
        ν = plant_ν,
        S_s = plant_S_s,
        conductivity_model,
        retention_model,
    )
    height = h_stem + h_leaf
    biomass =
        Canopy.PrescribedBiomassModel{FT}(; LAI, SAI, RAI, rooting_depth, height)
    energy = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy)
    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy_forcing = (; atmos, radiation, ground)
    canopy = Canopy.CanopyModel{FT}(
        canopy_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        prognostic_land_components,
        radiative_transfer,
        photosynthesis,
        conductance,
        hydraulics,
        energy,
        biomass,
    )

    # Integrated land model
    land = SoilCanopyModel{FT}(soilco2, soil, canopy)

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
    output_vars = ["sco2_ppm"]
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

    return simulation
end

# ── 7. Forward Model G ───────────────────────────────────────────────────────
function G(α_sx_val, Ea_sx_val, kM_sx_val, kM_o2_val)
    simulation = run_model(α_sx_val, Ea_sx_val, kM_sx_val, kM_o2_val)

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

# ── 8. EKI Calibration ───────────────────────────────────────────────────────
println("\n=== Setting up EKI calibration ===")

# 4 constrained Gaussian priors
priors = [
    PD.constrained_gaussian("α_sx", 23835.0, 12000.0, 1000.0, 200000.0),
    PD.constrained_gaussian("Ea_sx", 61000.0, 10000.0, 40000.0, 80000.0),
    PD.constrained_gaussian("kM_sx", 0.005, 0.003, 1e-5, 0.1),
    PD.constrained_gaussian("kM_o2", 0.004, 0.002, 1e-5, 0.1),
]
prior = PD.combine_distributions(priors)

# Initial ensemble
initial_ensemble = EKP.construct_initial_ensemble(rng, prior, ensemble_size)

# Create EKP object
ensemble_kalman_process = EKP.EnsembleKalmanProcess(
    initial_ensemble,
    y_obs,
    noise_cov,
    EKP.Inversion();
    scheduler = EKP.DataMisfitController(
        terminate_at = Inf,
        on_terminate = "continue",
    ),
    rng,
)

# EKI loop
Logging.with_logger(SimpleLogger(devnull, Logging.Error)) do
    for i in 1:N_iterations
        println("\n--- EKI Iteration $i / $N_iterations ---")
        params_i = EKP.get_ϕ_final(prior, ensemble_kalman_process)
        println("  Parameter means: ", round.(mean(params_i, dims = 2)[:, 1], sigdigits = 4))

        G_ens = hcat(
            [G(params_i[:, j]...) for j in 1:ensemble_size]...,
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

# ── 10. Visualization ────────────────────────────────────────────────────────
savedir = joinpath(@__DIR__, "calibrate_sco2_cper_output")
mkpath(savedir)

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
    title = "Soil CO₂ (ppm) at 6 cm: first vs last EKI iteration (n=$n_ens, iters 1 vs $last_iter)",
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
  kM_sx=$(round(best_params[3],sigdigits=4)), kM_o2=$(round(best_params[4],sigdigits=4))")

  # Plot best member vs obs
  fig = Figure(size=(1000, 500))
  ax = Axis(fig[1,1]; xlabel="Day index", ylabel="Soil CO₂ (ppm)", title="Best ensemble member vs NEON obs (502)")
  lines!(ax, 1:n_obs, y_obs; color=:black, linewidth=2, label="NEON obs")
  lines!(ax, 1:n_obs, last_G_ensemble[:, best_idx]; color=:blue, linewidth=1.5, label="Best member ($(best_idx))")
  axislegend(ax; position=:rt)
  save("sco2_best_vs_obs.png", fig)
