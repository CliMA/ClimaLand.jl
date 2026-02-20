"""
Calibrate canopy and DAMM soil CO₂ parameters against daily NEE, latent heat (LH),
and sensible heat (SH) observations at the DK-Sor (Soroe, Denmark) beech forest site
using Unscented Kalman Inversion (UKI).

Photosynthesis: PModel
Soil moisture stress: PiecewiseMoistureStressModel

12 calibrated parameters:
  9 canopy: moisture_stress_c, pmodel_cstar, pmodel_β, leaf_Cd,
            canopy_z_0m_coeff, canopy_z_0b_coeff, canopy_d_coeff,
            canopy_K_lw, canopy_emissivity
  3 soilco2 (DAMM): α_sx, kM_sx, kM_o2

Calibration window: summer only (Jun–Aug) to target the active growing season.
The full simulation (with 60-day spinup) still runs Jan–Dec.

Forcing: NetCDF meteorological data (PLUMBER2/CalLMIP format)
Observations: Daily aggregated NEE (gC/m²/d), Qle (W/m²), Qh (W/m²)
LAI: Copernicus LAI from the met NetCDF

Usage:
    julia --project=.buildkite experiments/integrated/fluxnet/calibrate_dk_sor.jl
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

import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using CairoMakie
CairoMakie.activate!()
using Statistics
using Logging
import Random
using NCDatasets
using TOML
using LinearAlgebra

# ── 2. Configuration ─────────────────────────────────────────────────────────
const FT = Float64
site_ID = "DK-Sor"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)
climaland_dir = pkgdir(ClimaLand)

# Data paths (local NetCDF files)
met_nc_path = joinpath(climaland_dir, "DK_Sor", "DK-Sor_1997-2014_FLUXNET2015_Met.nc")
flux_nc_path = joinpath(climaland_dir, "DK_Sor", "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")

spinup_days = 60
dt = Float64(450)  # 7.5 minutes

# Simulation period
sim_start_year = 2004
sim_end_year = 2005

# Site metadata
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

# Start with spinup before the calibration year
start_date = DateTime(sim_start_year, 1, 1) - Day(spinup_days)
stop_date = DateTime(sim_end_year + 1, 1, 1)  # end of sim_end_year
spinup_date = start_date + Day(spinup_days)

# Summer calibration windows (Jun 1 – Aug 31 of each year)
summer_ranges = [
    (Date(yr, 6, 1), Date(yr, 8, 31))
    for yr in sim_start_year:sim_end_year
]

println("Simulation: $start_date to $stop_date")
println("Spinup until: $spinup_date")
println("Calibration windows: $(["$s to $e" for (s,e) in summer_ranges])")

# UKI settings
rng_seed = 1234
rng = Random.MersenneTwister(rng_seed)
N_iterations = 10

# ── 3. One-time Setup (domain, forcing, LAI) ─────────────────────────────────
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

# Build forcing from NetCDF
toml_dict_base = LP.create_toml_dict(FT)
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_netcdf(
    met_nc_path,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict_base,
    FT,
)

# LAI from the met NetCDF (Copernicus LAI)
met_ds = NCDataset(met_nc_path, "r")
lai_data = Float64.(coalesce.(met_ds["LAI"][1, 1, :], NaN))
lai_times = met_ds["time"][:]
close(met_ds)

# Convert to seconds since start_date and build TVI
lai_seconds = [Float64(Second(t - Dates.Hour(time_offset) - start_date).value) for t in lai_times]
valid_lai = .!isnan.(lai_data)
LAI = TimeVaryingInput(lai_seconds[valid_lai], lai_data[valid_lai])

# ── 4. Load & Process Flux Observations (summer only) ────────────────────────
flux_ds = NCDataset(flux_nc_path, "r")
flux_times_dt = flux_ds["time"][:]

# Read flux variables (fill_value = 1e20 → missing → NaN)
nee_raw = Float64.(coalesce.(flux_ds["NEE_daily"][:], NaN))
qle_raw = Float64.(coalesce.(flux_ds["Qle_daily"][:], NaN))
qh_raw = Float64.(coalesce.(flux_ds["Qh_daily"][:], NaN))

# Uncertainty columns
nee_uc_raw = Float64.(coalesce.(flux_ds["NEE_uc_daily"][:], NaN))
qle_uc_raw = Float64.(coalesce.(flux_ds["Qle_uc_daily"][:], NaN))
qh_uc_raw = Float64.(coalesce.(flux_ds["Qh_uc_daily"][:], NaN))
close(flux_ds)

flux_dates = Date.(flux_times_dt)

# Filter: summer windows + valid data for all three fluxes
is_summer = map(flux_dates) do d
    any(s <= d <= e for (s, e) in summer_ranges)
end
valid_mask = is_summer .&
             .!isnan.(nee_raw) .& .!isnan.(qle_raw) .& .!isnan.(qh_raw) .&
             (abs.(nee_raw) .< 1e10) .& (abs.(qle_raw) .< 1e10) .& (abs.(qh_raw) .< 1e10)

obs_dates = flux_dates[valid_mask]
nee_obs = nee_raw[valid_mask]
qle_obs = qle_raw[valid_mask]
qh_obs = qh_raw[valid_mask]
nee_uc = nee_uc_raw[valid_mask]
qle_uc = qle_uc_raw[valid_mask]
qh_uc = qh_uc_raw[valid_mask]

n_obs = length(obs_dates)
println("Summer observation days: $n_obs (NEE+Qle+Qh)")

# Stack observations: [NEE_1,...,NEE_n, Qle_1,...,Qle_n, Qh_1,...,Qh_n]
y_obs = vcat(nee_obs, qle_obs, qh_obs)

# Noise covariance from uncertainty columns (diagonal)
nee_var = mean(filter(!isnan, nee_uc .^ 2))
qle_var = mean(filter(!isnan, qle_uc .^ 2))
qh_var = mean(filter(!isnan, qh_uc .^ 2))
noise_diag = vcat(
    fill(nee_var, n_obs),
    fill(qle_var, n_obs),
    fill(qh_var, n_obs),
)
noise_cov = Diagonal(noise_diag)
println("Noise variances - NEE: $(round(nee_var, sigdigits=3)) (gC/m²/d)², Qle: $(round(qle_var, sigdigits=3)) (W/m²)², Qh: $(round(qh_var, sigdigits=3)) (W/m²)²")

# ── 5. Model Function ────────────────────────────────────────────────────────
function run_model(params_vec)
    # Unpack parameters
    moisture_stress_c = params_vec[1]
    pmodel_cstar = params_vec[2]
    pmodel_beta = params_vec[3]
    leaf_Cd = params_vec[4]
    z_0m_coeff = params_vec[5]
    z_0b_coeff = params_vec[6]
    d_coeff = params_vec[7]
    K_lw = params_vec[8]
    emissivity = params_vec[9]
    alpha_sx = params_vec[10]
    kM_sx = params_vec[11]
    kM_o2 = params_vec[12]

    # Write TOML overrides
    toml_content = """
    [moisture_stress_c]
    value = $moisture_stress_c
    type = "float"

    [pmodel_cstar]
    value = $pmodel_cstar
    type = "float"

    ["pmodel_β"]
    value = $pmodel_beta
    type = "float"

    [leaf_Cd]
    value = $leaf_Cd
    type = "float"

    [canopy_z_0m_coeff]
    value = $z_0m_coeff
    type = "float"

    [canopy_z_0b_coeff]
    value = $z_0b_coeff
    type = "float"

    [canopy_d_coeff]
    value = $d_coeff
    type = "float"

    [canopy_K_lw]
    value = $K_lw
    type = "float"

    [canopy_emissivity]
    value = $emissivity
    type = "float"

    [soilCO2_pre_exponential_factor]
    value = $alpha_sx
    type = "float"

    [michaelis_constant]
    value = $kM_sx
    type = "float"

    [O2_michaelis_constant]
    value = $kM_o2
    type = "float"
    """
    tmp_toml_path = tempname() * ".toml"
    open(tmp_toml_path, "w") do io
        write(io, toml_content)
    end

    toml_dict = LP.create_toml_dict(FT; override_files = [tmp_toml_path])

    # Build canopy with PModel + PiecewiseMoistureStress
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

    forcing_nt = (;
        atmos = atmos,
        radiation = radiation,
        ground = ClimaLand.PrognosticGroundConditions{FT}(),
    )

    canopy = Canopy.CanopyModel{FT}(
        canopy_domain,
        forcing_nt,
        LAI,
        toml_dict;
        prognostic_land_components,
        photosynthesis = Canopy.PModel{FT}(canopy_domain, toml_dict),
        conductance = Canopy.PModelConductance{FT}(toml_dict),
        soil_moisture_stress = Canopy.PiecewiseMoistureStressModel{FT}(
            land_domain, toml_dict,
        ),
    )

    land = LandModel{FT}(
        (; atmos, radiation),
        LAI,
        toml_dict,
        land_domain,
        dt;
        prognostic_land_components,
        canopy,
    )

    function custom_set_ic!(Y, p, t, model)
        earth_param_set = ClimaLand.get_earth_param_set(model.soil)
        evaluate!(p.drivers.T, atmos.T, t)

        (; θ_r, ν, ρc_ds) = model.soil.parameters
        @. Y.soil.ϑ_l = θ_r + (ν - θ_r) / 2
        Y.soil.θ_i .= FT(0.0)
        ρc_s =
            ClimaLand.Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Y.soil.θ_i,
                ρc_ds,
                earth_param_set,
            )
        Y.soil.ρe_int .=
            ClimaLand.Soil.volumetric_internal_energy.(
                Y.soil.θ_i,
                ρc_s,
                p.drivers.T,
                earth_param_set,
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

        # SoilCO2 IC
        if !isnothing(model.soilco2)
            Y.soilco2.CO2 .= FT(0.000412)
            Y.soilco2.O2_f .= FT(0.21)
            # Prescribed SOC profile: 15 kgC/m³ at surface, 0.5 kgC/m³ at 1m depth
            SOC_top = FT(15.0)
            SOC_bot = FT(0.5)
            τ_soc = FT(1.0 / log(SOC_top / SOC_bot))
            z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
            @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
        end
    end

    # Diagnostics — use :short which includes nee, lhf, shf for LandModel
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    diags = ClimaLand.default_diagnostics(
        land,
        start_date;
        output_writer = output_writer,
        output_vars = :short,
        reduction_period = :daily,
    )

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

    rm(tmp_toml_path; force = true)
    return simulation
end

# ── 6. Forward Model G ───────────────────────────────────────────────────────
first_G_call = true
function G(params_vec)
    global first_G_call
    local simulation
    try
        simulation = run_model(params_vec)
    catch e
        println("  ERROR in run_model: ", sprint(showerror, e))
        return zeros(3 * n_obs)
    end

    # Extract daily diagnostics (verbose on first call for debugging)
    nee_result = extract_daily_diag(simulation, "nee_1d_average", obs_dates; verbose = first_G_call)
    qle_result = extract_daily_diag(simulation, "lhf_1d_average", obs_dates; verbose = first_G_call)
    qh_result = extract_daily_diag(simulation, "shf_1d_average", obs_dates; verbose = first_G_call)
    first_G_call = false

    # Convert NEE from model units (mol CO₂/m²/s) to gC/m²/d
    nee_gC = nee_result .* 12.0 .* 86400.0

    # Stack: [NEE, Qle, Qh]
    result = vcat(nee_gC, qle_result, qh_result)

    # Replace NaN/Inf with 0 to prevent UKI from crashing
    n_bad = count(x -> isnan(x) || isinf(x), result)
    if n_bad > 0
        println("  WARNING: $n_bad / $(length(result)) NaN/Inf values in G output, replacing with 0")
        replace!(x -> isnan(x) || isinf(x) ? 0.0 : x, result)
    end
    return result
end

function extract_daily_diag(simulation, diag_name, target_dates; verbose = false)
    # Search all diagnostic writers for the requested variable
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

    (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
        writer,
        diag_name,
    )

    model_dates_dt = times isa Vector{DateTime} ? times : date.(times)
    model_dates = Date.(model_dates_dt)

    if verbose
        println("  $diag_name: $(length(model_dates)) model dates, range $(first(model_dates)) to $(last(model_dates))")
        println("  Target dates: $(length(target_dates)), range $(first(target_dates)) to $(last(target_dates))")
        n_match = count(d -> d in Set(model_dates), target_dates)
        println("  Matching dates: $n_match / $(length(target_dates))")
    end

    model_dict = Dict{Date, Float64}()
    for (d, v) in zip(model_dates, data)
        model_dict[d] = Float64(v)
    end

    result = Float64[]
    for d in target_dates
        push!(result, get(model_dict, d, NaN))
    end
    return result
end

# ── 7. UKI Calibration ───────────────────────────────────────────────────────
println("\n=== Setting up UKI calibration (12 params, summer only) ===")

# 12 constrained Gaussian priors
priors = [
    # Canopy parameters
    PD.constrained_gaussian("moisture_stress_c", 0.5, 0.3, 0.01, 5.0),
    PD.constrained_gaussian("pmodel_cstar", 0.43, 0.15, 0.05, 2.0),
    PD.constrained_gaussian("pmodel_β", 51.0, 20.0, 5.0, 500.0),
    PD.constrained_gaussian("leaf_Cd", 0.1, 0.05, 0.005, 1.0),
    PD.constrained_gaussian("canopy_z_0m_coeff", 0.05, 0.03, 0.001, 0.3),
    PD.constrained_gaussian("canopy_z_0b_coeff", 0.001, 0.0005, 1e-5, 0.01),
    PD.constrained_gaussian("canopy_d_coeff", 0.1, 0.05, 0.001, 0.95),
    PD.constrained_gaussian("canopy_K_lw", 0.85, 0.25, 0.1, 2.0),
    PD.constrained_gaussian("canopy_emissivity", 0.97, 0.02, 0.9, 1.0),
    # DAMM soil CO₂ parameters
    PD.constrained_gaussian("α_sx", 25000.0, 10000.0, 1000.0, 200000.0),
    PD.constrained_gaussian("kM_sx", 0.01, 0.005, 1e-4, 0.1),
    PD.constrained_gaussian("kM_o2", 0.01, 0.005, 1e-4, 0.1),
]
prior = PD.combine_distributions(priors)

ensemble_kalman_process = EKP.EnsembleKalmanProcess(
    y_obs,
    noise_cov,
    EKP.Unscented(prior; α_reg = 1.0, update_freq = 1);
    rng,
)

# UKI loop
n_ens = EKP.get_N_ens(ensemble_kalman_process)
println("Ensemble size: $n_ens (for $(length(priors)) parameters)")

for i in 1:N_iterations
    println("\n--- UKI Iteration $i / $N_iterations ---")
    params_i = EKP.get_ϕ_final(prior, ensemble_kalman_process)
    println("  Parameter means: ", round.(mean(params_i, dims = 2)[:, 1], sigdigits = 4))

    G_ens = hcat(
        [G(params_i[:, j]) for j in 1:n_ens]...,
    )
    EKP.update_ensemble!(ensemble_kalman_process, G_ens)
    println("  Ensemble updated.")
end

# ── 8. Results ────────────────────────────────────────────────────────────────
final_params = EKP.get_ϕ_mean_final(prior, ensemble_kalman_process)
param_names = [
    "moisture_stress_c", "pmodel_cstar", "pmodel_β", "leaf_Cd",
    "canopy_z_0m_coeff", "canopy_z_0b_coeff", "canopy_d_coeff",
    "canopy_K_lw", "canopy_emissivity",
    "α_sx", "kM_sx", "kM_o2",
]
println("\n=== Calibration Complete ===")
println("Final parameter means:")
for (name, val) in zip(param_names, final_params)
    println("  $name = $(round(val, sigdigits = 5))")
end

# ── 9. Visualization ─────────────────────────────────────────────────────────
savedir = joinpath(@__DIR__, "calibrate_dk_sor_output")
mkpath(savedir)

# Parameter convergence + error plot
dim_size = sum(length.(EKP.batch(prior)))
fig1 = CairoMakie.Figure(size = ((dim_size + 1) * 300, 400))
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
CairoMakie.save(joinpath(savedir, "dk_sor_param_convergence.png"), fig1)
println("Saved: dk_sor_param_convergence.png")

# First vs last iteration: model vs obs for all 3 fluxes
first_G_ensemble = EKP.get_g(ensemble_kalman_process, 1)
last_iter = EKP.get_N_iterations(ensemble_kalman_process)
last_G_ensemble = EKP.get_g(ensemble_kalman_process, last_iter)

fig2 = CairoMakie.Figure(size = (1400, 1000))
flux_labels = ["NEE (gC/m²/d)", "Latent Heat (W/m²)", "Sensible Heat (W/m²)"]
for (fi, flux_label) in enumerate(flux_labels)
    idx_start = (fi - 1) * n_obs + 1
    idx_end = fi * n_obs
    obs_slice = y_obs[idx_start:idx_end]

    ax = Axis(
        fig2[fi, 1];
        title = "$flux_label: first vs last UKI iteration (summer only)",
        xlabel = "Day index (Jun–Aug)",
        ylabel = flux_label,
    )

    for j in 1:size(first_G_ensemble, 2)
        lines!(ax, 1:n_obs, first_G_ensemble[idx_start:idx_end, j];
               color = (:red, 0.3), linewidth = 0.5)
    end
    for j in 1:size(last_G_ensemble, 2)
        lines!(ax, 1:n_obs, last_G_ensemble[idx_start:idx_end, j];
               color = (:blue, 0.3), linewidth = 0.5)
    end
    lines!(ax, 1:n_obs, obs_slice;
           color = (:black, 0.8), linewidth = 2)

    if fi == 1
        axislegend(
            ax,
            [
                LineElement(color = :red, linewidth = 2),
                LineElement(color = :blue, linewidth = 2),
                LineElement(color = :black, linewidth = 3),
            ],
            ["First ensemble", "Last ensemble", "Observations"];
            position = :rt,
            framevisible = false,
        )
    end
end

CairoMakie.resize_to_layout!(fig2)
CairoMakie.save(joinpath(savedir, "dk_sor_G_first_and_last.png"), fig2)
println("Saved: dk_sor_G_first_and_last.png")
println("\nDone.")
