# # Observation-Based Site-Level CES Tutorial (Calibrate–Emulate–Sample)
#
# This tutorial extends the observation-based site calibration experiment to
# the full [CalibrateEmulateSample.jl](https://github.com/CliMA/CalibrateEmulateSample.jl)
# (CES) workflow. We calibrate two parameters (`Vcmax25`, `g1`) against real
# FLUXNET LHF observations at US-MOz, then use the EKI trajectory to train a GP
# emulator, and finally run MCMC over the surrogate to get the joint posterior.
#
# Using real observations (rather than synthetic perfect-model data) means the
# posterior captures genuine model-data discrepancy. The full posterior from
# MCMC reveals parameter correlations, multi-modality, and uncertainty that a
# single EKI MAP estimate cannot show. Both parameters are calibrated within
# prior ranges where LHF is sensitive to each: Vcmax25 ∈ [0, 2e-3] keeps the
# ensemble in a regime where canopy transpiration is significant, and g1
# modulates stomatal aperture directly.
#
# ## Overview
#
# The tutorial covers:
# 1. Setting up a land surface model for FLUXNET site US-MOz
# 2. Obtaining FLUXNET LHF observations
# 3. **Calibrate**: EKI for Vcmax25 and g1
# 4. **Emulate**: GP surrogate trained on the EKI trajectory
# 5. **Sample**: RWMH-MCMC joint posterior over both parameters
# 6. Visualizing the 2D posterior, marginals, and posterior predictive
#
#nb # ## Prerequisites
#nb #
#nb # First, ensure you have the required packages installed:
#nb ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
#nb using Pkg
#nb required_pkgs = [
#nb     "ClimaLand",
#nb     "ClimaDiagnostics",
#nb     "CairoMakie",
#nb     "EnsembleKalmanProcesses",
#nb     "CalibrateEmulateSample",
#nb     "Random",
#nb     "Logging",
#nb ]
#nb Pkg.add(required_pkgs)
#nb plotting_pkgs = ["ClimaAnalysis", "GeoMakie", "Printf", "Statistics"]
#nb Pkg.add(plotting_pkgs)
#nb Pkg.instantiate()
#nb Pkg.precompile()

# ## Setup and Imports

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Canopy
using ClimaLand.Simulations
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Parameters as LP
using ClimaDiagnostics: ClimaDiagnostics
import ClimaUtilities.TimeManager: date
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using CalibrateEmulateSample.Emulators
using CalibrateEmulateSample.MarkovChainMonteCarlo
using CalibrateEmulateSample.Utilities
using LinearAlgebra: LinearAlgebra
using CairoMakie
CairoMakie.activate!()
using Statistics
using Logging
using Random: Random
using Dates

# ## Configuration and Site Setup

rng_seed = 1234
rng = Random.MersenneTwister(rng_seed)
const FT = Float32

toml_dict = LP.create_toml_dict(FT)
site_ID = "US-MOz"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID)

(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

# Two months of data (May–Jun 2010)
start_date = DateTime(2010, 5, 1)
stop_date = DateTime(2010, 7, 1)
Δt = 450.0

# ## Domain and Forcing Setup
#
# Domain, forcing, and LAI are defined once outside the forward model function
# because each EKI iteration only varies (Vcmax25, g1), not the forcing.

zmin = FT(-2)
zmax = FT(0)
domain = Column(; zlim = (zmin, zmax), nelements = 10, longlat = (long, lat))

# ## Model Setup

function model(Vcmax25, g1)
    Vcmax25 = FT(Vcmax25)
    g1 = FT(g1)
    forcing = FluxnetSimulations.prescribed_forcing_fluxnet(
        site_ID,
        lat,
        long,
        time_offset,
        atmos_h,
        start_date,
        toml_dict,
        FT,
    )
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        domain.space.surface,
        start_date,
        stop_date,
    )
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
    retention_parameters = (;
        ν = FT(0.55),
        θ_r = FT(0.04),
        K_sat = FT(4e-7),
        hydrology_cm = ClimaLand.Soil.vanGenuchten{FT}(;
            α = FT(0.05),
            n = FT(2.0),
        ),
    )
    composition_parameters =
        (; ν_ss_om = FT(0.1), ν_ss_quartz = FT(0.1), ν_ss_gravel = FT(0.0))
    soil = ClimaLand.Soil.EnergyHydrology{FT}(
        domain,
        forcing,
        toml_dict;
        prognostic_land_components,
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
        runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(),
        retention_parameters,
        composition_parameters,
        S_s = FT(1e-2),
    )
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    canopy_forcing = (;
        atmos = forcing.atmos,
        radiation = forcing.radiation,
        ground = ClimaLand.PrognosticGroundConditions{FT}(),
    )
    photosyn_defaults =
        Canopy.clm_photosynthesis_parameters(surface_domain.space.surface)
    photosynthesis = Canopy.FarquharModel{FT}(
        surface_domain,
        toml_dict;
        photosynthesis_parameters = (;
            fractional_c3 = photosyn_defaults.fractional_c3,
            Vcmax25,
        ),
    )
    conductance =
        Canopy.MedlynConductanceModel{FT}(surface_domain, toml_dict; g1)
    canopy = Canopy.CanopyModel{FT}(
        surface_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        photosynthesis,
        conductance,
        prognostic_land_components,
    )
    land_model = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
        soil,
        canopy,
    )
    set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        time_offset,
        land_model,
    )
    output_vars = ["lhf"]
    diagnostics = ClimaLand.default_diagnostics(
        land_model,
        start_date;
        output_writer = ClimaDiagnostics.Writers.DictWriter(),
        output_vars,
        reduction_period = :halfhourly,
    )
    simulation = Simulations.LandSimulation(
        start_date,
        stop_date,
        Δt,
        land_model;
        set_ic!,
        updateat = Second(Δt),
        user_callbacks = (),
        diagnostics,
    )
    solve!(simulation)
    return simulation
end

function get_lhf(simulation)
    return ClimaLand.Diagnostics.diagnostic_as_vectors(
        simulation.diagnostics[1].output_writer,
        "lhf_30m_average",
    )
end

function get_diurnal_average(var, start_date, spinup_date, stop_date)
    (times, data) = var
    model_dates = if times isa Vector{DateTime}
        times
    else
        date.(times)
    end
    spinup_idx = findfirst(spinup_date .<= model_dates)
    stop_idx = findlast(model_dates .< stop_date)
    model_dates = model_dates[spinup_idx:stop_idx]
    data = data[spinup_idx:stop_idx]
    hour_of_day = Hour.(model_dates)
    mean_by_hour = [mean(data[hour_of_day .== Hour(i)]) for i in 0:23]
    return mean_by_hour
end

function G(Vcmax25, g1)
    simulation = model(Vcmax25, g1)
    lhf = get_lhf(simulation)
    observation =
        Float64.(
            get_diurnal_average(
                lhf,
                simulation.start_date,
                simulation.start_date + Day(20),
                stop_date,
            ),
        )
    return observation
end

# ## Step 1 — Calibrate (Ensemble Kalman Inversion)
#
# Load real FLUXNET observations at US-MOz and run EKI to calibrate Vcmax25
# and g1. With real observations there is inherent model-data discrepancy:
# EKI finds the parameter combination that minimises the LHF residual over the
# two-month window.

dataset = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
observations = get_diurnal_average(
    (dataset.UTC_datetime, dataset.lhf),
    start_date,
    start_date + Day(20),
    stop_date,
)

# Observation noise: flat 5% variance covariance (24 × 24)
## noise_var is defined once and used for both EKI and the GP emulator likelihood.
noise_var = 0.05
noise_covariance = noise_var * EKP.I

# Priors: Vcmax25 ∈ [0, 2e-3], g1 ∈ [0, 1000]
#
# The Vcmax25 prior is centred at 1e-3 mol m⁻² s⁻¹ so that canopy transpiration
# contributes meaningfully to LHF and the forward map is sensitive to both
# parameters. A prior centred at 1e-4 mol m⁻² s⁻¹ would put the ensemble in a
# soil-evaporation-dominated regime where LHF is nearly flat in Vcmax25,
# making it unidentifiable from LHF observations alone.
priors = [
    PD.constrained_gaussian("Vcmax25", 1e-3, 5e-4, 0, 2e-3),
    PD.constrained_gaussian("g1", 150, 90, 0, 1000),
]
prior = PD.combine_distributions(priors)

ensemble_size = 10
N_iterations = 4

initial_ensemble = EKP.construct_initial_ensemble(rng, prior, ensemble_size)

ensemble_kalman_process = EKP.EnsembleKalmanProcess(
    initial_ensemble,
    observations,
    noise_covariance,
    EKP.Inversion();
    scheduler = EKP.DataMisfitController(;
        terminate_at = Inf,
        on_terminate = "continue",
    ),
    rng,
)

Logging.with_logger(SimpleLogger(devnull, Logging.Error)) do
    for i in 1:N_iterations
        println("EKI iteration $i / $N_iterations")
        params_i = EKP.get_ϕ_final(prior, ensemble_kalman_process)
        G_ens = hcat([G(params_i[:, j]...) for j in 1:ensemble_size]...)
        EKP.update_ensemble!(ensemble_kalman_process, G_ens)
    end
end

eki_map = EKP.get_ϕ_mean_final(prior, ensemble_kalman_process)
@info "EKI MAP: Vcmax25=$(round(eki_map[1]; sigdigits=4))  g1=$(round(eki_map[2]; sigdigits=4))"

# Calibration diagnostic plots
dim_size = sum(length.(EKP.batch(prior)))
fig = CairoMakie.Figure(; size = ((dim_size + 1) * 500, 500))
for i in 1:dim_size
    EKP.Visualize.plot_ϕ_over_iters(
        fig[1, i],
        ensemble_kalman_process,
        prior,
        i,
    )
end
EKP.Visualize.plot_error_over_iters(
    fig[1, dim_size + 1],
    ensemble_kalman_process,
)
CairoMakie.save("obs_ces_params_and_error.png", fig)
# ![](obs_ces_params_and_error.png)

fig2 = CairoMakie.Figure(; size = (900, 400))
first_G = EKP.get_g(ensemble_kalman_process, 1)
last_iter = EKP.get_N_iterations(ensemble_kalman_process)
last_G = EKP.get_g(ensemble_kalman_process, last_iter)
ax2 = Axis(
    fig2[1, 1];
    title = "G ensemble: first vs last EKI iteration",
    xlabel = "Hour of day",
    ylabel = "LHF [W m⁻²]",
)
for g in eachcol(first_G)
    lines!(ax2, 0:23, g; color = (:red, 0.4), linewidth = 1.5)
end
for g in eachcol(last_G)
    lines!(ax2, 0:23, g; color = (:blue, 0.4), linewidth = 1.5)
end
lines!(ax2, 0:23, observations; color = :black, linewidth = 3)
axislegend(
    ax2,
    [
        LineElement(; color = :red, linewidth = 2),
        LineElement(; color = :blue, linewidth = 2),
        LineElement(; color = :black, linewidth = 3),
    ],
    ["First ensemble", "Last ensemble", "FLUXNET observations"];
    position = :rb,
    framevisible = false,
)
CairoMakie.resize_to_layout!(fig2)
CairoMakie.save("obs_ces_G_first_last.png", fig2)
# ![](obs_ces_G_first_last.png)

# ## Step 2 — Emulate (Gaussian Process Surrogate)
#
# Train a GP emulator on all (θ, G(θ)) pairs collected across EKI iterations.
# With 2 parameters and 24 observation dimensions, CES decorrelates the output
# space (via PCA) before fitting one GP per principal component.

n_iter_used = EKP.get_N_iterations(ensemble_kalman_process)
@info "Building GP emulator from $n_iter_used iterations ($(n_iter_used * ensemble_size) training points)..."

input_output_pairs =
    Utilities.get_training_points(ensemble_kalman_process, n_iter_used)

n_obs = length(observations)
## Same noise_var used in EKI ensures a consistent observation error model.
obs_noise_cov = noise_var * Matrix(LinearAlgebra.I, n_obs, n_obs)

gppackage = GPJL()
gauss_proc = GaussianProcess(gppackage; noise_learn = false)
emulator = Emulator(gauss_proc, input_output_pairs; obs_noise_cov)
optimize_hyperparameters!(emulator)
@info "GP emulator trained."

# ## Step 3 — Sample (MCMC Joint Posterior)
#
# RWMH-MCMC over the trained GP surrogate produces samples from the joint
# posterior p(Vcmax25, g1 | LHF_FLUXNET). The joint posterior can reveal
# parameter correlations that are invisible in the marginals.

init_params = EKP.get_u_mean_final(ensemble_kalman_process)

mcmc = MCMCWrapper(
    RWMHSampling(),
    Float64.(observations),
    prior,
    emulator;
    init_params,
)

new_step =
    optimize_stepsize(mcmc; init_stepsize = 0.1, N = 2_000, discard_initial = 0)
@info "MCMC step size: $new_step"

n_samples = 50_000
discard_initial = 2_000
chain = MarkovChainMonteCarlo.sample(
    mcmc,
    n_samples;
    stepsize = new_step,
    discard_initial,
)
posterior = MarkovChainMonteCarlo.get_posterior(mcmc, chain)

constrained_posterior = Emulators.transform_unconstrained_to_constrained(
    prior,
    MarkovChainMonteCarlo.get_distribution(posterior),
)

post_Vcmax25 = vec(constrained_posterior["Vcmax25"])
post_g1 = vec(constrained_posterior["g1"])

@info "Posterior Vcmax25: mean=$(round(mean(post_Vcmax25); sigdigits=4))  std=$(round(std(post_Vcmax25); sigdigits=3))"
@info "Posterior g1:      mean=$(round(mean(post_g1); sigdigits=4))  std=$(round(std(post_g1); sigdigits=3))"
@info "95% CI Vcmax25: $(round.(quantile(post_Vcmax25, [0.025, 0.975]); sigdigits=4))"
@info "95% CI g1:      $(round.(quantile(post_g1, [0.025, 0.975]); sigdigits=4))"

# ## UQ Results: Marginal Posteriors vs Priors
#
# Each panel shows the 1D marginal posterior for one parameter, overlaid with
# prior samples and the EKI MAP estimate. The posteriors are narrower than the
# priors, showing that the FLUXNET observations are informative. Note that with
# only two months of data, the posteriors may still be broad — this is physically
# correct and precisely what UQ is meant to reveal.

rng_prior = Random.MersenneTwister(99)
prior_ens = EKP.construct_initial_ensemble(rng_prior, prior, 5_000)
prior_constrained = EKP.transform_unconstrained_to_constrained(prior, prior_ens)

fig3 = CairoMakie.Figure(; size = (1100, 450))

ax_v = Axis(
    fig3[1, 1];
    title = "Marginal posterior: Vcmax25",
    xlabel = "Vcmax25 [mol m⁻² s⁻¹]",
    ylabel = "Density",
)
density!(ax_v, prior_constrained[1, :]; label = "Prior", color = (:grey, 0.4))
density!(
    ax_v,
    post_Vcmax25;
    label = "Posterior (CES-MCMC)",
    color = (:steelblue, 0.6),
)
vlines!(
    ax_v,
    [eki_map[1]];
    color = :orange,
    linewidth = 2.5,
    linestyle = :dash,
    label = "EKI MAP",
)
axislegend(ax_v; position = :rt, framevisible = false)

ax_g = Axis(
    fig3[1, 2];
    title = "Marginal posterior: g1",
    xlabel = "g1 [Pa^0.5]",
    ylabel = "Density",
)
density!(ax_g, prior_constrained[2, :]; label = "Prior", color = (:grey, 0.4))
density!(ax_g, post_g1; label = "Posterior (CES-MCMC)", color = (:coral, 0.6))
vlines!(
    ax_g,
    [eki_map[2]];
    color = :orange,
    linewidth = 2.5,
    linestyle = :dash,
    label = "EKI MAP",
)
axislegend(ax_g; position = :rt, framevisible = false)

CairoMakie.resize_to_layout!(fig3)
CairoMakie.save("obs_ces_marginal_posteriors.png", fig3)
# ![](obs_ces_marginal_posteriors.png)

# ## Joint Posterior (Vcmax25 vs g1)
#
# The scatter plot of joint posterior samples reveals any correlation between
# the two parameters. A negative correlation would mean that higher Vcmax25
# can be compensated by lower g1 and still produce similar LHF — a fundamental
# equifinality issue in land surface modelling that EKI cannot show but UQ can.

fig4 = CairoMakie.Figure(; size = (600, 500))
ax4 = Axis(
    fig4[1, 1];
    title = "Joint posterior: Vcmax25 vs g1",
    xlabel = "Vcmax25 [mol m⁻² s⁻¹]",
    ylabel = "g1 [Pa^0.5]",
)
scatter!(
    ax4,
    post_Vcmax25[1:5:end],
    post_g1[1:5:end];
    color = (:steelblue, 0.15),
    markersize = 4,
    label = "Posterior samples",
)
scatter!(
    ax4,
    [eki_map[1]],
    [eki_map[2]];
    color = :orange,
    markersize = 14,
    marker = :star5,
    label = "EKI MAP",
)
axislegend(ax4; position = :rt, framevisible = false)
CairoMakie.resize_to_layout!(fig4)
CairoMakie.save("obs_ces_joint_posterior.png", fig4)
# ![](obs_ces_joint_posterior.png)

# ## Posterior Predictive Check
#
# We evaluate the GP *emulator* (not the land model) at random posterior samples
# and compare the predicted diurnal LHF against FLUXNET. Each call to
# `Emulators.predict` takes microseconds, making this check cheap. The spread of
# curves directly visualises the observation-constrained uncertainty in the model
# output — this is the key advantage of CES over a single EKI MAP estimate.

n_pp = 50
## unconstrained MCMC samples for the emulator (Dict keyed by parameter name)
post_params_unc_dict = MarkovChainMonteCarlo.get_distribution(posterior)
param_names_vec = EKP.get_name(prior)
unc_samples = vcat([post_params_unc_dict[nm] for nm in param_names_vec]...)  # (2 × n_samples)
pp_indices = rand(rng, 1:size(unc_samples, 2), n_pp)

fig5 = CairoMakie.Figure(; size = (900, 400))
ax5 = Axis(
    fig5[1, 1];
    title = "GP emulator posterior predictive vs FLUXNET",
    xlabel = "Hour of day",
    ylabel = "LHF [W m⁻²]",
)
for idx in pp_indices
    θ = unc_samples[:, idx]
    pred_mean, _ = Emulators.predict(emulator, reshape(θ, :, 1))
    lines!(
        ax5,
        0:23,
        vec(pred_mean);
        color = (:steelblue, 0.2),
        linewidth = 1.0,
    )
end
lines!(
    ax5,
    0:23,
    observations;
    color = :black,
    linewidth = 3,
    label = "FLUXNET observations",
)
axislegend(ax5; position = :rb, framevisible = false)
CairoMakie.resize_to_layout!(fig5)
CairoMakie.save("obs_ces_predictive.png", fig5)
# ![](obs_ces_predictive.png)
#
# The posterior predictive spread quantifies how much uncertainty in Vcmax25 and
# g1 translates to uncertainty in predicted LHF — the ultimate output-space UQ.
