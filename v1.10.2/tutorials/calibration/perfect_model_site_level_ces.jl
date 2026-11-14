# # Perfect Model Site-Level CES Tutorial (Calibrate–Emulate–Sample)
#
# This tutorial extends the perfect model calibration experiment to the full
# [CalibrateEmulateSample.jl](https://github.com/CliMA/CalibrateEmulateSample.jl)
# (CES) workflow. After EKI calibration converges to a MAP estimate, we train a
# Gaussian Process emulator on the ensemble trajectory, then run Markov Chain
# Monte Carlo (MCMC) using the model emulator to obtain the full posterior
# uncertainty distribution.
#
# The CES workflow is motivated by the fact that land model runs are expensive.
# EKI calibration finds the MAP with O(ensemble size × iterations) model evaluations.
# Full Bayesian UQ via direct MCMC commonly requires in excess of O(10⁴) and up to
# O(10⁷) evaluations depending on the Markov chain mixing properties. CES solves this by:
#   1. **Calibrate**: EKI explores parameter space efficiently and finds an optimal set of parameters.
#   2. **Emulate**: Train a GP on the EKI trajectory in parameter space (one time expense; reusable).
#   3. **Sample**: MCMC carried out using the GP in place of the model — each evaluation costs microseconds.
#
# In this tutorial we calibrate `Vcmax25` against synthetic LHF at US-MOz.
#
# ## Overview
#
# The tutorial covers:
# 1. Setting up a land surface model for a FLUXNET site (US-MOz)
# 2. Creating a synthetic observation dataset
# 3. **Calibrate**: Ensemble Kalman Inversion (EKI)
# 4. **Emulate**: Gaussian Process surrogate trained on the EKI trajectory
# 5. **Sample**: RWMH-MCMC over the GP posterior
# 6. Visualizing the posterior vs prior vs true value
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
#
# We fix the random seed for reproducibility, choose single-precision floats, and
# load the default ClimaLand parameter set. The simulation targets the US-MOz
# FLUXNET site; its location and flux-tower height are looked up from the bundled
# site metadata.

rng_seed = 1234;
rng = Random.MersenneTwister(rng_seed);
const FT = Float32;

toml_dict = LP.create_toml_dict(FT);
site_ID = "US-MOz";
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID);

(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val));
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val));

(start_date, stop_date) =
    FluxnetSimulations.get_data_dates(site_ID, time_offset);
stop_date = DateTime(2010, 4, 1, 6, 30);
Δt = 450.0;

# ## Domain and Forcing Setup
#
# The model runs on a single soil column 2 m deep with 10 vertical elements,
# located at the site's longitude and latitude.

zmin = FT(-2);
zmax = FT(0);
domain = Column(; zlim = (zmin, zmax), nelements = 10, longlat = (long, lat));

# ## Model Setup

function model(Vcmax25)
    ## Forcing is created inside the function so each model call is self-contained.
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
    Vcmax25 = FT(Vcmax25)
    ground = ClimaLand.PrognosticGroundConditions{FT}()
    prognostic_land_components = (:canopy, :snow, :soil)
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    canopy_forcing = (; forcing.atmos, forcing.radiation, ground)
    photosyn_defaults =
        Canopy.clm_photosynthesis_parameters(canopy_domain.space.surface)
    photosynthesis = Canopy.FarquharModel{FT}(
        canopy_domain,
        toml_dict;
        photosynthesis_parameters = (;
            fractional_c3 = photosyn_defaults.fractional_c3,
            Vcmax25,
        ),
    )
    conductance = Canopy.MedlynConductanceModel{FT}(
        canopy_domain,
        toml_dict;
        g1 = FT(141),
    )
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        canopy_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        photosynthesis,
        prognostic_land_components,
        conductance,
    )
    land_model = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
        canopy,
    )
    set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        time_offset,
        land_model,
    )
    output_vars = ["shf", "lhf"]
    diagnostics = ClimaLand.default_diagnostics(
        land_model,
        start_date;
        output_writer = ClimaDiagnostics.Writers.DictWriter(),
        output_vars,
        reduction_period = :hourly,
    )
    simulation = Simulations.LandSimulation(
        start_date,
        stop_date,
        Δt,
        land_model;
        set_ic!,
        user_callbacks = (),
        diagnostics,
    )
    solve!(simulation)
    return simulation
end

function get_lhf(simulation)
    return ClimaLand.Diagnostics.diagnostic_as_vectors(
        simulation.diagnostics[1].output_writer,
        "lhf_1h_average",
    )
end

function get_diurnal_average(var, start_date, spinup_date)
    (times, data) = var
    model_dates = if times isa Vector{DateTime}
        times
    else
        Second.(getproperty.(times, :counter)) .+ start_date
    end
    spinup_idx = findfirst(spinup_date .<= model_dates)
    model_dates = model_dates[spinup_idx:end]
    data = data[spinup_idx:end]
    hour_of_day = Hour.(model_dates)
    mean_by_hour = [mean(data[hour_of_day .== Hour(i)]) for i in 0:23]
    return mean_by_hour
end

function G(Vcmax25)
    simulation = model(Vcmax25)
    lhf = get_lhf(simulation)
    observation =
        Float64.(
            get_diurnal_average(
                lhf,
                simulation.start_date,
                simulation.start_date + Day(20),
            ),
        )
    return observation
end

# ## Step 1 — Calibrate (Ensemble Kalman Inversion)
#
# Generate synthetic "truth" from the model at a known true parameter value,
# then run EKI to recover it from noisy observations.
#
# The true value is set to 1.5e-4 mol m⁻² s⁻¹ (150 μmol m⁻² s⁻¹) — one standard
# deviation above the prior mean of 1e-4. This is within the physically realistic
# range of `Vcmax25` for C3 vegetation, and a regime in which the forward map G
# (diurnal LHF) responds strongly to `Vcmax25`, so EKI can distinguish the true
# parameter from the prior center.

true_Vcmax25 = 1.5e-4;
observations = G(true_Vcmax25);

# The noise variance is used for both EKI and the GP emulator likelihood so that
# both steps see a consistent observation error model. Here
# `σ = √noise_var ≈ 5 W m⁻²` (about 7% of the ~70 W m⁻² synthetic peak LHF),
# large enough to produce visible ensemble spread in the plots while still
# allowing EKI to recover the true value.
noise_var = 25.0;
noise_covariance = noise_var * EKP.I;

# We place a constrained Gaussian prior on `Vcmax25` bounded to `[0, 2e-4]`, then
# draw the initial ensemble and run three EKI iterations.
prior = PD.constrained_gaussian("Vcmax25", 1e-4, 5e-5, 0, 2e-4);

ensemble_size = 10;
N_iterations = 3;

initial_ensemble = EKP.construct_initial_ensemble(rng, prior, ensemble_size);

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
);

Logging.with_logger(SimpleLogger(devnull, Logging.Error)) do
    for i in 1:N_iterations
        println("EKI iteration $i / $N_iterations")
        params_i = EKP.get_ϕ_final(prior, ensemble_kalman_process)
        G_ens = hcat([G(params_i[:, j]...) for j in 1:ensemble_size]...)
        EKP.update_ensemble!(ensemble_kalman_process, G_ens)
    end
end

eki_map = EKP.get_ϕ_mean_final(prior, ensemble_kalman_process)[1];
@info "EKI MAP estimate: Vcmax25 = $eki_map  (true: $true_Vcmax25)"

# We now plot two calibration diagnostics: the parameter convergence and error
# decay over iterations, and the forward-model ensemble compared against the
# synthetic observations.
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
CairoMakie.save("perfect_model_ces_params_and_error.png", fig);
# ![](perfect_model_ces_params_and_error.png)

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
lines!(ax2, 0:23, observations; color = :black, linewidth = 3)
for g in eachcol(first_G)
    lines!(ax2, 0:23, g; color = (:red, 0.4), linewidth = 1.5)
end
for g in eachcol(last_G)
    lines!(ax2, 0:23, g; color = (:blue, 0.4), linewidth = 1.5)
end
axislegend(
    ax2,
    [
        LineElement(; color = :red, linewidth = 2),
        LineElement(; color = :blue, linewidth = 2),
        LineElement(; color = :black, linewidth = 3),
    ],
    ["First ensemble", "Last ensemble", "Truth"];
    position = :rb,
    framevisible = false,
)
CairoMakie.resize_to_layout!(fig2)
CairoMakie.save("perfect_model_ces_G_first_last.png", fig2);
# ![](perfect_model_ces_G_first_last.png)

# ## Step 2 — Emulate (Gaussian Process Surrogate)
#
# The EKI calibration produced (parameter, G(parameter)) pairs for all ensemble
# members across all iterations. We now train a GP emulator on this dataset.
# The surrogate maps from parameter space (1D) to observation space (24D,
# diurnal LHF) and is orders of magnitude cheaper to evaluate than the land model.

n_iter_used = EKP.get_N_iterations(ensemble_kalman_process);
@info "Building GP emulator from $n_iter_used EKI iterations ($(n_iter_used * ensemble_size) training points)..."

# First we extract all `(θ_unconstrained, G(θ))` pairs from the EKP object; these
# are the training data for the emulator.
input_output_pairs =
    Utilities.get_training_points(ensemble_kalman_process, n_iter_used);

# Next we build, configure, and optimize the GP emulator. `GPJL` uses
# GaussianProcesses.jl (squared-exponential kernel by default). The decorrelator
# inside `Emulator` projects the 24-dimensional output onto its principal
# components so that each GP only needs to model a one-dimensional output.
# `encoder_kwargs` bundles the noise covariance and prior covariance from the
# calibration step — the recommended way to configure the `Emulator` encoder.
encoder_kwargs = Utilities.encoder_kwargs_from(ensemble_kalman_process, prior);
gppackage = GPJL();
gauss_proc = GaussianProcess(gppackage; noise_learn = false);
emulator = Emulator(gauss_proc, input_output_pairs; encoder_kwargs);
Logging.with_logger(SimpleLogger(devnull, Logging.Error)) do
    optimize_hyperparameters!(emulator)
end
@info "GP emulator trained."

# ### Emulator Validation
#
# Before running MCMC, verify that the GP accurately reproduces the model
# outputs on the EKI training data. Red curves are the actual G(θ) values;
# blue curves are the emulator's predictions at the same inputs. Close overlap
# confirms the surrogate is fit for purpose before sampling begins.
train_in = EKP.DataContainers.get_inputs(input_output_pairs);
train_out = EKP.DataContainers.get_outputs(input_output_pairs);
emul_pred_train, _ = Emulators.predict(emulator, train_in);
n_train = size(train_in, 2);
val_idx = rand(rng, 1:n_train, min(10, n_train));

fig_val = CairoMakie.Figure(; size = (900, 400))
ax_val = Axis(
    fig_val[1, 1];
    title = "Emulator validation: GP predictions vs EKI training data",
    xlabel = "Hour of day",
    ylabel = "LHF [W m⁻²]",
)
lines!(
    ax_val,
    0:23,
    observations;
    color = :black,
    linewidth = 3,
    label = "Truth",
)
for j in val_idx
    lines!(ax_val, 0:23, train_out[:, j]; color = (:red, 0.5), linewidth = 1.5)
end
for j in val_idx
    lines!(
        ax_val,
        0:23,
        emul_pred_train[:, j];
        color = (:steelblue, 0.5),
        linewidth = 1.5,
    )
end
axislegend(
    ax_val,
    [
        LineElement(; color = :black, linewidth = 3),
        LineElement(; color = :red, linewidth = 2),
        LineElement(; color = :steelblue, linewidth = 2),
    ],
    ["Truth", "Model G(θ)", "Emulator GP(θ)"];
    position = :rt,
    framevisible = false,
)
CairoMakie.resize_to_layout!(fig_val)
CairoMakie.save("perfect_model_ces_emulator_validation.png", fig_val);
# ![](perfect_model_ces_emulator_validation.png)

# ## Step 3 — Sample (MCMC Posterior)
#
# Run Random Walk Metropolis-Hastings MCMC over the trained GP emulator.
# Each likelihood evaluation calls the GP instead of the land model —
# enabling 50 000 samples in seconds rather than years.

# We initialize the MCMC chain at the EKI posterior mean (in unconstrained
# space) and wrap it in an `MCMCWrapper` tied to the trained emulator.
init_params = EKP.get_u_mean_final(ensemble_kalman_process);

mcmc = MCMCWrapper(
    RWMHSampling(),
    Float64.(observations),
    prior,
    emulator;
    init_params,
);

# We first tune the step size to target a ~0.2 acceptance rate, then draw the
# samples. `init_stepsize` is set near the expected value so the doubling search
# converges in few steps. Both calls are wrapped to suppress their progress bars,
# which otherwise clutter the rendered tutorial.
init_stepsize = 1.0;
n_samples = 50_000;
discard_initial = 2_000;
new_step = redirect_stderr(devnull) do
    Logging.with_logger(SimpleLogger(devnull, Logging.Error)) do
        optimize_stepsize(mcmc; init_stepsize, N = 2_000, discard_initial = 0)
    end
end;
@info "MCMC step size: init=$init_stepsize  learned=$new_step"
chain = redirect_stderr(devnull) do
    Logging.with_logger(SimpleLogger(devnull, Logging.Error)) do
        MarkovChainMonteCarlo.sample(
            mcmc,
            n_samples;
            stepsize = new_step,
            discard_initial,
        )
    end
end;
posterior = MarkovChainMonteCarlo.get_posterior(mcmc, chain);

# Finally we transform the samples from unconstrained to constrained (physical)
# space and summarize the posterior.
constrained_posterior = Emulators.transform_unconstrained_to_constrained(
    prior,
    MarkovChainMonteCarlo.get_distribution(posterior),
);

post_Vcmax25 = vec(constrained_posterior["Vcmax25"]);
post_mean = mean(post_Vcmax25);
post_std = std(post_Vcmax25);
ci_lo, ci_hi = quantile(post_Vcmax25, [0.025, 0.975]);

@info "Posterior Vcmax25: mean=$(round(post_mean; sigdigits=4))  std=$(round(post_std; sigdigits=3))"
@info "95% CI: [$(round(ci_lo; sigdigits=4)), $(round(ci_hi; sigdigits=4))]"
@info "EKI MAP: $(round(eki_map; sigdigits=4))  true: $true_Vcmax25"

# ## UQ Results: Posterior vs Prior vs True Value
#
# Here we show the key results from this tutorial. First, we have obtained a probability distribution for the parameter we calibrated, which gives us a quantified uncertainty in its value. Second, the posterior is much narrower than the prior, indicating that the data provided useful constraints on this parameter. Additionally, the posterior is centered near the true value, demonstrating successful uncertainty quantification.
# The EKI MAP gives a single point estimate; MCMC over the emulator gives
# the full posterior credible interval.

fig3 = CairoMakie.Figure(; size = (700, 450))
ax3 = Axis(
    fig3[1, 1];
    title = "Posterior vs Prior: Vcmax25 (perfect model)",
    xlabel = "Vcmax25 [mol m⁻² s⁻¹]",
    ylabel = "Density",
)

## Sample the prior for comparison.
rng_prior = Random.MersenneTwister(99)
prior_ens = EKP.construct_initial_ensemble(rng_prior, prior, 5_000)
prior_constrained = EKP.transform_unconstrained_to_constrained(prior, prior_ens)

density!(ax3, vec(prior_constrained); label = "Prior", color = (:grey, 0.4))
density!(
    ax3,
    post_Vcmax25;
    label = "Posterior (CES-MCMC)",
    color = (:steelblue, 0.6),
)
vlines!(
    ax3,
    [true_Vcmax25];
    color = :red,
    linewidth = 2.5,
    label = "True value",
)
vlines!(
    ax3,
    [eki_map];
    color = :orange,
    linewidth = 2.5,
    linestyle = :dash,
    label = "EKI MAP",
)
axislegend(ax3; position = :rt, framevisible = false)
CairoMakie.resize_to_layout!(fig3)
CairoMakie.save("perfect_model_ces_posterior.png", fig3);
# ![](perfect_model_ces_posterior.png)

# ## Posterior Predictive Check
#
# The Emulator Validation section above checked that the emulator reproduces the
# true model output (in an L² sense) on the EKI training points. That is necessary
# but not sufficient: a small output error can still bias the inferred posterior.
# The posterior predictive check closes this gap by pushing posterior samples back
# through the TRUE forward model `G` and asking whether the observations are
# plausible under them.
#
# Each evaluation of `G` is a full land-model run, so we do not recompute this
# during the docs build. The figure below was generated offline with 50 posterior
# and 50 prior samples: it shows the posterior predictive spread (blue) collapsing
# around the truth (black) relative to the much wider prior predictive (gray) —
# the constraint that calibration buys us. To regenerate it (optionally with more
# samples), run:
#
# ```julia
# n_pp = 50
# post_unc = MarkovChainMonteCarlo.get_distribution(posterior)["Vcmax25"]
# pp_idx = rand(rng, 1:size(post_unc, 2), n_pp)
# post_constrained =
#     EKP.transform_unconstrained_to_constrained(prior, post_unc[:, pp_idx])
# prior_unc = EKP.construct_initial_ensemble(rng, prior, n_pp)
# prior_constrained =
#     EKP.transform_unconstrained_to_constrained(prior, prior_unc)
# fig = CairoMakie.Figure(; size = (900, 400))
# ax = Axis(fig[1, 1]; xlabel = "Hour of day", ylabel = "LHF [W m⁻²]")
# for j in 1:n_pp  # prior predictive
#     lines!(ax, 0:23, G(prior_constrained[1, j]); color = (:gray, 0.2))
# end
# for j in 1:n_pp  # posterior predictive
#     lines!(ax, 0:23, G(post_constrained[1, j]); color = (:steelblue, 0.35))
# end
# lines!(ax, 0:23, observations; color = :black, linewidth = 3)
# ```

pp_img = "perfect_model_ces_predictive_offline.png";
pp_src = joinpath("..", "..", "tutorials", "calibration", pp_img);
isfile(pp_src) && cp(pp_src, pp_img; force = true);
# ![](perfect_model_ces_predictive_offline.png)
