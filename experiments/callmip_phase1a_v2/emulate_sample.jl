"""
CalibrateEmulateSample (CES v1) uncertainty quantification for DK-Sor Phase 1a.

FIX 2 (minibatch / dimension mismatch):
  EKI was run with an ObservationSeries + RandomFixedSizeMinibatcher.  Each EKI
  iteration produced a different minibatch, so the observation vector length
  changes between iterations.  This makes the (u, G) pairs from different
  iterations incompatible as GP training data — the output dimension varies.

  SOLUTION: after EKI converges, we re-run ALL N_ens final-iteration members on
  a single FIXED window covering the full calibration record (1997–2012, ~5844
  days).  This yields a consistent (33-parameter, 3×5844)-output training matrix
  that the GP can learn.  We then load this fixed-window output as the GP input.

  The extra compute cost is ~N_ens × 1 forward-model call (same as 1 EKI iteration).

FIX 4 (deprecated CES API):
  Uses the current CES v1 Emulator constructor with encoder_schedule, not the
  deprecated obs_noise_cov / normalize_inputs / decorrelate kwargs.

Workflow:
  1. Load EKI output (final iteration parameters)
  2. Re-run all final-iteration members on full 1997-2012 window (if not cached)
  3. Load the fixed-window output matrix as GP training data
  4. Hold out 5 members for validation (must pass before MCMC)
  5. Train GP emulator
  6. Validate emulator: scatter plot + per-dim RMSE
  7. Run pCN-MCMC (100,000 steps)
  8. Save posterior, plot prior vs posterior marginals

Usage:
    julia --project=experiments/callmip_phase1a_v2 \\
          experiments/callmip_phase1a_v2/emulate_sample.jl
"""

using JLD2, TOML
using Plots
using Distributions, LinearAlgebra, Statistics, Random

using CalibrateEmulateSample
using CalibrateEmulateSample.Utilities
using CalibrateEmulateSample.ParameterDistributions
using CalibrateEmulateSample.Emulators
using CalibrateEmulateSample.MarkovChainMonteCarlo
using CalibrateEmulateSample.DataContainers
import CalibrateEmulateSample.EnsembleKalmanProcesses as EKP
import CalibrateEmulateSample.EnsembleKalmanProcesses.ParameterDistributions as PD

import ClimaCalibrate
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Canopy
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities
using Dates
using NCDatasets
using Statistics

rng = Random.MersenneTwister(60732)

# ── Settings ─────────────────────────────────────────────────────────────────
retain_var_in  = 0.99    # fraction of input variance retained after PCA
retain_var_out = 0.95    # fraction of output variance retained after SVD
chain_length   = 100_000  # MCMC chain length
init_stepsize  = 0.1
N_holdout      = 5        # members held out from GP training for validation

# ── Paths ─────────────────────────────────────────────────────────────────────
const climaland_dir  = abspath(joinpath(@__DIR__, "..", ".."))
const exp_dir        = @__DIR__
const cal_output_dir = joinpath(exp_dir, "output_eki")
const obs_filepath   = joinpath(exp_dir, "observations.jld2")
const fixedwin_dir   = joinpath(exp_dir, "output_fixed_window")   # re-run output
const ces_out_dir    = joinpath(exp_dir, "output_ces")
isdir(ces_out_dir)    || mkpath(ces_out_dir)
isdir(fixedwin_dir)   || mkpath(fixedwin_dir)

# ── Priors ────────────────────────────────────────────────────────────────────
include(joinpath(exp_dir, "priors.jl"))
prior, priors_vec = build_dk_sor_priors()
param_names = [only(PD.get_name(d)) for d in priors_vec]
@info "$(length(param_names))-parameter prior loaded"

# ── Load final-iteration EKP ──────────────────────────────────────────────────
function find_latest_ekp_path(cal_out_dir)
    iter = -1
    while isfile(joinpath(cal_out_dir,
                          "iteration_$(lpad(iter+1, 3, '0'))", "eki_file.jld2"))
        iter += 1
    end
    iter == -1 && error("No eki_file.jld2 found in $cal_out_dir")
    path = joinpath(cal_out_dir, "iteration_$(lpad(iter, 3, '0'))", "eki_file.jld2")
    @info "Loading EKP from iteration $iter: $path"
    return path, iter
end

ekp_path_str, final_iter = find_latest_ekp_path(cal_output_dir)
ekp = JLD2.load_object(ekp_path_str)
N_ens = EKP.get_N_ens(ekp)
@info "EKP loaded: ensemble size = $N_ens, iterations = $(length(EKP.get_u(ekp))-1)"

# ── Re-run final iteration on full 1997–2012 fixed window ─────────────────────
#
# WHY: The EKI used minibatching, so each iteration's G matrix has different
# column structure (different window).  The GP needs consistent (u, G) pairs
# across all training points.  We re-run the final-iteration parameter ensemble
# on a single FIXED window (full 1997-2012) to produce a consistent G matrix.
#
# COST: N_ens × 1 forward-model evaluation (~same as 1 EKI iteration).

const DT   = Float64(900)
const met_nc_path = joinpath(climaland_dir, "DK_Sor",
                             "DK-Sor_1997-2014_FLUXNET2015_Met.nc")

include(joinpath(exp_dir, "model_interface.jl"))

"""Run all final-iteration members on the full 1997-2012 window, save G matrix."""
function run_fixed_window_ensemble(ekp, final_iter)
    g_cache = joinpath(fixedwin_dir, "G_fixed_window.jld2")
    if isfile(g_cache)
        @info "Loading cached fixed-window G matrix from $g_cache"
        return JLD2.load(g_cache, "G_matrix"), JLD2.load(g_cache, "u_matrix")
    end

    @info "Re-running $N_ens final-iteration members on full 1997-2012 window…"

    # Final-iteration parameter ensemble (unconstrained space)
    # get_u(ekp, i) returns the parameter ensemble that was used to produce G_i
    u_final = EKP.get_u(ekp, final_iter)   # (n_params × N_ens)

    # Load obs to get window_dates for a single year to determine output length
    obs_data = JLD2.load(obs_filepath)

    # Full window: 1997-01-01 → 2012-12-31 (excluding leap days for consistency)
    full_dates_raw = collect(Date(1997,1,1):Day(1):Date(2012,12,31))
    full_dates = [d for d in full_dates_raw if !(month(d)==2 && day(d)==29)]
    n_days = length(full_dates)
    g_len  = 3 * n_days   # [NEE..., LE..., H...]

    G_matrix = zeros(g_len, N_ens)

    # We need to run each member with parameters from u_final, col m
    # Convert unconstrained → constrained → TOML → model
    for m in 1:N_ens
        @info "  Fixed-window run: member $m / $N_ens"
        # constrained params
        u_col    = u_final[:, m]
        ϕ_col    = EKP.ParameterDistributions.transform_unconstrained_to_constrained(
            prior, u_col,
        )

        # Write TOML
        member_dir = joinpath(fixedwin_dir, "member_$(lpad(m, 3, '0'))")
        isdir(member_dir) || mkpath(member_dir)
        param_toml = joinpath(member_dir, "parameters.toml")
        open(param_toml, "w") do io
            for (name, val) in zip(param_names, ϕ_col)
                used_in = name in ("soilCO2_reference_rate", "michaelis_constant",
                                   "O2_michaelis_constant", "soilCO2_activation_energy") ?
                          "[\"Land\"]" : "[\"getindex\"]"
                println(io, "[\"$name\"]"); println(io, "value = $(Float64(val))")
                println(io, "type  = \"float\""); println(io, "used_in = $used_in")
                println(io)
            end
        end

        toml_dict = LP.create_toml_dict(Float64; override_files = [param_toml])

        # Spinup: 2 years before 1997 → start 1995-01-01
        sim_start = DateTime(1995, 1, 1)
        sim_stop  = DateTime(2013, 1, 1)

        try
            land, forcing, ν, θ_r = build_dk_sor_model(
                Float64, sim_start, sim_stop, toml_dict, met_nc_path,
            )
            set_ic! = make_dk_sor_ic(forcing.atmos, ν, θ_r)

            output_writer = ClimaDiagnostics.Writers.DictWriter()
            diags = ClimaLand.default_diagnostics(
                land, sim_start;
                output_writer,
                output_vars      = :short,
                reduction_period = :daily,
            )
            sim = LandSimulation(
                sim_start, sim_stop, DT, land;
                set_ic!, updateat = Second(DT), user_callbacks = (), diagnostics = diags,
            )
            solve!(sim)

            nee_dates, nee_vals = extract_daily_diag(sim, "nee_1d_average")
            _,         qle_vals = extract_daily_diag(sim, "lhf_1d_average")
            _,         qh_vals  = extract_daily_diag(sim, "shf_1d_average")

            nd = Dict(zip(nee_dates, nee_vals))
            qd = Dict(zip(nee_dates, qle_vals))
            hd = Dict(zip(nee_dates, qh_vals))

            nee_out = [get(nd, d, 0.0) * 12.0 * 86400.0 for d in full_dates]
            qle_out = [get(qd, d, 0.0) for d in full_dates]
            qh_out  = [get(hd, d, 0.0) for d in full_dates]

            G_matrix[1:n_days,             m] = nee_out
            G_matrix[(n_days+1):(2*n_days), m] = qle_out
            G_matrix[(2*n_days+1):(3*n_days), m] = qh_out
        catch e
            @error "Member $m fixed-window run failed" exception = e
            G_matrix[:, m] .= NaN
        end
    end

    JLD2.jldsave(g_cache; G_matrix, u_matrix = u_final)
    @info "Fixed-window G matrix saved → $g_cache"
    return G_matrix, u_final
end

G_fixed, u_fixed = run_fixed_window_ensemble(ekp, final_iter)

# ── Build GP training data ────────────────────────────────────────────────────
# Hold out N_holdout members for validation
all_members  = 1:N_ens
test_members = collect(1:N_holdout)              # first N members as holdout
train_members = collect((N_holdout+1):N_ens)

train_pairs = PairedDataContainer(
    u_fixed[:, train_members],
    G_fixed[:, train_members],
)
test_pairs  = PairedDataContainer(
    u_fixed[:, test_members],
    G_fixed[:, test_members],
)

# ── Prune NaN / extreme columns ───────────────────────────────────────────────
function prune_bad_columns(pairs::PairedDataContainer; threshold::Real = 1_000)
    out      = get_outputs(pairs)
    bad_nan  = findall(vec(sum(isnan.(out); dims=1)) .> 0)
    bad_ext  = findall(vec(maximum(abs.(out); dims=1)) .> threshold)
    bad      = union(bad_nan, bad_ext)
    isempty(bad_nan) || @warn "Pruning $(length(bad_nan)) NaN columns"
    isempty(bad_ext) || @warn "Pruning $(length(bad_ext)) extreme columns"
    keep = setdiff(1:size(out,2), bad)
    PairedDataContainer(get_inputs(pairs)[:, keep], out[:, keep])
end

train_pairs = prune_bad_columns(train_pairs)
test_pairs  = prune_bad_columns(test_pairs)

@info "After pruning: $(size(get_inputs(train_pairs),2)) train / " *
      "$(size(get_inputs(test_pairs),2)) test members"

# ── Load full-record obs noise covariance for MCMC ────────────────────────────
obs_data   = JLD2.load(obs_filepath)
noise_full = obs_data["noise_full"]   # Diagonal covariance (1997-2012)
Γ = Matrix(noise_full)

# ── Build GP emulator (CES v1 encoder_schedule API) ──────────────────────────
emulator_path = joinpath(ces_out_dir, "emulator_fixed_window.jld2")

if isfile(emulator_path)
    @info "Loading cached emulator from $emulator_path"
    emulator = JLD2.load(emulator_path, "emulator")
else
    @info "Training GP emulator (train size = $(size(get_inputs(train_pairs),2)))…"

    gppackage = GPJL()
    mlt       = GaussianProcess(gppackage; noise_learn = false)

    # CES v1 API: encoder_schedule replaces deprecated kwargs
    #   decorrelate_sample_cov   : PCA on EKP input sample covariance
    #   decorrelate_structure_mat: SVD on observation noise structure matrix
    encoder_schedule = [
        (decorrelate_sample_cov(retain_var = retain_var_in),     "in"),
        (decorrelate_structure_mat(retain_var = retain_var_out), "out"),
    ]
    emulator = Emulator(
        mlt,
        train_pairs;
        encoder_schedule = deepcopy(encoder_schedule),
        encoder_kwargs   = (; obs_noise_cov = Γ),
        verbose          = true,
    )
    optimize_hyperparameters!(emulator)
    JLD2.save(emulator_path, "emulator", emulator)
    @info "Emulator saved → $emulator_path"
end

# ── Emulator validation (Ollie's requirement — non-negotiable) ────────────────
y_mean_train, _ = Emulators.predict(emulator, get_inputs(train_pairs);
                                    transform_to_real = true)
y_mean_test,  _ = Emulators.predict(emulator, get_inputs(test_pairs);
                                    transform_to_real = true)

y_true_train = get_outputs(train_pairs)
y_true_test  = get_outputs(test_pairs)

mse_train = mean((y_true_train .- y_mean_train) .^ 2)
mse_test  = mean((y_true_test  .- y_mean_test)  .^ 2)

# Correlation coefficient on held-out set
r_train = cor(vec(y_true_train), vec(y_mean_train))
r_test  = cor(vec(y_true_test),  vec(y_mean_test))

println("\n╔══════════════════════════════════════════╗")
println("║        Emulator validation               ║")
println("╠══════════════════════════════════════════╣")
println("║  GP MSE (train): $(lpad(round(mse_train, sigdigits=4), 16))      ║")
println("║  GP MSE (test) : $(lpad(round(mse_test,  sigdigits=4), 16))      ║")
println("║  GP r   (train): $(lpad(round(r_train,   sigdigits=4), 16))      ║")
println("║  GP r   (test) : $(lpad(round(r_test,    sigdigits=4), 16))      ║")
println("╚══════════════════════════════════════════╝")

if r_test < 0.6
    @warn "Emulator test correlation < 0.6 (r = $(round(r_test; sigdigits=3))). " *
          "Consider adding more EKI iterations or checking for model blow-ups. " *
          "Proceeding with MCMC but posterior quality may be limited."
end

# Scatter plot: GP predicted vs actual (latent output space)
p_scatter = scatter(
    vec(y_true_train), vec(y_mean_train);
    label = "train", alpha = 0.3, markersize = 2, markerstrokewidth = 0,
    xlabel = "True (latent output)", ylabel = "GP predicted (latent output)",
    title = "Emulator: predicted vs actual",
)
scatter!(vec(y_true_test), vec(y_mean_test);
         label = "test", alpha = 0.5, markersize = 3, markerstrokewidth = 0)
lo = min(minimum(y_true_train), minimum(y_true_test))
hi = max(maximum(y_true_train), maximum(y_true_test))
plot!(p_scatter, [lo, hi], [lo, hi]; lc = :black, ls = :dash, label = "1:1")
savefig(p_scatter, joinpath(ces_out_dir, "emulator_validation_scatter.png"))
@info "Emulator scatter → $(joinpath(ces_out_dir, "emulator_validation_scatter.png"))"

# ── MCMC posterior sampling ───────────────────────────────────────────────────
# NOTE: MCMCWrapper accepts the full observation series (ObservationSeries) or
# a raw y_obs vector + noise covariance.  We use the full-record fixed-window
# observation vector (y_obs_full) directly, constructing a synthetic single-
# observation EKP for consistency.
y_obs_full = obs_data["y_obs_full"]

# Build a minimal ObservationSeries covering the full record for MCMC
obs_for_mcmc = EKP.Observation(
    Dict(
        "samples"     => y_obs_full,
        "covariances" => Γ,
        "names"       => "full_1997_2012",
    ),
)

# init_params: EKI optimal (final iteration mean)
u_final_mean = vec(mean(u_fixed; dims=2))

mcmc = MCMCWrapper(
    pCNMHSampling(),
    obs_for_mcmc,
    prior,
    emulator;
    init_params = u_final_mean,
)

@info "Optimising MCMC step size…"
new_step = try
    optimize_stepsize(rng, mcmc; init_stepsize, N = 2000, discard_initial = 0)
catch e
    fallback = init_stepsize / 10
    @warn "optimize_stepsize failed ($e). Using fallback = $fallback"
    fallback
end
@info "Step size: $new_step"

chain = MarkovChainMonteCarlo.sample(
    rng, mcmc, chain_length;
    stepsize = new_step, discard_initial = 2_000,
)
display(chain)

posterior = MarkovChainMonteCarlo.get_posterior(mcmc, chain)

# ── Back-transform ────────────────────────────────────────────────────────────
posterior_samples = reduce(vcat,
    [get_distribution(posterior)[name] for name in get_name(posterior)])
constrained_posterior   = PD.transform_unconstrained_to_constrained(
    prior, posterior_samples,
)
constrained_ekp_optimal = PD.transform_unconstrained_to_constrained(
    prior, u_final_mean,
)

# ── Save ──────────────────────────────────────────────────────────────────────
posterior_path = joinpath(ces_out_dir, "posterior_fixed_window.jld2")
JLD2.save(
    posterior_path,
    "posterior",               posterior,
    "constrained_posterior",   constrained_posterior,
    "constrained_ekp_optimal", constrained_ekp_optimal,
    "param_names",             param_names,
)
@info "Posterior saved → $posterior_path"

# ── Summary ───────────────────────────────────────────────────────────────────
post_mean = vec(mean(constrained_posterior; dims=2))
post_std  = vec(std(constrained_posterior; dims=2))

println("\n╔══════════════════════════════════════════════════════════════════╗")
println("║   Posterior summary (constrained / physical space)               ║")
println("╠══════════════════════════════╦═══════════════╦══════════════════╣")
println("║ parameter                    ║  EKI optimal  ║  post μ ± σ      ║")
println("╠══════════════════════════════╬═══════════════╬══════════════════╣")
for (i, name) in enumerate(param_names)
    println("║ $(rpad(name, 29))║ $(lpad(round(constrained_ekp_optimal[i]; sigdigits=4), 13))  ║ " *
            "$(rpad("$(round(post_mean[i]; sigdigits=4)) ± $(round(post_std[i]; sigdigits=3))", 16))║")
end
println("╚══════════════════════════════╩═══════════════╩══════════════════╝")

# Flag pmodel_β collapse (known issue from old PR)
β_idx = findfirst(==("pmodel_β_c3"), param_names)
if !isnothing(β_idx) && post_mean[β_idx] < 30.0
    @warn "pmodel_β_c3 posterior mean ($(round(post_mean[β_idx]; sigdigits=3))) is below 30." *
          " This may indicate prior misspecification in this high-dimensional problem." *
          " Per Ollie: not a failure — use EKI-optimal for CalLMIP submission, document here."
end

# ── Plot prior vs posterior ───────────────────────────────────────────────────
p_post = plot(prior; fill = :lightgray, label = "prior")
plot!(posterior; fill = :darkblue, alpha = 0.5, label = "posterior")
for (i, sp) in enumerate(p_post.subplots)
    vline!(sp, [constrained_ekp_optimal[i]]; lc = :red, ls = :dash, lw = 2, label = "EKI")
    title!(sp, param_names[i]; titlefontsize = 7)
end
plot!(p_post; size = (1800, 1400), margin = 4Plots.mm)
savefig(p_post, joinpath(ces_out_dir, "posterior_marginals.png"))
@info "Posterior marginals → $(joinpath(ces_out_dir, "posterior_marginals.png"))"

@info "=== CES pipeline complete. Next: run run_posterior_ensemble.jl ==="
