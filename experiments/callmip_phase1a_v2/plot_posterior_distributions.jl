"""
Plot prior vs posterior marginal distributions for all 16 calibrated parameters.

Usage:
    julia --project=experiments/callmip_phase1a_v2 \\
          experiments/callmip_phase1a_v2/plot_posterior_distributions.jl
"""

import JLD2
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using Plots
using Statistics

const exp_dir    = @__DIR__
const ces_out_dir = joinpath(exp_dir, "output_ces")
const fig_dir    = joinpath(exp_dir, "figures")
isdir(fig_dir) || mkpath(fig_dir)

include(joinpath(exp_dir, "priors.jl"))
prior, priors_vec = build_dk_sor_priors()
param_names = [only(PD.get_name(d)) for d in priors_vec]
n_params = length(param_names)

# ── Load posterior MCMC samples ──────────────────────────────────────────────
post_file = joinpath(ces_out_dir, "posterior_fixed_window.jld2")
isfile(post_file) || error("Posterior not found: $post_file\nRun emulate_sample.jl first.")
constrained_posterior = JLD2.load(post_file, "constrained_posterior")   # (n_params × n_samples)

# ── Draw prior samples ────────────────────────────────────────────────────────
using Random
rng = Random.MersenneTwister(7)
n_samples = size(constrained_posterior, 2)
u_prior = PD.sample(rng, prior, n_samples)   # unconstrained (n_params × n_samples)
constrained_prior = EKP.transform_unconstrained_to_constrained(prior, u_prior)

# ── Panel plot ────────────────────────────────────────────────────────────────
n_cols = 4
n_rows = ceil(Int, n_params / n_cols)

figs = []
for i in 1:n_params
    pr_vals = constrained_prior[i, :]
    po_vals = constrained_posterior[i, :]

    lo = min(quantile(pr_vals, 0.01), quantile(po_vals, 0.01))
    hi = max(quantile(pr_vals, 0.99), quantile(po_vals, 0.99))
    bins = range(lo, hi; length = 40)

    fig = histogram(pr_vals;
                    bins = bins, alpha = 0.45, label = "prior",
                    color = :dodgerblue, normalize = :pdf)
    histogram!(fig, po_vals;
               bins = bins, alpha = 0.65, label = "posterior",
               color = :orangered, normalize = :pdf)
    title!(fig, param_names[i]; titlefontsize = 8)
    xlabel!(fig, ""; guidefontsize = 6)
    push!(figs, fig)
end

big_fig = plot(figs...; layout = (n_rows, n_cols), size = (1200, 900),
               plot_title = "DK-Sor CalLMIP Phase 1a — Prior vs Posterior")
fpath = joinpath(fig_dir, "03_posterior_marginals.png")
savefig(big_fig, fpath)
@info "Saved $fpath"

# ── Print summary statistics ──────────────────────────────────────────────────
println("\nParameter summary: prior mean → posterior mean (posterior std)")
for i in 1:n_params
    pr_m = mean(constrained_prior[i, :])
    po_m = mean(constrained_posterior[i, :])
    po_s = std(constrained_posterior[i, :])
    @printf "  %-40s  prior=%9.4g  →  post=%9.4g ± %9.4g\n" param_names[i] pr_m po_m po_s
end
using Printf
