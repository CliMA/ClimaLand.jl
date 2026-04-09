"""
Posterior parameter distribution plots for CalLMIP Phase 1a — DK-Sor.

Generates two figures that characterise the CES posterior over the 16
calibrated parameters, without requiring forward model runs:

  output_evaluation/posterior_marginals.png
      — 4×4 grid of marginal histograms (prior in grey, posterior in blue)
        with prior mean (red dashed) and EKI optimal (green solid) marked.

  output_evaluation/posterior_correlation.png
      — 16×16 Pearson correlation heatmap of the posterior samples.

Requires: experiments/callmip_uq_dk_sor/output_posterior_uq/posterior_its1to8.jld2
Usage:
    julia --project=experiments/callmip_uq_dk_sor \\
          experiments/callmip_uq_dk_sor/plot_posterior_distributions.jl
"""

println("[1/4] Loading packages..."); flush(stdout)
using JLD2, Statistics, LinearAlgebra
using Plots, Printf
println("[2/4] Packages loaded."); flush(stdout)

# ── Paths ──────────────────────────────────────────────────────────────────────
const exp_dir      = @__DIR__
const posterior_dir = joinpath(exp_dir, "output_posterior_uq")
const out_dir       = joinpath(exp_dir, "output_evaluation")
isdir(out_dir) || mkdir(out_dir)

# ── Load CES posterior ──────────────────────────────────────────────────────────
println("[3/4] Loading posterior samples..."); flush(stdout)
ces_files = filter(f -> startswith(f, "posterior_its") && endswith(f, ".jld2"),
                   readdir(posterior_dir))
isempty(ces_files) && error("No posterior JLD2 found in $posterior_dir")
d = JLD2.load(joinpath(posterior_dir, last(sort(ces_files))))

param_names    = String.(d["param_names"])           # (16,)
post_samples   = Float64.(d["constrained_posterior"]) # (16, 100000)
ekp_optimal    = Float64.(d["constrained_ekp_optimal"])
n_params, n_samples = size(post_samples)

# Prior means & sigmas (from priors.jl — constrained_gaussian μ, σ arguments)
prior_means  = [0.5, 0.43, 20.0, 0.1, 0.10, 0.001, 0.65, 0.85, 0.97,
                25000.0, 0.01, 0.01, 61000.0, 1.0, 0.1, 2500.0]
prior_sigmas = [0.3, 0.15, 8.0, 0.05, 0.04, 0.0005, 0.12, 0.25, 0.02,
                10000.0, 0.005, 0.005, 15000.0, 0.5, 0.05, 1500.0]
@assert length(prior_means) == n_params
@assert length(prior_sigmas) == n_params

# Short display labels
labels = ["moisture_stress_c", "pmodel_cstar", "pmodel_β", "leaf_Cd",
          "z_0m_coeff", "z_0b_coeff", "d_coeff", "K_lw",
          "emissivity", "soilCO2_preexp", "Km_CO2", "Km_O2",
          "soilCO2_Ea", "root_N_ratio", "stem_N_ratio", "ac_canopy"]

println("[4/4] Generating figures..."); flush(stdout)
gr()

# ── Figure 1: Marginal distributions (4 × 4 grid) ─────────────────────────────
nrows, ncols = 4, 4
plt_grid = []
for i in 1:n_params
    samps = post_samples[i, :]
    μ, σ  = prior_means[i], prior_sigmas[i]
    opt   = ekp_optimal[i]

    # Clip to [μ−4σ, μ+4σ] for sensible x-axis (prior range)
    xmin = max(minimum(samps), μ - 4σ)
    xmax = min(maximum(samps), μ + 4σ)

    p = histogram(samps;
        bins = range(xmin, xmax; length = 51),
        normalize = :pdf,
        color = :steelblue, alpha = 0.7, label = "",
        title  = labels[i],
        titlefontsize = 8,
        xlabel = "",
        ylabel = i in (1, 5, 9, 13) ? "PDF" : "",
        xtickfontsize = 6, ytickfontsize = 6,
        xlims = (xmin, xmax),
    )
    vline!(p, [μ];   color = :red,   lw = 1.5, ls = :dash,  label = "prior μ")
    vline!(p, [opt]; color = :green, lw = 1.5, ls = :solid, label = "EKI opt")
    push!(plt_grid, p)
end

fig1 = plot(plt_grid...;
    layout = (nrows, ncols),
    size   = (1600, 1000),
    plot_title = "DK-Sor: Posterior Marginal Distributions (blue) vs Prior Mean (red) vs EKI Optimal (green)",
    plot_titlefontsize = 10,
    margin = 3Plots.mm,
)
out1 = joinpath(out_dir, "posterior_marginals.png")
savefig(fig1, out1)
@info "Saved posterior_marginals.png"

# ── Figure 2: Posterior correlation heatmap ────────────────────────────────────
# Subsample for memory efficiency
idx = 1:10:n_samples
S   = post_samples[:, idx]          # (16, 10000)
R   = cor(S')                       # (16, 16) Pearson correlation

# Threshold weak correlations for readability
hmap_vals = R

fig2 = heatmap(
    labels, labels, hmap_vals;
    color    = :RdBu,
    clims    = (-1.0, 1.0),
    aspect_ratio = :equal,
    size     = (900, 820),
    title    = "DK-Sor: Posterior Parameter Correlations",
    titlefontsize = 12,
    xrotation = 45,
    xtickfontsize = 8,
    ytickfontsize = 8,
    colorbar_title = "Pearson r",
)
# Annotate cells with |r| > 0.3
for i in 1:n_params, j in 1:n_params
    rij = round(R[i,j], digits=2)
    if abs(rij) > 0.3 && i != j
        annotate!(fig2, j, i, Plots.text(@sprintf("%.2f", rij), 6, :center,
                                          rij > 0 ? :darkred : :darkblue))
    end
end
out2 = joinpath(out_dir, "posterior_correlation.png")
savefig(fig2, out2)
@info "Saved posterior_correlation.png"

println()
println("=" ^ 60)
println("  Posterior distribution plots saved to output_evaluation/")
println("    posterior_marginals.png")
println("    posterior_correlation.png")
println("=" ^ 60)
