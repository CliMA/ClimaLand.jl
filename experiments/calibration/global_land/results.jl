import EnsembleKalmanProcesses as EKP
using JLD2

path_to_results = "calibration_output_utki_imposeprior_nocos_evenlargervar/iteration_003/eki_file.jld2"

uki = JLD2.load_object(path_to_results)

errors = uki.error
normalized_errors = errors ./ errors[1] .* 100

names = ["pc", "sc", "h_leaf", "a", "K_sat_plant", "α_leaf_scaler", "τ_leaf_scaler", "α_soil_dry_scaler", "α_snow", "α_soil_scale",]
n_iterations =3
params = EKP.get_ϕ_mean.(Ref(prior), Ref(uki), collect(1:n_iterations))
initial_params = params[1]
calibrated_params = params[end]
relative_change = calibrated_params ./ initial_params .* 100

EKP.get_ϕ(prior, uki, 1)
EKP.get_ϕ(prior, uki, 2)

# Other things to do, to look at how calibration improved model - data

# Plot data vs season, and model vs season before and after calibration, for all locations, and all variables, maybe for a subset of locations for visibility (2 or 3)

# Global plots: data vs model, actual values, anomalies, RMSE, seasonality
# (use CI leaderboard figures, and ClimaLand leaderboard)



# call this from e.g. examples/Sinusoid (so plots are in the pkg)
using Plots
using EnsembleKalmanProcesses
using EnsembleKalmanProcesses.ParameterDistributions

using DelimitedFiles
using LinearAlgebra, Statistics, Random

prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -Inf, 0);
prior_sc = EKP.constrained_gaussian("sc", 5e-6, 2e-6, 0, Inf);
prior_K_sat_plant = EKP.constrained_gaussian("K_sat_plant", 3e-8, 2e-8, 0.0, Inf);
prior_a = EKP.constrained_gaussian("a", 0.00196, 0.00049, 0.0, 0.00588);
prior_h_leaf = EKP.constrained_gaussian("h_leaf", 4.0, 2.0, 0.5, 9.0);
prior_α_snow = EKP.constrained_gaussian("α_snow", 0.6, 0.2, 0.0, 1.0); # only applies physical bounds
prior_α_soildry_scale = EKP.constrained_gaussian("α_soil_dry_scaler", 1.15, 0.05, 1.0, 1.3); # capped at 1 from below to ensure alpha_wet < alpha_dry, from above to ensure that α_soil_scaler * α_soil_dry_scaler <~1.6, which is when dry albedos start to exceed 1
prior_α_soil_scale = EKP.constrained_gaussian("α_soil_scaler", 1, 0.15, 0.7, 1.3); # new parameter
prior_τ_leaf_scale = EKP.constrained_gaussian("τ_leaf_scaler", 1, 0.15, 0.5, 1.5);
prior_α_leaf_scale = EKP.constrained_gaussian("α_leaf_scaler", 1, 0.15, 0.5, 1.5);
prior = EKP.combine_distributions([
    prior_pc,
    prior_sc,
    prior_h_leaf,
    prior_a,
    prior_K_sat_plant,
    prior_α_leaf_scale,
    prior_τ_leaf_scale,
    prior_α_soildry_scale,
    prior_α_snow,
    prior_α_soil_scale,
]);

init = params[1]
final = params[end]

p = plot(prior, xtickfontsize=8, dpi=300)
# square updates blue, regular magenta
# uki updates solid, utki dashed
for (i,sp) in enumerate(p.subplots)
    vline!(sp, [init[i,1]], lc="black", lw=4)
    vline!(sp, [init[i,j] for j in 2:size(init,2)], lc="black", lw=4, alpha=0.3, linestyle=:dash)
    vline!(sp, [final[i,1]], lc="magenta", lw=4)
    vline!(sp, [final[i,j] for j in 2:size(init,2)], lc="magenta", lw=4, alpha=0.3, linestyle=:dash)
end

savefig(p, "uki_climaland_iteration.png")


p = plot(prior, xtickfontsize=8, dpi=300)
init = EKP.get_ϕ(prior, uki, 1)
final = EKP.get_ϕ(prior, uki, 3)
# square updates blue, regular magenta
# uki updates solid, utki dashed
for (i,sp) in enumerate(p.subplots)
    vline!(sp, [init[i,1]], lc="black", lw=4)
    vline!(sp, [maximum(init[i,2:end])], lc="black", lw=4, alpha=0.3, linestyle=:dash)
    vline!(sp, [minimum(init[i,2:end])], lc="black", lw=4, alpha=0.3, linestyle=:dash)
    vline!(sp, [final[i,1]], lc="magenta", lw=4)
    vline!(sp, [maximum(final[i,2:end])], lc="magenta", lw=4, alpha=0.3, linestyle=:dash)
    vline!(sp, [minimum(final[i,2:end])], lc="magenta", lw=4, alpha=0.3, linestyle=:dash)
end

savefig(p, "uki_climaland_iteration_1_3.png")

p = plot(prior, xtickfontsize=8, dpi=300)
init = EKP.get_ϕ(prior, uki, 1)
final = EKP.get_ϕ(prior, uki, 4)
# square updates blue, regular magenta
# uki updates solid, utki dashed
for (i,sp) in enumerate(p.subplots)
    vline!(sp, [init[i,1]], lc="black", lw=4)
    vline!(sp, [maximum(init[i,2:end])], lc="black", lw=4, alpha=0.3, linestyle=:dash)
    vline!(sp, [minimum(init[i,2:end])], lc="black", lw=4, alpha=0.3, linestyle=:dash)
    vline!(sp, [final[i,1]], lc="magenta", lw=4)
    vline!(sp, [maximum(final[i,2:end])], lc="magenta", lw=4, alpha=0.3, linestyle=:dash)
    vline!(sp, [minimum(final[i,2:end])], lc="magenta", lw=4, alpha=0.3, linestyle=:dash)
end

savefig(p, "NEWPLOT.png")
display(p)


p = plot(prior, xtickfontsize=8, dpi=300, constrained = false)
init = EKP.get_u( uki, 1)
final = EKP.get_u( uki, 2)
# square updates blue, regular magenta
# uki updates solid, utki dashed
for (i,sp) in enumerate(p.subplots)
    vline!(sp, [init[i,1]], lc="black", lw=4)
    vline!(sp, [maximum(init[i,2:end])], lc="black", lw=4, alpha=0.3, linestyle=:dash)
    vline!(sp, [minimum(init[i,2:end])], lc="black", lw=4, alpha=0.3, linestyle=:dash)
    vline!(sp, [final[i,1]], lc="magenta", lw=4)
    vline!(sp, [maximum(final[i,2:end])], lc="magenta", lw=4, alpha=0.3, linestyle=:dash)
    vline!(sp, [minimum(final[i,2:end])], lc="magenta", lw=4, alpha=0.3, linestyle=:dash)
end

savefig(p, "NEWPLOT3.png")
