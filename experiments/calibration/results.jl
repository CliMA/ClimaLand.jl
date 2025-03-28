import EnsembleKalmanProcesses as EKP
using JLD2

using Plots
using DelimitedFiles
using LinearAlgebra, Statistics, Random

folder_name = "calibration_output_utki_TEMP"

iterations = readdir(folder_name)

# split this and explain steps
final_iteration = parse(Int, last(sort(first.(filter(!isnothing, match.(r"^iteration_(\d*)$", iterations))))))

iteration_string = "iteration_"*lpad(final_iteration, 3, "0")
path_to_results = joinpath(folder_name, iteration_string, "eki_file.jld2")
uki = JLD2.load_object(path_to_results);
n_ensembles = EKP.get_N_ens(uki)
errors = uki.error
normalized_errors = errors ./ errors[1] .* 100

g_mean_first = EKP.get_g_mean(uki, 1) # is 0 not there?
g_mean_last = EKP.get_g_mean(uki, final_iteration) # why is there 4 and not 5 ???
g_all = [EKP.get_g_mean(uki, i) for i in 1:final_iteration]

training_locations = JLD2.load_object("locations.jld2")

y = EKP.get_obs(uki) # Ollie, how do I get y for minibatcher (year) 1, 2, 3... epoch 1?

obs_series = EKP.get_observation_series(uki)
y_obs = obs_series.observations
y_all = [EKP.get_obs(y_obs[i]) for i in 1:final_iteration]
[y_all[i] |> mean for i in 1:final_iteration]

# get LHF
"""
    n: 1=lhf, 2=shf, 3=swu, 4=lwu
"""
function slicevar(v, n)
    idx = vcat([(i+(n-1)*4:i+(n-1)*4+3) for i in 1:16:length(v)]...)
    return v[idx]
end
iteration_n = 1
LHF_y_n = slicevar(y_all[iteration_n], 1)
SHF_y_n = slicevar(y_all[iteration_n], 2)
SWU_y_n = slicevar(y_all[iteration_n], 3)
LWU_y_n = slicevar(y_all[iteration_n], 4)

LHF_g_n = slicevar(g_all[iteration_n], 1)
SHF_g_n = slicevar(g_all[iteration_n], 2)
SWU_g_n = slicevar(g_all[iteration_n], 3)
LWU_g_n = slicevar(g_all[iteration_n], 4)

function RMSE(x, z)
    return sqrt(mean((x.-z).^2))
end

function error_abs(x, z)
    return mean(abs.(x.-z))
end

stds = [RMSE(y, g[:,i]) for i in 1:size(g)[2]]
error_abs_all = [error_abs(y, g[:,i]) for i in 1:size(g)[2]]

#sqrt_errors_squared = zeros(length(y), n_ensembles)
#[sqrt_errors_squared[:, i] = (y-g[:,i]).^2 for i in 1:n_ensembles]
#stds = zeros(n_ensembles)
#[stds[i] = sqrt.(mean(sqrt_errors_squared[i])) for i in 1:n_ensembles]
params = EKP.get_ϕ.(Ref(prior), Ref(uki), collect(1:final_iteration))

#julia> params = EKP.get_ϕ.(Ref(prior), Ref(uki), collect(1:2))
#2-element Vector{Matrix{Float64}}:
# [0.6197082648086368 0.8140179617290034 0.3776083422584018]
# [0.1861335919481284 0.1912956574662672 0.1810796330207397]

# Plotting

# first season
lons = map(x -> x[1], training_locations)
lats = map(x -> x[2], training_locations)
y_winter = y[1:4:end]
y_spring = y[2:4:end]
y_summer = y[3:4:end]
y_fall = y[4:4:end]
#g = g[:, 3]
g_winter = g_mean_last[1:4:end]
g_spring = g_mean_last[2:4:end]
g_summer = g_mean_last[3:4:end]
g_fall = g_mean_last[4:4:end]

using Plots

# Scatter plot colored by z
p = Plots.scatter(lons, lats, marker_z = y_summer .- g_summer, legend = false, colorbar = true,
        xlabel = "lons", ylabel = "lats", title = "Scatter plot colored by g_summer")
display(p)




seasons_mean_g = [mean(season) for season in [g_winter, g_spring, g_summer, g_fall]]
seasons_mean_y = [mean(season) for season in [y_winter, y_spring, y_summer, y_fall]]
p = Plots.scatter(seasons_mean_g, label = "g")
Plots.scatter!(p, seasons_mean_y, label = "y")
display(p)

function global_av(x, lat)



    end











#using UnicodePlots
#A = reshape(y_winter, length(lons), length(lats))'

using UnicodePlots, StatsBase

# Define bin edges
n_lat_bins = 40
n_lon_bins = 80

lat_edges = range(minimum(lats), stop=maximum(lats), length=n_lat_bins + 1)
lon_edges = range(minimum(lons), stop=maximum(lons), length=n_lon_bins + 1)

# Initialize matrix and counters
A = zeros(n_lat_bins, n_lon_bins)
counts = zeros(Int, n_lat_bins, n_lon_bins)

# Bin the data
for i in 1:length(y_winter)
    lat_idx = searchsortedfirst(lat_edges, lats[i]) - 1
    lon_idx = searchsortedfirst(lon_edges, lons[i]) - 1

    if lat_idx in 1:n_lat_bins && lon_idx in 1:n_lon_bins
        A[lat_idx, lon_idx] += y_winter[i] - g_winter[i]
        counts[lat_idx, lon_idx] += 1
    end
end

# Average the values in each bin
for i in 1:n_lat_bins, j in 1:n_lon_bins
    if counts[i, j] > 0
        A[i, j] /= counts[i, j]
    else
        A[i, j] = NaN  # or 0 or some other fill value
    end
end

# Optional: replace NaNs with 0 for plotting
A .= coalesce.(A, 0.0)

# Plot
UnicodePlots.heatmap(A, xlabel="Longitude", ylabel="Latitude")
Plots.heatmap(A, xlabel="Longitude", ylabel="Latitude", c = :bluesreds)








## Priors
#prior_pc = EKP.constrained_gaussian("pc", -2e6, 1e6, -Inf, 0);
#prior_sc = EKP.constrained_gaussian("sc", 5e-6, 2e-6, 0, Inf);
prior_K_sat_plant =
    EKP.constrained_gaussian("K_sat_plant", 3e-8, 2e-8, 0.0, Inf);
prior_a = EKP.constrained_gaussian("a", 0.00196, 0.00049, 0.0, 0.00588);
#prior_h_leaf = EKP.constrained_gaussian("h_leaf", 4.0, 2.0, 0.5, 9.0);
prior_α_snow = EKP.constrained_gaussian("α_snow", 0.6, 0.2, 0.0, 1.0); # only applies physical bounds
prior_α_soildry_scale =
    EKP.constrained_gaussian("α_soil_dry_scaler", 1.15, 0.05, 1.0, 1.3); # capped at 1 from below to ensure alpha_wet < alpha_dry, from above to ensure that α_soil_scaler * α_soil_dry_scaler <~1.6, which is when dry albedos start to exceed 1
prior_α_soil_scale =
    EKP.constrained_gaussian("α_soil_scaler", 1, 0.15, 0.7, 1.3); # new parameter
prior_τ_leaf_scale =
    EKP.constrained_gaussian("τ_leaf_scaler", 1, 0.15, 0.5, 1.5);
prior_α_leaf_scale =
    EKP.constrained_gaussian("α_leaf_scaler", 1, 0.15, 0.5, 1.5);
prior = EKP.combine_distributions([
#    prior_pc,
#    prior_sc,
#    prior_h_leaf,
    prior_a,
    prior_K_sat_plant,
    prior_α_leaf_scale,
    prior_τ_leaf_scale,
    prior_α_soildry_scale,
    prior_α_snow,
    prior_α_soil_scale,
]);



names = [
#    "pc",
#    "sc",
#    "h_leaf",
#    "a",
#    "K_sat_plant",
    "α_leaf_scaler",
    "τ_leaf_scaler",
    "α_soil_dry_scaler",
    "α_snow",
    "α_soil_scale",
]

params_mean = EKP.get_ϕ_mean.(Ref(prior), Ref(uki), collect(1:final_iteration))
params = EKP.get_ϕ.(Ref(prior), Ref(uki), collect(1:1))
initial_params = params_mean[1]
calibrated_params = params_mean[end]
relative_change = calibrated_params ./ initial_params .* 100
Dict(zip(names, initial_params))
Dict(zip(names, calibrated_params))
Dict(zip(names, relative_change))


EKP.get_ϕ(prior, uki, 1)
EKP.get_ϕ(prior, uki, 2)

# Other things to do, to look at how calibration improved model - data

# Plot data vs season, and model vs season before and after calibration, for all locations, and all variables, maybe for a subset of locations for visibility (2 or 3)

# Global plots: data vs model, actual values, anomalies, RMSE, seasonality
# (use CI leaderboard figures, and ClimaLand leaderboard)

p = Plots.plot(prior, xtickfontsize=9, dpi=300)
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

savefig(p, "NEWPLOT.png")
display(p)

include("experiments/calibration/global_plots.jl")
outdir = "calibration_output_utki_TEMP/iteration_004/member_003/global_diagnostics/output_active/"
make_global_plots(outdir; root_save_path = pwd(), plot_net_radiation = false)
make_seasonal_plots(outdir; root_save_path = pwd(), plot_net_radiation = false)

# Leaderboard
include(joinpath(pkgdir(ClimaLand), "experiments", "long_runs", "leaderboard", "leaderboard.jl"))
compute_leaderboard(pwd(), outdir)


