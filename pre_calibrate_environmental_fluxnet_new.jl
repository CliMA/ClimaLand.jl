using Statistics
using StatsBase
using Dates
using Flux
using CairoMakie
using DelimitedFiles
using NCDatasets
using ClimaLand.Artifacts
using ClimaAnalysis
using StatsBase
include("models.jl")
(data, columns) = readdlm("all_data.csv", ',', header = true)
diffuse_fraction = readdlm("diffuse_fraction.txt", ',', header = false)

FT = Float32
albedo =
    FT.(data[:, findfirst(columns .== "sw_u")[2]]) ./
    FT.(data[:, findfirst(columns .== "sw_d")[2]]);
mask =
    albedo .< 0.5 .&&
    FT.(data[:, findfirst(columns .== "lai")[2]]) .> 2.5 .&&
    FT.(data[:, findfirst(columns .== "sw_d")[2]]) .> 10 .&&
    FT.(data[:, findfirst(columns .== "sw_u")[2]]) .> 0
μ = FT.(data[mask, findfirst(columns .== "cosZ")[2]]);
LAI = FT.(data[mask, findfirst(columns .== "lai")[2]]);
α_PAR_g = FT.(data[mask, findfirst(columns .== "a_par_soil")[2]]);
α_NIR_g = FT.(data[mask, findfirst(columns .== "a_nir_soil")[2]]);
SW_u = FT.(data[mask, findfirst(columns .== "sw_u")[2]]);
SW_d = FT.(data[mask, findfirst(columns .== "sw_d")[2]]);
lat = FT.(data[mask, findfirst(columns .== "lat")[2]]);
long = FT.(data[mask, findfirst(columns .== "long")[2]]);
site_ID = data[mask, findfirst(columns .== "siteID")[2]];
f_diff = FT.(diffuse_fraction[mask, 1])
function nearest_neighbor(x, x_nc, y, y_nc, v_nc)
    dist_x = abs.(x .- x_nc)
    dist_y = abs.(y .- y_nc)
    values = v_nc[dist_x .< 1, dist_y .< 1]
    mask = .~isnan.(values)
    return mean(values[mask])
end
N = sum(mask)
# Read in some other fields from net cdf that could be useful predictors
data_path = Artifacts.optimal_lai_initial_conditions_path()
ds = NCDataset(data_path)
gsl_env = ds["gsl"][:, :]
gs_vpd_env = ds["vpd_gs"][:, :]
ann_p_env = ds["precip_annual"][:, :]
ann_a0_env = ds["a0_annual"][:, :]
f0_env = ds["f0"][:, :]
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
base_lon = x_nc
base_lat = y_nc
gsl = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, gsl_env) for k in 1:N]
gs_vpd =
    [nearest_neighbor(long[k], x_nc, lat[k], y_nc, gs_vpd_env) for k in 1:N]
ann_p = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, ann_p_env) for k in 1:N]
f0 = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, f0_env) for k in 1:N]
ann_a0 = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, ann_a0_env) for k in 1:N]
close(ds)


data_path = Artifacts.clm_tree_share_path()
tree_share_var = ClimaAnalysis.OutputVar(data_path, "tree_share")
ts_env = ClimaAnalysis.shift_longitude(tree_share_var, -180.0, 180.0).data
x_nc = ClimaAnalysis.shift_longitude(tree_share_var, -180.0, 180.0).dims["lon"]
y_nc = ClimaAnalysis.shift_longitude(tree_share_var, -180.0, 180.0).dims["lat"]
ts = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, ts_env) for k in 1:N]

topo_data_path = "topographic_data.nc"
ds = NCDataset(topo_data_path)
elevation_topo = ds["mean_height"][:, :]
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
elevation =
    [nearest_neighbor(long[k], x_nc, lat[k], y_nc, elevation_topo) for k in 1:N]
close(ds)

vcmax25_data_path = "vcmax25_1M_average.nc"
ds = NCDataset(vcmax25_data_path)
vc = maximum(ds["vcmax25"][13:end,:, :], dims = 1)[1,:,:];
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
vcmax25 = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, vc) for k in 1:N]
close(ds)

tair_m_data_path = "tair_1M_average.nc"
ds = NCDataset(tair_m_data_path)
tair_min = minimum(ds["tair"][13:end,:, :], dims = 1)[1,:,:];
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
min_tair = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, tair_min) for k in 1:N]
close(ds)

tair_data_path = "../ClimaArtifacts/optimal_lai_inputs/tair_gs.nc"
ds = NCDataset(tair_data_path)
tair_env = ds["tair"][:, :]
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
gs_tair = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, tair_env) for k in 1:N]
close(ds)

m = 10
x = zeros(FT, N, m);
x[:, 1] .= ann_p ./ 365
x[:, 2] .= gsl ./ 365
x[:, 3] .= log.(gs_vpd)
x[:, 4] .= gs_tair
x[:, 5] .= log.(elevation .+ 1000)
x[:,6] .= log.(vcmax25 .+ 1e-7)
x[:,7] .= f0
x[:,8] .= min_tair ./ gs_tair
x[:,9] .= ann_a0
x[:,10] .= ts


x = permutedims(x, (2, 1))
x̄ = mean(x, dims = (2))
σ = std(x, dims = (2))
x̂ = @. Float32.((x-x̄)/σ)
μ = reshape(μ, (1, N))
LAI = reshape(LAI, (1, N))
α_PAR_g = reshape(α_PAR_g, (1, N))
α_NIR_g = reshape(α_NIR_g, (1, N))
SW_d = reshape(SW_d, (1, N))
SW_u = reshape(SW_u, (1, N))
f_diff = reshape(f_diff, (1,N))
 


m2 = m
mRT = LandRT2Stream(
    randn(Float32, (m2, m)),
    randn(Float32, (m2, m)),
    randn(Float32, (m2, 1)),
    randn(Float32, (m2, 1)),
    randn(Float32, (1, m2)),
    randn(Float32, (1, m2)),
    randn(Float32, (1, 1)),
    randn(Float32, (1, 1)),
    randn(Float32, (1, 1)),
    randn(Float32, (1, 1)),
    randn(Float32, (1, 1)),
    randn(Float32, (1, 1)),
)
opt = Flux.setup(Adam(FT(0.003), FT.((0.9, 0.99)), FT(1e-8)), mRT)
obs = SW_u ./ SW_d
loss(
    model,
    x,
    obs;
    μ = μ,
    LAI = LAI,
    α_PAR_g = α_PAR_g,
    α_NIR_g = α_NIR_g,
    SW_d = SW_d,
    f_diff = f_diff
) = Flux.mse(model(x, LAI, μ, α_PAR_g, α_NIR_g, SW_d, f_diff), obs)

@show loss(mRT, x̂, obs)
train_loss_values = []
for i in 1:1200
    Flux.train!(loss, mRT, [(x̂, obs)], opt)
    if i % 10 == 0
        L = loss(mRT, x̂, obs)
        push!(train_loss_values, L)
        @show L
    end
end

pred = mRT(x̂, LAI, μ, α_PAR_g, α_NIR_g, SW_d, f_diff);
@show compute_r2(obs[:], pred[:])
@show sqrt(loss(mRT, x̂, obs))
@show mean(pred[:] .- obs[:])
apar = get_α_PAR(mRT, x̂)
anir =get_α_NIR(mRT, x̂)
fig = Figure()
ax = Axis(fig[1,1], xlabel = "Leaf albedo", ylabel = "Frequency")
hist!(ax, apar[:], bins = 40, label = "PAR")
hist!(ax, anir[:], bins = 40, label = "NIR")
axislegend(ax; position = :rt)



fig = Figure()
ax = Axis(fig[1,1], xlabel = "Observed α", ylabel = "Predicted α")
h = fit(Histogram, (obs[:], pred[:]), (0:0.01:0.3, 0:0.01:0.3))
hm = heatmap!(ax, h.edges[1], h.edges[2], h.weights; colormap = :viridis)
Colorbar(fig[1, 2], hm, label = "Count")
lines!(ax, obs[:], obs[:], color = :red, linestyle = :dash)
xlims!(ax, (0, 0.3))
ylims!(ax, (0, 0.3))
save("obs_v_pred2s.png", fig)
fig = Figure()
ax = Axis(fig[1,1], xlabel = "Predicted", ylabel = "Residual")
h = fit(Histogram, (pred[:], pred[:] .- obs[:]), (0:0.01:0.3, -0.3:0.01:0.3))
hm = heatmap!(ax, h.edges[1], h.edges[2], h.weights; colormap = :viridis)
Colorbar(fig[1, 2], hm, label = "Count")
lines!(ax, pred[:], zeros(length(pred[:])), color = :red, linestyle = :dash)
xlims!(ax, (0, 0.3))
ylims!(ax, (-0.15, 0.15))
save("pred_v_res2s.png", fig)

function dependence_plot(y, y_name)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Observed α", ylabel =y_name)
    h = fit(Histogram, (obs[:], y), (0:0.01:0.3, minimum(y):(maximum(y) - minimum(y))/40:maximum(y)))
    hm = heatmap!(ax, h.edges[1], h.edges[2], h.weights; colormap = :viridis)
    Colorbar(fig[1, 2], hm, label = "Count")
    xlims!(ax, (0, 0.3))
    ylims!(ax, (minimum(y), maximum(y)))
    ax = Axis(fig[2,1], xlabel = "Predicted α", ylabel = "Vcmax25")
    h = fit(Histogram, (pred[:], y), (0:0.01:0.3, minimum(y):(maximum(y) - minimum(y))/40:maximum(y)))
    hm = heatmap!(ax, h.edges[1], h.edges[2], h.weights; colormap = :viridis)
    Colorbar(fig[2, 2], hm, label = "Count")
    xlims!(ax, (0, 0.3))
    ylims!(ax, (minimum(y), maximum(y)))
    save("$(y_name)_dep.png", fig)
end
    



####
#Two stream

clm_data_path = joinpath(Artifacts.clm_data_folder_path(; lowres = false),"vegetation_properties_map.nc")
χ_clm_var = ClimaAnalysis.OutputVar(clm_data_path, "xl")
χ_clm = ClimaAnalysis.shift_longitude(χ_clm_var, -180.0, 180.0).data
α_PAR_clm_var = ClimaAnalysis.OutputVar(clm_data_path, "rholvis")
α_PAR_clm = ClimaAnalysis.shift_longitude(α_PAR_clm_var, -180.0, 180.0).data;
α_NIR_clm_var = ClimaAnalysis.OutputVar(clm_data_path, "rholnir")
α_NIR_clm = ClimaAnalysis.shift_longitude(α_NIR_clm_var, -180.0, 180.0).data
τ_PAR_clm_var = ClimaAnalysis.OutputVar(clm_data_path, "taulvis")
τ_PAR_clm = ClimaAnalysis.shift_longitude(τ_PAR_clm_var, -180.0, 180.0).data
τ_NIR_clm_var = ClimaAnalysis.OutputVar(clm_data_path, "taulnir")
τ_NIR_clm = ClimaAnalysis.shift_longitude(τ_NIR_clm_var, -180.0, 180.0).data

x_clm = ClimaAnalysis.shift_longitude(τ_NIR_clm_var, -180.0, 180.0).dims["lon"]
y_clm = ClimaAnalysis.shift_longitude(τ_NIR_clm_var, -180.0, 180.0).dims["lat"]
χ_clm_sites = [nearest_neighbor(long[k], x_clm, lat[k], y_clm, χ_clm) for k in 1:N]
α_PAR_clm_sites = [nearest_neighbor(long[k], x_clm, lat[k], y_clm, α_PAR_clm) for k in 1:N]
α_NIR_clm_sites = [nearest_neighbor(long[k], x_clm, lat[k], y_clm, α_NIR_clm) for k in 1:N]
τ_PAR_clm_sites = [nearest_neighbor(long[k], x_clm, lat[k], y_clm, τ_PAR_clm) for k in 1:N]
τ_NIR_clm_sites = [nearest_neighbor(long[k], x_clm, lat[k], y_clm, τ_NIR_clm) for k in 1:N]
close(ds)
χ_clm_sites[isnan.(χ_clm_sites)] .= 0.25
α_PAR_clm_sites[isnan.(α_PAR_clm_sites)] .= 0.1
α_NIR_clm_sites[isnan.(α_NIR_clm_sites)] .= 0.45
τ_PAR_clm_sites[isnan.(τ_PAR_clm_sites)] .= 0.05
τ_NIR_clm_sites[isnan.(τ_NIR_clm_sites)] .= 0.25
modis_ci_artifact_path = Artifacts.modis_ci_data_folder_path()
ds = NCDataset(joinpath(modis_ci_artifact_path, "He_et_al_2012_1x1.nc"))
Ω = ds["ci"][:,:]
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
Ω_sites = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, Ω) for k in 1:N]
close(ds)

clm_par = ClimaLand.Canopy.canopy_sw_rt_two_stream.(ClimaLand.Canopy.CLMGFunction.(χ_clm_sites), Ω_sites, UInt64(20), α_PAR_clm_sites, τ_PAR_clm_sites, LAI[:], μ[:], α_PAR_g[:], f_diff[:])
clm_nir = ClimaLand.Canopy.canopy_sw_rt_two_stream.(ClimaLand.Canopy.CLMGFunction.(χ_clm_sites), Ω_sites, UInt64(20), α_NIR_clm_sites, τ_NIR_clm_sites, LAI[:], μ[:], α_NIR_g[:], f_diff[:])
SW_u_clm = ([f.refl for f in clm_par] .+ [f.refl for f in clm_nir]) ./ 2 .* SW_d[:]
α_clm = SW_u_clm ./ SW_d[:]
@show compute_r2(obs[:], α_clm)
@show sqrt(mean((α_clm .- obs[:]).^2))
@show mean(α_clm[:] .- obs[:])
fig = Figure()
ax = Axis(fig[1,1], xlabel = "Observed α", ylabel = "Predicted α")
h = fit(Histogram, (obs[:], α_clm[:]), (0:0.01:0.3, 0:0.01:0.3))
hm = heatmap!(ax, h.edges[1], h.edges[2], h.weights; colormap = :viridis)
Colorbar(fig[1, 2], hm, label = "Count")
lines!(ax, obs[:], obs[:], color = :red, linestyle = :dash)
xlims!(ax, (0, 0.3))
ylims!(ax, (0, 0.3))
save("obs_v_pred_clm2s.png", fig)
fig = Figure()
ax = Axis(fig[1,1], xlabel = "Predicted", ylabel = "Residual")
h = fit(Histogram, (α_clm[:], α_clm[:] .- obs[:]), (0:0.01:0.3, -0.3:0.01:0.3))
hm = heatmap!(ax, h.edges[1], h.edges[2], h.weights; colormap = :viridis)
Colorbar(fig[1, 2], hm, label = "Count")
lines!(ax, α_clm[:], zeros(N), color = :red, linestyle = :dash)
xlims!(ax, (0, 0.3))
ylims!(ax, (-0.15, 0.15))
save("pred_v_res_clm2s.png", fig)
