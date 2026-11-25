using ClimaAnalysis
data_path = Artifacts.optimal_lai_initial_conditions_path()
ds = NCDataset(data_path)
gsl_env = ds["gsl"][:, :]
gs_vpd_env = ds["vpd_gs"][:, :]
ann_p_env = ds["precip_annual"][:, :]
ann_a0_env = ds["a0_annual"][:, :]
f0_env = ds["f0"][:, :]
base_lon = ds["lon"][:]
base_lat = ds["lat"][:]
close(ds)
vcmax25_data_path = "vcmax25_1M_average.nc"
ds = NCDataset(vcmax25_data_path)
vc = maximum(ds["vcmax25"][13:end,:, :], dims = 1)[1,:,:];
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
@assert x_nc == base_lon
@assert y_nc == base_lat
close(ds)

tair_m_data_path = "tair_1M_average.nc"
ds = NCDataset(tair_m_data_path)
tair_min = minimum(ds["tair"][13:end,:, :], dims = 1)[1,:,:];
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
@assert x_nc == base_lon
@assert y_nc == base_lat
close(ds)

tair_data_path = "../ClimaArtifacts/optimal_lai_inputs/tair_gs.nc"
ds = NCDataset(tair_data_path)
tair_env = ds["tair"][:, :]
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
@assert x_nc == base_lon
@assert y_nc == base_lat
close(ds)

data_path = Artifacts.clm_tree_share_path()
tree_share_var = ClimaAnalysis.OutputVar(data_path, "tree_share")
ts_env = ClimaAnalysis.shift_longitude(tree_share_var, -180.0, 180.0).data

topo_data_path = "topographic_data.nc"
elev_var0 = ClimaAnalysis.OutputVar(topo_data_path, "mean_height")
elev_var = ClimaAnalysis.resampled_as(elev_var0; nan_threshold = 0.99, lon = base_lon, lat = base_lat)
elevation_topo = elev_var.data

x = zeros(length(base_lon), length(base_lat), 10);
x[:,:,1] .= (ann_p_env ./ 365 .- x̄[1])./ σ[1];
x[:,:,2] .= (gsl_env ./ 365 .- x̄[2])./ σ[2];
x[:,:,3] .= (log.(gs_vpd_env) .- x̄[3])./ σ[3];
x[:,:,4] .= (tair_env .- x̄[4])./ σ[4];
x[:,:,5] .= (log.(max.(elevation_topo .+1000, 1.0)) .- x̄[5])./ σ[5]; # not on correct grid
x[:,:,6] .= (log.(vc .+ 1e-7) .- x̄[6])./ σ[6];
x[:,:,7] .= (f0_env .- x̄[7])./ σ[7];
x[:,:,8] .= (tair_min ./ tair_env .- x̄[8])./ σ[8];
x[:,:,9] .= (ann_a0_env .- x̄[9])./ σ[9];
x[:,:,10] .= (ts_env .- x̄[10])./ σ[10];


path = Artifacts.modis_max_lai_data_path()
max_lai_var0 = ClimaAnalysis.OutputVar(path, "lai")
max_lai_var = ClimaAnalysis.resampled_as(max_lai_var0; lon = base_lon, lat = base_lat)
max_lai = max_lai_var.data

no_veg = max_lai .< 1
x[no_veg,:] .= NaN

function evaluate_nn(m, x)
    fPAR(x; M = m) = get_α_PAR(m, x)
    fNIR(x; M = m) = get_α_NIR(m, x)

    α_PAR = mapslices(fPAR, x, dims = (3));
    α_NIR = mapslices(fNIR, x, dims = (3));
    return α_PAR, α_NIR
end
non_nan_mean(x) = mean(x[.~ isnan.(x)])


(α_PAR_global, α_NIR_global) = evaluate_nn(mRT, x);
fig = Figure()
ax = Axis(fig[1,1], title = "αPAR")
hm = heatmap!(ax, α_PAR_global[:,:,1], colorrange = (0,1))
Colorbar(fig[1, 2], hm)
axb = Axis(fig[1,3], xlabel = "Mean", ylabel = "Latitude")
lines!(axb, mapslices(non_nan_mean, α_PAR_global[:,:,1], dims = 1)[:], base_lat)
ax2 = Axis(fig[2,1], title = "αNIR")
hm= heatmap!(ax2, α_NIR_global[:,:,1], colorrange = (0,1))
Colorbar(fig[2, 2], hm)
ax2b = Axis(fig[2,3], xlabel = "Mean", ylabel = "Latitude")
lines!(ax2b, mapslices(non_nan_mean, α_NIR_global[:,:,1], dims = 1)[:], base_lat)
save("global_results.png", fig)

α_PAR_clm[α_PAR_clm .==0] .= NaN
α_NIR_clm[α_NIR_clm .==0] .= NaN

fig = Figure()
ax = Axis(fig[1,1], title = "αPAR")
heatmap!(ax, α_PAR_clm, colorrange = (0,1))
Colorbar(fig[1, 2], hm)
axb = Axis(fig[1,3], xlabel = "Mean", ylabel = "Latitude")
lines!(axb, mapslices(non_nan_mean, α_PAR_clm, dims = 1)[:], y_clm)
ax2 = Axis(fig[2,1], title = "αNIR")
heatmap!(ax2, α_NIR_clm, colorrange = (0,1))
Colorbar(fig[2, 2], hm)
ax2b = Axis(fig[2,3], xlabel = "Mean", ylabel = "Latitude")
lines!(ax2b, mapslices(non_nan_mean, α_NIR_clm, dims = 1)[:], y_clm)
save("global_results_clm.png", fig)
