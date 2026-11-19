using Statistics
using Dates
using Flux
using CairoMakie
using DelimitedFiles
using NCDatasets
using ClimaLand.Artifacts
(data, columns) = readdlm("all_data.csv", ',', header = true)
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
T̄ = FT.(data[mask, findfirst(columns .== "temp")[2]]);# NaNs
P̄ = FT.(data[mask, findfirst(columns .== "precip")[2]]); # Nans
LAI = FT.(data[mask, findfirst(columns .== "lai")[2]]);
α_PAR_g = FT.(data[mask, findfirst(columns .== "a_par_soil")[2]]);
α_NIR_g = FT.(data[mask, findfirst(columns .== "a_nir_soil")[2]]);
SW_u = FT.(data[mask, findfirst(columns .== "sw_u")[2]]);
SW_d = FT.(data[mask, findfirst(columns .== "sw_d")[2]]);
lat = FT.(data[mask, findfirst(columns .== "lat")[2]]);
long = FT.(data[mask, findfirst(columns .== "long")[2]]);
site_ID = data[mask, findfirst(columns .== "siteID")[2]];
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
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
gsl = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, gsl_env) for k in 1:N]
gs_vpd =
    [nearest_neighbor(long[k], x_nc, lat[k], y_nc, gs_vpd_env) for k in 1:N]
ann_p = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, ann_p_env) for k in 1:N]
close(ds)

topo_data_path = "topographic_data.nc"
ds = NCDataset(topo_data_path)
μ_topo = ds["mean_slope"][:, :]
elevation_topo = ds["mean_height"][:, :]
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
slope = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, μ_topo) for k in 1:N]
elevation =
    [nearest_neighbor(long[k], x_nc, lat[k], y_nc, elevation_topo) for k in 1:N]
close(ds)

vcmax25_data_path = "vcmax25_1M_average.nc"
ds = NCDataset(vcmax25_data_path)
vc = maximum(ds["vcmax25"][13:end,:, :], dims = 1)[1,:,:];
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
vcmax25 = log10.([nearest_neighbor(long[k], x_nc, lat[k], y_nc, vc) for k in 1:N])
close(ds)

tair_data_path = "../ClimaArtifacts/optimal_lai_inputs/tair_gs.nc"
ds = NCDataset(tair_data_path)
tair_env = ds["tair"][:, :]
x_nc = ds["lon"][:]
y_nc = ds["lat"][:]
gs_tair = [nearest_neighbor(long[k], x_nc, lat[k], y_nc, tair_env) for k in 1:N]
close(ds)

#annual precip, gs_vpd, gsl, elevation, gs_tair; separate: μ, LAI, α_PAR_g, α_NIR_g
m = 6
x = zeros(FT, N, m);
x[:, 1] .= ann_p ./ 365
x[:, 2] .= gsl
x[:, 3] .= log.(gs_vpd)
x[:, 4] .= gs_tair
x[:, 5] .= log.(elevation .+ 100)
x[:,6] .= vcmax25
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

    
function compute_r2(y, y_hat)
    ss_res = sum((y .- y_hat) .^ 2)
    ss_tot = sum((y .- mean(y)) .^ 2)
    return 1.0 - (ss_res / ss_tot)
end

struct LandRT2
    a_PAR_w1::Any
    a_NIR_w1::Any
    a_PAR_b1::Any
    a_NIR_b1::Any
    a_PAR_w2::Any
    a_NIR_w2::Any
    a_PAR_b2::Any
    a_NIR_b2::Any
end
Flux.@layer LandRT2
function reflected(α_c, α_g, T)
    reflected_downwards_pass = α_c * (1 - T)
    upwelling_from_soil = α_g * T
    reflected_upwards_pass = upwelling_from_soil * (1 -T) * α_c
    transmitted_upwards_pass = upwelling_from_soil * T
    upwelling_from_land = reflected_downwards_pass + reflected_upwards_pass + transmitted_upwards_pass
    return upwelling_from_land
end

function (m::LandRT2)(x̂, LAI, μ, α_PAR_g, α_NIR_g, SW_d)
    α_PAR = sigmoid(m.a_PAR_w2 * sigmoid.(m.a_PAR_w1 * x̂ .+ m.a_PAR_b1) .+ m.a_PAR_b2)
    α_NIR = sigmoid(m.a_NIR_w2 * sigmoid.(m.a_NIR_w1 * x̂ .+ m.a_NIR_b1) .+ m.a_NIR_b2)
    G = 1/2
    K = min.(G ./ μ, 100.0f0)
    transmitted_fraction = @. exp(-K * LAI)
    upwelling_from_land_par = @. reflected(α_PAR, α_PAR_g, transmitted_fraction)
    upwelling_from_land_nir = @. reflected(α_NIR, α_NIR_g, transmitted_fraction)
    return @. (upwelling_from_land_nir + upwelling_from_land_par)/2
end
m2 = m
mRT = LandRT2(
    randn(Float32, (m2, m)),
    randn(Float32, (m2, m)),
    randn(Float32, (m2, 1)),
    randn(Float32, (m2, 1)),
    randn(Float32, (1, m2)),
    randn(Float32, (1, m2)),
    randn(Float32, (1, 1)),
    randn(Float32, (1, 1)),
)

opt = Flux.setup(Adam(FT(0.005), FT.((0.9, 0.99)), FT(1e-8)), mRT)
loss(
    model,
    x,
    obs;
    μ = μ,
    LAI = LAI,
    α_PAR_g = α_PAR_g,
    α_NIR_g = α_NIR_g,
    SW_d = SW_d,
) = Flux.mse(model(x, LAI, μ, α_PAR_g, α_NIR_g, SW_d), obs)
obs = SW_u ./ SW_d
@show loss(mRT, x̂, obs)

for i in 1:2000
    Flux.train!(loss, mRT, [(x̂, obs)], opt)
    if i % 10 == 0
        @show loss(mRT, x̂, obs)
    end
end
    
pred = mRT(x̂, LAI, μ, α_PAR_g, α_NIR_g, SW_d);

@show compute_r2(obs, pred)
@show sqrt(loss(mRT, x̂, obs))
@show mean(pred[:] .- obs[:])

sites = unique(site_ID)
errors = []
lais = []
for site in sites
    z = site_ID .== site .&& obs[:] .> 10
    push!(errors, median(((pred .- obs) ./ (obs .+ eps(FT)))[z]))
    push!(lais, maximum(LAI[:][z]))
end
