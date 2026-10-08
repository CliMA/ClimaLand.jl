using Statistics
using NCDatasets
using Insolation
using Interpolations
import ClimaLand.Parameters as LP
using ClimaLand:Artifacts
using Dates
using Flux
using CairoMakie
using ClimaAnalysis
using Functors

FT = Float32
toml_dict = LP.create_toml_dict(FT)
insol_params = LP.LandParameters(toml_dict).insol_params
data_path = Artifacts.optimal_lai_initial_conditions_path()
lai_data_path = Artifacts.modis_max_lai_data_path()
topo_data_path = "topographic_data.nc"
tair_data_path = "../ClimaArtifacts/optimal_lai_inputs/tair_gs.nc"
#Read in base data and regrid the others to that resolution
ds = NCDataset(lai_data_path)
maxLAI = ds["lai"][:,:]
lat = ds["lat"][:]
lon = ds["lon"][:]
close(ds)


ds = NCDataset(data_path)
gsl_env = ds["gsl"][:,:]
gs_vpd_env = ds["vpd_gs"][:,:]
ann_p_env = ds["precip_annual"][:,:]
ds_tair = NCDataset(tair_data_path)
tair_env = ds_tair["tair"][:,:]
lat_env = ds["lat"][:]
lon_env = ds["lon"][:]
close(ds)
close(ds_tair)

ds = NCDataset(topo_data_path)
μ_topo = ds["mean_slope"][:,:]
elevation_topo = ds["norm_stdz"][:,:]
lat_topo = ds["lat"][:]
lon_topo = ds["lon"][:]
close(ds)

# Below I get 2001 data for LAI and albedo
albedo_path = "/Users/katherinedeck/Downloads/albedo.nc"
α_var = ClimaAnalysis.Var.read_var(albedo_path; short_name = "albedo")
α_var = shift_longitude(α_var, -180.0, 180.0)
α = α_var.data[:,:,11:22];
times = α_var.dims["time"][11:22];
@assert α_var.dims["lon"][:] == lon
@assert α_var.dims["lat"][:] == lat
#Lai times do not coinicide exactly with albedo times. can refine later
ds = NCDataset(Artifacts.modis_lai_single_year_path(;year = Dates.year(DateTime(2001))))
lai = ds["lai"][:,:,:];
lai_time = ds["time"][:];
@assert ds["lat"][:] == lat
@assert ds["lon"][:] == lon
close(ds)

gridded_longitude = repeat(Float32.(lon),1,180);
gridded_latitude = transpose(repeat(Float32.(lat),1,360));
μ_itp = Interpolations.linear_interpolation((lon_topo,lat_topo), μ_topo;extrapolation_bc = Interpolations.Flat())
μ = μ_itp.(gridded_longitude, gridded_latitude)
elev_itp = Interpolations.linear_interpolation((lon_topo,lat_topo), elevation_topo;extrapolation_bc = Interpolations.Flat())
elevation = elev_itp.(gridded_longitude, gridded_latitude)

gsl_itp = Interpolations.linear_interpolation((lon_env,lat_env), gsl_env;extrapolation_bc = Interpolations.Flat())
gsl = gsl_itp.(gridded_longitude, gridded_latitude)
gs_vpd_itp = Interpolations.linear_interpolation((lon_env,lat_env), gs_vpd_env;extrapolation_bc = Interpolations.Flat())
gs_vpd = gs_vpd_itp.(gridded_longitude, gridded_latitude)
gs_tair_itp = Interpolations.linear_interpolation((lon_env,lat_env), tair_env;extrapolation_bc = Interpolations.Flat())
gs_tair = gs_tair_itp.(gridded_longitude, gridded_latitude)
ann_p_itp = Interpolations.linear_interpolation((lon_env,lat_env), ann_p_env;extrapolation_bc = Interpolations.Flat())
ann_p = ann_p_itp.(gridded_longitude, gridded_latitude)


zenith_only = (args...) -> Insolation.insolation(args...).μ

# Restrict model to areas with LAI > 1 make sure no ocean points accidentally included
id = 6
mask = mask = (.~isnan.(gsl)) .&& (lai[:,:,id] .> 2) .&& (α[:,:,id] .< 0.7) .&& (gridded_latitude .>-50)
# Normalize static data
x = permutedims(cat(maxLAI[mask], gsl[mask], gs_tair[mask], log10.(gs_vpd[mask]), elevation[mask], log10.(ann_p[mask]), log10.(μ[mask]), dims = 2),(2,1));
x̄ = mean(x, dims = (2))
σ = std(x, dims = (2))
x̂ = @. Float32.((x-x̄)/σ)
xLAI = Float32.(reshape(lai[mask,id], (1,sum(mask)))) ;
obs = Float32.(reshape(α[mask, id], (1, sum(mask))))

# Create data batches
date = times[id]
correct_format_date = DateTime(year(date)) + Month(month(date)-1) + Day(day(date)-1)
hod = correct_format_date .+ Hour.(0:2:22)
cosθs = zeros(360,180, length(hod));
for i in 1:length(hod)
    cosθs[:,:,i] .= zenith_only.(hod[i], gridded_latitude, gridded_longitude, insol_params)
end
cosθs = cosθs[mask, :];
xcos = reshape(cosθs, (1,sum(mask),length(hod)));

struct LandRT2
    GW
    Gb
    Aw
    Ac
    Ab
end
Flux.@layer LandRT2

function (m::LandRT2)(x̂,xcos, xLAI, α_soil)
    χ = sigmoid(m.GW * x̂ .+ m.Gb) .- Float32(0.4) # restricted to -0.4, 0.6
    ϕ1 = Float32(0.5) .- Float32(0.633) .*χ  .+ Float32(0.33) .* χ.^2
    ϕ2 = Float32(0.877) .*( 1 .- 2 .* ϕ1)
    
    G = ϕ1 .+ ϕ2 .* max.(xcos, Float32(0.001))
    K = min.(G ./ max.(xcos, Float32(0.01)), 1e3)
    α =  sigmoid(m.Aw * x̂ .+ m.Ab .+ m.Ac .* max.(xcos, Float32(0.001)))
    transmitted_fraction = @. exp(-K * xLAI);
    reflected_downwards_pass = @. α * (1 - transmitted_fraction);
    upwelling_from_soil = @. α_soil * transmitted_fraction;
    reflected_upwards_pass =
        @. upwelling_from_soil * (1 - transmitted_fraction) * α;
    transmitted_upwards_pass = @. upwelling_from_soil * transmitted_fraction;
    upwelling_from_land = @. reflected_downwards_pass + reflected_upwards_pass + transmitted_upwards_pass;
    return mean(upwelling_from_land, dims = 3)[:,:,1]
end

mRT = LandRT2(zeros(Float32, (1,7)), zeros(Float32, 1),  zeros(Float32, (1,7)),zeros(Float32, 1) .+0.01, zeros(Float32, 1) .+0.5)
opt = Flux.setup(Adam(FT(0.005),FT.((0.9,0.99)),FT(1e-8)), mRT)
loss(model, x, obs; xcos=xcos, xLAI=xLAI,α_soil = α_soil) = Flux.mse(model(x, xcos, xLAI,α_soil), obs)
@show loss(mRT, x̂, obs)

for i in 1:500
    Flux.train!(loss, mRT, [(x̂,obs)], opt)
    if i %10 == 0
        @show loss(mRT, x̂, obs)
    end 
end


#####


struct LandRTNN
    layers::NamedTuple
end
Functors.@functor LandRTNN
#Flux.trainable(m::LandRTNN) = (Flux.trainable(m.χ)..., Flux.trainable(m.α)...)


#
#Step1: Take Input and compute G, compute α in parallel
α_leaf_function = Chain(Dense(7 => 14), Dense(15 =>1, sigmoid))
G_function = Chain(Parallel(-, Dense(7 => 14), Dense(14 =>1, sigmoid), Float32(0.4)),
                   Parallel(+, ϕ1, 
step1 = Parallel(land_albedo, G_function, α_leaf_function)
    
function (m::LandRTNN)(x̂,xcos, xLAI, α_soil)
    χ = mm.layers.χ(x̂) .- Float32(0.4) # restricted to -0.4, 0.6
    ϕ1 = Float32(0.5) .- Float32(0.633) .*χ  .+ Float32(0.33) .* χ.^2
    ϕ2 = Float32(0.877) .*( 1 .- 2 .* ϕ1)
    
    G = ϕ1 .+ ϕ2 .* max.(xcos, Float32(0.001))
    K = min.(G ./ max.(xcos, Float32(0.01)), 1e3)
    α = mm.layers.α(x̂)
    transmitted_fraction = @. exp(-K * xLAI);
    reflected_downwards_pass = @. α * (1 - transmitted_fraction);
    upwelling_from_soil = @. α_soil * transmitted_fraction;
    reflected_upwards_pass =
        @. upwelling_from_soil * (1 - transmitted_fraction) * α;
    transmitted_upwards_pass = @. upwelling_from_soil * transmitted_fraction;
    upwelling_from_land = @. reflected_downwards_pass + reflected_upwards_pass + transmitted_upwards_pass;
    return mean(upwelling_from_land, dims = 3)[:,:,1]
end
c1 = Chain(Dense(7 => 14), Dense(14 => 1, sigmoid));
c2 = Chain(Dense(7 => 14), Dense(14 => 1, sigmoid));
layers = (; χ = c1, α = c2);
mm = LandRTNN(layers);

opt = Flux.setup(Adam(FT(0.005),FT.((0.9,0.99)),FT(1e-8)), mm)
loss(model, x, obs; xcos=xcos, xLAI=xLAI,α_soil = α_soil) = Flux.mse(model(x, xcos, xLAI,α_soil), obs)
@show loss(mm, x̂, obs)
for i in 1:500
    Flux.train!(loss, mm, [(x̂,obs)], opt)
    if i %10 == 0
        @show loss(mm, x̂, obs)
    end 
end
