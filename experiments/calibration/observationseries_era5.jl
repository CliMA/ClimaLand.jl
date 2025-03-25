using ClimaAnalysis
using ClimaLand
using ClimaLand.Artifacts
using ClimaCore
using Statistics
import EnsembleKalmanProcesses as EKP
using LinearAlgebra
FT = Float64

era5_path = joinpath(
    ClimaLand.Artifacts.era5_surface_data_1979_2024_path(),
    "era5_monthly_averages_surface_single_level_197901-202410.nc",
)

# Load data
lhf = OutputVar(era5_path, "mslhf")
#shf = OutputVar(era5_path, "msshf")
swu_temp = OutputVar(era5_path, "msuwswrf")

land_mask_var = ClimaAnalysis.apply_landmask(lhf)
lhf.data[isnan.(land_mask_var.data)]

simdir = SimDir(
    "calibration_output_utki_sample/iteration_000/member_001/global_diagnostics/output_active/",
)
lhf_out = get(simdir, "lhf")
swu_out = get(simdir, "swu")

lhf_on_diagnostic_grid = ClimaAnalysis.resampled_as(lhf, lhf_out)
land_mask_var =
    slice(ClimaAnalysis.apply_landmask(lhf_on_diagnostic_grid), time = 0)

training_locations = []
for (i, lon) in enumerate(ClimaAnalysis.longitudes(land_mask_var))
    for (j, lat) in enumerate(ClimaAnalysis.latitudes(land_mask_var))
        if isnan(land_mask_var.data[i, j]) && -60 <= lat <= 60
            push!(training_locations, (lon, lat))
        end
    end
end

swu = ClimaAnalysis.resampled_as(swu_temp, swu_out, dim_names = "longitude")

# Get slices (time series for each location) and variance
#lhf_all_seasons = []
#shf_all_seasons = []
swu_all_seasons = []
#lhf_seasonal_vars = []
#shf_seasonal_vars = []
swu_seasonal_vars = []
#lhf_seasonal_means = []
#shf_seasonal_means = []
swu_seasonal_means = []
for (lon, lat) in training_locations
    # Initialize temporary storage for each location
#    lhf_var = Vector{Vector{FT}}(undef, 12)
#    shf_var = Vector{Vector{FT}}(undef, 12)
    swu_var = Vector{Vector{FT}}(undef, 12)
#    lhf_slice = slice(lhf, longitude = lon, latitude = lat) # 1 slice has 550 data (~ 45years*12months)
#    shf_slice = slice(shf, longitude = lon, latitude = lat)
    swu_slice = slice(swu, longitude = lon, latitude = lat)
#    lhf_seasonal = [
#        mean(lhf_slice.data[i:(i + 2)]) for
#        i in 1:3:(length(lhf_slice.data) - 2)
#    ]
#    shf_seasonal = [
#        mean(shf_slice.data[i:(i + 2)]) for
#        i in 1:3:(length(shf_slice.data) - 2)
#    ]
    swu_seasonal = [
        mean(swu_slice.data[i:(i + 2)]) for
        i in 1:3:(length(swu_slice.data) - 2)
    ]
#    lhf_var = [FT(var(lhf_seasonal[i:4:end])) for i in 1:4]
#    shf_var = [FT(var(shf_seasonal[i:4:end])) for i in 1:4]
    swu_var = [FT(var(swu_seasonal[i:4:end])) for i in 1:4]
#    lhf_seasonal_mean = [FT(mean(.-lhf_seasonal[i:4:end])) for i in 1:4]
#    shf_seasonal_mean = [FT(mean(.-shf_seasonal[i:4:end])) for i in 1:4]
    swu_seasonal_mean = [FT(mean(swu_seasonal[i:4:end])) for i in 1:4]
#    push!(lhf_all_seasons, FT.(.-lhf_seasonal))
#    push!(shf_all_seasons, FT.(.-shf_seasonal))
    push!(swu_all_seasons, FT.(swu_seasonal))
#    push!(lhf_seasonal_vars, lhf_var)
#    push!(shf_seasonal_vars, shf_var)
    push!(swu_seasonal_vars, swu_var)
#    push!(lhf_seasonal_means, lhf_seasonal_mean)
#    push!(shf_seasonal_means, shf_seasonal_mean)
    push!(swu_seasonal_means, swu_seasonal_mean)
end

#lhf_q = quantile(vcat(lhf_seasonal_vars...), [0.05, 0.99])
#shf_q = quantile(vcat(shf_seasonal_vars...), [0.05, 0.99])
swu_q = quantile(vcat(swu_seasonal_vars...), [0.2, 0.99])

if any([
#    lhf_q[2] / lhf_q[1] > 1e8,
#    shf_q[2] / shf_q[1] > 1e8,
    swu_q[2] / swu_q[1] > 1e8,
])
    throw(
        ArgumentError(
            "quantiles of noise variance exceed 1e8, may cause unstable update \n investigate and adjust quantiles.",
        ),
    )
end

for i in 1:length(swu_seasonal_vars)
#    lhf_seasonal_vars[i] = min.(max.(lhf_seasonal_vars[i], lhf_q[1]), lhf_q[2])
#    shf_seasonal_vars[i] = min.(max.(shf_seasonal_vars[i], shf_q[1]), shf_q[2])
    swu_seasonal_vars[i] = min.(max.(swu_seasonal_vars[i], swu_q[1]), swu_q[2])
#    lhf_seasonal_vars[i] = lhf_seasonal_vars[i] .+ lhf_seasonal_means[i].^2 .* 0.01
#    shf_seasonal_vars[i] = shf_seasonal_vars[i] .+ shf_seasonal_means[i].^2 .* 0.01
    swu_seasonal_vars[i] = swu_seasonal_vars[i] .+ swu_seasonal_means[i].^2 .* 0.01
end

obs_y = []
n_samples = 10
for y in 31:(31 + n_samples) # 31 to start in 2009, see below
    obs_y_temp = []
    for (i, (lon, lat)) in enumerate(training_locations)
#        lhf_obs = EKP.Observation(
#            Dict(
#                "samples" =>
#                    lhf_all_seasons[i][(1 + (4 * (y - 1))):(4 + (4 * (y - 1)))],
#                "covariances" => Diagonal(lhf_seasonal_vars[i]),
#                "names" => "lhf_$(lon)_$(lat)_$(y)",
#            ),
#        )
#        shf_obs = EKP.Observation(
#            Dict(
#                "samples" =>
#                    shf_all_seasons[i][(1 + (4 * (y - 1))):(4 + (4 * (y - 1)))],
#                "covariances" => Diagonal(shf_seasonal_vars[i]),
#                "names" => "shf_$(lon)_$(lat)_$(y)",
#            ),
#        )
        swu_obs = EKP.Observation(
            Dict(
                "samples" =>
                    swu_all_seasons[i][(1 + (4 * (y - 1))):(4 + (4 * (y - 1)))],
                "covariances" => Diagonal(swu_seasonal_vars[i]),
                "names" => "swu_$(lon)_$(lat)_$(y)",
            ),
        )
        push!(obs_y_temp, EKP.combine_observations([swu_obs])) # ([lhf_obs, shf_obs, swu_obs]))
    end
    push!(obs_y, EKP.combine_observations(obs_y_temp))
end

m_size = 1 # 1 year
given_batches =
    [collect(((i - 1) * m_size + 1):(i * m_size)) for i in 1:n_samples]
minibatcher = EKP.FixedMinibatcher(given_batches)
o_names = ["$i" for i in 2009:(2009 + n_samples - 1)]

observationseries = EKP.ObservationSeries(obs_y, minibatcher, o_names)
