#using NCDatasets
using ClimaAnalysis
using ClimaLand
using ClimaLand.Artifacts
using ClimaCore
using Statistics
import EnsembleKalmanProcesses as EKP
using LinearAlgebra
#using Dates

#era5_ds = Dataset(
#    joinpath(
#        ClimaLand.Artifacts.era5_surface_data2008_path(),
#        "era5_monthly_surface_fluxes_200801-200812.nc",
#    ),
#)

era5_path = joinpath(
                     ClimaLand.Artifacts.era5_surface_data_1979_2024_path(),
                     "era5_monthly_averages_surface_single_level_197901-202410.nc",
                    )

# Load data
lhf = OutputVar(era5_path, "mslhf")
shf = OutputVar(era5_path, "msshf")
swu = OutputVar(era5_path, "msuwswrf")

land_mask_var = ClimaAnalysis.apply_landmask(lhf)
lhf.data[isnan.(land_mask_var.data)]

simdir = SimDir("calibration_output_utki_sample/iteration_000/member_001/global_diagnostics/output_active/")
lhf_out = get(simdir, "lhf")

lhf_on_diagnostic_grid = ClimaAnalysis.resampled_as(lhf, lhf_out)
land_mask_var = slice(ClimaAnalysis.apply_landmask(lhf_on_diagnostic_grid), time = 0)

training_locations = []
for (i, lon) in enumerate(ClimaAnalysis.longitudes(land_mask_var))
    for (j, lat) in enumerate(ClimaAnalysis.latitudes(land_mask_var))
        if isnan(land_mask_var.data[i,j]) && -60 <= lat <= 60
            push!(training_locations, (lon, lat))
        end
    end
end

# length of training_locations is 7053

# lhf_vector is (lon, lat)
# ^ this is new training_locations
# need to run a new forward_model first, just for 1 or 2 months, to retrieve the grid

#lhf_era5 = lhf_on_diagnostic_grid.data[isnan.(land_mask_var.data)]


# Get all land coordinates:
#domain = ClimaLand.global_domain(FT; nelements = nelements)
#surface_space = domain.space.surface
#mask = ClimaLand.landsea_mask(surface_space)
#coords = ClimaCore.Fields.coordinate_field(surface_space) # Field
#lat = coords.lat # Field
#lon = coords.long # Field
#lat = parent(lat)[parent(mask) .== 1]
#lon = parent(lon)[parent(mask) .== 1]
##training_locations = SVector{length(lon)}([(lon[i], lat[i]) for i in eachindex(lon)]...)
##training_locations = [SVector(lon[i], lat[i]) for i in eachindex(lon)]
##training_locations = collect(zip(lon, lat))
#training_locations = collect(zip(Array(lon), Array(lat)))
#
## Load data
#lhf = OutputVar(era5_path, "mslhf")
#shf = OutputVar(era5_path, "msshf")
#swu = OutputVar(era5_path, "msuwswrf")

# Get slices (time series for each location) and variance
lhf_slices = []
shf_slices = []
swu_slices = []
lhf_vars = []
shf_vars = []
swu_vars = []
for (lon, lat) in training_locations
        # Initialize temporary storage for each location
    lhf_var = Vector{Vector{FT}}(undef, 12)
    shf_var = Vector{Vector{FT}}(undef, 12)
    swu_var = Vector{Vector{FT}}(undef, 12)
    lhf_slice = slice(lhf, longitude = lon, latitude = lat) # 1 slice has 550 data (~ 45years*12months)
    shf_slice = slice(shf, longitude = lon, latitude = lat)
    swu_slice = slice(swu, longitude = lon, latitude = lat)
    lhf_seasonal = [mean(lhf_slice.data[i:i+2]) for i in 1:3:(length(lhf_slice.data)-2)]
    shf_seasonal = [mean(shf_slice.data[i:i+2]) for i in 1:3:(length(shf_slice.data)-2)]
    swu_seasonal = [mean(swu_slice.data[i:i+2]) for i in 1:3:(length(swu_slice.data)-2)]
    lhf_var = [FT(var(lhf_slice.data[i:4:end])) for i in 1:4]
    shf_var = [FT(var(shf_slice.data[i:4:end])) for i in 1:4]
    swu_var = [FT(var(swu_slice.data[i:4:end])) for i in 1:4]
    push!(lhf_slices, FT.(.-lhf_seasonal))
    push!(shf_slices, FT.(.-shf_seasonal))
    push!(swu_slices, FT.(swu_seasonal))
    push!(lhf_vars, lhf_var)
    push!(shf_vars, shf_var)
    push!(swu_vars, swu_var)
end

lhf_q = quantile(vcat(lhf_vars...), [0.05, 0.99])
shf_q = quantile(vcat(shf_vars...), [0.05, 0.99])
swu_q = quantile(vcat(swu_vars...), [0.2, 0.99])

if any([lhf_q[2]/lhf_q[1] > 1e8,shf_q[2]/shf_q[1] > 1e8,swu_q[2]/swu_q[1] > 1e8])
    throw(ArgumentError("quantiles of noise variance exceed 1e8, may cause unstable update \n investigate and adjust quantiles."))
end

for i = 1:length(lhf_vars)
    lhf_vars[i] = min.(max.(lhf_vars[i], lhf_q[1]), lhf_q[2])
    shf_vars[i] = min.(max.(shf_vars[i], shf_q[1]), shf_q[2])
    swu_vars[i] = min.(max.(swu_vars[i], swu_q[1]), swu_q[2])
end

obs_y = []
n_samples = 10
for y in 31:31+n_samples # 31 to start in 2009, see below
    obs_y_temp = []
    for (i, (lon, lat)) in enumerate(training_locations)
        lhf_obs = EKP.Observation(
                                  Dict(
                                       "samples" => lhf_slices[i][1+(4*(y-1)):4+(4*(y-1))],
                                       "covariances" => Diagonal(lhf_vars[i]),
                                       "names" => "lhf_$(lon)_$(lat)_$(y)",
                                      ),
                                 )
        shf_obs = EKP.Observation(
                                  Dict(
                                       "samples" => shf_slices[i][1+(4*(y-1)):4+(4*(y-1))],
                                       "covariances" => Diagonal(shf_vars[i]),
                                       "names" => "shf_$(lon)_$(lat)_$(y)",
                                      ),
                                 )
        swu_obs = EKP.Observation(
                                  Dict(
                                       "samples" => swu_slices[i][1+(4*(y-1)):4+(4*(y-1))],
                                       "covariances" => Diagonal(swu_vars[i]),
                                       "names" => "swu_$(lon)_$(lat)_$(y)",
                                      ),
                                 )
        push!(obs_y_temp, EKP.combine_observations([lhf_obs, shf_obs, swu_obs]))
    end
    push!(obs_y, EKP.combine_observations(obs_y_temp))
end

m_size = 1 # 1 year
given_batches = [collect(((i - 1) * m_size + 1):(i * m_size)) for i in 1:n_samples]
minibatcher = EKP.FixedMinibatcher(given_batches)
o_names = ["$i" for i in 2009:2009+n_samples-1]

observationseries =  EKP.ObservationSeries(obs_y, minibatcher, o_names)

#= we need to start in 2009, so at year = 31
# Model starts in 2008 but discards first year
for i = 1:45
    println(i, " ", 1978+i)
end
=#


#=


"""
    era5_monthly_stdevs(varname::String)
Returns the year over year standard deviation of a given variable from ERA5.
i.e, a 12 element vector of the standard deviation of the variable for each
month taken across all data years. Assumes you're looking for a variable
contained in the era5_monthly_averages_surface_single_level_197901-202410.nc
data file.
"""
function era5_monthly_stdevs(varname::String, lat, lon)
    # Read the data for this variable from all data years

    # Create a static array to store the stddevs in
    stdevs = zeros(FT, 12)

    # Find the nearest latitude and longitude to the target location in the data
    lat_ind = argmin(abs.(era5_all_ds["latitude"][:] .- lat))
    lon_ind = argmin(abs.(era5_all_ds["longitude"][:] .- lon))

    for i in 1:12
        monthly_data = era5_all_ds[varname][
            lon_ind,
            lat_ind,
            findall((x) -> month(x) == i, era5_all_ds["time"]),
        ]
        stdevs[i] = std(monthly_data)
    end

    return stdevs
end

# Make the ERA5 target
ERA5_target = zeros(FT, l_obs*3)
ERA5_std = zeros(FT, l_obs*3)
close_location = (x, y) -> abs(x - y) < 7e-1
for i in 1:length(training_locations)
    current_loc = training_locations[i]
    lat = current_loc[2]
    lon = current_loc[1]
#for (lon, lat) in training_locations[1:5] # crashes if too many, not sure why. MVector?
    # Fetch slices of lhf and shf era5 data from the era5 dataset
    lat_ind, lon_ind =
        findall((x) -> close_location(x, lat), era5_ds["latitude"][:])[1],
        findall((x) -> close_location(x, lon + 180), era5_ds["longitude"][:])[1]
    lhf_loc = vec(era5_ds["mslhf"][lon_ind, lat_ind, :][:, end, :])
    shf_loc = vec(era5_ds["msshf"][lon_ind, lat_ind, :][:, end, :])
    lwu_loc = vec(era5_ds["msuwlwrf"][lon_ind, lat_ind, :][:, end, :])
    swu_loc = vec(era5_ds["msuwswrf"][lon_ind, lat_ind, :][:, end, :])
    lhf_ERA5_std = era5_monthly_stdevs("mslhf", lat, lon)
    shf_ERA5_std = era5_monthly_stdevs("msshf", lat, lon)
    lwu_ERA5_std = era5_monthly_stdevs("msuwlwrf", lat, lon)
    swu_ERA5_std = era5_monthly_stdevs("msuwswrf", lat, lon)

    # Add the observations to the target
    loc_id = (i-1) * 12 * 4
    ERA5_target[loc_id+1:loc_id+12] .= .-lhf_loc[1:12]
    ERA5_target[loc_id+12+1:loc_id+24] .= .-shf_loc[1:12]
    ERA5_target[loc_id+24+1:loc_id+36] .= lwu_loc[1:12]
    ERA5_target[loc_id+36+1:loc_id+48] .= swu_loc[1:12]

    ERA5_std[loc_id+1:loc_id+12] .= lhf_ERA5_std[1:12]
    ERA5_std[loc_id+12+1:loc_id+24] .= shf_ERA5_std[1:12]
    ERA5_std[loc_id+24+1:loc_id+36] .= lwu_ERA5_std[1:12]
    ERA5_std[loc_id+36+1:loc_id+48] .= swu_ERA5_std[1:12]
end

#full_obs_era5 = EKP.combine_observations(ERA5_target)
#observations = EKP.get_obs(full_obs_era5)
observations = mean(reshape(ERA5_target, 3, :), dims=1)[:]
noise_era5 = ERA5_std
noise_era5 = mean(reshape(noise_era5, 3, :), dims=1)[:]

=#

#noise_era5 = Vector(reduce(vcat, noise_era5)) # not sure why MVector

# to do, use climaanalysis instead of netcdf,
# make sure same grid is used
# and timestamps align
#=

filepath =joinpath(
           ClimaLand.Artifacts.era5_surface_data2008_path(),
           "era5_monthly_surface_fluxes_200801-200812.nc",
       )
data = OutputVar(filepath)
# lat and lon needs to be defined the same for model and era5
# this can be done with the function below (but it needs to be fixed)
# center_longitude!(data, 180)
slice(data, latitude = 0, longitude = 45)


# to ensure observation (era5) and model align in time, use the same approach as
# https://github.com/CliMA/ClimaLand.jl/blob/main/experiments/long_runs/leaderboard/data_sources.jl#L95

=#
