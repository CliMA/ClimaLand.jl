using NCDatasets

era5_ds = Dataset(
    joinpath(
        ClimaLand.Artifacts.era5_surface_data2008_path(),
        "era5_monthly_surface_fluxes_200801-200812.nc",
    ),
)

era5_all_ds = Dataset(
                      joinpath(
                               ClimaLand.Artifacts.era5_surface_data_1979_2024_path(),
                               "era5_monthly_averages_surface_single_level_197901-202410.nc",
                              ),
                     )
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
ERA5_target = zeros(FT, l_obs)
ERA5_std = zeros(FT, l_obs)
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
    loc_id = (i-1) * 10 * 4
    ERA5_target[loc_id+1:loc_id+10] .= lhf_loc[3:12]
    ERA5_target[loc_id+10+1:loc_id+20] .= shf_loc[3:12]
    ERA5_target[loc_id+20+1:loc_id+30] .= lwu_loc[3:12]
    ERA5_target[loc_id+30+1:loc_id+40] .= swu_loc[3:12]

    ERA5_std[loc_id+1:loc_id+10] .= lhf_ERA5_std[3:12]
    ERA5_std[loc_id+10+1:loc_id+20] .= shf_ERA5_std[3:12]
    ERA5_std[loc_id+20+1:loc_id+30] .= lwu_ERA5_std[3:12]
    ERA5_std[loc_id+30+1:loc_id+40] .= swu_ERA5_std[3:12]
end

#full_obs_era5 = EKP.combine_observations(ERA5_target)
#observations = EKP.get_obs(full_obs_era5)
observations = ERA5_target
noise_era5 = ERA5_std
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

