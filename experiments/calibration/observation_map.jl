using ClimaAnalysis
using Statistics
import EnsembleKalmanProcesses as EKP

"""
    process_member_data(simdir, iteration, m, training_locations)

Process data from a simulation directory for a specific ensemble member.

# Arguments
- `simdir`: Simulation directory containing the data
- `iteration`: Current iteration of the calibration
- `m`: Ensemble member index
- `training_locations`: Coordinates of locations for the loss function

# Returns
- Array of processed observations including seasonal averages of lhf, shf, swu, and lwu
  at all training locations
"""
function process_member_data(simdir, iteration, m, training_locations)
    variables = ["lhf", "shf", "swu", "lwu"]
    day_in_seconds = 86400
    first_december = 11 * 30 * day_in_seconds
    year_in_seconds = 360 * day_in_seconds
    last_december = first_december + year_in_seconds - 15*day_in_seconds
    # TO DO: update ClimaAnalysis to get Date dimension to make this better
    data = Dict(var => window(shift_to_start_of_previous_month(get(simdir; short_name = var)), "time", left = first_december, right = last_december) for var in variables)

    obs = []

    for (lon, lat) in training_locations
        for varrr in varrriables
            loc_data = ClimaAnalysis.slice(data[varrr], lon = lon, lat = lat)
            seasonal_data = split_by_season_across_time(loc_data)
            # returns a vector of CA varrrs split by season, for each year
            # It starts in DJF (1 element)
            # It has a length of 9 (4 seasons * 2 years + 1 because 1st DJF has 1 month, last DJF has 2)
            seasonal = [mean(seasonal_data[i].data) for i in 1:length(seasonal_data)]


# OLD code below, with manual indexing instead of ClimaAnalysis function split_by_season_across_time

#            loc_data = ClimaAnalysis.slice(data[varrr], lon = lon, lat = lat).data;
#            seasonal =
#                [mean(loc_data[i:(i + 2)]) for i in 1:3:(length(loc_data) - 2)][5:8]
            push!(obs, seasonal)
        end
    end

    return vcat(obs...)
end

function CAL.observation_map(iteration)
    single_member_dims = (l_obs,)
    G_ensemble = Array{Float64}(undef, single_member_dims..., ensemble_size)

    for m in 1:ensemble_size
        member_path = path_to_ensemble_member(caldir, iteration, m)
        simdir_path =
            joinpath(member_path, "global_diagnostics", "output_active")
        if isdir(simdir_path)
            simdir = SimDir(simdir_path)
            G_ensemble[:, m] .=
                process_member_data(simdir, iteration, m, training_locations)
        else
            G_ensemble[:, m] .= NaN
        end
    end

    return G_ensemble
end
