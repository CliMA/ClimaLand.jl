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
    data = Dict(var => get(simdir; short_name = var) for var in variables)
    obs = []

    for (lon, lat) in training_locations
        for var in variables
            loc_data = ClimaAnalysis.slice(data[var], lon = lon, lat = lat).data
            seasonal =
                [mean(loc_data[i:(i + 2)]) for i in 1:3:(length(loc_data) - 2)][5:8]
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
