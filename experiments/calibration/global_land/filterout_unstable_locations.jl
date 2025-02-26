caldir = "calibration_output"

ensemble_size = 6
n_iterations = 2
n_locations = 500 # number of random locations
nelements = (50, 10) # resolution for model and era5
l_obs = 10*4*n_locations # 10 months * 4 variables * n_locations

using ClimaLand
dir = pkgdir(ClimaLand)

include(joinpath(dir, "experiments/calibration/shared/rand_locations.jl"))
FT = Float64
using ClimaUtilities
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
using ClimaCore
include(joinpath(dir, "experiments/calibration/shared/make_training_locations.jl"))

using ClimaAnalysis
import EnsembleKalmanProcesses as EKP

import ClimaCalibrate:
    forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL

using Statistics
#nonan_training_locations = SVector(deleteat!(collect(training_locations), rows_to_remove)...)
#rows_to_remove = [41, 59, 70] # iteration 0, most in member 5. member 1 has no nan.
#rows_to_remove = [66] # iteration 1, member 5

# explain things here - tutorial style
function make_rows_to_remove(simdir)

    lhf = get(simdir; short_name = "lhf")
    shf = get(simdir; short_name = "shf")
    lwu = get(simdir; short_name = "lwu")
    swu = get(simdir; short_name = "swu")

    rows_to_remove=[]
    # Initialize an empty list to store observations
    obs_list = []
    i=1
    # Loop over each location
    for (lon, lat) in training_locations # training_locations
        # Slice lhf and shf at the given longitude and latitude
        lhf_loc = ClimaAnalysis.slice(lhf, lon = lon, lat = lat)
        shf_loc = ClimaAnalysis.slice(shf, lon = lon, lat = lat)
        lwu_loc = ClimaAnalysis.slice(lwu, lon = lon, lat = lat)
        swu_loc = ClimaAnalysis.slice(swu, lon = lon, lat = lat)

        # Create Observation objects for lhf and shf
        lhf_obs = EKP.Observation(
            Dict(
                 "samples" => lhf_loc.data[3:12],
                "covariances" => cov(lhf_loc.data[3:12]) * EKP.I,
                "names" => "lhf_$(lon)_$(lat)",
            ),
        )
        shf_obs = EKP.Observation(
            Dict(
                "samples" => shf_loc.data[3:12],
                "covariances" => cov(shf_loc.data[3:12]) * EKP.I,
                "names" => "shf_$(lon)_$(lat)",
            ),
        )

        lwu_obs = EKP.Observation(
            Dict(
                "samples" => lwu_loc.data[3:12],
                "covariances" => cov(lwu_loc.data[3:12]) * EKP.I,
                "names" => "lwu_$(lon)_$(lat)",
            ),
        )

        swu_obs = EKP.Observation(
            Dict(
                "samples" => swu_loc.data[3:12],
                "covariances" => cov(swu_loc.data[3:12]) * EKP.I,
                "names" => "swu_$(lon)_$(lat)",
            ),
        )
        # Add the observations to the list
        #if all(x -> all(!isnan, x), [lhf_loc.data, shf_loc.data, lwu_loc.data, swu_loc.data])
            push!(obs_list, lhf_obs)
            push!(obs_list, shf_obs)
            push!(obs_list, lwu_obs)
            push!(obs_list, swu_obs)
        #end

        if all(x -> all(!isnan, x), [lhf_loc.data, shf_loc.data, lwu_loc.data, swu_loc.data])
            println("no nan", ", m = ", m, ", i = ", i)
        else
            println("some nans", ", m = ", m, ", i = ", i)
            push!(rows_to_remove, i)
        end
        i+=1
    end

    # Combine all observations into a single observation
    full_obs = EKP.combine_observations(obs_list)
    obs = EKP.get_obs(full_obs)

    return rows_to_remove
end

function remove_unstable_locations!(iteration, training_locations)
    #single_member_dims = (l_obs,)
    #G_ensemble = Array{Float64}(undef, single_member_dims..., ensemble_size)

#    rows_to_remove = []
    for m in 1:ensemble_size
        member_path = path_to_ensemble_member(caldir, iteration, m)
        simdir_path =
            joinpath(member_path, "global_diagnostics", "output_active")
        if isdir(simdir_path)
            simdir = SimDir(simdir_path)
            rows_to_remove = make_rows_to_remove(simdir)
            println(m)
            training_locations = SVector(deleteat!(collect(training_locations), rows_to_remove)...)
        end
    end

#    rows_to_remove = unique(rows_to_remove)
    # now in G_ensemble, remove the locations (12 rows for each index of rows_to_remove)
    # then remove those for observation_data too
    return training_locations
end

for i = 0:n_iterations-1
    training_locations = remove_unstable_locations!(i, training_locations)
end

using StaticArrays

# Define the output file name
filename = "stable_training_locations.jl"

# Open the file and write the SVector definition
open(filename, "w") do io
    println(io, "using StaticArrays\n")
    println(io, "const training_locations = SVector{$(length(training_locations)), Tuple{Float64, Float64}}(")

    # Write each tuple with a comma at the end
    for elem in training_locations
        println(io, "    $(elem),")
    end

    println(io, ")")  # Close the SVector
end

println("Vector written to $filename")

