using ClimaAnalysis
using Statistics
import EnsembleKalmanProcesses as EKP

function process_member_data(simdir, iteration, m)

    lhf = get(simdir; short_name = "lhf")
    shf = get(simdir; short_name = "shf")
    swu = get(simdir; short_name = "swu")
    lwu = get(simdir; short_name = "lwu")

    # Initialize an empty list to store observations
    obs = []
    # Loop over each location
    for (lon, lat) in training_locations # training_locations
        lhf_loc = ClimaAnalysis.slice(lhf, lon = lon, lat = lat)
        shf_loc = ClimaAnalysis.slice(shf, lon = lon, lat = lat)
        swu_loc = ClimaAnalysis.slice(swu, lon = lon, lat = lat)
        lwu_loc = ClimaAnalysis.slice(lwu, lon = lon, lat = lat)

        lhf_seasonal = [
            mean(lhf_loc.data[i:(i + 2)]) for
            i in 1:3:(length(lhf_loc.data) - 2)
        ][5:8] # 5:8 for seasons of 2nd year
        shf_seasonal = [
            mean(shf_loc.data[i:(i + 2)]) for
            i in 1:3:(length(shf_loc.data) - 2)
        ][5:8]
        swu_seasonal = [
            mean(swu_loc.data[i:(i + 2)]) for
            i in 1:3:(length(swu_loc.data) - 2)
        ][5:8]
        lwu_seasonal = [
            mean(lwu_loc.data[i:(i + 2)]) for
            i in 1:3:(length(lwu_loc.data) - 2)
        ][5:8]

        push!(obs, lhf_seasonal)
        push!(obs, shf_seasonal)
        push!(obs, swu_seasonal)
        push!(obs, lwu_seasonal)
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
            G_ensemble[:, m] .= process_member_data(simdir, iteration, m)
        else
            G_ensemble[:, m] .= NaN
        end
    end

    return G_ensemble
end
