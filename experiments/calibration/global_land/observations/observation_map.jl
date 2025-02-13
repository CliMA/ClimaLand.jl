using ClimaAnalysis
import EnsembleKalmanProcesses as EKP

# explain things here - tutorial style
function process_member_data(simdir)

    lhf = get(simdir; short_name = "lhf")
    shf = get(simdir; short_name = "shf")
    lwu = get(simdir; short_name = "lwu")
    swu = get(simdir; short_name = "swu")

    # Initialize an empty list to store observations
    obs_list = []
    # Loop over each location
    for (lon, lat) in locations
        # Slice lhf and shf at the given longitude and latitude
        lhf_loc = ClimaAnalysis.slice(lhf, lon = lon, lat = lat)
        shf_loc = ClimaAnalysis.slice(shf, lon = lon, lat = lat)
        lwu_loc = ClimaAnalysis.slice(lwu, lon = lon, lat = lat)
        swu_loc = ClimaAnalysis.slice(swu, lon = lon, lat = lat)

        # Create Observation objects for lhf and shf
        lhf_obs = EKP.Observation(
            Dict(
                "samples" => lhf_loc.data,
                "covariances" => cov(lhf_loc.data) * EKP.I,
                "names" => "lhf_$(lon)_$(lat)",
            ),
        )
        shf_obs = EKP.Observation(
            Dict(
                "samples" => shf_loc.data,
                "covariances" => cov(shf_loc.data) * EKP.I,
                "names" => "shf_$(lon)_$(lat)",
            ),
        )

        lwu_obs = EKP.Observation(
            Dict(
                "samples" => lwu_loc.data,
                "covariances" => cov(lwu_loc.data) * EKP.I,
                "names" => "lwu_$(lon)_$(lat)",
            ),
        )

        swu_obs = EKP.Observation(
            Dict(
                "samples" => swu_loc.data,
                "covariances" => cov(swu_loc.data) * EKP.I,
                "names" => "swu_$(lon)_$(lat)",
            ),
        )
        # Add the observations to the list
        push!(obs_list, lhf_obs)
        push!(obs_list, shf_obs)
        push!(obs_list, lwu_obs)
        push!(obs_list, swu_obs)
    end

    # Combine all observations into a single observation
    full_obs = EKP.combine_observations(obs_list)
    obs = EKP.get_obs(full_obs)
    return obs
end

function CAL.observation_map(iteration)
    single_member_dims = (1200,) # 4 variables * 12 months * 25 locations
    G_ensemble = Array{Float64}(undef, single_member_dims..., ensemble_size)

    for m in 1:ensemble_size
        member_path = path_to_ensemble_member(caldir, iteration, m)
        simdir_path =
            joinpath(member_path, "global_diagnostics", "output_active")
        if isdir(simdir_path)
            simdir = SimDir(simdir_path)
            G_ensemble[:, m] .= process_member_data(simdir)
        else
            G_ensemble[:, m] .= NaN
        end
    end
    return G_ensemble
end
