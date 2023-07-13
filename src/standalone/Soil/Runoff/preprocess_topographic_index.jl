using NCDatasets
using JLD2
#using Distributions
include("distributions.jl")
function main(; resolution, dist_name)

    data = NCDataset(
        "/Users/katherinedeck/Downloads/6b0c4358-2bf3-4924-aa8f-793d468b92be(1)/ga2.nc",
    )
    lat = data["lat"][:]
    lon = data["lon"][:]

    (lat_min, lat_max) = extrema(lat)
    (lon_min, lon_max) = extrema(lon)

    lat_count = Int(ceil((lat_max - lat_min) / resolution)) + 1
    lon_count = Int(ceil((lon_max - lon_min) / resolution)) + 1
    if dist_name == "inv_gamma"
        dist_type = Distribtuions.InverInvGammaDistribution()
    elseif dist_name == "gamma"
        dist_type = GammaDistribution()
    elseif dist_name == "percentile"
        dist_type = PercentileDistribution()
    elseif dist_name == "lognormal"
        dist_type = LogNormalDistribution()
    elseif dist_name == "frechet"
        dist_type = FrechetDistribution()
    end

    parameters = zeros((lat_count, lon_count, 7 + n_params(dist_type)))

    for lat_id in 1:1:lat_count
        @info lat_id / lat_count
        for lon_id in 1:1:lon_count
            lat_indices =
                (lat .>= lat_min + resolution * (lat_id - 1)) .&
                (lat .< lat_min + resolution * lat_id)
            lon_indices =
                (lon .>= lon_min + resolution * (lon_id - 1)) .&
                (lon .< lon_min + resolution * lon_id)
            ϕ = data["Band1"][lon_indices, lat_indices][:]
            present = .~(typeof.(ϕ) .<: Missing) .& (ϕ .> 0)
            present_count = sum(present)
            zeros = .~(typeof.(ϕ) .<: Missing) .& (ϕ .== 0)
            

            if present_count > 100
                x̄, var, params = fit(dist_type, Float32.(ϕ[present]))
                if ~isnan(params[1])
                    
                    parameters[lat_id, lon_id, 1] = x̄
                    parameters[lat_id, lon_id, 2] = var
                    parameters[lat_id, lon_id, 3] = mean(lat[lat_indices])
                    parameters[lat_id, lon_id, 4] = mean(lon[lon_indices])
                    parameters[lat_id, lon_id, 5] =
                        log_likelihood(dist_type, ϕ[present], params)
                    parameters[lat_id, lon_id, 6] = present_count
                    parameters[lat_id, lon_id, 7] = sum(zeros)
                    
                    parameters[lat_id, lon_id, 8:(7 + n_params(dist_type))] .=
                        params
                end
            end
        end
    end

    save("./$(dist_name)_$resolution.jld2", "parameters", parameters)
    #=
    grid_lat = zeros(lat_count)
    grid_lon = zeros(lon_count)
    for lat_id in 1:1:lat_count
        lat_indices = (lat .>= lat_min + resolution * (lat_id-1)) .& (lat .< lat_min + resolution * lat_id)
        grid_lat[lat_id] = mean(lat[lat_indices])
    end

    for lon_id in 1:1:lon_count
        lon_indices = (lon .>= lon_min + resolution * (lon_id-1)) .& (lon .< lon_min + resolution * lon_id)
        grid_lon[lon_id] = mean(lon[lon_indices])
    end

    Plots.heatmap(grid_lon[1:end-1], grid_lat[1:end-1], parameters[1:end-1,1:end-1,3])
    Plots.heatmap(grid_lon[1:end-1], grid_lat[1:end-1], parameters[1:end-1,1:end-1,4], clims = (0,50))
    Plots.savefig("alpha_inv_gamma.png")
    Plots.heatmap(grid_lon[1:end-1], grid_lat[1:end-1], parameters[1:end-1,1:end-1,5], clims = (0,300))
    Plots.savefig("beta_inv_gamma.png")
    =#
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(resolution = parse(Float64, ARGS[1]), dist_name = ARGS[2])
end
