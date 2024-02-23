using NCDatasets
using Statistics
function main(; resolution)

    data = NCDataset(
        "/Users/katherinedeck/Downloads/6b0c4358-2bf3-4924-aa8f-793d468b92be(1)/ga2.nc",
    )
    lat = data["lat"][:]
    lon = data["lon"][:]

    (lat_min, lat_max) = extrema(lat)
    (lon_min, lon_max) = extrema(lon)

    lat_count = Int(ceil((lat_max - lat_min) / resolution)) + 1
    lon_count = Int(ceil((lon_max - lon_min) / resolution)) + 1

    parameters = zeros((lat_count, lon_count, 4))

    for lat_id in 1:1:lat_count
        @info lat_id / lat_count
        for lon_id in 1:1:lon_count
            lat_indices =
                Array(1:1:length(lat))[(lat .>= lat_min + resolution * (lat_id - 1)) .& (lat .< lat_min + resolution * lat_id)]
            lon_indices =
                Array(1:1:length(lon))[(lon .>= lon_min + resolution * (lon_id - 1)) .& (lon .< lon_min + resolution * lon_id)]
            ϕ = data["Band1"][lon_indices, lat_indices][:]
            present = .~(typeof.(ϕ) .<: Missing)
            present_count = sum(present)
            zero_count = sum(.~(typeof.(ϕ) .<: Missing) .& (ϕ .== 0))
            #present_nonzero = .~(typeof.(ϕ) .<: Missing) .& (ϕ .> 0)
            if present_count / prod(size(ϕ)) > 0.5
                # Land
                parameters[lat_id, lon_id, 1] = zero_count ./ present_count
                parameters[lat_id, lon_id, 2] = mean(ϕ[present])
                fmax = sum(ϕ[present] .> mean(ϕ[present])) ./ sum(present)
                parameters[lat_id, lon_id, 3] = fmax
                parameters[lat_id, lon_id, 4] = 1.0
            else
                nothing # all set to zero
            end

        end
    end
    ds = NCDataset("means_ll_$resolution.nc", "c")
    defDim(ds, "lon", lon_count)
    defDim(ds, "lat", lat_count)
    ds.attrib["title"] = "Topographic Index Data"

    la = defVar(ds, "lat", Float32, ("lat",))
    la[:] =
        (0.5:1:(lat_count - 0.5)) ./ lat_count .* (lat_max - lat_min) .+ lat_min

    lo = defVar(ds, "lon", Float32, ("lon",))
    lo[:] =
        (0.5:1:(lon_count - 0.5)) ./ lon_count .* (lon_max - lon_min) .+ lon_min

    mean_ϕ = defVar(ds, "ϕ_mean", Float32, ("lat", "lon"))
    mean_ϕ[:, :] = parameters[:, :, 2]
    f0 = defVar(ds, "fraction_zero", Float32, ("lat", "lon"))
    f0[:, :] = parameters[:, :, 1]
    fmax = defVar(ds, "fmax", Float32, ("lat", "lon"))
    fmax[:, :] = parameters[:, :, 3]
    mask = defVar(ds, "landsea_mask", Float32, ("lat", "lon"))
    mask[:, :] = parameters[:, :, 4]
    mean_ϕ.attrib["units"] = "unitless"
    f0.attrib["units"] = "unitless"
    fmax.attrib["units"] = "unitless"
    mask.attrib["units"] = "unitless"
    la.attrib["units"] = "degrees_north"
    la.attrib["standard_name"] = "latitude"
    lo.attrib["standard_name"] = "longitude"
    lo.attrib["units"] = "degrees_east"
    mean_ϕ.attrib["comments"] = "The mean topographic index"
    fmax.attrib["comments"] = "The maximum saturated fraction from TOPMODEL"
    f0.attrib["comments"] = "The fraction of zero values in cell"
    mask.attrib["comments"] = "1.0 if cell is > 50% land"
    close(ds)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(resolution = parse(Float64, ARGS[1]))
end
