# Writes a state saved by `run_trendy.jl` (every field of `Y`, HDF5 from
# `ClimaLand.save_checkpoint`) as a NetCDF file on the longitude-latitude grid of the 1°
# box domain: subsurface fields on (lon, lat, z), surface fields on (lon, lat), and fields
# with several components on (lon, lat, component). Every node of the domain is written,
# including those the land mask leaves out, so the file holds the state exactly. Each
# variable keeps its path in `Y` in the attribute `climaland_field`. The file is then read
# back and compared with the state.
#
# Usage (CPU): julia --project=.buildkite state_to_netcdf.jl <state.hdf5> <out.nc> [<title>]
import ClimaComms
import ClimaCore
using ClimaLand
using NCDatasets
using Dates

const FT = Float64
const HDF5 = ClimaCore.InputOutput.HDF5
state_path, out_path = ARGS[1:2]
title = length(ARGS) >= 3 ? join(ARGS[3:end], " ") : "ClimaLand state"

context = ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())
domain = ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))
surface = domain.space.surface
coords = ClimaCore.Fields.coordinate_field(surface)
lat_nodes = vec(Array(parent(coords.lat)))
lon_nodes = vec(Array(parent(coords.long)))
z_layers = Array(parent(ClimaCore.Fields.coordinate_field(domain.space.subsurface).z))[:, 1, 1, 1, 1]
mask = ClimaCore.Spaces.get_mask(surface)
land = vec(Array(parent(mask.is_active)))
lon = sort(unique(lon_nodes))
lat = sort(unique(lat_nodes))
ilon = [searchsortedfirst(lon, x) for x in lon_nodes]
ilat = [searchsortedfirst(lat, x) for x in lat_nodes]
nn = length(lon_nodes)
@assert length(lon) * length(lat) == nn "the nodes do not form a full longitude-latitude grid"
@assert all(lon[ilon] .== lon_nodes) && all(lat[ilat] .== lat_nodes)

ascii_name(path) = foldl((s, r) -> replace(s, r), ("ϑ" => "vartheta", "θ" => "theta", "ρ" => "rho", "∫" => "int_", "." => "_"); init = path)

fields = Pair{String, Array{FT}}[]
HDF5.h5open(state_path, "r") do file
    function walk(group, prefix)
        for key in keys(group)
            obj = group[key]
            path = isempty(prefix) ? key : "$prefix.$key"
            obj isa HDF5.Group ? walk(obj, path) : push!(fields, path => read(obj))
        end
    end
    walk(file["fields/Y"], "")
    global epoch = HDF5.read_attribute(file, "epoch")
    global time_s = HDF5.read_attribute(file, "time")
end
println(length(fields), " fields of Y; land nodes ", count(land), " of ", nn)

# parent layout (levels, i, j, components, nodes) -> (lon, lat, levels) or (lon, lat, components)
function to_grid(data)
    nlev, ncomp = size(data, 1), size(data, 4)
    @assert size(data, 2) == size(data, 3) == 1 && size(data, 5) == nn
    k = nlev > 1 ? nlev : ncomp
    out = zeros(FT, length(lon), length(lat), k)
    for h in 1:nn, m in 1:k
        out[ilon[h], ilat[h], m] = nlev > 1 ? data[m, 1, 1, 1, h] : data[1, 1, 1, m, h]
    end
    return nlev > 1 ? (out, :z) : ncomp > 1 ? (out, :component) : (out[:, :, 1], :surface)
end

isfile(out_path) && rm(out_path)
NCDataset(out_path, "c") do ds
    ds.attrib["title"] = title
    ds.attrib["source"] = "ClimaLand state $(basename(state_path)), written by state_to_netcdf.jl"
    ds.attrib["epoch"] = string(epoch)
    ds.attrib["time_seconds"] = Float64(time_s)
    ds.attrib["created"] = string(Dates.today())
    defDim(ds, "lon", length(lon))
    defDim(ds, "lat", length(lat))
    defDim(ds, "z", length(z_layers))
    v = defVar(ds, "lon", FT, ("lon",); attrib = ["units" => "degrees_east", "standard_name" => "longitude"])
    v[:] = lon
    v = defVar(ds, "lat", FT, ("lat",); attrib = ["units" => "degrees_north", "standard_name" => "latitude"])
    v[:] = lat
    v = defVar(ds, "z", FT, ("z",); attrib = ["units" => "m", "long_name" => "height of the layer centre (negative below the surface)"])
    v[:] = z_layers
    m = defVar(ds, "land_mask", Int8, ("lon", "lat"); attrib = ["long_name" => "1 where the land mask of the domain is active"])
    m[:, :] = to_grid(reshape(Float64.(land), 1, 1, 1, 1, nn))[1] .|> Int8
    for (path, data) in fields
        grid, kind = to_grid(data)
        any(!isfinite, grid) && @warn "$path holds non-finite values" count(!isfinite, grid)
        dims = kind == :z ? ("lon", "lat", "z") : kind == :component ? ("lon", "lat", "component_$(size(grid, 3))") : ("lon", "lat")
        kind == :component && !haskey(ds.dim, dims[3]) && defDim(ds, dims[3], size(grid, 3))
        var = defVar(ds, ascii_name(path), FT, dims; deflatelevel = 4, attrib = ["climaland_field" => path])
        var[:] = grid
    end
end

# read back and compare node by node, bit for bit (NaN matches NaN)
same(x, y) = x == y || (isnan(x) && isnan(y))
mismatches = 0
NCDataset(out_path, "r") do ds
    for (path, data) in fields
        grid = Array(ds[ascii_name(path)])
        for h in 1:nn
            a = ndims(grid) == 3 ? grid[ilon[h], ilat[h], :] : [grid[ilon[h], ilat[h]]]
            b = size(data, 1) > 1 ? data[:, 1, 1, 1, h] : data[1, 1, 1, :, h]
            global mismatches += count(!, same.(a, b))
        end
    end
end
println("wrote ", out_path, " (", round(filesize(out_path) / 2^20; digits = 1), " MiB); values differing after reading back: ", mismatches)
mismatches == 0 || error("the NetCDF file does not reproduce the state")
