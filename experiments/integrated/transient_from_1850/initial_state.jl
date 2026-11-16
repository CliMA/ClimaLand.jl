import ClimaCore
import NCDatasets
const HDF5 = ClimaCore.InputOutput.HDF5

"""
    state_fields(Y)

The fields of the state `Y` as `path => field` pairs, with paths such as `"soil.ϑ_l"`.
Fields with several components (such as `canopy.photosynthesis.acclimated`) are one
entry.
"""
function state_fields(Y, prefix = "")
    out = Pair{String, Any}[]
    for name in propertynames(Y)
        x = getproperty(Y, name)
        path = isempty(prefix) ? string(name) : "$prefix.$name"
        if x isa ClimaCore.Fields.Field
            push!(out, path => x)
        else
            append!(out, state_fields(x, path))
        end
    end
    return out
end

"""The index of `x` in the sorted `grid`, or 0 if `x` is not one of its points."""
function grid_index(grid, x)
    i = searchsortedfirst(grid, x - 1e-9)
    return i <= length(grid) && abs(grid[i] - x) < 1e-9 ? i : 0
end

"""
    set_state_from_netcdf!(Y, path)

Sets every field of `Y` from the NetCDF file `path` holding a full state (the
`climaland_1850_spunup_ic` artifact: one variable per field, found by its attribute
`climaland_field`, on the longitude-latitude grid of the 1° box domain). The values are
copied node by node, so the state is reproduced exactly; the domain of `Y` must have
the file's grid and layers.
"""
function set_state_from_netcdf!(Y, path)
    surface_space = axes(Y.canopy.energy.T)
    coords = ClimaCore.Fields.coordinate_field(surface_space)
    lon_nodes = vec(Array(parent(coords.long)))
    lat_nodes = vec(Array(parent(coords.lat)))
    z_field = ClimaCore.Fields.coordinate_field(axes(Y.soil.ϑ_l)).z
    z_nodes = Array(parent(z_field))[:, 1, 1, 1, 1]
    NCDatasets.NCDataset(path) do ds
        lon, lat, z = Array(ds["lon"]), Array(ds["lat"]), Array(ds["z"])
        ilon = grid_index.(Ref(lon), lon_nodes)
        ilat = grid_index.(Ref(lat), lat_nodes)
        same_layers =
            length(z) == length(z_nodes) && all(abs.(z .- z_nodes) .< 1e-9)
        if !(all(>(0), ilon) && all(>(0), ilat) && same_layers)
            error("$path is not on the grid of this domain")
        end
        variables = Dict(
            ds[name].attrib["climaland_field"] => name for
            name in keys(ds) if haskey(ds[name].attrib, "climaland_field")
        )
        for (field_path, field) in state_fields(Y)
            haskey(variables, field_path) ||
                error("$path has no variable for $field_path")
            grid = Array(ds[variables[field_path]])
            data = Array(parent(field))
            nlev, ncomp = size(data, 1), size(data, 4)
            for h in eachindex(ilon), m in 1:max(nlev, ncomp)
                value =
                    ndims(grid) == 3 ? grid[ilon[h], ilat[h], m] :
                    grid[ilon[h], ilat[h]]
                nlev > 1 ? (data[m, 1, 1, 1, h] = value) :
                (data[1, 1, 1, m, h] = value)
            end
            copyto!(parent(field), data)
        end
    end
    return Y
end

"""
    set_state_from_checkpoint!(Y, path)

Sets every field of `Y` from a state saved by `ClimaLand.save_checkpoint` (HDF5), by
copying its datasets into the parent arrays of `Y`. Restarts and branches (S1 from S2 at
1901, S2E from S2 at 1979) use it: ClimaCore cannot rebuild the meshes of
longitude-latitude box domains from the file, so the state is built from the model and
filled.
"""
function set_state_from_checkpoint!(Y, path)
    HDF5.h5open(path, "r") do file
        copy_group!(Y, file["fields/Y"], "Y")
    end
    return Y
end

function copy_group!(target, group, label)
    for key in keys(group)
        obj = group[key]
        sub = getproperty(target, Symbol(key))
        if obj isa HDF5.Group
            copy_group!(sub, obj, "$label.$key")
        else
            data = read(obj)
            dest = parent(sub)
            size(data) == size(dest) || error(
                "$label.$key: size $(size(data)) in the file, $(size(dest)) in the state",
            )
            copyto!(dest, data)
        end
    end
end

"""
    set_initial_state!(Y, path)

Sets `Y` from `path`: a NetCDF state (`.nc`) or a checkpoint (`.hdf5`).
"""
set_initial_state!(Y, path) =
    endswith(path, ".nc") ? set_state_from_netcdf!(Y, path) :
    set_state_from_checkpoint!(Y, path)
