import ClimaCore
const HDF5 = ClimaCore.InputOutput.HDF5

"""
    read_into!(target, path, name)

Copies the datasets of the field `name` written by `ClimaCore.InputOutput` in the HDF5
file `path` (the state `"Y"` of `ClimaLand.save_checkpoint`, or the sums `"acc"` of
`save_spinup`) into the parent arrays of the matching fields of `target`. The spaces of
`target` are kept: ClimaCore cannot rebuild the meshes of longitude-latitude box domains
from a file ("Invalid coord type LatPoint").
"""
function read_into!(target, path, name)
    HDF5.h5open(path, "r") do file
        copy_group!(target, file["fields/$name"], name)
    end
    return target
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
            size(data) == size(dest) ||
                error("$label.$key: size $(size(data)) in the file, $(size(dest)) in the target")
            copyto!(dest, data)
        end
    end
end

"""
    read_root_attribute(path, name)

The attribute `name` of the root group of the HDF5 file `path`.
"""
read_root_attribute(path, name) =
    HDF5.h5open(file -> HDF5.read_attribute(file, name), path, "r")

"""
    read_dataset(path, name, chain...)

The array of the dataset `fields/<name>/<chain...>` of the HDF5 file `path`, in the
layout of the parent array of the field.
"""
read_dataset(path, name, chain...) =
    HDF5.h5open(path, "r") do file
        read(file[join(("fields", name, string.(chain)...), "/")])
    end
