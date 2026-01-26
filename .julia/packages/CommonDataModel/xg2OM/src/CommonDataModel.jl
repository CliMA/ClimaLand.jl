module CommonDataModel

using Base.Broadcast: Broadcasted, BroadcastStyle, DefaultArrayStyle

using CFTime
using Dates
using Printf
using Preferences
using DataStructures
import Base:
    LogicalIndex,
    checkbounds,
    close,
    collect,
    display,
    filter,
    getindex,
    in,
    isopen,
    iterate,
    ndims,
    reduce,
    show,
    size,
    write

import Statistics
import Statistics:
    maximum,
    mean,
    median,
    minimum,
    std,
    sum,
    var




include("CatArrays.jl")
include("types.jl")
include("dataset.jl")
include("variable.jl")
include("cfvariable.jl")
include("attribute.jl")
include("dimension.jl")
include("cfconventions.jl")
include("multifile.jl")
include("defer.jl")
include("subvariable.jl")
include("select.jl")
include("aggregation.jl")
include("groupby.jl")
include("rolling.jl")
include("memory_dataset.jl")

end # module CommonDataModel

#  LocalWords:  AbstractDataset NetCDF GRIB ds AbstractVariable
#  LocalWords:  varname dimnames iterable attribnames attrib dataset
#  LocalWords:  groupnames AbstractArray
