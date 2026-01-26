# rolling means (also called running means) and other reductions


struct OverlappingGroupMapping{T,TRC} <: AbstractGroupMapping
    coord::Vector{T}
    rolling_classes::TRC
    groupindex_to_dataindex::Vector{Vector{Int}}
    dataindex_to_groupindex::Vector{Vector{Int}}
end



# multiple classes for overlapping groups
struct MClass{T <: NTuple}
    class::T
end

Base.length(mc::MClass) = 1
Base.getindex(mc::MClass,i) = mc
Base.in(item,collection::MClass) = item in collection.class

function OverlappingGroupMapping(coord::AbstractVector{T},rolling_classes) where T
    # k is data index
    # ku is the group index

    coord_to_dataindex = Dict{T,Int}()

    for k = 1:length(coord)
        coord_to_dataindex[coord[k]] = k
    end

    groupindex_to_dataindex = [Int[] for i in 1:length(rolling_classes)]
    dataindex_to_groupindex = [Int[] for i in 1:length(coord)]

    for ku = 1:length(rolling_classes)
        for class in rolling_classes[ku].class
            k = get(coord_to_dataindex,class,nothing)

            if !isnothing(k)
                push!(groupindex_to_dataindex[ku],k)
                push!(dataindex_to_groupindex[k],ku)
            end
        end
    end

    return OverlappingGroupMapping(coord,rolling_classes,groupindex_to_dataindex,dataindex_to_groupindex)
end

dataindices(gmap::OverlappingGroupMapping,ku) = gmap.groupindex_to_dataindex[ku]
groupindices(gmap::OverlappingGroupMapping,k) = gmap.dataindex_to_groupindex[k]

ngroups(gmap::OverlappingGroupMapping) = length(gmap.groupindex_to_dataindex)
ndata(gmap::OverlappingGroupMapping) = length(gmap.dataindex_to_groupindex)

function groupsubset(gmap::OverlappingGroupMapping,kus)
    OverlappingGroupMapping(gmap.coord,gmap.rolling_classes[kus])
end

groupsubset(gmap::OverlappingGroupMapping,kus::Colon) = gmap

grouplabel(gmap::OverlappingGroupMapping,ku) = gmap.rolling_classes[ku].class

#--------------
"""
    gv = CommonDataModel.rolling(v::AbstractVariable,:coordname => n)
    gv = CommonDataModel.rolling(v::AbstractVariable,:coordname => rolling_classes)

Create a grouped variable `gv` whose elements composed by all elements in `v`
grouped by a rolling window of length `n` along the coordinate variable `coordname`.
One can also specify a vector classes (`rolling_classes`) with as many elements
as they are groups and the elements coorespond to the values of the coordinate.
Unlike `CommonDataModel.groupby`, the groups defined by `rolling` can overlap.

The grouped variable `gv` and be reduced using the functions `sum` `mean`,
`median`, `var` or `std`, for example `gr = mean(gv)`.
The result `gr` is a lazy structure representing the
outcome of these operations performed over the grouped dimension. Only when the
result `gr` is indexed the actually values are computed.

Broadcasting for `gv` is overloaded. Broadcasting over all elements of
`gv` means that a mapping function is to be applied to all elements of `gv`
before a possible the reduction.

This operations is also called "running mean" when using `mean` as reduction
function.

Example:

```julia
using NCDatasets, Dates
using CommonDataModel: rolling

# create same test data

time = DateTime(2000,1,1):Day(1):DateTime(2009,12,31);  # 10 years
data = rand(Float32.(-9:99),360,180,length(time));
fname = "test_file.nc"
ds = NCDataset(fname,"c");
defVar(ds,"time",time,("time",));
defVar(ds,"data",data,("lon","lat","time"));

# running weekly mean

gv = rolling(ds["data"],:time => 7)
length(gv)
# output 3653 as a mean will be compute for all time instance (including a
# partial mean and the beginning and the end)
size(gv[1])
# output 360 x 180 x 4: the first time slice will only be a mean of 4 values
size(gv[4])
# output 360 x 180 x 7: data from the first 7 days

# compute basic statistics

using Statistics
weekly_running_mean = mean(gv);
size(weekly_running_mean)
# 360 x 180 x 3653 array with the running weekly mean

# get a regular julia array
weekly_running_mean_array = weekly_running_mean[:,:,:];
typeof(weekly_running_mean_array)
# Array{Float32, 3}

# computing a centred 3 monthly mean taking into account that month do not have the
# same length

rolling_classes = [(t - Dates.Month(1)):Day(1):(t + Dates.Month(1)) for t in time]
extrema(length.(rolling_classes))
# output (60, 63)
data_3monthly = mean(rolling(ds["data"],:time => rolling_classes))[:,:,:];

close(ds)
```

"""
function rolling(v,(coordname,rolling_classes)::Pair)
    coord = v[coordname][:]
    _rolling(v,coord,coordname,rolling_classes)
end

function _rolling(v,coord,coordname::SymbolOrString,rolling_classes::AbstractVector)
    rolling_classes_ = [MClass(tuple(r...)) for r in rolling_classes]

    groupmap = OverlappingGroupMapping(coord,rolling_classes_)

    group_fun = identity # still necessary?
    dim = findfirst(==(Symbol(coordname)),Symbol.(dimnames(v)))
    map_fun = identity
    return GroupedVariable(v,coordname,group_fun,groupmap,dim,map_fun)
end

function _rolling(v,coord,coordname::SymbolOrString,nn::Integer)

    # if nn = 7, n0:n1 = -3:3
    # if nn = 8, n0:n1 = -4:3
    # for even nn where are biased towards lower indices as in xarray
    # for odd nn where are perfectly centred

    n0 = nn ÷ 2
    n1 = nn - n0 - 1

    rolling_classes = [coord[max(n-n0,1):min(n+n1,length(coord))]
                       for n = 1:length(coord)]

    return _rolling(v,coord,coordname,rolling_classes)
end
