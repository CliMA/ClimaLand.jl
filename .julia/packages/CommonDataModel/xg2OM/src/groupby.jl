#
# types
#

# mappings between groups and source data

abstract type AbstractGroupMapping; end

# non-overlapping groups
struct GroupMapping{TClass,TUC} <: AbstractGroupMapping where {TClass <: AbstractVector{T}, TUC <: Union{T,AbstractVector{T}}}  where T
    class::TClass
    unique_class::TUC
end

struct GroupedDataset{TDS<:AbstractDataset,TF,TGM,TM}
    ds::TDS # dataset
    coordname::Symbol
    group_fun::TF # mapping function
    groupmap::TGM
    map_fun::TM
end

struct ReducedGroupedDataset{TDS,TF,TGM,TM,TRF} <: AbstractDataset
    ds::TDS # dataset
    coordname::Symbol
    group_fun::TF # mapping function
    groupmap::TGM
    map_fun::TM
    reduce_fun::TRF
end

struct GroupedVariable{TV,TF,TGM,TM,TG} <: AbstractVector{TG} where TV <: AbstractArray{T,N} where {T,N}
    v::TV # dataset
    coordname::Symbol
    # mapping function to create groups when applied to coordinate
    group_fun::TF
    groupmap::TGM
    dim::Int
    map_fun::TM
end

# reduce_fun is e.g. sum, mean, var,...
# map_fun is the mapping function applied before reduction
struct ReducedGroupedVariable{T,N,TGV,TF}  <: AbstractVariable{T,N}
    gv::TGV
    reduce_fun::TF
    _attrib::OrderedDict{String,Any}
end


function dataindices(gmap::GroupMapping,ku)
    class_ku = gmap.unique_class[ku]
    return findall(==(class_ku),gmap.class)
end

function groupindices(gmap::GroupMapping,k)
    # Types like Dates.Month behave different than normal scalars, e.g.
    # julia> length(Dates.Month(1))
    # ERROR: MethodError: no method matching length(::Month)

    if gmap.unique_class isa DatePeriod
        if gmap.unique_class == gmap.class[k]
            return (1,)
        else
            return ()
        end
    end

    for ku = 1:length(gmap.unique_class)
        if gmap.unique_class[ku] == gmap.class[k]
            return (ku,)
        end
    end
    return ()
end

grouplabel(gmap::GroupMapping,ku) = gmap.unique_class[ku]


# Types like Dates.Month behave different than normal scalars, e.g.
# julia> length(Dates.Month(1))
# ERROR: MethodError: no method matching length(::Month)

_length(x) = length(x)
_length(x::DatePeriod) = 1

ngroups(gmap::GroupMapping) = _length(gmap.unique_class)
ndata(gmap::GroupMapping) = length(gmap.class)


function groupsubset(gmap::GroupMapping,kus)
    return GroupMapping(gmap.class,gmap.unique_class[kus])
end

groupsubset(gmap::GroupMapping,kus::Colon) = gmap

#--------------

# helper function

_indices_helper(j,ku,i,val) = ()
_indices_helper(j,ku,i,val,ind1::Integer,indices...) = _indices_helper(j,ku,i+1,val,indices...)
_indices_helper(j,ku,i,val,ind1,indices...) = ((i == j ? ku : val), _indices_helper(j,ku,i+1,val,indices...)...)


# _indices(j,ku,val,indices) produces a tuple with val for every dimension
# except `ku` for the the `j` dimension for an array A after subsetting it as
# `A[indices]`.

_indices(j,ku,val,indices) = _indices_helper(j,ku,1,val,indices...)
_dest_indices(j,ku,indices) = _indices_helper(j,ku,1,:,indices...)


@inline size_getindex(array,indexes...) = _size_getindex(array,(),1,indexes...)
@inline _size_getindex(array,sh,n,i::Integer,indexes...) = _size_getindex(array,sh,                   n+1,indexes...)
@inline _size_getindex(array::AbstractArray,sh,n,i::Colon,  indexes...) = _size_getindex(array,(sh...,size(array,n)),n+1,indexes...)
@inline _size_getindex(sz::Tuple,sh,n,i::Colon,  indexes...) = _size_getindex(sz,(sh...,sz[n]),n+1,indexes...)
@inline _size_getindex(array,sh,n,i,         indexes...) = _size_getindex(array,(sh...,length(i)),    n+1,indexes...)
@inline _size_getindex(array,sh,n) = sh

#
# methods with ReducedGroupedDataset as main argument
#

Base.keys(gds::ReducedGroupedDataset) = keys(gds.ds)

function variable(gds::ReducedGroupedDataset,varname::SymbolOrString)
    v = variable(gds.ds,varname)

    dim = findfirst(==(gds.coordname),Symbol.(dimnames(v)))
    if isnothing(dim)
        return v
    else
        gv = GroupedVariable(
            v,
            gds.coordname,
            gds.group_fun,
            gds.groupmap,
            dim,
            gds.map_fun)
        return ReducedGroupedVariable(gv,gds.reduce_fun)
    end
end


attribnames(gds::ReducedGroupedDataset) = attribnames(gds.ds)
attrib(gds::ReducedGroupedDataset,attribname::SymbolOrString) =
    attrib(gds.ds,attribname)

#
# methods with GroupedVariable as main argument
#

function Base.show(io::IO,::MIME"text/plain",gv::GroupedVariable)
    println(io,length(gv),"-element ",
            "variable grouped by '",gv.coordname,"'",
            " of array '",name(gv.v),"' ",
            join(string.(size(gv.v)),"×"),
            " (",
            join(dimnames(gv.v),"×"),
            ")",
            )
end

Base.show(io::IO,gv::GroupedVariable) = Base.show(io,MIME"text/plain",gv)
Base.ndims(gv::GroupedVariable) = 1
Base.size(gv::GroupedVariable) = (ngroups(gv.groupmap),)
Base.eltype(gv::GroupedVariable{TV,TF,TGM,TM,TG}) where {TV,TF,TGM,TM,TG} = TG

function Base.getindex(gv::GroupedVariable,k::Integer)
    class_k,indices = group(gv,k)
    return gv.map_fun(Array(selectdim(gv.v,gv.dim,indices)))
end

function group(gv::GroupedVariable,ku::Integer)
    # make generic
    class_ku = grouplabel(gv.groupmap,ku)
    k = dataindices(gv.groupmap,ku)
    return class_ku, k
end

function _mapreduce(map_fun,reduce_op,gv::GroupedVariable{TV},indices;
                 init = reduce(reduce_op,T[])) where TV <: AbstractArray{T,N} where {T,N}
    data = gv.v
    dim = findfirst(==(Symbol(gv.coordname)),Symbol.(dimnames(data)))
    group_fun = gv.group_fun

    groupmap = groupsubset(gv.groupmap,indices[dim])

    nclass = ngroups(groupmap)
    sz_all = ntuple(i -> (i == dim ? nclass : size(data,i) ),ndims(data))
    sz = size_getindex(sz_all,indices...)

    data_by_class = Array{T,length(sz)}(undef,sz)
    data_by_class .= init

    count = zeros(Int,nclass)
    for k = 1:size(data,dim)

        for ku in groupindices(groupmap,k)
            dest_ind = _dest_indices(dim,ku,indices)
            src_ind = ntuple(i -> (i == dim ? k : indices[i] ),ndims(data))
            data_by_class_ind = view(data_by_class,dest_ind...)

            data_by_class_ind .= reduce_op.(
                data_by_class_ind,
                map_fun(data[src_ind...]))
            count[ku] += 1
        end
    end

    return data_by_class,reshape(count,_indices(dim,length(count),1,indices))
end



function _mapreduce_aggregation(map_fun,ag,gv::GroupedVariable{TV},indices) where TV <: AbstractArray{T,N} where {T,N}
    data = gv.v
    dim = findfirst(==(Symbol(gv.coordname)),Symbol.(dimnames(data)))
    group_fun = gv.group_fun
    groupmap = groupsubset(gv.groupmap,indices[dim])

    nclass = ngroups(groupmap)
    sz_all = ntuple(i -> (i == dim ? nclass : size(data,i) ),ndims(data))
    sz = size_getindex(sz_all,indices...)

    data_by_class = fill(ag(T),sz)

    count = zeros(Int,nclass)
    for k = 1:size(data,dim)
        for ku in groupindices(groupmap,k)
            dest_ind = _dest_indices(dim,ku,indices)
            src_ind = ntuple(i -> (i == dim ? k : indices[i] ),ndims(data))
            data_by_class_ind = view(data_by_class,dest_ind...)
            std_data_ind = map_fun(data[src_ind...])
            data_by_class_ind .= update.(data_by_class_ind,std_data_ind)
        end
    end

    return result.(data_by_class)
end

function _reduce(args...; kwargs...)
    _mapreduce(identity,args...; kwargs...)
end

struct GroupedVariableStyle <: BroadcastStyle end
Base.BroadcastStyle(::Type{<:GroupedVariable}) = GroupedVariableStyle()

"""
    A = find_gv(T,As)

returns the first type T among the arguments.
"""
find_gv(T,bc::Base.Broadcast.Broadcasted) = find_gv(T,bc.args)
find_gv(T,args::Tuple) = find_gv(T,find_gv(T,args[1]), Base.tail(args))
find_gv(T,x) = x
find_gv(T,::Tuple{}) = nothing
find_gv(::Type{T},a::T, rest) where T = a
find_gv(T,::Any, rest) = find_gv(T,rest)

function Base.similar(bc::Broadcasted{GroupedVariableStyle}, ::Type{ElType})  where ElType
    A = find_gv(GroupedVariable,bc)
    return A
end

function Base.broadcasted(::GroupedVariableStyle,f::Function,A::GroupedVariable{TV,TF,TGM,TM,TG}) where {TV,TF,TGM,TM,TG}
    # TODO change output TG

    map_fun = ∘(f,A.map_fun)
    TM2 = typeof(map_fun)
    TG2 = TG

    ff = map_fun ∘ Array ∘ selectdim
    #TG = Base.return_types(selectdim,(TV,Int,Int,))[1]
    TG2 = Base.return_types(ff,(TV,Int,Int,))[1]

    GroupedVariable{TV,TF,TGM,TM2,TG2}(
        A.v,A.coordname,A.group_fun,A.groupmap,A.dim,map_fun)
end

function GroupedVariable(v::TV,coordname,group_fun::TF,groupmap,dim,map_fun::TM) where TV <: AbstractVariable where {TF,TM}
    TGM = typeof(groupmap)

    #TG = Base.return_types(selectdim,(TV,Int,Int,))[1]
    TG = Base.return_types(_array_selectdim,(TV,Int,Vector{Int}))[1]

    @debug "inferred types" TV TF TGM TM TG
#    groupmap = GroupMapping(class,unique_class)

    GroupedVariable{TV,TF,TGM,TM,TG}(
        v,Symbol(coordname),group_fun,groupmap,dim,map_fun)
end


"""
    gv = CommonDataModel.groupby(v::AbstractVariable,:coordname => group_fun)
    gv = CommonDataModel.groupby(v::AbstractVariable,"coordname" => group_fun)

Create a grouped variable `gv` whose elements composed by all elements in `v`
whose corresponding coordinate variable (with the name `coordname`) map to the
same value once the group function `group_fun` is applied to the coordinate.

The grouped variable `gv` and be reduced using the functions `sum` `mean`,
`median`, `var` or `std`, for example `gr = mean(gv)`.
The result `gr` is a lazy structure representing the
outcome of these operations performed over the grouped dimension. Only when the
result `gr` is indexed the actually values are computed.

Broadcasting for `gv` and `gr` is overloaded. Broadcasting over all elements of
`gv` means that a mapping function is to be applied to all elements of `gv`
before a possible the reduction.
Broadcasting over `gr`, for example `v .- gr` mean that `gr` is broadcasted over the
full size of `v` according to the grouping function.

Example:

```julia
using NCDatasets, Dates
using CommonDataModel: @groupby, groupby

# create same test data

time = DateTime(2000,1,1):Day(1):DateTime(2009,12,31);  # 10 years
data = rand(Float32.(-9:99),360,180,length(time));
fname = "test_file.nc"
ds = NCDataset(fname,"c");
defVar(ds,"time",time,("time",));
defVar(ds,"data",data,("lon","lat","time"));

# group over month

gv = @groupby(ds["data"],Dates.Month(time))
# or
# gv = groupby(ds["data"],:time => Dates.Month)
length(gv)
# output 12 as they are all 12 months in this dataset
size(gv[1])
# 360 x 180 x 310 array with all time slices corresponding to the 1st month

# the variable `gv` is equivalent to the following operation
time_month = Dates.Month.(ds[:time][:])
gv2 = [ds[:data][:,:,findall(time_month .== m)] for m in sort(unique(time_month))];

# compute basic statistics

using Statistics
monthly_mean = mean(gv);
size(monthly_mean)
# 360 x 180 x 12 array with the monthly mean

# get a regular julia array
monthly_mean_array = monthly_mean[:,:,:];
typeof(monthly_mean_array)
# Array{Float32, 3}

# substact from data the corresponding monthly mean
monthly_anomalies = data .- mean(gv);

close(ds)
```

"""
function groupby(v::AbstractVariable,(coordname,group_fun)::Pair{<:SymbolOrString,TF}) where TF
    # for NCDatasets 0.12
    c = v[String(coordname)][:]
    class = group_fun.(c)
    unique_class = sort(unique(class))
    dim = findfirst(==(Symbol(coordname)),Symbol.(dimnames(v)))
    map_fun = identity
    groupmap = GroupMapping(class,unique_class)
    return GroupedVariable(v,Symbol(coordname),group_fun,groupmap,dim,map_fun)
end


function groupby(ds::AbstractDataset,(coordname,group_fun)::Pair{<:SymbolOrString,TF}) where TF
    c = ds[String(coordname)][:]
    class = group_fun.(c)
    unique_class = sort(unique(class))
    map_fun = identity
    groupmap = GroupMapping(class,unique_class)
    return GroupedDataset(ds,Symbol(coordname),group_fun,groupmap,map_fun)
end

"""
    gv = CommonDataModel.@groupby(v,group_fun(coordname))


Create a grouped variable `gv` whose elements are composed by all elements in `v`
whose corresponding coordinate variable (with the name `coordname`) map to the
same value once the group function `group_fun` is applied to the coordinate variable.

See [`groupby`](@ref CommonDataModel.groupby) for more information.
"""
macro groupby(vsym,expression)
    (param, newsym),exp = scan_coordinate_name(expression)
    fun = :($newsym -> $exp)
    return :(groupby($(esc(vsym)),$(Meta.quot(param)) => $fun))
end

function ReducedGroupedVariable(gv::GroupedVariable,reduce_fun)
    T = eltype(gv.v)
    @debug "inference " T reduce_fun  Base.return_types(reduce_fun, (Vector{T},))
    N = ndims(gv.v)
    _attrib = OrderedDict(gv.v.attrib)
    ReducedGroupedVariable{T,N,typeof(gv),typeof(reduce_fun)}(gv,reduce_fun,_attrib)
end

function ReducedGroupedDataset(gds::GroupedDataset,reduce_fun)
    return ReducedGroupedDataset(
        gds.ds,
        gds.coordname,
        gds.group_fun,
        gds.groupmap,
        gds.map_fun,
        reduce_fun,
    )
end

"""
    gr = reduce(f,gv::GroupedVariable)

Reduce the grouped variable `gv` along grouped dimension using the function `f`.
The function `f` will be called as `f(x,dims=d)` where `x` array (an element
of `gv`) and `d` is an integer of the dimension overwhich one need to reduce
`x`.
"""
Base.reduce(f::Function,gv::GroupedVariable) = ReducedGroupedVariable(gv,f)
Base.reduce(f::typeof(hcat),gv::GroupedVariable) = ReducedGroupedVariable(gv,f)
Base.reduce(f::typeof(vcat),gv::GroupedVariable) = ReducedGroupedVariable(gv,f)

Base.reduce(f::Function,gds::GroupedDataset) = ReducedGroupedDataset(gds,f)

for fun in (:maximum, :mean, :median, :minimum, :std, :sum, :var)
    @eval $fun(gv::GroupedVariable) = reduce($fun,gv)
    @eval $fun(gds::GroupedDataset) = reduce($fun,gds)
end

# methods with ReducedGroupedVariable as main argument

Base.ndims(gr::ReducedGroupedVariable) = ndims(gr.gv.v)
Base.size(gr::ReducedGroupedVariable) = ntuple(ndims(gr)) do i
    if i == gr.gv.dim
        length(gr.gv)
    else
        size(gr.gv.v,i)
    end
end

dimnames(gr::ReducedGroupedVariable) = dimnames(gr.gv.v)
name(gr::ReducedGroupedVariable) = name(gr.gv.v)


attribnames(gr::ReducedGroupedVariable) = collect(keys(gr._attrib))
attrib(gr::ReducedGroupedVariable,attribname::SymbolOrString) = gr._attrib[attribname]
defAttrib(gr::ReducedGroupedVariable,attribname::SymbolOrString,value) =
    gr._attrib[attribname] = value

struct ReducedGroupedVariableStyle <: BroadcastStyle end

Base.BroadcastStyle(::Type{<:ReducedGroupedVariable}) = ReducedGroupedVariableStyle()
Base.BroadcastStyle(::DefaultArrayStyle,::ReducedGroupedVariableStyle) = ReducedGroupedVariableStyle()
Base.BroadcastStyle(::ReducedGroupedVariableStyle,::DefaultArrayStyle) = ReducedGroupedVariableStyle()

function Base.similar(bc::Broadcasted{ReducedGroupedVariableStyle}, ::Type{ElType})  where ElType
    # Scan the inputs for the ReducedGroupedVariable:
    A = find_gv(ReducedGroupedVariable,bc)
    return similar(A.gv.v)
end

# _array_selectdim_indices(ind,dim,i,sz...)
# returns a tuple (:,:,:,i,:,:,:) where the i is at the dim position
# in total there are as many indices as elements in the tuple sz
# (typically the size of the array)

_array_selectdim_indices(ind,dim,i,sz1,rest...) = _array_selectdim_indices((ind...,(length(ind) == dim-1 ? i : (:))),dim,i,rest...)
_array_selectdim_indices(ind,dim,i) = ind


# indices_B is not type-stable as dim is not know at compile type
# but if i is a range (e.g. 1:2), then the type-unstability does not propagate
function _array_selectdim(B,dim,i)
    indices_B = _array_selectdim_indices((),dim,i,size(B)...)
    return B[indices_B...]
end


_broadcasted_array_selectdim(A::ReducedGroupedVariable,dim,indices,k) = _array_selectdim(A,dim,k:k)
_broadcasted_array_selectdim(A,dim,indices,k) = _array_selectdim(A,dim,indices)

function broadcasted_gvr!(C,f,A,B)
    gr = find_gv(ReducedGroupedVariable,(A,B))
    gv = gr.gv
    dim = gr.gv.dim

    for k = 1:length(gv)
        class_k, indices = group(gv,k)

        selectdim(C,dim,indices) .= broadcast(
            f,
            _broadcasted_array_selectdim(A,dim,indices,k),
            _broadcasted_array_selectdim(B,dim,indices,k))
    end

    return C
end


Base.broadcasted(::ReducedGroupedVariableStyle,f::Function,A,B::ReducedGroupedVariable) =
    broadcasted_gvr!(similar(A),f,A,B)
Base.broadcasted(::ReducedGroupedVariableStyle,f::Function,A::ReducedGroupedVariable,B) =
    broadcasted_gvr!(similar(B),f,A,B)

function Base.broadcasted(::ReducedGroupedVariableStyle,f::Function,A::ReducedGroupedVariable,B::ReducedGroupedVariable)
    # undecided what to do
    # method needs to be there to avoid ambiguities
    error("unimplemented");
end

function Base.Array(gr::ReducedGroupedVariable)
    gr[ntuple(i -> Colon(),ndims(gr))...]
end

function Base.getindex(gr::ReducedGroupedVariable{T,N,TGV,typeof(sum)},indices::TIndices...) where {T,N,TGV}
    data,count = _mapreduce(gr.gv.map_fun,+,gr.gv,indices)
    data
end

function Base.getindex(gr::ReducedGroupedVariable{T,N,TGV,typeof(mean)},indices::TIndices...) where {T,N,TGV}
    data,count = _mapreduce(gr.gv.map_fun,+,gr.gv,indices)
    data ./ count
end


function Base.getindex(gr::ReducedGroupedVariable{T,N,TGV,TF},indices::TIndices...) where TF <: Union{typeof(var),typeof(maximum),typeof(minimum)} where {T,N,TGV}

    return _mapreduce_aggregation(
        gr.gv.map_fun,aggregator(TF),gr.gv,indices);
end


function Base.getindex(gr::ReducedGroupedVariable{T,N,TGV,typeof(std)},indices::TIndices...) where {T,N,TGV}

    return sqrt.(_mapreduce_aggregation(
        gr.gv.map_fun,VarianceWelfordAggregation,gr.gv,indices))
end


_dim_after_getindex(dim,ind::TIndices,other...) = _dim_after_getindex(dim+1,other...)
_dim_after_getindex(dim,ind::Integer,other...) = _dim_after_getindex(dim,other...)
_dim_after_getindex(dim) = dim

function Base.getindex(gr::ReducedGroupedVariable{T},indices::TIndices...) where T
    gv = gr.gv
    sz = size_getindex(gr,indices...)
    data_by_class = Array{T}(undef,sz)

    # after indexing some dimensions are not longer present
    cdim = _dim_after_getindex(0,indices[1:(gv.dim-1)]...) + 1

    indices_source = ntuple(ndims(gr)) do i
        if i == gv.dim
            (:)
        else
            indices[i]
        end
    end

    for (kl,ku) in enumerate(to_indices(gr,indices)[gv.dim])
        dest_ind = _dest_indices(gv.dim,kl,indices)
        data = gv[ku]
        data_by_class[dest_ind...] = gr.reduce_fun(gv.map_fun(data[indices_source...]),dims=cdim)
    end

    return data_by_class
end


function dataset(gr::ReducedGroupedVariable)
    gv = gr.gv
    ds = dataset(gv.v)

    return ReducedGroupedDataset(
        ds,gv.coordname,gv.group_fun,
        gv.groupmap,
        gv.map_fun,
        gr.reduce_fun,
    )
end
