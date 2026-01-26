
const SymbolOrString = Union{Symbol, AbstractString}

const TIndices = Union{<:Integer,AbstractVector{<:Integer},Colon}

"""

`AbstractDataset` is a collection of multidimensional variables (for example a
NetCDF or GRIB file)

A data set `ds` of a type derived from `AbstractDataset` should implemented at minimum:

* `Base.key(ds)`: return a list of variable names as strings
* `variable(ds,varname::String)`: return an array-like data structure (derived from `AbstractVariable`) of the variables corresponding to `varname`. This array-like data structure should follow the CF semantics.
* `dimnames(ds)`: should be an iterable with all dimension names in the data set  `ds`
* `dim(ds,name)`: dimension value corresponding to name

Optionally a data set can have attributes and groups:

* `attribnames(ds)`: should be an iterable with all attribute names
* `attrib(ds,name)`: attribute value corresponding to name
* `groupnames(ds)`: should be an iterable with all group names
* `group(ds,name)`: group corresponding to the name

For a writable dataset, one should also implement:
* `defDim`: define a dimension
* `defAttrib`: define a attribute
* `defVar`: define a variable
* `defGroup`: define a group
"""
abstract type AbstractDataset
end


"""
`AbstractVariable{T,N}` is a subclass of `AbstractArray{T, N}`. A variable `v` of a type derived from `AbstractVariable` should implement:

* `name(v)`: should be the name of variable within the data set
* `dimnames(v)`: should be a iterable data structure with all dimension names of the variable `v`
* `dataset(v)`: the parent dataset containing `v`
* `Base.size(v)`: the size of the variable
* `Base.getindex(v,indices...)`: get the data of `v` at the provided indices

Optionally a variable can have attributes:

* `attribnames(v)`: should be an iterable with all attribute names
* `attrib(v,name)`: attribute value corresponding to name

For a writable dataset, one should also implement:
* `defAttrib`: define a attribute
* `Base.setindex!(v,data,indices...)`: set the data in `v` at the provided indices

"""
abstract type AbstractVariable{T,N} <: AbstractArray{T, N}
end


"""
Variable (with applied transformations following the CF convention)
attrib can have different attributes as the parent variables
(used in GRIBDatasets to map from grib attributes to CF attributes)
"""
struct CFVariable{T,N,TV,TA,TSA}  <: AbstractVariable{T, N}
    # this var is generally a `Variable` type
    var::TV
    # Dict-like object for all attributes
    attrib::TA
    # a named tuple with fill value, scale factor, offset,...
    # immutable for type-stability
    _storage_attrib::TSA
end


"""
A collection of attributes with a Dict-like interface dispatching to
`attribnames`, `attrib`, `defAttrib` for `keys`, `getindex` and `setindex!`
respectively.
"""
struct Attributes{TDS<:Union{AbstractDataset,AbstractVariable}} <: AbstractDict{SymbolOrString,Any}
    ds::TDS
end

"""
A collection of dimensions with a Dict-like interface dispatching to
`dimnames`, `dim`, `defDim` for `keys`, `getindex` and `setindex!`
respectively.
"""
struct Dimensions{TDS<:AbstractDataset} <: AbstractDict{SymbolOrString,Any}
    ds::TDS
end


"""
A collection of groups with a Dict-like interface dispatching to
`groupnames` and `group` for `keys` and `getindex` respectively.
"""
struct Groups{TDS<:AbstractDataset} <: AbstractDict{SymbolOrString,Any}
    ds::TDS
end



struct CFStdName
    name::Symbol
end

# Multi-file related type definitions

struct MFVariable{T,N,M,TA,TDS} <: AbstractVariable{T,N}
    ds::TDS
    var::CatArrays.CatArray{T,N,M,TA}
    dimnames::NTuple{N,String}
    varname::String
end

struct MFCFVariable{T,N,M,TA,TV,TDS} <: AbstractVariable{T,N}
    ds::TDS
    cfvar::CatArrays.CatArray{T,N,M,TA}
    var::TV
    dimnames::NTuple{N,String}
    varname::String
end

struct MFDataset{T,N,S<:AbstractString} <: AbstractDataset where T <: AbstractDataset
    ds::Array{T,N}
    aggdim::S
    isnewdim::Bool
    constvars::Vector{Symbol}
    _boundsmap::Dict{String,String}
end

# DeferDataset are Dataset which are open only when there are accessed and
# closed directly after. This is necessary to work with a large number
# of files (e.g. more than 1000).

struct Resource
    filename::String
    mode::String
    args::Dict{Symbol,Any}
    metadata::OrderedDict
end

struct DeferDataset{TDS} <: AbstractDataset
    r::Resource
    groupname::String
    data::OrderedDict
    _boundsmap::Dict{String,String}
end

struct DeferVariable{T,N,TDS} <: AbstractVariable{T,N}
    r::Resource
    varname::String
    data::OrderedDict
end

# view of subsets

struct SubVariable{T,N,TA,TI,TAttrib,TV} <: AbstractVariable{T,N}
    parent::TA
    indices::TI
    attrib::TAttrib
    # unpacked variable
    var::TV
end

struct SubDataset{TD,TI,TA,TG}  <: AbstractDataset
    ds::TD
    indices::TI
    attrib::TA
    group::TG
end


const Iterable = Union{Attributes,Dimensions,Groups,AbstractDataset}


# memory dataset

mutable struct ResizableArray{T,N} <: AbstractArray{T,N}
    A::AbstractArray{T,N}
    fillvalue::T
end

struct MemoryVariable{T,N,TP,TA <: AbstractArray{T,N}} <: AbstractVariable{T,N}
    parent_dataset::TP
    name::String
    dimnames::NTuple{N,String}
    data::TA
    _attrib::OrderedDict{String,Any}
end

struct MemoryDataset{TP <: Union{Nothing,AbstractDataset},Tmasingvalue} <: AbstractDataset
    parent_dataset::TP
    name::String # "/" for root group
    dimensions::OrderedDict{String,Int}
    variables::OrderedDict{String,MemoryVariable}
    _attrib::OrderedDict{String,Any}
    unlimited::Vector{String}
    _group::OrderedDict{String,Any}
    maskingvalue::Tmasingvalue
end

