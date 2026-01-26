module DataAPI

"""
    defaultarray(T, N)

For a given element type `T` and number of dimensions `N`, return the appropriate array
type.

The default definition returns `Array{T, N}`. This function is useful for custom types that
have a more efficient vectorized representation (usually using SOA optimizations).

This generic function is owned by DataAPI.jl itself, which is the sole provider of the
default definition.
"""
function defaultarray end
defaultarray(::Type{T}, N) where {T} = Array{T, N}

"""
    refarray(A::AbstractArray)

For a given array `A`, potentially return an optimized "ref array" representation of the
original array, which can allow for faster comparison and sorting.

The default definition just returns the input array. This function is useful for custom
array types which already store a "hashed"-like representation of elements where testing
equality or permuting elements in place can be much faster than the original scalar value,
like pooled arrays.

This generic function is owned by DataAPI.jl itself, which is the sole provider of the
default definition.
"""
function refarray end
refarray(A::AbstractArray) = A

"""
    refvalue(A, x)

For the *original* array `A`, and a "ref value" `x` taken from `refarray(A)`, return the
appropriate *original* value. `refvalue(A, refarray(A)[I...])` must be equal to `A[I...]`.

By default, `refvalue(A, x)` returns `x` (since `refarray(A)` returns `A` by default).
This allows recovering an original array element after operating on the "ref values".

This generic function is owned by DataAPI.jl itself, which is the sole provider of the
default definition.
"""
function refvalue end
refvalue(A::AbstractArray, x) = x

"""
    refpool(A)

Whenever available, return an indexable object `pool` such that, given the *original* array `A` and
a "ref value" `x` taken from `refarray(A)`, `pool[x]` is the appropriate *original* value. Return
`nothing` if such object is not available.

By default, `refpool(A)` returns `nothing`.

If `refpool(A)` is not `nothing`, then `refpool(A)[refarray(A)[I...]]`
must be equal to (according to `isequal`) and of the same type as `A[I...]`,
and the object returned by `refpool(A)` must implement the iteration and
indexing interfaces as well as the `length`, `eachindex`, `keys`, `values`, `pairs`,
`firstindex`, `lastindex`, and `eltype` functions
in accordance with the `AbstractArray` interface.

This generic function is owned by DataAPI.jl itself, which is the sole provider of the
default definition.
"""
function refpool end
refpool(A::AbstractArray) = nothing

"""
    invrefpool(A)

Whenever available, return an indexable object such that given an array `A`
for which `refpool(A)` is not `nothing`:

* for any valid index `x` into `refpool(A)`, `invrefpool(A)[refpool(A)[x]]` is equal to `x`
  (according to `isequal`) and of the same type as `x`;
* for any valid index `ix` into `invrefpool(A)` , `refpool(A)[invrefpool(A)[ix]]` is equal to `ix`
  (according to `isequal`) and of the same type as `ix`.

Additionally it is required that for `invrefpool(A)` the following methods are defined:

* `Base.haskey`: allowing to check if `ix` is a valid index into it.
* `Base.get`: allowing to get a value from it or a passed default value if it is not present.

By default, `invrefpool(A)` returns `nothing`.

If `invrefpool(A)` is not `nothing`, then `refpool(A)` also must not be `nothing`.

This generic function is owned by DataAPI.jl itself, which is the sole provider of the
default definition.
"""
function invrefpool end
invrefpool(A::AbstractArray) = nothing

"""
    describe(io::IO, x)

For an object `x`, print descriptive statistics to `io`.

This generic function is owned by StatsBase.jl, which is the sole provider of the default
definition.
"""
function describe end

# Sentinel type needed to make `levels` inferrable
struct _Default end

"""
    levels(x; skipmissing=true)

Return a vector of unique values which occur or could occur in collection `x`.
`missing` values are skipped unless `skipmissing=false` is passed.

Values are returned in the preferred order for the collection,
with the result of [`sort`](@ref) as a default.
If the collection is not sortable then the order of levels is unspecified.

Contrary to [`unique`](@ref), this function may return values which do not
actually occur in the data, and does not preserve their order of appearance in `x`.
"""
@inline levels(x; skipmissing::Union{Bool, _Default}=_Default()) =
    skipmissing isa _Default || skipmissing ?
        _levels_skipmissing(x) : _levels_missing(x)

# The `which` check is here for backward compatibility:
# if a type implements a custom `levels` method but does not support
# keyword arguments, `levels(x, skipmissing=true/false)` will dispatch
# to the fallback methods here, and we take care of calling that custom method
function _levels_skipmissing(x)
    if which(DataAPI.levels, Tuple{typeof(x)}) === which(DataAPI.levels, Tuple{Any})
        T = Base.nonmissingtype(eltype(x))
        u = unique(x)
        # unique returns its input with copying for ranges
        # (and possibly for other types guaranteed to hold unique values)
        nmu = (u isa AbstractRange || u === x || Base.mightalias(u, x)) ?
            filter(!ismissing, u) : filter!(!ismissing, u)
        levs = convert(AbstractArray{T}, nmu)
        try
            sort!(levs)
        catch
        end
        return levs
    else
        return levels(x)
    end
end

function _levels_missing(x)
    if which(DataAPI.levels, Tuple{typeof(x)}) === which(DataAPI.levels, Tuple{Any})
        u = convert(AbstractArray{eltype(x)}, unique(x))
        # unique returns its input with copying for ranges
        # (and possibly for other types guaranteed to hold unique values)
        levs = (x isa AbstractRange || u === x || Base.mightalias(u, x)) ?
            Base.copymutable(u) : u
        try
            sort!(levs)
        catch
        end
        return levs
    # This is a suboptimal fallback since it does a second pass over the data
    elseif any(ismissing, x)
        return [levels(x); missing]
    else
        return convert(AbstractArray{eltype(x)}, levels(x))
    end
end

"""
    Between(first, last)

Select the columns between `first` and `last` (including both) from a table.
"""
struct Between{T1 <: Union{Int, Symbol}, T2 <: Union{Int, Symbol}}
    first::T1
    last::T2
end

Between(x::AbstractString, y::AbstractString) = Between(Symbol(x), Symbol(y))
Between(x::Union{Int, Symbol}, y::AbstractString) = Between(x, Symbol(y))
Between(x::AbstractString, y::Union{Int, Symbol}) = Between(Symbol(x), y)

"""
    All()

Select all columns.
"""
struct All{T<:Tuple}
    cols::T
    function All(args...)
        if !isempty(args)
            Base.depwarn("All(args...) is deprecated, use Cols(args...) instead", :All)
        end
        return new{typeof(args)}(args)
    end
end

"""
    Cols(cols...; operator=union)
    Cols(f::Function; operator=union)

Select columns matching specifications in `cols`. If `cols == ()`, select no columns.

If the only positional argument is a `Function` `f` then select the columns whose
names passed to the `f` predicate as strings return `true`.

When multiple `cols` selectors are passed then the sets of columns selected by them
are passed to `operator` as positional arguments.
`operator` should be a set operation function, like `union`, `intersect`, `setdiff`, and `symdiff`
defined in Base Julia. By default `operator=union` in which case all columns matching
at least one selector are returned.
"""
struct Cols{T<:Tuple}
    cols::T
    operator
    Cols(args...; operator=union) = new{typeof(args)}(args, operator)
end

"""
    BroadcastedSelector(selector)

Wrapper type around a `Between`, `All` or `Cols` indicating that
an operation should be applied to each column included by the wrapped selector.

# Examples
```jldoctest
julia> using DataAPI

julia> DataAPI.Between(:a, :e) .=> sin
DataAPI.BroadcastedSelector{DataAPI.Between{Symbol, Symbol}}(DataAPI.Between{Symbol, Symbol}(:a, :e)) => sin

julia> DataAPI.Cols(r"x") .=> [sum, prod]
2-element Vector{Pair{DataAPI.BroadcastedSelector{DataAPI.Cols{Tuple{Regex}}}, _A} where _A}:
 DataAPI.BroadcastedSelector{DataAPI.Cols{Tuple{Regex}}}(DataAPI.Cols{Tuple{Regex}}((r"x",))) => sum
 DataAPI.BroadcastedSelector{DataAPI.Cols{Tuple{Regex}}}(DataAPI.Cols{Tuple{Regex}}((r"x",))) => prod
```
"""
struct BroadcastedSelector{T}
    sel::T
    BroadcastedSelector(sel) = new{typeof(sel)}(sel)
end

Base.Broadcast.broadcastable(x::Between) = Ref(BroadcastedSelector(x))
Base.Broadcast.broadcastable(x::All) = Ref(BroadcastedSelector(x))
Base.Broadcast.broadcastable(x::Cols) = Ref(BroadcastedSelector(x))

"""
    unwrap(x)

For a given scalar argument `x`, potentially "unwrap" it to return the base wrapped value.
Useful as a generic API for wrapper types when the original value is needed.

The default definition just returns `x` itself, i.e. no unwrapping is performned.

This generic function is owned by DataAPI.jl itself, which is the sole provider of the
default definition.
"""
function unwrap end
unwrap(x) = x

# The database-style join methods for tabular data type.
# The common interface is `*join(x, y; ...)` and use the keyword arguments
# for the join criteria.
# See the design of DataFrames.jl also.
function innerjoin end
function outerjoin end
function rightjoin end
function leftjoin end
function semijoin end
function antijoin end
function crossjoin end

"""
    nrow(t)

Return the number of rows of table `t`.
"""
function nrow end

"""
    ncol(t)

Return the number of columns of table `t`.
"""
function ncol end

"""
    allcombinations(sink, ...)

Create table from all combinations of values in passed arguments
using a `sink` function to materialize the table.
"""
function allcombinations end

const STYLE_INFO = """
One of the uses of the metadata `style` is decision
how the metadata should be propagated when `x` is transformed. This interface
defines the `:default` style that indicates that metadata should not be propagated
under any operations (it is only preserved when a copy of the source table is
performed). All types supporting metadata allow at least this style.
"""

const COL_INFO = """
`col` must have a type that is supported by table `x` for column indexing.
Following the Tables.jl contract `Symbol` and `Int` are always allowed.
Throw an error if `col` is not a column of `x`.
"""

"""
    metadatasupport(T::Type)

Return a `NamedTuple{(:read, :write), Tuple{Bool, Bool}}` indicating whether
values of type `T` support metadata.

The `read` field indicates whether reading metadata with the [`metadata`](@ref)
and [`metadatakeys`]](@ref) functions is supported.

The `write` field indicates whether modifying metadata with the [`metadata!`](@ref),
[`deletemetadata!`](@ref), and [`emptymetadata!`](@ref) functions is supported.
"""
metadatasupport(::Type) = (read=false, write=false)

"""
    colmetadatasupport(T::Type)

Return a `NamedTuple{(:read, :write), Tuple{Bool, Bool}}` indicating whether
values of type `T` support column metadata.

The `read` field indicates whether reading metadata with the [`colmetadata`](@ref)
and [`colmetadatakeys`](@ref) functions is supported.

The `write` field indicates whether modifying metadata with the [`colmetadata!`](@ref),
[`deletecolmetadata!`](@ref), and [`emptycolmetadata!`](@ref) functions is supported.
"""
colmetadatasupport(::Type) = (read=false, write=false)

"""
    metadata(x, key::AbstractString, [default]; style::Bool=false)

Return metadata value associated with object `x` for key `key`. Throw an error
if `x` does not support reading metadata or does not have a mapping for `key`.

If `style=true` return a tuple of metadata value and metadata style. Metadata
style is an additional information about the kind of metadata that is stored for
the `key`.

$STYLE_INFO

If `default` is passed then return it if reading metadata is supported but
mapping for `key` is missing. If `style=true` return `(default, :default)`.
"""
function metadata end

"""
    metadata(x; style::Bool=false)

Return a dictionary mapping all metadata keys to metadata values associated
with object `x`. Throw an error if `x` does not support reading metadata.

If `style=true` values are tuples of metadata value and metadata style. Metadata
style is an additional information about the kind of metadata that is stored for
the `key`.

$STYLE_INFO

The returned dictionary may be freshly allocated on each call to `metadata` and
is considered to be owned by `x` so it must not be modified.
"""
function metadata(x::T; style::Bool=false) where {T}
    if !metadatasupport(T).read
        throw(ArgumentError("Objects of type $T do not support reading metadata"))
    end
    return Dict(key => metadata(x, key, style=style) for key in metadatakeys(x))
end

"""
    metadatakeys(x)

Return an iterator of metadata keys for which `metadata(x, key)` returns a
metadata value.
Throw an error if `x` does not support reading metadata.
"""
function metadatakeys end

"""
    metadata!(x, key::AbstractString, value; style::Symbol=:default)

Set metadata for object `x` for key `key` to have value `value`
and style `style` (`:default` by default) and return `x`.
Throw an error if `x` does not support setting metadata.

$STYLE_INFO
"""
function metadata! end

"""
    deletemetadata!(x, key::AbstractString)

Delete metadata for object `x` for key `key` and return `x`
(if metadata for `key` is not present do not perform any action).
Throw an error if `x` does not support metadata deletion.
"""
function deletemetadata! end

"""
    emptymetadata!(x)

Delete all metadata for object `x`.
Throw an error if `x` does not support metadata deletion.
"""
function emptymetadata! end

"""
    colmetadata(x, col, key::AbstractString, [default]; style::Bool=false)

Return metadata value associated with table `x` for column `col` and key `key`.
Throw an error if `x` does not support reading metadata for column `col` or `x`
supports reading metadata, but does not have a mapping for column `col` for `key`.

$COL_INFO

If `style=true` return a tuple of metadata value and metadata style. Metadata
style is an additional information about the kind of metadata that is stored for
the `key`.

$STYLE_INFO

If `default` is passed then return it if `x` supports reading metadata and has
column `col` but mapping for `key` is missing.
If `style=true` return `(default, :default)`.
"""
function colmetadata end

"""
    colmetadata(x, [col]; style::Bool=false)

If `col` is not passed return a dictionary mapping columns represented as
`Symbol` that have associated metadata to dictionaries mapping all
metadata keys to metadata values associated with table `x` for a given column.

If `col` is passed return a dictionary mapping all column metadata keys to metadata values
associated with column `col` of table `x`. Throw an error if `x` does not
support reading metadata for column `col` or column `col` is not present in `x`.

If `style=true` values are tuples of metadata value and metadata style. Metadata
style is an additional information about the kind of metadata that is stored for
the `key`.

$STYLE_INFO

The returned dictionary may be freshly allocated on each call to `colmetadata`
and is considered to be owned by `x` so it must not be modified.
"""
function colmetadata(x::T, col; style::Bool=false) where {T}
    if !colmetadatasupport(T).read
        throw(ArgumentError("Objects of type $T do not support reading column metadata"))
    end
    return Dict(key => colmetadata(x, col, key, style=style) for key in colmetadatakeys(x, col))
end

function colmetadata(x::T; style::Bool=false) where {T}
    if !colmetadatasupport(T).read
        throw(ArgumentError("Objects of type $T do not support reading column metadata"))
    end
    return Dict(col => Dict(key => colmetadata(x, col, key, style=style) for key in keys)
                for (col, keys) in colmetadatakeys(x))
end

"""
    colmetadatakeys(x, [col])

If `col` is passed return an iterator of metadata keys for which
`metadata(x, col, key)` returns a metadata value. Throw an error if `x` does not
support reading column metadata or if `col` is not a column of `x`.

`col` must have a type that is supported by table `x` for column indexing.
Following the Tables.jl contract `Symbol` and `Int` are always allowed.

If `col` is not passed return an iterator of `col => colmetadatakeys(x, col)`
pairs for all columns that have metadata, where `col` are `Symbol`.
If `x` does not support column metadata return `()`.
"""
function colmetadatakeys end

"""
    colmetadata!(x, col, key::AbstractString, value; style::Symbol=:default)

Set metadata for table `x` for column `col` for key `key` to have value `value`
and style `style` (`:default` by default) and return `x`.
Throw an error if `x` does not support setting metadata for column `col`.

$COL_INFO

$STYLE_INFO
"""
function colmetadata! end

"""
    deletecolmetadata!(x, col, key::AbstractString)

Delete metadata for table `x` for column `col` for key `key` and return `x`
(if metadata for `key` is not present do not perform any action).
Throw an error if `x` does not support metadata deletion for column `col`.
"""
function deletecolmetadata! end

"""
    emptycolmetadata!(x, [col])

Delete all metadata for table `x` for column `col`.
If `col` is not passed delete all column level metadata for table `x`.
Throw an error if `x` does not support metadata deletion for column `col`.
"""
function emptycolmetadata! end

"""
    rownumber(row)

Return the row number of `row` in the source table.
"""
function rownumber end

# To avoid type piracy and method ambiguities, implementations of `groupby` 
# must restrict the first argument to a type defined in the same package. 
function groupby end

end # module
