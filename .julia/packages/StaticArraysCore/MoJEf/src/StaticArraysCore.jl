module StaticArraysCore

export SArray, SMatrix, SVector
export MArray, MMatrix, MVector
export SizedArray, SizedMatrix, SizedVector
export FieldArray, FieldMatrix, FieldVector
export Size

"""
    abstract type StaticArray{S, T, N} <: AbstractArray{T, N}
    StaticScalar{T}     = StaticArray{Tuple{}, T, 0}
    StaticVector{N,T}   = StaticArray{Tuple{N}, T, 1}
    StaticMatrix{N,M,T} = StaticArray{Tuple{N,M}, T, 2}

`StaticArray`s are Julia arrays with fixed, known size.

## Dev docs

They must define the following methods:
 - Constructors that accept a flat tuple of data.
 - `getindex()` with an integer (linear indexing) (preferably `@inline` with `@boundscheck`).
 - `Tuple()`, returning the data in a flat Tuple.

It may be useful to implement:

- `similar_type(::Type{MyStaticArray}, ::Type{NewElType}, ::Size{NewSize})`, returning a
  type (or type constructor) that accepts a flat tuple of data.

For mutable containers you may also need to define the following:

 - `setindex!` for a single element (linear indexing).
 - `similar(::Type{MyStaticArray}, ::Type{NewElType}, ::Size{NewSize})`.
 - In some cases, a zero-parameter constructor, `MyStaticArray{...}()` for uninitialized data
   is assumed to exist.

See also
[`SVector`](@ref StaticArraysCore.SVector),
[`SMatrix`](@ref StaticArraysCore.SMatrix),
[`SArray`](@ref StaticArraysCore.SArray),
[`MVector`](@ref), [`MMatrix`](@ref), [`MArray`](@ref),
[`SizedArray`](@ref),
[`FieldVector`](@ref StaticArraysCore.FieldVector),
[`FieldMatrix`](@ref), and [`FieldArray`](@ref)
"""
abstract type StaticArray{S <: Tuple, T, N} <: AbstractArray{T, N} end
const StaticScalar{T} = StaticArray{Tuple{}, T, 0}
const StaticVector{N, T} = StaticArray{Tuple{N}, T, 1}
const StaticMatrix{N, M, T} = StaticArray{Tuple{N, M}, T, 2}
const StaticVecOrMat{T} = Union{StaticVector{<:Any, T}, StaticMatrix{<:Any, <:Any, T}}

# The ::Tuple variants exist to make sure that anything that calls with a tuple
# instead of a Tuple gets through to the constructor, so the user gets a nice
# error message
tuple_length(T::Tuple) = length(T)
tuple_prod(T::Tuple) = prod(T)
tuple_minimum(T::Tuple) = minimum(T)

tuple_svec(::Type{T}) where {T<:Tuple} = T.parameters

# Julia `Base` provides a function for this purpose since v1.1: `fieldtypes`.
# We don't use it yet, though, because it wrecks type inference with Julia
# v1.6.
Base.@pure tuple_tuple(::Type{T}) where {T<:Tuple} = (tuple_svec(T)...,)

tuple_length(::Type{<:NTuple{n, Any}}) where {n} = n
tuple_prod(::Type{T}) where {T<:Tuple} = mapreduce(Int, *, tuple_tuple(T), init = 1)
tuple_minimum(::Type{Tuple{}}) = 0
tuple_minimum(::Type{T}) where {T<:Tuple{Any,Vararg}} = mapreduce(Int, min, tuple_tuple(T), init = typemax(Int))

"""
    size_to_tuple(::Type{S}) where S<:Tuple

Converts a size given by `Tuple{N, M, ...}` into a tuple `(N, M, ...)`.
"""
size_to_tuple(::Type{T}) where {T<:Tuple} = tuple_tuple(T)

# Something doesn't match up type wise
@generated function check_array_parameters(::Type{Size}, ::Type{T}, ::Type{Val{N}}, ::Type{Val{L}}) where {Size,T,N,L}
    if !all(x->isa(x, Int), Size.parameters)
        return :(throw(ArgumentError("Static Array parameter Size must be a tuple of Ints (e.g. `SArray{Tuple{3,3}}` or `SMatrix{3,3}`).")))
    end

    if L != tuple_prod(Size) || L < 0 || tuple_minimum(Size) < 0 || tuple_length(Size) != N
        return :(throw(ArgumentError("Size mismatch in Static Array parameters. Got size $Size, dimension $N and length $L.")))
    end

    return nothing
end

# Cast any Tuple to an TupleN{T}
@generated function convert_ntuple(::Type{T}, d::NTuple{N,Any}) where {N,T}
    exprs = ntuple(i -> :(convert(T, d[$i])), Val(N))
    return quote
        Base.@_inline_meta
        $(Expr(:tuple, exprs...))
    end
end


"""
    SArray{S, T, N, L}(x::NTuple{L})
    SArray{S, T, N, L}(x1, x2, x3, ...)

Construct a statically-sized array `SArray`. Since this type is immutable, the data must be
provided upon construction and cannot be mutated later. The `S` parameter is a Tuple-type
specifying the dimensions, or size, of the array - such as `Tuple{3,4,5}` for a 3×4×5-sized
array. The `N` parameter is the dimension of the array; the `L` parameter is the `length`
of the array and is always equal to `prod(S)`. Constructors may drop the `L`, `N` and `T`
parameters if they are inferrable from the input (e.g. `L` is always inferrable from `S`).

    SArray{S}(a::Array)

Construct a statically-sized array of dimensions `S` (expressed as a `Tuple{...}`) using
the data from `a`. The `S` parameter is mandatory since the size of `a` is unknown to the
compiler (the element type may optionally also be specified).
"""
struct SArray{S <: Tuple, T, N, L} <: StaticArray{S, T, N}
    data::NTuple{L,T}

    function SArray{S, T, N, L}(x::NTuple{L,T}) where {S<:Tuple, T, N, L}
        check_array_parameters(S, T, Val{N}, Val{L})
        new{S, T, N, L}(x)
    end

    function SArray{S, T, N, L}(x::NTuple{L,Any}) where {S<:Tuple, T, N, L}
        check_array_parameters(S, T, Val{N}, Val{L})
        new{S, T, N, L}(convert_ntuple(T, x))
    end
end

@inline SArray{S,T,N}(x::Tuple) where {S<:Tuple,T,N} = SArray{S,T,N,tuple_prod(S)}(x)

"""
    SVector{S, T}(x::NTuple{S, T})
    SVector{S, T}(x1, x2, x3, ...)

Construct a statically-sized vector `SVector`. Since this type is immutable,
the data must be provided upon construction and cannot be mutated later.
Constructors may drop the `T` and `S` parameters if they are inferrable from the
input (e.g. `SVector(1,2,3)` constructs an `SVector{3, Int}`).

    SVector{S}(vec::Vector)

Construct a statically-sized vector of length `S` using the data from `vec`.
The parameter `S` is mandatory since the length of `vec` is unknown to the
compiler (the element type may optionally also be specified).
"""
const SVector{S, T} = SArray{Tuple{S}, T, 1, S}

"""
    SMatrix{S1, S2, T, L}(x::NTuple{L, T})
    SMatrix{S1, S2, T, L}(x1, x2, x3, ...)

Construct a statically-sized matrix `SMatrix`. Since this type is immutable,
the data must be provided upon construction and cannot be mutated later. The
`L` parameter is the `length` of the array and is always equal to `S1 * S2`.
Constructors may drop the `L`, `T` and even `S2` parameters if they are inferrable
from the input (e.g. `L` is always inferrable from `S1` and `S2`).

    SMatrix{S1, S2}(mat::Matrix)

Construct a statically-sized matrix of dimensions `S1 × S2` using the data from
`mat`. The parameters `S1` and `S2` are mandatory since the size of `mat` is
unknown to the compiler (the element type may optionally also be specified).
"""
const SMatrix{S1, S2, T, L} = SArray{Tuple{S1, S2}, T, 2, L}

# MArray

"""
    MArray{S, T, N, L}(undef)
    MArray{S, T, N, L}(x::NTuple{L})
    MArray{S, T, N, L}(x1, x2, x3, ...)

Construct a statically-sized, mutable array `MArray`. The data may optionally be
provided upon construction and can be mutated later. The `S` parameter is a Tuple-type
specifying the dimensions, or size, of the array - such as `Tuple{3,4,5}` for a 3×4×5-sized
array. The `N` parameter is the dimension of the array; the `L` parameter is the `length`
of the array and is always equal to `prod(S)`. Constructors may drop the `L`, `N` and `T`
parameters if they are inferrable from the input (e.g. `L` is always inferrable from `S`).

    MArray{S}(a::Array)

Construct a statically-sized, mutable array of dimensions `S` (expressed as a `Tuple{...}`)
using the data from `a`. The `S` parameter is mandatory since the size of `a` is unknown to
the compiler (the element type may optionally also be specified).
"""
mutable struct MArray{S <: Tuple, T, N, L} <: StaticArray{S, T, N}
    data::NTuple{L,T}

    function MArray{S,T,N,L}(x::NTuple{L,T}) where {S<:Tuple,T,N,L}
        check_array_parameters(S, T, Val{N}, Val{L})
        new{S,T,N,L}(x)
    end

    function MArray{S,T,N,L}(x::NTuple{L,Any}) where {S<:Tuple,T,N,L}
        check_array_parameters(S, T, Val{N}, Val{L})
        new{S,T,N,L}(convert_ntuple(T, x))
    end

    function MArray{S,T,N,L}(::UndefInitializer) where {S<:Tuple,T,N,L}
        check_array_parameters(S, T, Val{N}, Val{L})
        new{S,T,N,L}()
    end
end

@inline MArray{S,T,N}(x::Tuple) where {S<:Tuple,T,N} = MArray{S,T,N,tuple_prod(S)}(x)

"""
    MVector{S,T}(undef)
    MVector{S,T}(x::NTuple{S, T})
    MVector{S,T}(x1, x2, x3, ...)

Construct a statically-sized, mutable vector `MVector`. Data may optionally be
provided upon construction, and can be mutated later. Constructors may drop the
`T` and `S` parameters if they are inferrable from the input (e.g.
`MVector(1,2,3)` constructs an `MVector{3, Int}`).

    MVector{S}(vec::Vector)

Construct a statically-sized, mutable vector of length `S` using the data from
`vec`. The parameter `S` is mandatory since the length of `vec` is unknown to the
compiler (the element type may optionally also be specified).
"""
const MVector{S, T} = MArray{Tuple{S}, T, 1, S}

"""
    MMatrix{S1, S2, T, L}(undef)
    MMatrix{S1, S2, T, L}(x::NTuple{L, T})
    MMatrix{S1, S2, T, L}(x1, x2, x3, ...)

Construct a statically-sized, mutable matrix `MMatrix`. The data may optionally
be provided upon construction and can be mutated later. The `L` parameter is the
`length` of the array and is always equal to `S1 * S2`. Constructors may drop
the `L`, `T` and even `S2` parameters if they are inferrable from the input
(e.g. `L` is always inferrable from `S1` and `S2`).

    MMatrix{S1, S2}(mat::Matrix)

Construct a statically-sized, mutable matrix of dimensions `S1 × S2` using the data from
`mat`. The parameters `S1` and `S2` are mandatory since the size of `mat` is
unknown to the compiler (the element type may optionally also be specified).
"""
const MMatrix{S1, S2, T, L} = MArray{Tuple{S1, S2}, T, 2, L}


# SizedArray

"""
    SizedArray{Tuple{dims...}}(array)

Wraps an `AbstractArray` with a static size, so to take advantage of the (faster)
methods defined by StaticArrays.jl. The size is checked once upon
construction to determine if the number of elements (`length`) match, but the
array may be reshaped.

The aliases `SizedVector{N}` and `SizedMatrix{N,M}` are provided as more
convenient names for one and two dimensional `SizedArray`s. For example, to
wrap a 2×3 array `a` in a `SizedArray`, use `SizedMatrix{2,3}(a)`.
"""
struct SizedArray{S<:Tuple,T,N,M,TData<:AbstractArray{T,M}} <: StaticArray{S,T,N}
    data::TData

    function SizedArray{S,T,N,M,TData}(a::TData) where {S<:Tuple,T,N,M,TData<:AbstractArray{T,M}}
        Base.require_one_based_indexing(a)
        if size(a) != size_to_tuple(S) && size(a) != (tuple_prod(S),)
            @static if VERSION >= v"1.8"
                throw(DimensionMismatch(lazy"Dimensions $(size(a)) don't match static size $S"))
            else
                throw(DimensionMismatch("Dimensions $(size(a)) don't match static size $S"))
            end
        end
        return new{S,T,N,M,TData}(a)
    end

    function SizedArray{S,T,N,1,TData}(::UndefInitializer) where {S<:Tuple,T,N,TData<:AbstractArray{T,1}}
        return new{S,T,N,1,TData}(TData(undef, tuple_prod(S)))
    end
    function SizedArray{S,T,N,N,TData}(::UndefInitializer) where {S<:Tuple,T,N,TData<:AbstractArray{T,N}}
        return new{S,T,N,N,TData}(TData(undef, size_to_tuple(S)...))
    end
end


"""
    SizedVector{S, T} = SizedArray{Tuple{S}, T, 1, 1}

Wraps a one-dimensional `AbstractArray` with static length `S` and element type `T`,
leveraging the performance optimizations of StaticArrays.jl.

For detailed usage and functionality, refer to the documentation of `SizedArray`.
"""
const SizedVector{S,T} = SizedArray{Tuple{S},T,1,1}

"""
    SizedMatrix{S1,S2,T} = SizedArray{Tuple{S1,S2},T,2}

Wraps a two-dimensional `AbstractArray` with static dimensions `S1` by `S2` and element type `T`,
leveraging the performance optimizations of StaticArrays.jl.

For detailed usage and functionality, refer to the documentation of `SizedArray`.
"""
const SizedMatrix{S1,S2,T} = SizedArray{Tuple{S1,S2},T,2}

# FieldArray

"""
    abstract type FieldArray{N, T, D} <: StaticArray{N, T, D}

Inheriting from this type will make it easy to create your own rank-D tensor types. A `FieldArray`
will automatically define `getindex` and `setindex!` appropriately. An immutable
`FieldArray` will be as performant as an `SArray` of similar length and element type,
while a mutable `FieldArray` will behave similarly to an `MArray`.

Note that you must define the fields of any `FieldArray` subtype in column major order. If you
want to use an alternative ordering you will need to pay special attention in providing your
own definitions of `getindex`, `setindex!` and tuple conversion.

If you define a `FieldArray` which is parametric on the element type you should
consider defining `similar_type` as in the `FieldVector` example.


# Example

```julia
struct Stiffness <: FieldArray{Tuple{2,2,2,2}, Float64, 4}
    xxxx::Float64
    yxxx::Float64
    xyxx::Float64
    yyxx::Float64
    xxyx::Float64
    yxyx::Float64
    xyyx::Float64
    yyyx::Float64
    xxxy::Float64
    yxxy::Float64
    xyxy::Float64
    yyxy::Float64
    xxyy::Float64
    yxyy::Float64
    xyyy::Float64
    yyyy::Float64
end
```
"""
abstract type FieldArray{N, T, D} <: StaticArray{N, T, D} end

"""
    abstract type FieldMatrix{N1, N2, T} <: FieldArray{Tuple{N1, N2}, 2}

Inheriting from this type will make it easy to create your own rank-two tensor types. A `FieldMatrix`
will automatically define `getindex` and `setindex!` appropriately. An immutable
`FieldMatrix` will be as performant as an `SMatrix` of similar length and element type,
while a mutable `FieldMatrix` will behave similarly to an `MMatrix`.

Note that the fields of any subtype of `FieldMatrix` must be defined in column
major order unless you are willing to implement your own `getindex`.

If you define a `FieldMatrix` which is parametric on the element type you
should consider defining `similar_type` as in the `FieldVector` example.

# Example

```julia
struct Stress <: FieldMatrix{3, 3, Float64}
    xx::Float64
    yx::Float64
    zx::Float64
    xy::Float64
    yy::Float64
    zy::Float64
    xz::Float64
    yz::Float64
    zz::Float64
end
```

Note that the fields of any subtype of `FieldMatrix` must be defined in column major order.
This means that formatting of constructors for literal `FieldMatrix` can be confusing. For example

```julia-repl
julia> sigma = Stress(1.0, 2.0, 3.0,
                      4.0, 5.0, 6.0,
                      7.0, 8.0, 9.0)

3×3 Stress with indices SOneTo(3)×SOneTo(3):
 1.0  4.0  7.0
 2.0  5.0  8.0
 3.0  6.0  9.0
```

will give you the transpose of what the multi-argument formatting suggests. For clarity,
you may consider using the alternative

```julia-repl
julia> sigma = Stress(SA[1.0 2.0 3.0;
                         4.0 5.0 6.0;
                         7.0 8.0 9.0])

3×3 Stress with indices SOneTo(3)×SOneTo(3):
 1.0  2.0  3.0
 4.0  5.0  6.0
 7.0  8.0  9.0
```
"""
abstract type FieldMatrix{N1, N2, T} <: FieldArray{Tuple{N1, N2}, T, 2} end

"""
    abstract type FieldVector{N, T} <: FieldArray{Tuple{N}, 1}

Inheriting from this type will make it easy to create your own vector types. A `FieldVector`
will automatically define `getindex` and `setindex!` appropriately. An immutable
`FieldVector` will be as performant as an `SVector` of similar length and element type,
while a mutable `FieldVector` will behave similarly to an `MVector`.

If you define a `FieldVector` which is parametric on the element type you
should consider defining `similar_type` to preserve your array type through
array operations as in the example below.

# Example

```julia
struct Vec3D{T} <: FieldVector{3, T}
    x::T
    y::T
    z::T
end

StaticArrays.similar_type(::Type{<:Vec3D}, ::Type{T}, s::Size{(3,)}) where {T} = Vec3D{T}
```
"""
abstract type FieldVector{N, T} <: FieldArray{Tuple{N}, T, 1} end

# Add a new BroadcastStyle for StaticArrays, derived from AbstractArrayStyle
# A constructor that changes the style parameter N (array dimension) is also required
struct StaticArrayStyle{N} <: Base.Broadcast.AbstractArrayStyle{N} end
StaticArrayStyle{M}(::Val{N}) where {M,N} = StaticArrayStyle{N}()

"""
    similar_type(static_array)
    similar_type(static_array, T)
    similar_type(array, ::Size)
    similar_type(array, T, ::Size)

Returns a constructor for a statically-sized array similar to the input array
(or type) `static_array`/`array`, optionally with different element type `T` or size
`Size`. If the input `array` is not a `StaticArray` then the `Size` is mandatory.

This differs from `similar()` in that the resulting array type may not be
mutable (or define `setindex!()`), and therefore the returned type may need to
be *constructed* with its data.

Note that the (optional) size *must* be specified as a static `Size` object (so the compiler
can infer the result statically).

New types should define the signature `similar_type(::Type{A},::Type{T},::Size{S}) where {A<:MyType,T,S}`
if they wish to overload the default behavior.
"""
function similar_type end

"""
    Dynamic()

Used to signify that a dimension of an array is not known statically.
"""
struct Dynamic end

const StaticDimension = Union{Int, Dynamic}

"""
    Size(dims::Int...)

`Size` is used extensively throughout the `StaticArrays` API to describe _compile-time_
knowledge of the size of an array. The dimensions are stored as a type parameter and are
statically propagated by the compiler, resulting in efficient, type-inferrable code. For
example, to create a static matrix of zeros, use `A = zeros(SMatrix{3,3})`. The static
size of `A` can be obtained by `Size(A)`. (rather than `size(zeros(3,3))`, which returns
`Base.Tuple{2,Int}`).

Note that if dimensions are not known statically (e.g., for standard `Array`s),
[`Dynamic()`](@ref) should be used instead of an `Int`.

    Size(a::AbstractArray)
    Size(::Type{T<:AbstractArray})

The `Size` constructor can be used to extract static dimension information from a given
array. For example:

```jldoctest
julia> Size(zeros(SMatrix{3, 4}))
Size(3, 4)

julia> Size(zeros(3, 4))
Size(StaticArraysCore.Dynamic(), StaticArraysCore.Dynamic())
```

This has multiple uses, including "trait"-based dispatch on the size of a statically-sized
array. For example:

```julia
det(x::StaticMatrix) = _det(Size(x), x)
_det(::Size{(1,1)}, x::StaticMatrix) = x[1,1]
_det(::Size{(2,2)}, x::StaticMatrix) = x[1,1]*x[2,2] - x[1,2]*x[2,1]
# and other definitions as necessary
```

"""
struct Size{S}
    function Size{S}() where {S}
        new{S::Tuple{Vararg{StaticDimension}}}()
    end
end

Base.@pure Size(s::Tuple{Vararg{StaticDimension}}) = Size{s}()
Base.@pure Size(s::StaticDimension...) = Size{s}()
Size() = Size(())
Size(::Type{T}) where {T<:Tuple} = Size{tuple_tuple(T)}()

Base.show(io::IO, ::Size{S}) where {S} = print(io, "Size", S)

function missing_size_error(::Type{SA}) where SA
    error("""
        The size of type `$SA` is not known.

        If you were trying to construct (or `convert` to) a `StaticArray` you
        may need to add the size explicitly as a type parameter so its size is
        inferrable to the Julia compiler (or performance would be terrible). For
        example, you might try

            m = zeros(3,3)
            SMatrix(m)            # this error
            SMatrix{3,3}(m)       # correct - size is inferrable
            SArray{Tuple{3,3}}(m) # correct, note Tuple{3,3}
        """)
end

Size(a::T) where {T<:AbstractArray} = Size(T)
Size(::Type{SA}) where {SA <: StaticArray} = missing_size_error(SA)
Size(::Type{SA}) where {SA <: StaticArray{S}} where {S<:Tuple} = @isdefined(S) ? Size(S) : missing_size_error(SA)

Size(::Type{<:AbstractArray{<:Any, N}}) where {N} = Size(ntuple(_ -> Dynamic(), Val(N)))

end # module
