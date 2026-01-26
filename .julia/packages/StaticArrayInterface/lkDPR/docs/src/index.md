```@meta
CurrentModule = StaticArrayInterface
```

# StaticArrayInterface

Designs for new Base array interface primitives, used widely through scientific machine learning (SciML) and other organizations

## StaticArrayInterfaceCore

StaticArrayInterfaceCore is a smaller set of the StaticArrayInterface setup which defines the subset which has no compile time impact.
This for example includes simple functions like `StaticArrayInterfaceCore.zeromatrix` which have simple few dispatch definitions
and no dependency on other libraries such as Static.jl. Notably, Static.jl currently has issues with invalidations
(https://github.com/SciML/Static.jl/issues/52), and thus anything with static outputs are in the domain of StaticArrayInterface.jl
proprer.

## Extenion Packages

StaticArrayInterface.jl uses extension packages in order to add support for popular libraries to its interface functions. These packages are:

- StaticArrays.jl
- OffsetArrays.jl

## Inheriting Array Traits

Creating an array type with unique behavior in Julia is often accomplished by creating a lazy wrapper around previously defined array types 
(e.g. [composition by inheritance](https://en.wikipedia.org/wiki/Composition_over_inheritance)).
This allows the new array type to inherit functionality by redirecting methods to the parent array (e.g., `Base.size(x::Wrapper) = size(parent(x))`).
Generic design limits the need to define an excessive number of methods like this.
However, methods used to describe a type's traits often need to be explicitly defined for each trait method.
If the the underlying data and access to it are unchanged by it's wrapper the 
[`StaticArrayInterface.is_forwarding_wrapper`](@ref) trait can signal to other trait methods to access its parent data structure.
Supporting this for a new type only requires defines these methods:

```julia
StaticArrayInterface.is_forwarding_wrapper(::Type{<:NewType}) = true
StaticArrayInterface.parent_type(::Type{<:NewType}) = NewTypeParent
Base.parent(x::NewType) = x.parent
```

For those authoring new trait methods, this may change the default definition from `has_trait(::Type{T}) where {T} = false`, to:
```julia
function has_trait(::Type{T}) where {T}
    if is_forwarding_wrapper(T)
        return false
    else
        return has_trait(parent_type(T))
    end
end
```

Most traits in `StaticArrayInterface` are a variant on this pattern.
If the trait in question may be altered by a wrapper array, this pattern should be altered or may be inappropriate.

## Static Traits

The size along one or more dimensions of an array may be known at compile time.
`StaticArrayInterface.known_size` is useful for extracting this information from array types and `StaticArrayInterface.static_size` is useful for extracting this information from an instance of an array.
For example:

```julia
julia> a = ones(3)';

julia> StaticArrayInterface.static_size(a)
(static(1), 3)

julia> StaticArrayInterface.known_size(typeof(a))
(1, nothing)

```

This is useful for dispatching on known information about the size of an array:
```julia
fxn(x) = _fxn(StaticArrayInterface.static_size(x), x)
_fxn(sz::Tuple{StaticInt{S1},StaticInt{S2}}, x) where {S1,S2} = ...
_fxn(sz::Tuple{StaticInt{3},StaticInt{3}}, x) = ...
_fxn(sz::Tuple{Int,StaticInt{S2}}, x) where {S2} = ...
_fxn(sz::Tuple{StaticInt{S1},Int}, x) where {S1} = ...
_fxn(sz::Tuple{Int,Int}, x) = ...
```

Methods should avoid forcing conversion to static sizes when dynamic sizes could potentially be returned.
Fore example, `fxn(x) = _fxn(Static.static(StaticArrayInterface.static_size(x)), x)` would result in dynamic dispatch if `x` is an instance of `Matrix`.
Additionally, `StaticArrayInterface.static_size` should only be used outside of generated functions to avoid possible world age issues.

Generally, `StaticArrayInterface.static_size` uses the return of `known_size` to form a static value for those dimensions with known length and only queries dimensions corresponding to `nothing`.
For example, the previous example had a known size of `(1, nothing)`.
Therefore, `StaticArrayInterface.static_size` would have compile time information about the first dimension returned as `static(1)` and would only look up the size of the second dimension at run time.
This means the above example `StaticArrayInterface.static_size(a)` would lower to code similar to this at compile time: `Static.StaticInt(1), Base.arraysize(x, 1)`.
Generic support for `StaticArrayInterface.known_size` relies on calling `known_length` for each type returned from `axes_types`.
Therefore, the recommended approach for supporting static sizing in newly defined array types is defining a new `axes_types` method.

Static information related to subtypes of `AbstractRange` include `known_length`, `known_first`, `known_step`, and `known_last`.

## Dimensions

Methods such as `size(x, dim)` need to map `dim` to the dimensions of `x`.
Typically, `dim` is an `Int` with an invariant mapping to the dimensions of `x`.
Some methods accept `:` or a tuple of dimensions as an argument.
`StaticArrayInterface` also considers `StaticInt` a viable dimension argument.

[`StaticArrayInterface.to_dims`](@ref) helps ensure that `dim` is converted to a viable dimension mapping in a manner that helps with type stability.
For example, all `Integers` passed to `to_dims` are converted to `Int` (unless `dim` is a `StaticInt`).
This is also useful for arrays that uniquely label dimensions, in which case `to_dims` serves as a safe point of hooking into existing methods with dimension arguments.
`StaticArrayInterface` also defines native `Symbol` to `Int` and `StaticSymbol` to `StaticInt` mapping  for arrays defining [`StaticArrayInterface.dimnames`](@ref).

Methods requiring dimension specific arguments should use some variation of the following pattern.

```julia
f(x, dim) = f(x, StaticArrayInterface.to_dims(x, dim))
f(x, dim::Int) = ...
f(x, dim::StaticInt) = ...
```

If `x`'s first dimension is named `:dim_1` then calling `f(x, :dim_1)` would result in `f(x, 1)`.
If users knew they always wanted to call `f(x, 2)` then they could define `h(x) = f(x, static(2))`, ensuring `f` passes along that information while compiling.

New types defining dimension names can do something similar to:

```julia
using Static
using StaticArrayInterface

struct StaticDimnames{dnames} end  # where dnames::Tuple{Vararg{Symbol}}

StaticArrayInterface.known_dimnames(::Type{StaticDimnames{dnames}}) where {dnames} = dnames
StaticArrayInterface.dimnames(::StaticDimnames{dnames}) where {dnames} = static(dnames)

struct DynamicDimnames{N}
    dimnames::NTuple{N,Symbol}
end
StaticArrayInterface.known_dimnames(::Type{DynamicDimnames{N}}) where {N} = ntuple(_-> nothing, Val(N))
StaticArrayInterface.dimnames(x::DynamicDimnames) = getfield(x, :dimnames)

```

Notice that `DynamicDimnames` returns `nothing` instead of a symbol for each dimension.
This indicates dimension names are present for `DynamicDimnames` but that information is nothing at compile time.

Dimension names should be appropriately propagated between nested arrays using `StaticArrayInterface.to_parent_dims`. 
This allows types such as `SubArray` and `PermutedDimsArray` to work with named dimensions.
Similarly, other methods that return information corresponding to dimensions (e.g., `ArrayInterfce.size`, `StaticArrayInterface.static_axes`) use `to_parent_dims` to appropriately propagate parent information.

## Axes

Where Julia's currently documented [array interface]( https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array) requires defining `Base.size`, 
StaticArrayInterface instead requires defining [`StaticArrayInterface.static_axes`](@ref) and [`StaticArrayInterface.axes_types`](@ref).
`StaticArrayInterface.axes_types(::Type{T})` facilitates propagation of a number of traits known at compile time (e.g., `known_size`, `known_offsets`) 
and `StaticArrayInterface.static_axes(::AbstractArray)` replaces `Base.OneTo` with `StaticArrayInterface.OptionallyStaticUnitRange` in situations where static 
information would otherwise be lost.
`StaticArrayInterface.static_axes(::AbstractArray, dim)` utilizes `to_dims`, [as described elsewhere](#dimensions).

### Simple Wrappers

Let's say we have a new array type doesn't affect axes then this is as simple as:
```julia
Base.axes(x::SimpleWrapper) = StaticArrayInterface.static_axes(parent(x))
Base.axes(x::SimpleWrapper, dim) = StaticArrayInterface.static_axes(parent(x), dim)
StaticArrayInterface.axes_types(::Type{T}) where {T<:SimpleWrapper} = axes_types(parent_type(T))
```

To reiterate, `StaticArrayInterface.static_axes` improves on `Base.axes` for few Base array types but is otherwise identical.
Therefore, the first method simply ensures you don't have to define multiple parametric methods for your new type to preserve 
statically sized nested axes (e.g., `SimpleWrapper{T,N,<:Transpose{T,<:AbstractVector}}`).
This is otherwise identical to standard inheritance by composition.

### When to Discard Axis Information

Occasionally the parent array's axis information can't be preserved.
For example, we can't map axis information from the parent array of `Base.ReshapedArray`.
In this case we can simply build axes from the new size information.

```julia
StaticArrayInterface.axes_types(T::Type{<:ReshapedArray}) = NTuple{ndims(T),OneTo{Int}}
StaticArrayInterface.static_axes(A::ReshapedArray) = map(OneTo, size(A))
```

### New Axis Types

`OffsetArray` changes the first index for each axis.
It produces axes of type `IdOffsetRange`, which contains the value of the relative offset and the parent axis.

```julia
using StaticArrayInterface: axes_types, parent_type, to_dims
# Note that generating a `Tuple` type piecewise like may be type unstable and should be
# tested using `Test.@inferred`. It's often necessary to use generated function
# (`@generated`) or methods defined in Static.jl.
@generated function StaticArrayInterface.axes_types(::Type{A}) where {A<:OffsetArray}
    out = Expr(:curly, :Tuple)
    P = parent_type(A)
    for dim in 1:ndims(A)
        # offset relative to parent array
        O = relative_known_offsets(A, dim)
        if O === nothing  # offset is not known at compile time and is an `Int`
            push!(out.args, :(IdOffsetRange{Int, axes_types($P, $(static(dim)))}))
        else # offset is known, therefore it is a `StaticInt`
            push!(out.args, :(IdOffsetRange{StaticInt{$O}, axes_types($P, $(static(dim))}))
        end
    end
end
function Base.axes(A::OffsetArray)
    map(IdOffsetRange, StaticArrayInterface.static_axes(parent(A)), relative_offsets(A))
end
function Base.axes(A::OffsetArray, dim)
    d = to_dims(A, dim)
    IdOffsetRange(StaticArrayInterface.static_axes(parent(A), d), relative_offsets(A, d))
end
```

Defining these two methods ensures that other array types that wrap `OffsetArray` and appropriately define these methods 
propagate offsets independent of any dependency on `OffsetArray`.
It is entirely optional to define `StaticArrayInterface.static_size` for `OffsetArray` because the size can be 
derived from the axes.
However, in this particularly case we should also define
 `StaticArrayInterface.static_size(A::OffsetArray)  = StaticArrayInterface.static_size(parent(A))` because the relative offsets 
 attached to `OffsetArray` do not change the size but may hide static sizes if using a relative offset that is defined with an `Int`.

## Processing Indices (`static_to_indices`)

For most users, the only reason you should use `StaticArrayInterface.static_to_indices` over `Base.to_indices` is that it's 
faster and perhaps some of the more detailed benefits described in the [`to_indices`](@ref) doc string.
For those interested in how this is accomplished, the following steps (beginning with the `to_indices(A::AbstractArray, I::Tuple)`) 
are used to accomplish this:

1. The number of dimensions that each indexing argument in `I` corresponds to is determined using using the [`ndims_index`](@ref) and [`is_splat_index`](@ref) traits.
2. A non-allocating reference to each axis of `A` is created (`lazy_axes(A) -> axs`). These are aligned to each the index arguments using information from the first step. For example, if an index argument maps to a single dimension then it is paired with `axs[dim]`. In the case of multiple dimensions it is paired with `CartesianIndices(axs[dim_1], ... axs[dim_n])`. These pairs are further processed using `to_index(axis, I[n])`.
3. Tuples returned from `to_index` are flattened out so that there are no nested tuples returned from `to_indices`.

Entry points:

* `static_to_indices(::ArrayType, indices)` : dispatch on unique array type `ArrayType`
* `to_index(axis, ::IndexType)` : dispatch on a unique indexing type, `IndexType`. `StaticArrayInterface.ndims_index(::Type{IndexType})` should also be defined in this case.
* `to_index(S::IndexStyle, axis, index)` : The index style `S` that corresponds to `axis`. This is 

