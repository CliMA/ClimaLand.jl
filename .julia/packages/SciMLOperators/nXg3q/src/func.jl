#
"""
Matrix free operator given by a function

$(FIELDS)
"""
mutable struct FunctionOperator{iip, oop, mul5, T <: Number, F, Fa, Fi, Fai, Tr, U, P, Tt,
    C, iType, oType} <: AbstractSciMLOperator{T}
    """ Function with signature op(v, u, p, t) and (if isinplace) op(w, v, u, p, t) """
    op::F
    """ Adjoint operator"""
    op_adjoint::Fa
    """ Inverse operator"""
    op_inverse::Fi
    """ Adjoint inverse operator"""
    op_adjoint_inverse::Fai
    """ Traits """
    traits::Tr
    """ State """
    u::U
    """ Parameters """
    p::P
    """ Time """
    t::Tt
    """ Cache """
    cache::C
end

function FunctionOperator(op, op_adjoint, op_inverse, op_adjoint_inverse, traits, u, p, t,
        cache, ::Type{iType}, ::Type{oType}) where {iType, oType}
    iip = traits.isinplace
    oop = traits.outofplace
    mul5 = traits.has_mul5
    T = traits.T

    return FunctionOperator{iip, oop, mul5, T, typeof(op), typeof(op_adjoint),
        typeof(op_inverse), typeof(op_adjoint_inverse), typeof(traits), typeof(u), typeof(p),
        typeof(t), typeof(cache), iType, oType}(op, op_adjoint, op_inverse,
        op_adjoint_inverse, traits, u, p, t, cache)
end

function set_op(
        f::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C, iType,
            oType},
        op) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C, iType, oType}
    return FunctionOperator{
        iip, oop, mul5, T, typeof(op), Fa, Fi, Fai, Tr, U, P, Tt, C, iType,
        oType}(
        op, f.op_adjoint, f.op_inverse, f.op_adjoint_inverse, f.traits, f.u, f.p, f.t,
        f.cache)
end

function set_op_adjoint(
        f::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C,
            iType, oType},
        op_adjoint) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt,
        C, iType, oType}
    return FunctionOperator{
        iip, oop, mul5, T, F, typeof(op_adjoint), Fi, Fai, Tr, U, P, Tt,
        C, iType, oType}(f.op, op_adjoint, f.op_inverse, f.op_adjoint_inverse, f.traits,
        f.u, f.p, f.t, f.cache)
end

function set_op_inverse(
        f::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C,
            iType, oType},
        op_inverse) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt,
        C, iType, oType}
    return FunctionOperator{
        iip, oop, mul5, T, F, Fa, typeof(op_inverse), Fai, Tr, U, P, Tt,
        C, iType, oType}(f.op, f.op_adjoint, op_inverse, f.op_adjoint_inverse, f.traits,
        f.u, f.p, f.t, f.cache)
end

function set_op_adjoint_inverse(
        f::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr,
            U, P, Tt, C, iType, oType},
        op_adjoint_inverse) where {iip, oop, mul5, T, F, Fa,
        Fi, Fai, Tr, U, P, Tt, C, iType, oType}
    return FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, typeof(op_adjoint_inverse), Tr,
        U, P, Tt, C, iType, oType}(f.op, f.op_adjoint, f.op_inverse, op_adjoint_inverse,
        f.traits, f.u, f.p, f.t, f.cache)
end

function set_traits(
        f::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C,
            iType, oType},
        traits) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C,
        iType, oType}
    return FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, typeof(traits), U, P, Tt,
        C, iType, oType}(f.op, f.op_adjoint, f.op_inverse, f.op_adjoint_inverse, traits,
        f.u, f.p, f.t, f.cache)
end

function set_u(
        f::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C,
            iType, oType},
        u) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C, iType,
        oType}
    return FunctionOperator{
        iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, typeof(u), P, Tt, C, iType,
        oType}(
        f.op, f.op_adjoint, f.op_inverse, f.op_adjoint_inverse, f.traits, u, f.p, f.t,
        f.cache)
end

function set_p(
        f::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C,
            iType, oType},
        p) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C, iType,
        oType}
    return FunctionOperator{
        iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, typeof(p), Tt, C, iType,
        oType}(
        f.op, f.op_adjoint, f.op_inverse, f.op_adjoint_inverse, f.traits, f.u, p, f.t,
        f.cache)
end

function set_t(
        f::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C, iType,
            oType},
        t) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C, iType, oType}
    return FunctionOperator{
        iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, typeof(t), C, iType,
        oType}(
        f.op, f.op_adjoint, f.op_inverse, f.op_adjoint_inverse, f.traits, f.u, f.p, t,
        f.cache)
end

function set_cache(
        f::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C,
            iType, oType},
        cache) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C,
        iType, oType}
    return FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, typeof(cache),
        iType, oType}(f.op, f.op_adjoint, f.op_inverse, f.op_adjoint_inverse, f.traits,
        f.u, f.p, f.t, cache)
end

function input_eltype(::FunctionOperator{
        iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C,
        iType, oType
}) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C, iType, oType}
    return iType
end

function output_eltype(::FunctionOperator{
        iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C,
        iType, oType
}) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C, iType, oType}
    return oType
end

"""
$(SIGNATURES)

Wrap callable object `op` within an `AbstractSciMLOperator`.

## Mathematical Description

```julia
L = FunctionOperator(op, v, w; kwargs...)
```

where ``w = L(u,p,t)*v`` is done matrix-free given the function
definition ``w = op(v,u,p,t)``.

## Arguments

`op` is assumed to have signature

    op(v, u, p, t; <accepted_kwargs>) -> w

or

    op(w, v, u, p, t; <accepted_kwargs>) -> [modifies w]

and optionally

    op(w, v, u, p, t, α, β; <accepted_kwargs>) -> [modifies w]

where `u`, `v`, `w` are `AbstractArray`s, `p` is a parameter object, and
`t`, `α`, `β` are scalars. The first signature corresponds to applying
the operator with `Base.*`, and the latter two correspond to the
three-argument, and the five-argument `mul!` respectively.

`input` and `output` prototype `AbstractArray`s are required for
determining operator traits such as `eltype`, `size`, and for
preallocating cache. If `output` array is not provided, the output
is assumed to be of the same type and share as the input.

## Keyword Arguments

Keyword arguments are used to pass in the adjoint evaluation function,
`op_adjoint`, the inverse function, `op_inverse`, and the adjoint-inverse
function `adjoint_inverse`. All are assumed to have the same calling signature and
below traits.

## Traits

Keyword arguments are used to set operator traits, which are assumed to be
uniform across `op`, `op_adjoint`, `op_inverse`, `op_adjoint_inverse`.

  - `u` - Prototype of the state struct passed to the operator during evaluation, i.e. `L(u, p, t)`. `u` is set to `nothing` if no value is provided.
  - `p` - Prototype of parameter struct passed to the operator during evaluation, i.e. `L(u, p, t)`. `p` is set to `nothing` if no value is provided.
  - `t` - Protype of scalar time variable passed to the operator during evaluation. `t` is set to `zero(T)` if no value is provided.
  - `accepted_kwargs` - `Val` of a `Tuple` of `Symbol`s for zero-allocation kwarg filtering. Corresponds to the keyword arguments accepted by `op*`, and `update_coefficients[!]`. For example, if `op` accepts kwarg `scale`, as in `op(u, p, t; scale)`, then `accepted_kwargs = Val((:scale,))`. Plain tuples like `(:scale,)` are deprecated but still supported.
  - `T` - `eltype` of the operator. If no value is provided, the constructor inferrs the value from types of `input`, and `output`
  - `isinplace` - `true` if the operator can be used is a mutating way with in-place allocations. This trait is inferred if no value is provided.
  - `outofplace` - `true` if the operator can be used is a non-mutating way with in-place allocations. This trait is inferred if no value is provided.
  - `has_mul5` - `true` if the operator provides a five-argument `mul!` via the signature `op(v, u, p, t, α, β; <accepted_kwargs>)`. This trait is inferred if no value is provided.
  - `isconstant` - `true` if the operator is constant, and doesn't need to be updated via `update_coefficients[!]` during operator evaluation. Defaults to false.
  - `islinear` - `true` if the operator is linear. Defaults to `false`.
  - `isconvertible` - `true` a cheap `convert(AbstractMatrix, L.op)` method is available. Defaults to `false`.
  - `batch` - Boolean indicating if the input/output arrays comprise of batched column-vectors stacked in a matrix. If `true`, the input/output arrays must be `AbstractVecOrMat`s, and the length of the  second dimension (the batch dimension) must be the same. The batch dimension is not involved in size computation. For example, with `batch = true`, and `size(output), size(input) = (M, K), (N, K)`, the `FunctionOperator` size is set to `(M, N)`. If `batch = false`, which is the default, the `input`/`output` arrays may of any size so long as `ndims(input) == ndims(output)`, and the `size` of `FunctionOperator` is set to `(length(input), length(output))`.
  - `ifcache` - Allocate cache arrays in constructor. Defaults to `true`. Cache can be generated afterwards by calling `cache_operator(L, input, output)`
  - `cache` - Pregenerated cache arrays for in-place evaluations. Expected to be of type and shape `(similar(input), similar(output),)`. The constructor generates cache if no values are provided. Cache generation by the constructor can be disabled by setting the kwarg `ifcache = false`.
  - `opnorm` - The norm of `op`. Can be a `Number`, or function `opnorm(p::Integer)`. Defaults to `nothing`.
  - `issymmetric` - `true` if the operator is linear and symmetric. Defaults to `false`.
  - `ishermitian` - `true` if the operator is linear and hermitian. Defaults to `false`.
  - `isposdef` - `true` if the operator is linear and positive-definite. Defaults to `false`.
  - `kwargs` - Keyword arguments for cache initialization. If `accepted_kwargs` is provided, the corresponding keyword arguments must be passed.
"""
function FunctionOperator(op,
        input::AbstractArray,
        output::AbstractArray = input; op_adjoint = nothing,
        op_inverse = nothing,
        op_adjoint_inverse = nothing, u = nothing, p = nothing,
        t::Union{Number, Nothing} = nothing,
        accepted_kwargs::Union{Nothing, Val, NTuple{N, Symbol}} = nothing,

        # traits
        T::Union{Type{<:Number}, Nothing} = nothing,
        isinplace::Union{Nothing, Bool, Val} = nothing,
        outofplace::Union{Nothing, Bool, Val} = nothing,
        has_mul5::Union{Nothing, Bool, Val} = nothing,
        isconstant::Bool = false,
        islinear::Bool = false,
        isconvertible::Bool = false, batch::Bool = false,
        ifcache::Union{Bool, Val} = Val(true),
        cache::Union{Nothing, NTuple{2}} = nothing,

        # LinearAlgebra traits
        opnorm = nothing,
        issymmetric::Union{Bool, Val} = Val(false),
        ishermitian::Union{Bool, Val} = Val(false),
        isposdef::Bool = false,
        kwargs...) where {N}

    # establish types

    # store eltype of input/output for caching with ComposedOperator.
    _T = T === nothing ? promote_type(eltype(input), eltype(output)) : T
    _t = t === nothing ? zero(real(_T)) : t

    isinplace isa Val && (@assert _unwrap_val(isinplace) isa Bool)
    outofplace isa Val && (@assert _unwrap_val(outofplace) isa Bool)
    has_mul5 isa Val && (@assert _unwrap_val(has_mul5) isa Bool)
    issymmetric isa Val && (@assert _unwrap_val(issymmetric) isa Bool)
    ishermitian isa Val && (@assert _unwrap_val(ishermitian) isa Bool)

    @assert _T<:Number """The `eltype` of `FunctionOperator`, as well as
  the `input`/`output` arrays must be `<:Number`."""

    # establish sizes

    @assert ndims(output)==ndims(input) """`input`/`output` arrays,
  ($(typeof(input)), $(typeof(output))) provided to `FunctionOperator`
  do not have the same number of dimensions. Further, if `batch = true`,
  then both arrays must be `AbstractVector`s, or both must be
  `AbstractMatrix` types."""

    if batch
        if !isa(input, AbstractVecOrMat)
            msg = """`FunctionOperator` constructed with `batch = true` only
                accepts `AbstractVecOrMat` types with
                `size(L, 2) == size(u, 1)`."""
            throw(ArgumentError(msg))
        end

        if input isa AbstractMatrix
            # assume batches are 2nd dimension of `AbstractVecOrMat`
            if size(input, 2) != size(output, 2)
                msg = """ Batch size (length of second dimension) in `input`/
                    `output` arrays to `FunctionOperator` is not equal. Input
                    array, $(typeof(input)), has size $(size(input)), whereas
                    output array, $(typeof(output)), has size
                    $(size(output))."""
                throw(ArgumentError(msg))
            end
        end
    end

    sizes = size(input), size(output)

    _size = if batch
        # assume batches are 2nd dimension of `AbstractVecOrMat`
        (size(output, 1), size(input, 1))
    else
        (length(output), length(input))
    end

    # evaluation signatures

    _isinplace = if isinplace === nothing
        Val(hasmethod(op, typeof((output, input, u, p, _t))))
    elseif isinplace isa Bool
        Val(isinplace)
    else
        isinplace
    end

    _outofplace = if outofplace === nothing
        Val(hasmethod(op, typeof((input, u, p, _t))))
    elseif outofplace isa Bool
        Val(outofplace)
    else
        outofplace
    end

    if !_unwrap_val(_isinplace) & !_unwrap_val(_outofplace)
        @error """Please provide a function with signatures `op(v, u, p, t)` for
        applying the operator out-of-place, and/or the signature is
        `op(w, v, u, p, t)` for in-place application."""
    end

    _has_mul5 = if has_mul5 === nothing
        __and_val(__has_mul5(op, output, input, u, p, _t),
            __has_mul5(op_adjoint, input, output, u, p, _t),
            __has_mul5(op_inverse, output, input, u, p, _t),
            __has_mul5(op_adjoint_inverse, input, output, u, p, _t))
    elseif has_mul5 isa Bool
        Val(has_mul5)
    else
        has_mul5
    end

    # traits

    isreal = _T <: Real
    selfadjoint = _unwrap_val(ishermitian) | (isreal & _unwrap_val(issymmetric))
    adjointable = !(op_adjoint isa Nothing) | _unwrap_val(selfadjoint)
    invertible = !(op_inverse isa Nothing)

    if selfadjoint & (op_adjoint isa Nothing)
        _op_adjoint = op
    else
        _op_adjoint = op_adjoint
    end

    if selfadjoint & invertible & (op_adjoint_inverse isa Nothing)
        _op_adjoint_inverse = op_inverse
    else
        _op_adjoint_inverse = op_adjoint_inverse
    end

    if accepted_kwargs === nothing
        accepted_kwargs = Val(())
        kwargs = NamedTuple()
    else
        length(kwargs) != 0 ||
            throw(ArgumentError("No keyword arguments provided. When `accepted_kwargs` is provided, the corresponding keyword arguments must be passed for cache initialization."))
        kwargs = get_filtered_kwargs(kwargs, accepted_kwargs)
    end

    traits = (; islinear, isconvertible, isconstant, opnorm,
        issymmetric = _unwrap_val(issymmetric), ishermitian = _unwrap_val(ishermitian),
        isposdef, isinplace = _unwrap_val(_isinplace),
        outofplace = _unwrap_val(_outofplace), has_mul5 = _unwrap_val(_has_mul5),
        ifcache = _unwrap_val(ifcache), T = _T, batch, size = _size, sizes,
        accepted_kwargs, kwargs = kwargs)

    L = FunctionOperator{_unwrap_val(_isinplace), _unwrap_val(_outofplace),
        _unwrap_val(_has_mul5), _T, typeof(op), typeof(_op_adjoint), typeof(op_inverse),
        typeof(_op_adjoint_inverse), typeof(traits), typeof(u), typeof(p), typeof(_t), typeof(cache),
        eltype(input), eltype(output)}(op,
        _op_adjoint, op_inverse, _op_adjoint_inverse, traits, u, p, _t, cache)

    # create cache

    if _unwrap_val(ifcache) & (L.cache === nothing)
        L_cached = cache_operator(L, input)
    else
        L_cached = L
    end

    return L_cached
end

@inline __has_mul5(::Nothing, w, v, u, p, t) = Val(true)
@inline function __has_mul5(f::F, w, v, u, p, t) where {F}
    return Val(hasmethod(f, typeof((w, v, u, p, t, t, t))))
end
@inline __and_val(vs...) = mapreduce(_unwrap_val, *, vs)

function update_coefficients(L::FunctionOperator, u, p, t; kwargs...)

    # update u, p, t
    L = set_u(L, u)
    L = set_p(L, p)
    L = set_t(L, t)

    # filter and update kwargs
    filtered_kwargs = get_filtered_kwargs(kwargs, L.traits.accepted_kwargs)

    L = set_traits(L, merge(L.traits, (; kwargs = filtered_kwargs)))

    isconstant(L) && return L

    L = set_op(L, update_coefficients(L.op, u, p, t; filtered_kwargs...))
    L = set_op_adjoint(L, update_coefficients(L.op_adjoint, u, p, t; filtered_kwargs...))
    L = set_op_inverse(L, update_coefficients(L.op_inverse, u, p, t; filtered_kwargs...))
    L = set_op_adjoint_inverse(L,
        update_coefficients(L.op_adjoint_inverse, u, p, t; filtered_kwargs...))
end

function update_coefficients!(L::FunctionOperator, u, p, t; kwargs...)

    # update u, p, t
    L.u = u
    L.p = p
    L.t = t

    # filter and update kwargs
    filtered_kwargs = get_filtered_kwargs(kwargs, L.traits.accepted_kwargs)
    L.traits = merge(L.traits, (; kwargs = filtered_kwargs))

    isconstant(L) && return

    for op in getops(L)
        update_coefficients!(op, u, p, t; filtered_kwargs...)
    end

    nothing
end

function iscached(L::FunctionOperator)
    # L.traits.ifcache ? !isnothing(L.cache) : !L.traits.ifcache
    L.cache !== nothing
end

# Copy method to avoid aliasing
function Base.copy(L::FunctionOperator)
    FunctionOperator(
        L.op,
        L.op_adjoint,
        L.op_inverse,
        L.op_adjoint_inverse,
        L.traits,
        isdefined(L, :u) && L.u !== nothing ? copy(L.u) : nothing,
        isdefined(L, :p) && L.p !== nothing ? deepcopy(L.p) : nothing,
        L.t,
        L.cache === nothing ? nothing : deepcopy(L.cache),
        typeof(L).parameters[end - 1],  # iType
        typeof(L).parameters[end]      # oType
    )
end

# fix method amg bw AbstractArray, AbstractVecOrMat
cache_operator(L::FunctionOperator, u::AbstractArray) = _cache_operator(L, u)
cache_operator(L::FunctionOperator, u::AbstractVecOrMat) = _cache_operator(L, u)

function _cache_operator(L::FunctionOperator, u::AbstractArray)
    U = if L.traits.batch
        if !isa(u, AbstractVecOrMat)
            msg = """$L constructed with `batch = true` only accepts
                `AbstractVecOrMat` types with `size(L, 2) == size(u, 1)`."""
            throw(ArgumentError(msg))
        end

        if size(L, 2) != size(u, 1)
            msg = """Second dimension of $L of size $(size(L))
                is not consistent with first dimension of input array `u`
                of size $(size(u))."""
            throw(DimensionMismatch(msg))
        end

        M = size(L, 1)
        K = size(u, 2)
        size_out = u isa AbstractVector ? (M,) : (M, K)

        new_traits = merge(L.traits, (; sizes = (size(u), size_out)))
        L = set_traits(L, new_traits)

        u
    else
        if size(L, 2) != length(u)
            msg = """Length of input array, $(typeof(u)), of size $(size(u))
                not consistent with second dimension of $L of size
                $(size(L))."""
            throw(DimensionMismatch(msg))
        end

        reshape(u, L.traits.sizes[1])
    end

    L = cache_self(L, U)
    L = cache_internals(L, U)
    L
end

# fix method amg bw AbstractArray, AbstractVecOrMat
cache_self(L::FunctionOperator, v::AbstractArray) = _cache_self(L, v)
cache_self(L::FunctionOperator, v::AbstractVecOrMat) = _cache_self(L, v)

function _cache_self(
        L::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C,
            iType, oType},
        v::AbstractArray) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P,
        Tt, C, iType, oType}
    _v = similar(v, iType, L.traits.sizes[1])
    _w = similar(v, oType, L.traits.sizes[2])

    return set_cache(L, (_v, _w))
end

# fix method amg bw AbstractArray, AbstractVecOrMat
cache_internals(L::FunctionOperator, u::AbstractArray) = _cache_internals(L, u)
cache_internals(L::FunctionOperator, u::AbstractVecOrMat) = _cache_internals(L, u)

function _cache_internals(
        L::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt,
            C, iType, oType},
        u::AbstractArray) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr,
        U, P, Tt, C, iType, oType}
    newop = cache_operator(L.op, u)
    newop_adjoint = cache_operator(L.op_adjoint, u)
    newop_inverse = cache_operator(L.op_inverse, u)
    newop_adjoint_inverse = cache_operator(L.op_adjoint_inverse, u)

    return FunctionOperator{iip, oop, mul5, T, typeof(newop), typeof(newop_adjoint),
        typeof(newop_inverse), typeof(newop_adjoint_inverse), Tr, U, P, Tt, C, iType, oType}(
        newop, newop_adjoint, newop_inverse, newop_adjoint_inverse, L.traits, L.u, L.p, L.t,
        L.cache)
end

function Base.show(io::IO, L::FunctionOperator)
    M, N = size(L)
    print(io, "FunctionOperator($M × $N)")
end
Base.size(L::FunctionOperator) = L.traits.size

function Base.adjoint(L::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt,
        C, iType,
        oType
}) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C, iType, oType}
    (ishermitian(L) | (isreal(L) & issymmetric(L))) && return L

    has_adjoint(L) || return AdjointOperator(L)

    op = L.op_adjoint
    op_adjoint = L.op

    op_inverse = L.op_adjoint_inverse
    op_adjoint_inverse = L.op_inverse

    traits = merge(L.traits, (; size = reverse(size(L)), sizes = reverse(L.traits.sizes)))

    cache = iscached(L) ? reverse(L.cache) : nothing

    return FunctionOperator{iip, oop, mul5, T, typeof(op), typeof(op_adjoint),
        typeof(op_inverse), typeof(op_adjoint_inverse), typeof(traits), U, P, Tt,
        typeof(cache), oType, iType}(
        op, op_adjoint, op_inverse, op_adjoint_inverse, traits,
        L.u, L.p, L.t, cache)
end

function Base.inv(L::FunctionOperator{iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt,
        C, iType,
        oType
}) where {iip, oop, mul5, T, F, Fa, Fi, Fai, Tr, U, P, Tt, C, iType, oType}
    has_ldiv(L) || return InvertedOperator(L)

    op = L.op_inverse
    op_inverse = L.op

    op_adjoint = L.op_adjoint_inverse
    op_adjoint_inverse = L.op_adjoint

    opnorm = if L.traits.opnorm isa Number
        1 / L.traits.opnorm
    elseif L.traits.opnorm isa Nothing
        nothing
    else
        (p::Real) -> 1 / L.traits.opnorm(p)
    end
    traits = merge(L.traits,
        (; size = reverse(size(L)), sizes = reverse(L.traits.sizes), opnorm))

    cache = iscached(L) ? reverse(L.cache) : nothing

    return FunctionOperator{iip, oop, mul5, T, typeof(op), typeof(op_adjoint),
        typeof(op_inverse), typeof(op_adjoint_inverse), typeof(traits), U, P, Tt,
        typeof(cache), oType, iType}(
        op, op_adjoint, op_inverse, op_adjoint_inverse, traits,
        L.u, L.p, L.t, cache)
end

Base.convert(::Type{AbstractMatrix}, L::FunctionOperator) = convert(AbstractMatrix, L.op)

function Base.resize!(L::FunctionOperator, n::Integer)

    # input/output to `L` must be `AbstractVector`s
    if length(L.traits.sizes[1]) != 1
        msg = """`Base.resize!` is only supported by $L whose input/output
            arrays are `AbstractVector`s."""
        throw(MethodError(msg))
    end

    for op in getops(L)
        if hasmethod(resize!, typeof((op, n)))
            resize!(op, n)
        end
    end

    for v in L.cache
        resize!(v, n)
    end

    L.traits = (; L.traits..., size = (n, n), sizes = ((n,), (n,)))

    L
end

function LinearAlgebra.opnorm(L::FunctionOperator, p)
    L.traits.opnorm === nothing && error("""
      M.opnorm is nothing, please define opnorm as a function that takes one
      argument. E.g., `(p::Real) -> p == Inf ? 100 : error("only Inf norm is
      defined")`""")
    opn = L.traits.opnorm
    return opn isa Number ? opn : L.traits.opnorm(p)
end
LinearAlgebra.issymmetric(L::FunctionOperator) = L.traits.issymmetric
LinearAlgebra.ishermitian(L::FunctionOperator) = L.traits.ishermitian
LinearAlgebra.isposdef(L::FunctionOperator) = L.traits.isposdef

function getops(L::FunctionOperator)
    (;
        op = L.op,
        op_adjoint = L.op_adjoint,
        op_inverse = L.op_inverse,
        op_adjoint_inverse = L.op_adjoint_inverse)
end

islinear(L::FunctionOperator) = L.traits.islinear
isconvertible(L::FunctionOperator) = L.traits.isconvertible
isconstant(L::FunctionOperator) = L.traits.isconstant
has_adjoint(L::FunctionOperator) = !(L.op_adjoint isa Nothing)
has_mul(::FunctionOperator{iip}) where {iip} = true
has_mul!(::FunctionOperator{iip}) where {iip} = true
has_ldiv(L::FunctionOperator{iip}) where {iip} = !(L.op_inverse isa Nothing)
has_ldiv!(L::FunctionOperator{iip}) where {iip} = iip & !(L.op_inverse isa Nothing)

function _sizecheck(L::FunctionOperator, v, w)
    sizes = L.traits.sizes
    if L.traits.batch
        if !isnothing(v)
            if !isa(v, AbstractVecOrMat)
                msg = """$L constructed with `batch = true` only
                    accept input arrays that are `AbstractVecOrMat`s with
                    `size(L, 2) == size(v, 1)`. Received $(typeof(v))."""
                throw(ArgumentError(msg))
            end

            if size(L, 2) != size(v, 1)
                msg = """$L accepts input `AbstractVecOrMat`s of size
                    ($(size(L, 2)), K). Received array of size $(size(v))."""
                throw(DimensionMismatch(msg))
            end
        end # v

        if !isnothing(w)
            if !isa(w, AbstractVecOrMat)
                msg = """$L constructed with `batch = true` only
                    returns output arrays that are `AbstractVecOrMat`s with
                    `size(L, 1) == size(w, 1)`. Received $(typeof(w))."""
                throw(ArgumentError(msg))
            end

            if size(L, 1) != size(w, 1)
                msg = """$L accepts output `AbstractVecOrMat`s of size
                    ($(size(L, 1)), K). Received array of size $(size(w))."""
                throw(DimensionMismatch(msg))
            end
        end # w

        if !isnothing(v) & !isnothing(w)
            if size(v, 2) != size(w, 2)
                msg = """input array $v, and output array, $w, must have the
                    same batch size (i.e. length of second dimension). Got
                    $(size(v)), $(size(w)). If you encounter this error during
                    an in-place evaluation (`LinearAlgebra.mul!`, `ldiv!`),
                    ensure that the operator $L has been cached with an input
                    array of the correct size. Do so by calling
                    `L = cache_operator(L, v)`."""
                throw(DimensionMismatch(msg))
            end
        end # v, w

    else # !batch
        if !isnothing(v)
            if size(v) ∉ (sizes[1], tuple(size(L, 2)))
                msg = """$L received input array of size $(size(v)), but only
                    accepts input arrays of size $(sizes[1]), or vectors like
                    `vec(v)` of size $(tuple(prod(sizes[1])))."""
                throw(DimensionMismatch(msg))
            end
        end # v

        if !isnothing(w)
            if size(w) ∉ (sizes[2], tuple(size(L, 1)))
                msg = """$L received output array of size $(size(w)), but only
                    accepts output arrays of size $(sizes[2]), or vectors like
                    `vec(v)` of size $(tuple(prod(sizes[2])))"""
                throw(DimensionMismatch(msg))
            end
        end # w
    end # batch

    return
end

function _unvec(L::FunctionOperator, v, w)
    if L.traits.batch
        return v, w, false
    else
        sizes = L.traits.sizes

        # no need to vec since expected input/output are AbstractVectors
        if length(sizes[1]) == 1
            return v, w, false
        end

        vec_v = isnothing(v) ? false : size(v) != sizes[1]
        vec_w = isnothing(w) ? false : size(w) != sizes[2]

        if !isnothing(v) & !isnothing(w)
            if (vec_v & !vec_w) | (!vec_v & vec_w)
                msg = """Input / output to $L can either be of sizes
                    $(sizes[1]) / $(sizes[2]), or
                    $(tuple(prod(sizes[1]))) / $(tuple(prod(sizes[2]))). Got
                    $(size(v)), $(size(w))."""
                throw(DimensionMismatch(msg))
            end
        end

        V = vec_v ? reshape(v, sizes[1]) : v
        W = vec_w ? reshape(w, sizes[2]) : w
        vec_output = vec_v | vec_w

        return V, W, vec_output
    end
end

# operator application
function Base.:*(L::FunctionOperator{iip, true}, v::AbstractArray) where {iip}
    _sizecheck(L, v, nothing)
    V, _, vec_output = _unvec(L, v, nothing)

    W = L.op(V, L.u, L.p, L.t; L.traits.kwargs...)

    vec_output ? vec(W) : W
end

function Base.:*(L::FunctionOperator{iip, false}, v::AbstractArray) where {iip}
    _sizecheck(L, v, nothing)
    V, _, vec_output = _unvec(L, v, nothing)

    W, _ = L.cache
    W = copy(W)
    L.op(W, V, L.u, L.p, L.t; L.traits.kwargs...)

    vec_output ? vec(W) : W
end

function Base.:\(L::FunctionOperator{iip, true}, v::AbstractArray) where {iip}
    _sizecheck(L, nothing, v)
    _, V, vec_output = _unvec(L, nothing, v)

    W = L.op_inverse(V, L.u, L.p, L.t; L.traits.kwargs...)

    vec_output ? vec(W) : W
end

function Base.:\(L::FunctionOperator{iip, false}, v::AbstractArray) where {iip}
    _sizecheck(L, nothing, v)
    _, V, vec_output = _unvec(L, nothing, v)

    W, _ = L.cache
    W = copy(W)
    L.op_inverse(W, V, L.u, L.p, L.t; L.traits.kwargs...)

    vec_output ? vec(W) : W
end

function LinearAlgebra.mul!(w::AbstractArray, L::FunctionOperator{true}, v::AbstractArray)
    _sizecheck(L, v, w)
    V, W, vec_output = _unvec(L, v, w)
    L.op(W, V, L.u, L.p, L.t; L.traits.kwargs...)
    vec_output ? vec(W) : W
end

function LinearAlgebra.mul!(w::AbstractArray, L::FunctionOperator{false}, v::AbstractArray,
        args...)
    _sizecheck(L, v, w)
    V, W, vec_output = _unvec(L, v, w)
    W .= L.op(V, L.u, L.p, L.t; L.traits.kwargs...)
    vec_output ? vec(W) : W
end

function LinearAlgebra.mul!(w::AbstractArray, L::FunctionOperator{true, oop, false},
        v::AbstractArray, α, β) where {oop}
    _, Co = L.cache

    _sizecheck(L, v, w)
    V, W, _ = _unvec(L, v, w)

    copy!(Co, W)
    L.op(W, V, L.u, L.p, L.t; L.traits.kwargs...) # mul!(V, L, U)
    axpby!(β, Co, α, W)

    w
end

function LinearAlgebra.mul!(w::AbstractArray, L::FunctionOperator{true, oop, true},
        v::AbstractArray, α, β) where {oop}
    _sizecheck(L, v, w)
    V, W, _ = _unvec(L, v, w)

    L.op(W, V, L.u, L.p, L.t, α, β; L.traits.kwargs...)

    w
end

function LinearAlgebra.ldiv!(w::AbstractArray, L::FunctionOperator{true}, v::AbstractArray)
    _sizecheck(L, v, w)
    W, V, _ = _unvec(L, w, v)

    L.op_inverse(W, V, L.u, L.p, L.t; L.traits.kwargs...)

    w
end

function LinearAlgebra.ldiv!(L::FunctionOperator{true}, v::AbstractArray)
    W, _ = L.cache

    _sizecheck(L, nothing, v)
    V, _, vec_output = _unvec(L, v, nothing)

    copy!(W, V)
    L.op_inverse(W, V, L.u, L.p, L.t; L.traits.kwargs...) # ldiv!(U, L, V)

    V .= W
    vec_output ? vec(V) : V
end

function LinearAlgebra.ldiv!(w::AbstractArray, L::FunctionOperator{false}, v::AbstractArray)
    _sizecheck(L, v, w)
    W, V, _ = _unvec(L, w, v)

    W .= L.op_inverse(V, L.u, L.p, L.t; L.traits.kwargs...)

    w
end

function LinearAlgebra.ldiv!(L::FunctionOperator{false}, v::AbstractArray)
    _sizecheck(L, nothing, v)
    V, _, vec_output = _unvec(L, v, nothing)

    V .= L.op_inverse(V, L.u, L.p, L.t; L.traits.kwargs...) # ldiv!(W, L, V)

    vec_output ? vec(V) : V
end

# Out-of-place: v is action vector, u is update vector
function (L::FunctionOperator)(v::AbstractArray, u, p, t; kwargs...)
    L = update_coefficients(L, u, p, t; kwargs...)
    _sizecheck(L, v, nothing)
    V, _, vec_output = _unvec(L, v, nothing)

    # Apply the operator to action vector v after updating with u
    if L.traits.outofplace
        result = L.op(V, L.u, L.p, L.t; L.traits.kwargs...)
        return vec_output ? vec(result) : result
    else
        # For operators without out-of-place methods, use their in-place methods with a temporary
        Co = similar(V)
        L.op(Co, V, L.u, L.p, L.t; L.traits.kwargs...)
        return vec_output ? vec(Co) : Co
    end

    v
end

# In-place: w is destination, v is action vector, u is update vector
function (L::FunctionOperator)(w::AbstractArray, v::AbstractArray, u, p, t; kwargs...)
    update_coefficients!(L, u, p, t; kwargs...)

    # Check dimensions
    _sizecheck(L, v, w)
    V, W, _ = _unvec(L, v, w)

    # Apply the operator in-place to action vector v after updating with u
    if L.traits.isinplace
        L.op(W, V, L.u, L.p, L.t; L.traits.kwargs...)
    else
        # For operators without in-place methods, use their out-of-place methods
        result = L.op(V, L.u, L.p, L.t; L.traits.kwargs...)
        copyto!(W, result)
    end

    return w
end

# In-place with scaling: w = α*(L*v) + β*w
function (L::FunctionOperator)(w::AbstractArray, v::AbstractArray, u, p, t, α, β; kwargs...)
    update_coefficients!(L, u, p, t; kwargs...)

    # Check dimensions
    _sizecheck(L, v, w)
    V, W, _ = _unvec(L, v, w)

    # Apply the operator in-place to action vector v with scaling
    if L.traits.isinplace && L.traits.has_mul5
        # Direct 5-arg mul! if supported
        L.op(W, V, L.u, L.p, L.t, α, β; L.traits.kwargs...)
    elseif L.traits.isinplace
        # Use temporary for regular in-place
        temp = copy(W)
        L.op(W, V, L.u, L.p, L.t; L.traits.kwargs...)
        axpby!(β, temp, α, W)
    else
        # Out-of-place with scaling
        result = L.op(V, L.u, L.p, L.t; L.traits.kwargs...)
        axpby!(β, W, α, result)
    end

    return w
end
#
