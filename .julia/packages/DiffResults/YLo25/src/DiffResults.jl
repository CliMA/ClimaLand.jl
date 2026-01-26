module DiffResults

using StaticArraysCore: StaticArray, similar_type, Size

#########
# Types #
#########

abstract type DiffResult{O,V,D<:Tuple} end

struct ImmutableDiffResult{O,V,D<:Tuple} <: DiffResult{O,V,D}
    value::V
    derivs::D # ith element = ith-order derivative
    function ImmutableDiffResult(value::V, derivs::NTuple{O,Any}) where {O,V}
        return new{O,V,typeof(derivs)}(value, derivs)
    end
end

mutable struct MutableDiffResult{O,V,D<:Tuple} <: DiffResult{O,V,D}
    value::V
    derivs::D # ith element = ith-order derivative
    function MutableDiffResult(value::V, derivs::NTuple{O,Any}) where {O,V}
        return new{O,V,typeof(derivs)}(value, derivs)
    end
end

################
# Constructors #
################

"""
    DiffResult(value::Union{Number,AbstractArray}, derivs::Tuple{Vararg{Number}})
    DiffResult(value::Union{Number,AbstractArray}, derivs::Tuple{Vararg{AbstractArray}})

Return `r::DiffResult`, with output value storage provided by `value` and output derivative
storage provided by `derivs`.

In reality, `DiffResult` is an abstract supertype of two concrete types, `MutableDiffResult`
and `ImmutableDiffResult`. If all `value`/`derivs` are all `Number`s or `StaticArray`s,
then `r` will be immutable (i.e. `r::ImmutableDiffResult`). Otherwise, `r` will be mutable
(i.e. `r::MutableDiffResult`).

Note that `derivs` can be provide in splatted form, i.e. `DiffResult(value, derivs...)`.
"""
DiffResult

DiffResult(value::Number, derivs::Tuple{Vararg{Number}}) = ImmutableDiffResult(value, derivs)
DiffResult(value::Number, derivs::Tuple{Vararg{StaticArray}}) = ImmutableDiffResult(value, derivs)
DiffResult(value::StaticArray, derivs::Tuple{Vararg{StaticArray}}) = ImmutableDiffResult(value, derivs)
DiffResult(value::Number, derivs::Tuple{Vararg{AbstractArray}}) = MutableDiffResult(value, derivs)
DiffResult(value::AbstractArray, derivs::Tuple{Vararg{AbstractArray}}) = MutableDiffResult(value, derivs)
DiffResult(value::Union{Number,AbstractArray}, derivs::Union{Number,AbstractArray}...) = DiffResult(value, derivs)

"""
    GradientResult(x::AbstractArray)

Construct a `DiffResult` that can be used for gradient calculations where `x` is the input
to the target function.

Note that `GradientResult` allocates its own storage; `x` is only used for type and
shape information. If you want to allocate storage yourself, use the `DiffResult`
constructor instead.
"""
GradientResult(x::AbstractArray) = DiffResult(first(x), similar(x))
GradientResult(x::StaticArray) = DiffResult(first(x), x)

"""
    JacobianResult(x::AbstractArray)

Construct a `DiffResult` that can be used for Jacobian calculations where `x` is the input
to the target function. This method assumes that the target function's output dimension
equals its input dimension.

Note that `JacobianResult` allocates its own storage; `x` is only used for type and
shape information. If you want to allocate storage yourself, use the `DiffResult`
constructor instead.
"""
JacobianResult(x::AbstractArray) = DiffResult(similar(x), similar(x, length(x), length(x)))
JacobianResult(x::StaticArray) = DiffResult(x, zeros(similar_type(typeof(x), Size(length(x),length(x)))))

"""
    JacobianResult(y::AbstractArray, x::AbstractArray)

Construct a `DiffResult` that can be used for Jacobian calculations where `x` is the
input to the target function, and `y` is the output (e.g. when taking the Jacobian
of `f!(y, x)`).

Like the single argument version, `y` and `x` are only used for type and
shape information and are not stored in the returned `DiffResult`.
"""
JacobianResult(y::AbstractArray, x::AbstractArray) = DiffResult(similar(y), similar(y, length(y), length(x)))
JacobianResult(y::StaticArray, x::StaticArray) = DiffResult(y, zeros(similar_type(typeof(x), Size(length(y),length(x)))))

"""
    HessianResult(x::AbstractArray)

Construct a `DiffResult` that can be used for Hessian calculations where `x` is the
input to the target function.

Note that `HessianResult` allocates its own storage; `x` is only used for type and
shape information. If you want to allocate storage yourself, use the `DiffResult`
constructor instead.
"""
HessianResult(x::AbstractArray) = DiffResult(first(x), zeros(eltype(x), size(x)), similar(x, length(x), length(x)))
HessianResult(x::StaticArray) = DiffResult(first(x), x, zeros(similar_type(typeof(x), Size(length(x),length(x)))))

#############
# Interface #
#############

@generated function tuple_eltype(x::Tuple, ::Type{Val{i}}) where {i}
    return quote
        $(Expr(:meta, :inline))
        return $(x.parameters[i])
    end
end

@generated function tuple_setindex(x::NTuple{N,Any}, y, ::Type{Val{i}}) where {N,i}
    T = x.parameters[i]
    new_tuple = Expr(:tuple, [ifelse(i == n, :(convert($T, y)), :(x[$n])) for n in 1:N]...)
    return quote
        $(Expr(:meta, :inline))
        return $new_tuple
    end
end

Base.eltype(r::DiffResult) = eltype(typeof(r))

Base.eltype(::Type{D}) where {O,V,D<:DiffResult{O,V}} = eltype(V)

Base.:(==)(a::DiffResult, b::DiffResult) = a.value == b.value && a.derivs == b.derivs

Base.copy(r::DiffResult) = DiffResult(copy(r.value), map(copy, r.derivs))

# value/value! #
#--------------#

"""
    value(r::DiffResult)

Return the primal value stored in `r`.

Note that this method returns a reference, not a copy.
"""
value(r::DiffResult) = r.value

"""
    value!(r::DiffResult, x)

Return `s::DiffResult` with the same data as `r`, except for `value(s) == x`.

This function may or may not mutate `r`. If `r::ImmutableDiffResult`, a totally new
instance will be created and returned, whereas if `r::MutableDiffResult`, then `r` will be
mutated in-place and returned. Thus, this function should be called as `r = value!(r, x)`.
"""
value!(r::MutableDiffResult, x::Number) = (r.value = x; return r)
value!(r::MutableDiffResult, x::AbstractArray) = (copyto!(value(r), x); return r)
value!(r::ImmutableDiffResult{O,V}, x::Union{Number,AbstractArray}) where {O,V} = ImmutableDiffResult(convert(V, x), r.derivs)

"""
    value!(f, r::DiffResult, x)

Equivalent to `value!(r::DiffResult, map(f, x))`, but without the implied temporary
allocation (when possible).
"""
value!(f, r::MutableDiffResult, x::Number) = (r.value = f(x); return r)
value!(f, r::MutableDiffResult, x::AbstractArray) = (map!(f, value(r), x); return r)
value!(f, r::ImmutableDiffResult{O,V}, x::Number) where {O,V} = value!(r, convert(V, f(x)))
value!(f, r::ImmutableDiffResult{O,V}, x::AbstractArray) where {O,V} = value!(r, convert(V, map(f, x)))

# derivative/derivative! #
#------------------------#

"""
    derivative(r::DiffResult, ::Type{Val{i}} = Val{1})

Return the `ith` derivative stored in `r`, defaulting to the first derivative.

Note that this method returns a reference, not a copy.
"""
derivative(r::DiffResult, ::Type{Val{i}} = Val{1}) where {i} = r.derivs[i]

"""
    derivative!(r::DiffResult, x, ::Type{Val{i}} = Val{1})

Return `s::DiffResult` with the same data as `r`, except `derivative(s, Val{i}) == x`.

This function may or may not mutate `r`. If `r::ImmutableDiffResult`, a totally new
instance will be created and returned, whereas if `r::MutableDiffResult`, then `r` will be
mutated in-place and returned. Thus, this function should be called as
`r = derivative!(r, x, Val{i})`.
"""
function derivative!(r::MutableDiffResult, x::Number, ::Type{Val{i}} = Val{1}) where {i}
    r.derivs = tuple_setindex(r.derivs, x, Val{i})
    return r
end

function derivative!(r::MutableDiffResult, x::AbstractArray, ::Type{Val{i}} = Val{1}) where {i}
    copyto!(derivative(r, Val{i}), x)
    return r
end

function derivative!(r::ImmutableDiffResult, x::Union{Number,StaticArray}, ::Type{Val{i}} = Val{1}) where {i}
    return ImmutableDiffResult(value(r), tuple_setindex(r.derivs, x, Val{i}))
end

function derivative!(r::ImmutableDiffResult, x::AbstractArray, ::Type{Val{i}} = Val{1}) where {i}
    T = tuple_eltype(r.derivs, Val{i})
    return ImmutableDiffResult(value(r), tuple_setindex(r.derivs, T(x), Val{i}))
end

"""
    derivative!(f, r::DiffResult, x, ::Type{Val{i}} = Val{1})

Equivalent to `derivative!(r::DiffResult, map(f, x), Val{i})`, but without the implied
temporary allocation (when possible).
"""
function derivative!(f, r::MutableDiffResult, x::Number, ::Type{Val{i}} = Val{1}) where {i}
    r.derivs = tuple_setindex(r.derivs, f(x), Val{i})
    return r
end

function derivative!(f, r::MutableDiffResult, x::AbstractArray, ::Type{Val{i}} = Val{1}) where {i}
    map!(f, derivative(r, Val{i}), x)
    return r
end

function derivative!(f, r::ImmutableDiffResult, x::Number, ::Type{Val{i}} = Val{1}) where {i}
    return derivative!(r, f(x), Val{i})
end

function derivative!(f, r::ImmutableDiffResult, x::StaticArray, ::Type{Val{i}} = Val{1}) where {i}
    return derivative!(r, map(f, x), Val{i})
end

function derivative!(f, r::ImmutableDiffResult, x::AbstractArray, ::Type{Val{i}} = Val{1}) where {i}
    T = tuple_eltype(r.derivs, Val{i})
    return derivative!(r, map(f, T(x)), Val{i})
end

# special-cased methods #
#-----------------------#

"""
    gradient(r::DiffResult)

Return the gradient stored in `r`.

Equivalent to `derivative(r, Val{1})`.
"""
gradient(r::DiffResult) = derivative(r)

"""
    gradient!(r::DiffResult, x)

Return `s::DiffResult` with the same data as `r`, except `gradient(s) == x`.

Equivalent to `derivative!(r, x, Val{1})`; see `derivative!` docs for aliasing behavior.
"""
gradient!(r::DiffResult, x) = derivative!(r, x)

"""
    gradient!(f, r::DiffResult, x)

Equivalent to `gradient!(r::DiffResult, map(f, x))`, but without the implied temporary
allocation (when possible).

Equivalent to `derivative!(f, r, x, Val{1})`; see `derivative!` docs for aliasing behavior.
"""
gradient!(f, r::DiffResult, x) = derivative!(f, r, x)

"""
    jacobian(r::DiffResult)

Return the Jacobian stored in `r`.

Equivalent to `derivative(r, Val{1})`.
"""
jacobian(r::DiffResult) = derivative(r)

"""
    jacobian!(r::DiffResult, x)

Return `s::DiffResult` with the same data as `r`, except `jacobian(s) == x`.

Equivalent to `derivative!(r, x, Val{1})`; see `derivative!` docs for aliasing behavior.
"""
jacobian!(r::DiffResult, x) = derivative!(r, x)

"""
    jacobian!(f, r::DiffResult, x)

Equivalent to `jacobian!(r::DiffResult, map(f, x))`, but without the implied temporary
allocation (when possible).

Equivalent to `derivative!(f, r, x, Val{1})`; see `derivative!` docs for aliasing behavior.
"""
jacobian!(f, r::DiffResult, x) = derivative!(f, r, x)

"""
    hessian(r::DiffResult)

Return the Hessian stored in `r`.

Equivalent to `derivative(r, Val{2})`.
"""
hessian(r::DiffResult) = derivative(r, Val{2})

"""
    hessian!(r::DiffResult, x)

Return `s::DiffResult` with the same data as `r`, except `hessian(s) == x`.

Equivalent to `derivative!(r, x, Val{2})`; see `derivative!` docs for aliasing behavior.
"""
hessian!(r::DiffResult, x) = derivative!(r, x, Val{2})

"""
    hessian!(f, r::DiffResult, x)

Equivalent to `hessian!(r::DiffResult, map(f, x))`, but without the implied temporary
allocation (when possible).

Equivalent to `derivative!(f, r, x, Val{2})`; see `derivative!` docs for aliasing behavior.
"""
hessian!(f, r::DiffResult, x) = derivative!(f, r, x, Val{2})

###################
# Pretty Printing #
###################

Base.show(io::IO, r::ImmutableDiffResult) = print(io, "ImmutableDiffResult($(r.value), $(r.derivs))")

Base.show(io::IO, r::MutableDiffResult) = print(io, "MutableDiffResult($(r.value), $(r.derivs))")

end # module
