#

###
# operator update interface
###

"""
$SIGNATURES

The default update function for `AbstractSciMLOperator`s, a no-op that
leaves the operator state unchanged.
"""
DEFAULT_UPDATE_FUNC(A, u, p, t) = A

const UPDATE_COEFFS_WARNING = """
!!! warning
    The user-provided `update_func[!]` must not use `u` in
    its computation. Positional argument `(u, p, t)` to `update_func[!]` are
    passed down by `update_coefficients[!](L, u, p, t)`, where `u` is the
    input-vector to the composite `AbstractSciMLOperator`. For that reason,
    the values of `u`, or even shape, may not correspond to the input
    expected by `update_func[!]`. If an operator's state depends on its
    input vector, then it is, by definition, a nonlinear operator.
    We recommend sticking such nonlinearities in `FunctionOperator.`
    This topic is further discussed in
    (this issue)[https://github.com/SciML/SciMLOperators.jl/issues/159].
"""

"""
$SIGNATURES

Update the state of `L` based on `u`, input vector, `p` parameter object,
`t`, and keyword arguments. Internally, `update_coefficients` calls the
user-provided `update_func` method for every component operator in `L`
with the positional arguments `(u, p, t)` and keyword arguments
corresponding to the symbols provided to the operator via kwarg
`accepted_kwargs`.

This method is out-of-place, i.e. fully non-mutating and `Zygote`-compatible.

$(UPDATE_COEFFS_WARNING)

# Example

```
using SciMLOperators

mat_update_func = (A, u, p, t; scale = 1.0) -> p * p' * scale * t

M = MatrixOperator(zeros(4,4); update_func = mat_update_func,
                   accepted_kwargs = (:state,))

L = M + IdentityOperator(4)

u = rand(4)
p = rand(4)
t = 1.0

# Update the operator to `(u,p,t)` and apply it to `v`
L = update_coefficients(L, u, p, t; scale = 2.0)
result = L * v

# Or use the interface which separrates the update from the application
result = L(v, u, p, t; scale = 2.0)
```
"""
update_coefficients(L, u, p, t; kwargs...) = L

"""
$SIGNATURES

Update in-place the state of `L` based on `u`, input vector, `p` parameter
object, `t`, and keyword arguments. Internally, `update_coefficients!` calls
the user-provided mutating `update_func!` method for every component operator
in `L` with the positional arguments `(u, p, t)` and keyword arguments
corresponding to the symbols provided to the operator via kwarg
`accepted_kwargs`.

$(UPDATE_COEFFS_WARNING)

# Example

```
using SciMLOperators

_A = rand(4, 4)
mat_update_func! = (L, u, p, t; scale = 1.0) -> copy!(A, _A)

M = MatrixOperator(zeros(4,4); update_func! = mat_update_func!)

L = M + IdentityOperator(4)

u = rand(4)
p = rand(4)
t = 1.0

update_coefficients!(L, u, p, t)
L * v
```
"""
update_coefficients!(L, u, p, t; kwargs...) = nothing

# We cannot use @generate because we don't know the type structure of L in advance
function update_coefficients!(L::AbstractSciMLOperator, u, p, t; kwargs...)
    foreach(op -> update_coefficients!(op, u, p, t; kwargs...), getops(L))

    nothing
end

###
# operator evaluation interface
###

# Out-of-place: v is action vector, u is update vector
function (L::AbstractSciMLOperator)(v, u, p, t; kwargs...)
    update_coefficients(L, u, p, t; kwargs...) * v
end
# In-place: w is destination, v is action vector, u is update vector
function (L::AbstractSciMLOperator)(w, v, u, p, t; kwargs...)
    (update_coefficients!(L, u, p, t; kwargs...); mul!(w, L, v))
end
# In-place with scaling: w = α*(L*v) + β*w
function (L::AbstractSciMLOperator)(w, v, u, p, t, α, β; kwargs...)
    (update_coefficients!(L, u, p, t; kwargs...); mul!(w, L, v, α, β))
end

function (L::AbstractSciMLOperator)(w::Number, v::Number, u, p, t, args...; kwargs...)
    msg = """Nonallocating L(w, v, u, p, t) type methods are not available for
    subtypes of `Number`."""
    throw(ArgumentError(msg))
end

###
# operator caching interface
###

getops(L) = ()

"""
$SIGNATURES

Checks whether `L` has preallocated caches for inplace evaluations.
"""
function iscached(L::AbstractSciMLOperator)
    has_cache = hasfield(typeof(L), :cache) # TODO - confirm this is static
    isset = has_cache ? !isnothing(L.cache) : true

    return isset & all(iscached, getops(L))
end

"""
Check if `SciMLOperator` `L` has preallocated cache-arrays for in-place
computation.
"""
iscached(L) = true

iscached(::Union{# LinearAlgebra
    AbstractMatrix,
    UniformScaling,
    Factorization, # Base
    Number}) = true

"""
$SIGNATURES

Allocate caches for `L` for in-place evaluation with `u`-like input vectors.
"""
cache_operator(L, u) = L

function cache_operator(L::AbstractSciMLOperator, u::AbstractVecOrMat)
    L = cache_self(L, u)
    L = cache_internals(L, u)
    L
end

cache_self(L::AbstractSciMLOperator, ::AbstractVecOrMat) = L
cache_internals(L::AbstractSciMLOperator, ::AbstractVecOrMat) = L

###
# operator traits
###

Base.size(A::AbstractSciMLOperator, d::Integer) = d <= 2 ? size(A)[d] : 1
Base.eltype(::Type{AbstractSciMLOperator{T}}) where {T} = T
Base.eltype(::AbstractSciMLOperator{T}) where {T} = T

Base.oneunit(L::AbstractSciMLOperator) = one(L)
Base.oneunit(LType::Type{<:AbstractSciMLOperator}) = one(LType)

Base.iszero(::AbstractSciMLOperator) = false # TODO

"""
$SIGNATURES

Check if `adjoint(L)` is lazily defined.
"""
has_adjoint(L::AbstractSciMLOperator) = false # L', adjoint(L)
"""
$SIGNATURES

Check if `expmv!(w, L, v, t)`, equivalent to `mul!(w, exp(t * A), v)`, is
defined for `Number` `t`, and `AbstractArray`s `w, v` of appropriate sizes.
"""
has_expmv!(L::AbstractSciMLOperator) = false # expmv!(v, L, t, u)
"""
$SIGNATURES

Check if `expmv(L, v, t)`, equivalent to `exp(t * A) * v`, is defined for
`Number` `t`, and `AbstractArray` `u` of appropriate size.
"""
has_expmv(L::AbstractSciMLOperator) = false # v = exp(L, t, u)
"""
$SIGNATURES

Check if `exp(L)` is defined lazily defined.
"""
has_exp(L::AbstractSciMLOperator) = islinear(L)
"""
$SIGNATURES

Check if `L * v` is defined for `AbstractArray` `u` of appropriate size.
"""
has_mul(L::AbstractSciMLOperator) = true # du = L*u
"""
$SIGNATURES

Check if `mul!(w, L, v)` is defined for `AbstractArray`s `w, v` of
appropriate sizes.
"""
has_mul!(L::AbstractSciMLOperator) = true # mul!(du, L, u)
"""
$SIGNATURES

Check if `L \\ v` is defined for `AbstractArray` `v` of appropriate size.
"""
has_ldiv(L::AbstractSciMLOperator) = false # du = L\u
"""
$SIGNATURES

Check if `ldiv!(w, L, v)` is defined for `AbstractArray`s `w, v` of
appropriate sizes.
"""
has_ldiv!(L::AbstractSciMLOperator) = false # ldiv!(du, L, u)

### Extra standard assumptions

"""
$SIGNATURES

Checks if an `L`'s state is constant or needs to be updated by calling
`update_coefficients`.
"""
isconstant(::Union{# LinearAlgebra
    AbstractMatrix,
    UniformScaling,
    Factorization, # Base
    Number}) = true
isconstant(L::AbstractSciMLOperator) = all(isconstant, getops(L))
isconstant(L) = false

"""
    isconvertible(L) -> Bool

Checks if `L` can be cheaply converted to an `AbstractMatrix` via eager fusion.
"""
isconvertible(L::AbstractSciMLOperator) = all(isconvertible, getops(L))

function isconvertible(::Union{
        # LinearAlgebra
        AbstractMatrix,
        UniformScaling,
        Factorization,

        # Base
        Number,

        # SciMLOperators
        AbstractSciMLScalarOperator
})
    true
end

"""
    concretize(L) -> AbstractMatrix

    concretize(L) -> Number

Convert `SciMLOperator` to a concrete type via eager fusion. This method is a
no-op for types that are already concrete.
"""
function concretize(L::Union{# LinearAlgebra
        AbstractMatrix,
        Factorization, # SciMLOperators
        AbstractSciMLOperator
})
    convert(AbstractMatrix, L)
end

function concretize(L::Union{
        # LinearAlgebra
        UniformScaling,

        # Base
        Number,

        # SciMLOperators
        AbstractSciMLScalarOperator
})
    convert(Number, L)
end

"""
$SIGNATURES

Checks if `L` is a linear operator.
"""
islinear(::AbstractSciMLOperator) = false
islinear(L) = false

islinear(::Union{# LinearAlgebra
    AbstractMatrix,
    UniformScaling,
    Factorization, # Base
    Number
}) = true

has_mul(L) = false
has_mul(::Union{# LinearAlgebra
    AbstractVecOrMat,
    AbstractMatrix,
    UniformScaling, # Base
    Number
}) = true

has_mul!(L) = false
has_mul!(::Union{# LinearAlgebra
    AbstractVecOrMat,
    AbstractMatrix,
    UniformScaling, # Base
    Number
}) = true

has_ldiv(L) = false
has_ldiv(::Union{
    AbstractMatrix,
    Factorization,
    Number
}) = true

has_ldiv!(L) = false
has_ldiv!(::Union{
    Diagonal,
    Bidiagonal,
    Factorization
}) = true

has_adjoint(L) = islinear(L)
has_adjoint(::Union{# LinearAlgebra
    AbstractMatrix,
    UniformScaling,
    Factorization, # Base
    Number
}) = true

"""
Checks if `size(L, 1) == size(L, 2)`.
"""
issquare(L) = ndims(L) >= 2 && size(L, 1) == size(L, 2)
issquare(::AbstractVector) = false
issquare(::Union{# LinearAlgebra
    UniformScaling, # SciMLOperators
    AbstractSciMLScalarOperator, # Base
    Number
}) = true
issquare(A...) = @. (&)(issquare(A)...)

Base.length(L::AbstractSciMLOperator) = prod(size(L))
Base.ndims(L::AbstractSciMLOperator) = length(size(L))
Base.isreal(L::AbstractSciMLOperator{T}) where {T} = T <: Real
Base.Matrix(L::AbstractSciMLOperator) = Matrix(convert(AbstractMatrix, L))

LinearAlgebra.exp(L::AbstractSciMLOperator, t) = exp(t * L)
expmv(L::AbstractSciMLOperator, u, p, t) = exp(L, t) * u
expmv!(v, L::AbstractSciMLOperator, u, p, t) = mul!(v, exp(L, t), u)

###
# fallback implementations
###

function Base.conj(L::AbstractSciMLOperator)
    isreal(L) && return L
    if !isconvertible(L)
        @warn """using convert-based fallback for Base.conj"""
    end
    concretize(L) |> conj
end

function Base.:(==)(L1::AbstractSciMLOperator, L2::AbstractSciMLOperator)
    if !isconvertible(L1) || !isconvertible(L2)
        @warn """using convert-based fallback for Base.=="""
    end
    size(L1) != size(L2) && return false
    concretize(L1) == concretize(L2)
end

function Base.getindex(L::AbstractSciMLOperator, I::Vararg{Int, N}) where {N}
    if !isconvertible(L)
        @warn """using convert-based fallback for Base.getindex"""
    end
    concretize(L)[I...]
end

function Base.resize!(L::AbstractSciMLOperator, n::Integer)
    throw(MethodError(resize!, typeof.((L, n))))
end

LinearAlgebra.exp(L::AbstractSciMLOperator) = exp(Matrix(L))

function LinearAlgebra.opnorm(L::AbstractSciMLOperator, p::Real = 2)
    if !isconvertible(L)
        @warn """using convert-based fallback in LinearAlgebra.opnorm."""
    end
    opnorm(concretize(L), p)
end

for op in (:sum, :prod)
    @eval function Base.$op(L::AbstractSciMLOperator; kwargs...)
        if !isconvertible(L)
            @warn """using convert-based fallback in $($op)."""
        end
        $op(concretize(L); kwargs...)
    end
end

for pred in (:issymmetric,
    :ishermitian,
    :isposdef)
    @eval function LinearAlgebra.$pred(L::AbstractSciMLOperator)
        if !isconvertible(L)
            @warn """using convert-based fallback in $($pred)."""
        end
        $pred(concretize(L))
    end
end

function LinearAlgebra.mul!(v::AbstractArray, L::AbstractSciMLOperator, u::AbstractArray)
    if !isconvertible(L)
        @warn """using convert-based fallback in mul!."""
    end
    mul!(v, concretize(L), u)
end

function LinearAlgebra.mul!(v::AbstractArray,
        L::AbstractSciMLOperator,
        u::AbstractArray,
        α,
        β)
    if !isconvertible(L)
        @warn """using convert-based fallback in mul!."""
    end
    mul!(v, concretize(L), u, α, β)
end

#
