#
###
# AbstractSciMLScalarOperator interface
###

function (L::AbstractSciMLScalarOperator)(v::Number, u, p, t; kwargs...)
    L = update_coefficients(L, u, p, t; kwargs...)
    convert(Number, L) * v
end

SCALINGNUMBERTYPES = (:AbstractSciMLScalarOperator,
    :Number,
    :UniformScaling)

#=
The identity operator must be listed here
so that rules for combination with scalar
operators take precedence over rules for
combining with the identity operator when
the two are combined together.
=#
SCALINGCOMBINETYPES = (:AbstractSciMLOperator,
    :IdentityOperator)

Base.length(::AbstractSciMLScalarOperator) = 1
Base.size(::AbstractSciMLScalarOperator) = ()
Base.adjoint(α::AbstractSciMLScalarOperator) = conj(α)
Base.transpose(α::AbstractSciMLScalarOperator) = α

has_mul!(::AbstractSciMLScalarOperator) = true
isconcrete(::AbstractSciMLScalarOperator) = true
islinear(::AbstractSciMLScalarOperator) = true
has_adjoint(::AbstractSciMLScalarOperator) = true

Base.:*(α::AbstractSciMLScalarOperator, u::AbstractArray) = convert(Number, α) * u
Base.:\(α::AbstractSciMLScalarOperator, u::AbstractArray) = convert(Number, α) \ u

function LinearAlgebra.rmul!(u::AbstractArray, α::AbstractSciMLScalarOperator)
    rmul!(u, convert(Number, α))
end
function LinearAlgebra.lmul!(α::AbstractSciMLScalarOperator, u::AbstractArray)
    lmul!(convert(Number, α), u)
end
function LinearAlgebra.ldiv!(α::AbstractSciMLScalarOperator, u::AbstractArray)
    ldiv!(convert(Number, α), u)
end
function LinearAlgebra.ldiv!(v::AbstractArray,
        α::AbstractSciMLScalarOperator,
        u::AbstractArray)
    ldiv!(v, convert(Number, α), u)
end

function LinearAlgebra.mul!(v::AbstractArray,
        α::AbstractSciMLScalarOperator,
        u::AbstractArray)
    x = convert(Number, α)
    mul!(v, x, u)
end

function LinearAlgebra.mul!(v::AbstractArray,
        α::AbstractSciMLScalarOperator,
        u::AbstractArray,
        a::Union{Number, AbstractSciMLScalarOperator},
        b::Union{Number, AbstractSciMLScalarOperator})
    α = convert(Number, α)
    a = convert(Number, a)
    b = convert(Number, b)
    mul!(v, α, u, a, b)
end

function LinearAlgebra.axpy!(α::AbstractSciMLScalarOperator,
        x::AbstractArray,
        y::AbstractArray)
    α = convert(Number, α)
    axpy!(α, x, y)
end

function LinearAlgebra.axpby!(α::AbstractSciMLScalarOperator,
        x::AbstractArray,
        β::Number,
        y::AbstractArray)
    α = convert(Number, α)
    axpby!(α, x, β, y)
end

function LinearAlgebra.axpby!(α::Number,
        x::AbstractArray,
        β::AbstractSciMLScalarOperator,
        y::AbstractArray)
    β = convert(Number, β)
    axpby!(α, x, β, y)
end

function LinearAlgebra.axpby!(α::AbstractSciMLScalarOperator,
        x::AbstractArray,
        β::AbstractSciMLScalarOperator,
        y::AbstractArray)
    α = convert(Number, α)
    β = convert(Number, β)
    axpby!(α, x, β, y)
end

Base.:+(α::AbstractSciMLScalarOperator) = α

mutable struct ScalarOperator{T <: Number, F} <: AbstractSciMLScalarOperator{T}
    val::T
    update_func::F
end

"""
$SIGNATURES

Represents a linear scaling operator that may be applied to a `Number`,
or an `AbstractArray` subtype. Its state is updated by the user-provided
`update_func` during operator evaluation (`L([w,] v, u, p, t)`), or by
calls to `update_coefficients[!]`. Both recursively call the
update function, `update_func` which is assumed to have the signature:

    update_func(oldval::Number, u, p, t; <accepted kwargs>) -> newval

The set of keyword-arguments accepted by `update_func` must be provided
to `ScalarOperator` via the kwarg `accepted_kwargs` as a tuple of `Symbol`s.
`kwargs` cannot be passed down to `update_func` if `accepted_kwargs`
are not provided.

$(UPDATE_COEFFS_WARNING)

# Interface

Lazy scalar algebra is defined for `AbstractSciMLScalarOperator`s. The
interface supports lazy addition, subtraction, multiplication and division.

# Example

```
v = rand(4)
u = rand(4)
w = zeros(4)
p = nothing
t = 0.0

val_update = (a, u, p, t; scale = 0.0) -> scale
α = ScalarOperator(0.0; update_func = val_update, accepted_kwargs = (:scale,))
β = 2 * α + 3 / α

# Update β and evaluate with the new interface
result = β(v, u, p, t; scale = 1.0)

# In-place application
β(w, v, u, p, t; scale = 1.0)

# In-place with scaling
w_orig = copy(w)
α_val = 2.0
β_val = 0.5
β(w, v, u, p, t, α_val, β_val; scale = 1.0) # w = α_val*(β*v) + β_val*w
```
"""
function ScalarOperator(val;
        update_func = DEFAULT_UPDATE_FUNC,
        accepted_kwargs = nothing)
    update_func = preprocess_update_func(update_func, accepted_kwargs)
    ScalarOperator(val, update_func)
end

# constructors
Base.convert(T::Type{<:Number}, α::ScalarOperator) = convert(T, α.val)
Base.convert(::Type{ScalarOperator}, α::Number) = ScalarOperator(α)

ScalarOperator(α::AbstractSciMLScalarOperator) = α
ScalarOperator(λ::UniformScaling) = ScalarOperator(λ.λ)

# traits
Base.show(io::IO, α::ScalarOperator) = print(io, "ScalarOperator($(α.val))")
function Base.conj(α::ScalarOperator) # TODO - test
    val = conj(α.val)
    update_func = (oldval, u, p, t;
        kwargs...) -> α.update_func(oldval |> conj,
        u,
        p,
        t;
        kwargs...) |> conj
    ScalarOperator(val; update_func = update_func, accepted_kwargs = NoKwargFilter())
end

Base.one(::AbstractSciMLScalarOperator{T}) where {T} = ScalarOperator(one(T))
Base.zero(::AbstractSciMLScalarOperator{T}) where {T} = ScalarOperator(zero(T))

Base.one(::Type{<:AbstractSciMLScalarOperator}) = ScalarOperator(true)
Base.zero(::Type{<:AbstractSciMLScalarOperator}) = ScalarOperator(false)
Base.abs(α::ScalarOperator) = abs(α.val)

Base.iszero(α::ScalarOperator) = iszero(α.val)

getops(α::ScalarOperator) = (α.val,)
isconstant(α::ScalarOperator) = update_func_isconstant(α.update_func)
has_ldiv(α::ScalarOperator) = !iszero(α.val)
has_ldiv!(α::ScalarOperator) = has_ldiv(α)

function update_coefficients!(L::ScalarOperator, u, p, t; kwargs...)
    L.val = L.update_func(L.val, u, p, t; kwargs...)
    nothing
end

function SciMLOperators.update_coefficients(L::ScalarOperator, u, p, t; kwargs...)
    return ScalarOperator(L.update_func(L.val, u, p, t; kwargs...), L.update_func)
end

# Copy method to avoid aliasing
function Base.copy(L::ScalarOperator)
    ScalarOperator(L.val, L.update_func)
end

# Add ScalarOperator specific implementations for the new interface
function (α::ScalarOperator)(v::AbstractArray, u, p, t; kwargs...)
    α = update_coefficients(α, u, p, t; kwargs...)
    convert(Number, α) * v
end

function (α::ScalarOperator)(w::AbstractArray, v::AbstractArray, u, p, t; kwargs...)
    update_coefficients!(α, u, p, t; kwargs...)
    mul!(w, α, v)
end

function (α::ScalarOperator)(w::AbstractArray, v::AbstractArray, u, p, t, a, b; kwargs...)
    update_coefficients!(α, u, p, t; kwargs...)
    mul!(w, α, v, a, b)
end

"""
$TYPEDEF

Lazy addition of `AbstractSciMLScalarOperator`s
"""
struct AddedScalarOperator{T, O} <: AbstractSciMLScalarOperator{T}
    ops::O

    function AddedScalarOperator(ops::NTuple{N, AbstractSciMLScalarOperator}) where {N}
        @assert !isempty(ops)
        T = promote_type(eltype.(ops)...)
        new{T, typeof(ops)}(ops)
    end
end

# constructors
function AddedScalarOperator(ops::AbstractSciMLScalarOperator...)
    AddedScalarOperator(ops)
end

Base.:+(ops::AbstractSciMLScalarOperator...) = AddedScalarOperator(ops...)
function Base.:+(A::AddedScalarOperator, B::AddedScalarOperator)
    AddedScalarOperator(A.ops..., B.ops...)
end
function Base.:+(A::AbstractSciMLScalarOperator, B::AddedScalarOperator)
    AddedScalarOperator(A, B.ops...)
end
function Base.:+(A::AddedScalarOperator, B::AbstractSciMLScalarOperator)
    AddedScalarOperator(A.ops..., B)
end
function Base.:-(A::AbstractSciMLScalarOperator, B::AbstractSciMLScalarOperator)
    AddedScalarOperator(A, -B)
end

for op in (:-, :+)
    for T in SCALINGNUMBERTYPES[2:end]
        @eval Base.$op(α::AbstractSciMLScalarOperator, x::$T) = AddedScalarOperator(α,
            ScalarOperator($op(x)))
        @eval Base.$op(x::$T, α::AbstractSciMLScalarOperator) = AddedScalarOperator(
            ScalarOperator(x),
            $op(α))
    end
end

function Base.convert(T::Type{<:Number}, α::AddedScalarOperator)
    sum(convert.(T, α.ops))
end

function Base.show(io::IO, α::AddedScalarOperator)
    print(io, "(")
    show(io, α.ops[1])
    for i in 2:length(α.ops)
        print(io, " + ")
        show(io, α.ops[i])
    end
    print(io, ")")
end
Base.conj(L::AddedScalarOperator) = AddedScalarOperator(conj.(L.ops))

function update_coefficients(L::AddedScalarOperator, u, p, t; kwargs...)
    ops = ()
    for op in L.ops
        ops = (ops..., update_coefficients(op, u, p, t; kwargs...))
    end

    AddedScalarOperator(ops)
end

function update_coefficients!(L::AddedScalarOperator, u, p, t; kwargs...)
    for op in L.ops
        update_coefficients!(op, u, p, t; kwargs...)
    end
    nothing
end
function (α::AddedScalarOperator)(v::AbstractArray, u, p, t; kwargs...)
    α = update_coefficients(α, u, p, t; kwargs...)
    convert(Number, α) * v
end

function (α::AddedScalarOperator)(w::AbstractArray, v::AbstractArray, u, p, t; kwargs...)
    update_coefficients!(α, u, p, t; kwargs...)
    mul!(w, α, v)
end

function (α::AddedScalarOperator)(
        w::AbstractArray, v::AbstractArray, u, p, t, a, b; kwargs...)
    update_coefficients!(α, u, p, t; kwargs...)
    mul!(w, α, v, a, b)
end

getops(α::AddedScalarOperator) = α.ops

# Copy method to avoid aliasing
function Base.copy(L::AddedScalarOperator)
    AddedScalarOperator(map(copy, L.ops))
end

has_ldiv(α::AddedScalarOperator) = !iszero(convert(Number, α))
has_ldiv!(α::AddedScalarOperator) = has_ldiv(α)

"""
$TYPEDEF

Lazy multiplication of `AbstractSciMLScalarOperator`s
"""
struct ComposedScalarOperator{T, O} <: AbstractSciMLScalarOperator{T}
    ops::O

    function ComposedScalarOperator(ops::NTuple{N, AbstractSciMLScalarOperator}) where {N}
        @assert !isempty(ops)
        T = promote_type(eltype.(ops)...)
        new{T, typeof(ops)}(ops)
    end
end

# constructor
function ComposedScalarOperator(ops::AbstractSciMLScalarOperator...)
    ComposedScalarOperator(ops)
end

for op in (:*, :∘)
    @eval Base.$op(ops::AbstractSciMLScalarOperator...) = reduce($op, ops)
    @eval Base.$op(
        A::AbstractSciMLScalarOperator, B::AbstractSciMLScalarOperator) = ComposedScalarOperator(
        A,
        B)
    @eval Base.$op(A::ComposedScalarOperator,
        B::AbstractSciMLScalarOperator) = ComposedScalarOperator(
        A.ops...,
        B)
    @eval Base.$op(A::AbstractSciMLScalarOperator,
        B::ComposedScalarOperator) = ComposedScalarOperator(
        A,
        B.ops...)
    @eval Base.$op(A::ComposedScalarOperator,
        B::ComposedScalarOperator) = ComposedScalarOperator(
        A.ops...,
        B.ops...)

    for T in SCALINGNUMBERTYPES[2:end]
        @eval Base.$op(α::AbstractSciMLScalarOperator, x::$T) = ComposedScalarOperator(α,
            ScalarOperator(x))
        @eval Base.$op(x::$T, α::AbstractSciMLScalarOperator) = ComposedScalarOperator(
            ScalarOperator(x),
            α)
    end
end

# Different methods for constant ScalarOperators
for T in SCALINGNUMBERTYPES[2:end]
    @eval function Base.:*(α::ScalarOperator, x::$T)
        if isconstant(α)
            return ScalarOperator(concretize(α) * x)
        else
            return ComposedScalarOperator(α, ScalarOperator(x))
        end
    end
    @eval function Base.:*(x::$T, α::ScalarOperator)
        if isconstant(α)
            return ScalarOperator(concretize(α) * x)
        else
            return ComposedScalarOperator(ScalarOperator(x), α)
        end
    end
end

function Base.convert(T::Type{<:Number}, α::ComposedScalarOperator)
    iszero(α) && return zero(T)
    prod(convert.(T, α.ops))
end

function Base.show(io::IO, α::ComposedScalarOperator)
    print(io, "(")
    show(io, α.ops[1])
    for i in 2:length(α.ops)
        print(io, " * ")
        show(io, α.ops[i])
    end
    print(io, ")")
end
Base.conj(L::ComposedScalarOperator) = ComposedScalarOperator(conj.(L.ops))
Base.:-(α::AbstractSciMLScalarOperator{T}) where {T} = (-one(T)) * α

function update_coefficients(L::ComposedScalarOperator, u, p, t; kwargs...)
    ops = ()
    for op in L.ops
        ops = (ops..., update_coefficients(op, u, p, t; kwargs...))
    end

    ComposedScalarOperator(ops)
end
function update_coefficients!(L::ComposedScalarOperator, u, p, t; kwargs...)
    for op in L.ops
        update_coefficients!(op, u, p, t; kwargs...)
    end
    nothing
end
function (α::ComposedScalarOperator)(v::AbstractArray, u, p, t; kwargs...)
    α = update_coefficients(α, u, p, t; kwargs...)
    convert(Number, α) * v
end

function (α::ComposedScalarOperator)(w::AbstractArray, v::AbstractArray, u, p, t; kwargs...)
    update_coefficients!(α, u, p, t; kwargs...)
    mul!(w, α, v)
end

function (α::ComposedScalarOperator)(
        w::AbstractArray, v::AbstractArray, u, p, t, a, b; kwargs...)
    update_coefficients!(α, u, p, t; kwargs...)
    mul!(w, α, v, a, b)
end

getops(α::ComposedScalarOperator) = α.ops

# Copy method to avoid aliasing
function Base.copy(L::ComposedScalarOperator)
    ComposedScalarOperator(map(copy, L.ops))
end

has_ldiv(α::ComposedScalarOperator) = all(has_ldiv, α.ops)
has_ldiv!(α::ComposedScalarOperator) = all(has_ldiv!, α.ops)

"""
$TYPEDEF

Lazy inverse of `AbstractSciMLScalarOperator`s
"""
struct InvertedScalarOperator{T, λType} <: AbstractSciMLScalarOperator{T}
    λ::λType

    function InvertedScalarOperator(λ::AbstractSciMLScalarOperator{T}) where {T}
        new{T, typeof(λ)}(λ)
    end
end
#=
Keeping with the style, we avoid use of the generic InvertedOperator and instead
have a specialized type for this purpose that subtypes AbstractSciMLScalarOperator.
=#
Base.inv(L::AbstractSciMLScalarOperator) = InvertedScalarOperator(L)

for op in (:/,)
    for T in SCALINGNUMBERTYPES[2:end]
        @eval Base.$op(α::AbstractSciMLScalarOperator, x::$T) = α * inv(ScalarOperator(x))
        @eval Base.$op(x::$T, α::AbstractSciMLScalarOperator) = ScalarOperator(x) * inv(α)
    end

    @eval Base.$op(α::AbstractSciMLScalarOperator, β::AbstractSciMLScalarOperator) = α *
                                                                                     inv(β)
end

for op in (:\,)
    for T in SCALINGNUMBERTYPES[2:end]
        @eval Base.$op(α::AbstractSciMLScalarOperator, x::$T) = inv(α) * ScalarOperator(x)
        @eval Base.$op(x::$T, α::AbstractSciMLScalarOperator) = inv(ScalarOperator(x)) * α
    end

    @eval Base.$op(α::AbstractSciMLScalarOperator, β::AbstractSciMLScalarOperator) = inv(α) *
                                                                                     β
end

function Base.convert(T::Type{<:Number}, α::InvertedScalarOperator)
    inv(convert(Number, α.λ))
end

function Base.show(io::IO, α::InvertedScalarOperator)
    print(io, "1 / ")
    show(io, α.λ)
end
Base.conj(L::InvertedScalarOperator) = InvertedScalarOperator(conj(L.λ))

function update_coefficients(L::InvertedScalarOperator, u, p, t; kwargs...)
    InvertedScalarOperator(update_coefficients(L.λ, u, p, t; kwargs...))
end
function update_coefficients!(L::InvertedScalarOperator, u, p, t; kwargs...)
    update_coefficients!(L.λ, u, p, t; kwargs...)
    nothing
end
function (α::InvertedScalarOperator)(v::AbstractArray, u, p, t; kwargs...)
    α = update_coefficients(α, u, p, t; kwargs...)
    convert(Number, α) * v
end

function (α::InvertedScalarOperator)(w::AbstractArray, v::AbstractArray, u, p, t; kwargs...)
    update_coefficients!(α, u, p, t; kwargs...)
    mul!(w, α, v)
end

function (α::InvertedScalarOperator)(
        w::AbstractArray, v::AbstractArray, u, p, t, a, b; kwargs...)
    update_coefficients!(α, u, p, t; kwargs...)
    mul!(w, α, v, a, b)
end
getops(α::InvertedScalarOperator) = (α.λ,)

# Copy method to avoid aliasing
function Base.copy(L::InvertedScalarOperator)
    InvertedScalarOperator(copy(L.λ))
end

has_ldiv(α::InvertedScalarOperator) = has_mul(α.λ)
has_ldiv!(α::InvertedScalarOperator) = has_ldiv(α)
#
