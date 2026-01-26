#
TENSOR_PROD_DOC = """
Computes the lazy pairwise Kronecker product, or tensor product,
operator of `AbstractMatrix`, and `AbstractSciMLOperator` subtypes.
Calling `⊗(ops...)` is equivalent to `Base.kron(ops...)`. Fast
operator evaluation is performed without forming the full tensor
product operator.

```
TensorProductOperator(A, B) = A ⊗ B
TensorProductOperator(A, B, C) = A ⊗ B ⊗ C

(A ⊗ B)(v) = vec(B * reshape(v, M, N) * transpose(A))
```
where `M = size(B, 2)`, and `N = size(A, 2)`

# Example

```
using SciMLOperators, LinearAlgebra

# Create basic operators
A = rand(3, 3)
B = rand(4, 4)
A_op = MatrixOperator(A)
B_op = MatrixOperator(B)

# Create tensor product operator
T = A_op ⊗ B_op

# Apply to a vector using the new interface
v = rand(3*4)    # Action vector
u = rand(3*4)    # Update vector
p = nothing
t = 0.0

# Out-of-place application
result = T(v, u, p, t)

# For in-place operations, need to cache the operator first
T_cached = cache_operator(T, v)

# In-place application
w = zeros(size(T, 1))
T_cached(w, v, u, p, t)

# In-place with scaling
w_orig = copy(w)
α = 2.0
β = 0.5
T_cached(w, v, u, p, t, α, β) # w = α*(T*v) + β*w_orig
```
"""

"""
$SIGNATURES

$TENSOR_PROD_DOC
"""
struct TensorProductOperator{T, O, C} <: AbstractSciMLOperator{T}
    ops::O
    cache::C

    function TensorProductOperator(
            ops::NTuple{
                2,
                Union{AbstractMatrix, AbstractSciMLOperator}
            },
            cache::Union{Tuple, Nothing})
        T = promote_type(eltype.(ops)...)

        new{T,
            typeof(ops),
            typeof(cache)
        }(ops, cache)
    end
end

function TensorProductOperator(outer::Union{AbstractMatrix, AbstractSciMLOperator},
        inner::Union{AbstractMatrix, AbstractSciMLOperator};
        cache = nothing)
    outer = outer isa AbstractMatrix ? MatrixOperator(outer) : outer
    inner = inner isa AbstractMatrix ? MatrixOperator(inner) : inner

    TensorProductOperator((outer, inner), cache)
end

# constructors
TensorProductOperator(ops...) = reduce(TensorProductOperator, ops)
TensorProductOperator(op::AbstractSciMLOperator) = op
TensorProductOperator(op::AbstractMatrix) = MatrixOperator(op)
function TensorProductOperator(ii1::IdentityOperator, ii2::IdentityOperator)
    IdentityOperator(ii1.len * ii2.len)
end

"""
$SIGNATURES

$TENSOR_PROD_DOC
"""
⊗(ops::Union{AbstractMatrix, AbstractSciMLOperator}...) = TensorProductOperator(ops...)

"""
$SIGNATURES

Construct a lazy representation of the Kronecker product `A ⊗ B`. One of the
two factors can be an `AbstractMatrix`, which is then promoted to a
`MatrixOperator` automatically. To avoid fallback to the generic
[`Base.kron`](@ref), at least one of `A` and `B` must be an
`AbstractSciMLOperator`.
"""
Base.kron(A::AbstractSciMLOperator, B::AbstractSciMLOperator) = TensorProductOperator(A, B)
Base.kron(A::AbstractMatrix, B::AbstractSciMLOperator) = TensorProductOperator(A, B)
Base.kron(A::AbstractSciMLOperator, B::AbstractMatrix) = TensorProductOperator(A, B)

Base.kron(ops::AbstractSciMLOperator...) = TensorProductOperator(ops...)

function Base.convert(::Type{AbstractMatrix}, L::TensorProductOperator)
    kron(convert.(AbstractMatrix, L.ops)...)
end

#LinearAlgebra.opnorm(L::TensorProductOperator) = prod(opnorm, L.ops)

function Base.show(io::IO, L::TensorProductOperator)
    print(io, "(")
    show(io, L.ops[1])
    print(io, " ⊗ ")
    show(io, L.ops[2])
    print(io, ")")
end
Base.size(L::TensorProductOperator) = reduce(.*, size.(L.ops))

for op in (:adjoint,
    :transpose)
    @eval function Base.$op(L::TensorProductOperator)
        TensorProductOperator($op.(L.ops)...;
            cache = issquare(L.ops[2]) ? L.cache : nothing)
    end
end
function Base.conj(L::TensorProductOperator)
    TensorProductOperator(conj.(L.ops)...; cache = L.cache)
end

function update_coefficients(L::TensorProductOperator, u, p, t)
    ops = ()
    for op in L.ops
        ops = (ops..., update_coefficients(op, u, p, t))
    end

    @reset L.ops = ops
end

getops(L::TensorProductOperator) = L.ops

# Copy method to avoid aliasing
function Base.copy(L::TensorProductOperator)
    TensorProductOperator(
        map(copy, L.ops),
        L.cache === nothing ? nothing : deepcopy(L.cache)
    )
end

islinear(L::TensorProductOperator) = reduce(&, islinear.(L.ops))
isconvertible(::TensorProductOperator) = false
Base.iszero(L::TensorProductOperator) = reduce(|, iszero.(L.ops))
has_adjoint(L::TensorProductOperator) = reduce(&, has_adjoint.(L.ops))
has_mul(L::TensorProductOperator) = reduce(&, has_mul.(L.ops))
has_mul!(L::TensorProductOperator) = reduce(&, has_mul!.(L.ops))
has_ldiv(L::TensorProductOperator) = reduce(&, has_ldiv.(L.ops))
has_ldiv!(L::TensorProductOperator) = reduce(&, has_ldiv!.(L.ops))

factorize(L::TensorProductOperator) = TensorProductOperator(factorize.(L.ops)...)

# operator application
function Base.:*(L::TensorProductOperator, v::AbstractVecOrMat)
    outer, inner = L.ops

    _, ni = size(inner)
    _, no = size(outer)
    m, n = size(L)
    k = size(v, 2)

    U = reshape(v, (ni, no * k))
    C = stack([inner * _v for _v in eachcol(U)])

    V = outer_mul(L, v, C)

    v isa AbstractMatrix ? reshape(V, (m, k)) : reshape(V, (m,))
end

function Base.:\(L::TensorProductOperator, v::AbstractVecOrMat)
    outer, inner = L.ops

    mi, _ = size(inner)
    mo, _ = size(outer)
    m, n = size(L)
    k = size(v, 2)

    U = reshape(v, (mi, mo * k))
    C = inner \ U

    V = outer_div(L, v, C)

    v isa AbstractMatrix ? reshape(V, (n, k)) : reshape(V, (n,))
end

function cache_self(L::TensorProductOperator, v::AbstractVecOrMat)
    outer, inner = L.ops

    mi, ni = size(inner)
    mo, no = size(outer)
    k = size(v, 2)

    # 3 arg mul!
    c1 = lmul!(false, similar(v, (mi, no * k))) # c1 = inner * v
    c2 = lmul!(false, similar(v, (no, mi, k))) # permute (2, 1, 3)
    c3 = lmul!(false, similar(v, (mo, mi * k))) # c3 = outer * c2

    # 5 arg mul!
    c4 = lmul!(false, similar(v, (mo * mi, k))) # cache v in 5 arg mul!

    # 3 arg ldiv!
    if reduce(&, issquare.(L.ops))
        c5, c6, c7 = c1, c2, c3
    else
        c5 = lmul!(false, similar(v, (ni, mo * k))) # c5 = inner \ v
        c6 = lmul!(false, similar(v, (mo, ni, k))) # permute (2, 1, 3)
        c7 = lmul!(false, similar(v, (no, ni * k))) # c7 = outer \ c6
    end

    @reset L.cache = (c1, c2, c3, c4, c5, c6, c7)
    L
end

function cache_internals(L::TensorProductOperator, v::AbstractVecOrMat)
    if !iscached(L)
        L = cache_self(L, v)
    end

    outer, inner = L.ops

    mi, ni = size(inner)
    _, no = size(outer)
    k = size(v, 2)

    vinner = reshape(v, (ni, no * k))
    vouter = reshape(L.cache[2], (no, mi * k))

    @reset L.ops[2] = cache_operator(inner, vinner)
    @reset L.ops[1] = cache_operator(outer, vouter)
    L
end

function LinearAlgebra.mul!(w::AbstractVecOrMat,
        L::TensorProductOperator,
        v::AbstractVecOrMat)
    @assert iscached(L) """cache needs to be set up for operator of type
    $L. Set up cache by calling `cache_operator(L, u)`"""

    outer, inner = L.ops

    _, ni = size(inner)
    _, no = size(outer)
    k = size(v, 2)

    C1, C2, C3 = L.cache[1:3]
    U = reshape(v, (ni, no * k))

    #=
        v .= kron(B, A) * v
        V .= A * U * B'
    =#

    # C .= A * U
    mul!(C1, inner, U)

    # V .= U * B' <===> V' .= B * C'
    outer_mul!(w, L, v)

    w
end

function LinearAlgebra.mul!(w::AbstractVecOrMat,
        L::TensorProductOperator,
        v::AbstractVecOrMat,
        α,
        β)
    @assert iscached(L) """cache needs to be set up for operator of type
    $L. Set up cache by calling `cache_operator(L, u)`"""

    outer, inner = L.ops

    mi, ni = size(inner)
    mo, no = size(outer)
    k = size(v, 2)

    C1 = first(L.cache)
    U = reshape(v, (ni, no * k))

    """
        v .= α * kron(B, A) * u + β * v
        V .= α * (A * U * B') + β * v
    """

    # C .= A * U
    mul!(C1, inner, U)

    # V = α(C * B') + β(V)
    c = reshape(C1, (mi * no, k))
    outer_mul!(w, L, c, α, β)

    w
end

function LinearAlgebra.ldiv!(w::AbstractVecOrMat,
        L::TensorProductOperator,
        v::AbstractVecOrMat)
    @assert iscached(L) """cache needs to be set up for operator of type
    $L. Set up cache by calling `cache_operator(L, u)`"""

    outer, inner = L.ops

    mi, ni = size(inner)
    mo, no = size(outer)
    k = size(v, 2)

    C5 = L.cache[5]
    U = reshape(v, (mi, mo * k))

    """
        v .= kron(B, A) ldiv v
        V .= (A ldiv U) / B'
    """

    # C .= A \ U
    ldiv!(C5, inner, U)

    # V .= C / B' <==> V' .= B \ C'
    c = reshape(C5, (ni * mo, k))
    outer_div!(w, L, c)

    w
end

function LinearAlgebra.ldiv!(L::TensorProductOperator, v::AbstractVecOrMat)
    outer, inner = L.ops

    msg = "Two-argument ldiv! is only available for square operators"
    @assert issquare(L) msg
    @assert issquare(inner) msg
    @assert issquare(outer) msg

    @assert iscached(L) """cache needs to be set up for operator of type
    $L. Set up cache by calling `cache_operator(L, v)`"""

    mi = size(inner, 1)
    mo = size(outer, 1)
    k = size(v, 2)

    U = reshape(v, (mi, mo * k))

    """
        u .= kron(B, A) ldiv u
        U .= (A ldiv U) / B'
    """

    # U .= A \ U
    ldiv!(inner, U)

    # U .= U / B' <==> U' .= B \ U'
    outer_div!(L, v)

    v
end

# helper functions
const PERM = (2, 1, 3)

function outer_mul(L::TensorProductOperator, v::AbstractVecOrMat, C::AbstractVecOrMat)
    outer, inner = L.ops

    if outer isa IdentityOperator
        return C
    elseif outer isa ScaledOperator
        return outer.λ * outer_mul(outer.L, v, C)
    end

    k = size(v, 2)
    if k == 1
        return transpose(outer * transpose(C))
    end

    mi, _ = size(inner)
    mo, no = size(outer)
    #   m , n  = size(L)

    C = reshape(C, (mi, no, k))
    C = permutedims(C, PERM)
    C = reshape(C, (no, mi * k))

    V = outer * C
    V = reshape(V, (mo, mi, k))
    V = permutedims(V, PERM)

    V
end

function outer_mul!(w::AbstractVecOrMat, L::TensorProductOperator, v::AbstractVecOrMat)
    outer, inner = L.ops

    C1 = first(L.cache)

    if outer isa IdentityOperator
        copyto!(w, C1)
        return w
    elseif outer isa ScaledOperator
        outer_mul!(w, outer.L, v)
        lmul!(outer.λ, w)
        return w
    end

    mi, _ = size(inner)
    mo, no = size(outer)
    #   m , n  = size(L)
    k = size(v, 2)

    if k == 1
        W = reshape(w, (mi, mo))
        C1 = reshape(C1, (mi, no))
        mul!(transpose(W), outer, transpose(C1))
        return w
    end

    C2, C3 = L.cache[2:3]

    C1 = reshape(C1, (mi, no, k))
    permutedims!(C2, C1, PERM)
    C2 = reshape(C2, (no, mi * k))
    mul!(C3, outer, C2)
    C3 = reshape(C3, (mo, mi, k))
    W = reshape(w, (mi, mo, k))
    permutedims!(W, C3, PERM)

    w
end

function outer_mul!(w::AbstractVecOrMat, L::TensorProductOperator,
        v::AbstractVecOrMat, α, β)
    outer, inner = L.ops

    m, _ = size(L)
    k = size(v, 2)

    if outer isa IdentityOperator
        v = reshape(v, (m, k))
        axpby!(α, v, β, w)
        return w
    elseif outer isa ScaledOperator
        a = convert(Number, α * outer.λ)
        outer_mul!(w, outer.L, v, a, β)
        return w
    end

    mi, _ = size(inner)
    mo, no = size(outer)

    if k == 1
        W = reshape(w, (mi, mo))
        C = reshape(v, (mi, no))
        mul!(transpose(W), outer, transpose(C), α, β)
        return w
    end

    C2, C3, c4 = L.cache[2:4]

    C = reshape(v, (mi, no, k))
    permutedims!(C2, C, PERM)
    C2 = reshape(C2, (no, mi * k))
    mul!(C3, outer, C2)
    C3 = reshape(C3, (mo, mi, k))
    W = reshape(w, (mi, mo, k))
    copy!(c4, w)
    permutedims!(W, C3, PERM)
    axpby!(β, c4, α, w)

    w
end

function outer_div(L::TensorProductOperator, v::AbstractVecOrMat, C::AbstractVecOrMat)
    outer, inner = L.ops

    if outer isa IdentityOperator
        return C
    elseif outer isa ScaledOperator
        return outer.λ \ outer_div(outer.L, v, C)
    end

    k = size(v, 2)
    if k == 1
        return transpose(outer \ transpose(C))
    end

    _, ni = size(inner)
    mo, no = size(outer)

    C = reshape(C, (ni, mo, k))
    C = permutedims(C, PERM)
    C = reshape(C, (mo, ni * k))

    V = outer \ C
    V = reshape(V, (no, ni, k))
    V = permutedims(V, PERM)

    V
end

function outer_div!(v::AbstractVecOrMat, L::TensorProductOperator, c::AbstractVecOrMat)
    outer, inner = L.ops

    if outer isa IdentityOperator
        copyto!(v, c)
        return v
    elseif outer isa ScaledOperator
        outer_div!(v, outer.L, c)
        ldiv!(outer.λ, v)
        return v
    end

    mi, ni = size(inner)
    mo, no = size(outer)
    k = size(c, 2)

    if k == 1
        V = reshape(v, (ni, no))
        C = reshape(c, (ni, mo))
        ldiv!(transpose(V), outer, transpose(C))
        return v
    end

    C6, C7 = L.cache[6:7]

    C = reshape(c, (ni, mo, k))
    permutedims!(C6, C, PERM)
    C6 = reshape(C6, (mo, ni * k))
    ldiv!(C7, outer, C6)
    C7 = reshape(C7, (no, ni, k))
    V = reshape(v, (ni, no, k))
    permutedims!(V, C7, PERM)

    v
end

function outer_div!(L::TensorProductOperator, v::AbstractVecOrMat)
    outer, inner = L.ops

    if outer isa IdentityOperator
        return v
    elseif outer isa ScaledOperator
        outer_div!(outer.L, v)
        ldiv!(outer.λ, v)
        return u
    end

    _, ni = size(inner)
    _, no = size(outer)
    k = size(v, 2)

    U = reshape(v, (ni, no * k))

    if k == 1
        ldiv!(outer, transpose(U))
        return v
    end

    C = first(L.cache)

    U = reshape(U, (ni, no, k))
    C = reshape(C, (no, ni, k))
    permutedims!(C, U, PERM)
    C = reshape(C, (no, ni * k))
    ldiv!(outer, C)
    C = reshape(C, (no, ni, k))
    permutedims!(U, C, PERM)

    v
end

# Out-of-place: v is action vector, u is update vector
function (L::TensorProductOperator)(v::AbstractVecOrMat, u, p, t; kwargs...)
    L = update_coefficients(L, u, p, t; kwargs...)
    L * v
end

# In-place: w is destination, v is action vector, u is update vector
function (L::TensorProductOperator)(
        w::AbstractVecOrMat, v::AbstractVecOrMat, u, p, t; kwargs...)
    update_coefficients!(L, u, p, t; kwargs...)
    mul!(w, L, v)
    return w
end

# In-place with scaling: w = α*(L*v) + β*w
function (L::TensorProductOperator)(
        w::AbstractVecOrMat, v::AbstractVecOrMat, u, p, t, α, β; kwargs...)
    update_coefficients!(L, u, p, t; kwargs...)
    mul!(w, L, v, α, β)
    return w
end
#
