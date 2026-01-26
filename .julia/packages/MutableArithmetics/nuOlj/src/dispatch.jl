# Copyright (c) 2019 MutableArithmetics.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v.2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

# TODO: this file contains a large number of method specializations to intercept
# "externally owned" method calls by dispatching on type parameters (rather than
# outermost wrapper type). This is generally bad practice, but refactoring this
# code to use a different mechanism would be a lot of work. In the future, this
# interception code would be more easily/robustly replaced by using a tool like
# https://github.com/jrevels/Cassette.jl.

abstract type AbstractMutable end

function Base.sum(
    a::AbstractArray{T};
    dims = :,
    init = zero(promote_operation(+, T, T)),
) where {T<:AbstractMutable}
    if dims !== Colon()
        # We cannot use `mapreduce` with `add!!` instead of `Base.add_mul` like
        # `operate(sum, ...)` because the same instance given at `init` is used
        # at several places.
        return mapreduce(identity, Base.add_sum, a; dims, init)
    end
    return operate(sum, a; init)
end

# When doing `x'y` where the elements of `x` and/or `y` are arrays, redirecting
# to `dot(x, y)` is not equivalent to `x'y` as it will call dot recursively on
# the elements of `x` and `y`. See
# https://github.com/JuliaLang/julia/issues/35174
# For this reason, a `_dot_nonrecursive` function was added that does not
# recursively call `dot`:
# https://github.com/JuliaLang/julia/commit/eae3216416453b53631afa6c803591cf2c5ae5b3
#
# However, it does not exploit mutability and returns
# `zero(eltype(lhs)) * zero(eltype(rhs))` in case the arrays are empty which
# creates type instability for some types for which this type is not invariant
# under addition.

# TODO: LinearAlgebra should have a documented function so that we don't have to
# overload an internal function
function LinearAlgebra._dot_nonrecursive(
    lhs::AbstractArray{<:AbstractMutable},
    rhs::AbstractArray,
)
    return fused_map_reduce(add_mul, lhs, rhs)
end
function LinearAlgebra._dot_nonrecursive(
    lhs::AbstractArray,
    rhs::AbstractArray{<:AbstractMutable},
)
    return fused_map_reduce(add_mul, lhs, rhs)
end
function LinearAlgebra._dot_nonrecursive(
    lhs::AbstractArray{<:AbstractMutable},
    rhs::AbstractArray{<:AbstractMutable},
)
    return fused_map_reduce(add_mul, lhs, rhs)
end

for A in (LinearAlgebra.Symmetric, LinearAlgebra.Hermitian, AbstractArray)
    B = A{<:AbstractMutable}
    @eval begin
        LinearAlgebra.dot(x::$A, y::$B) = operate(LinearAlgebra.dot, x, y)
        LinearAlgebra.dot(x::$B, y::$A) = operate(LinearAlgebra.dot, x, y)
        LinearAlgebra.dot(x::$B, y::$B) = operate(LinearAlgebra.dot, x, y)
    end
end

# Special-case because the the base version wants to do
# fill!(::Array{AbstractVariableRef}, zero(GenericAffExpr{Float64,eltype(x)}))
_one_indexed(A) = all(x -> isa(x, Base.OneTo), axes(A))

function LinearAlgebra.diagm_container(
    size,
    kv::Pair{<:Integer,<:AbstractVector{<:AbstractMutable}}...,
)
    T = promote_type(map(x -> promote_type(eltype(x.second)), kv)...)
    U = promote_type(T, promote_operation(zero, T))
    return zeros(U, LinearAlgebra.diagm_size(size, kv...)...)
end

function LinearAlgebra.diagm(x::AbstractVector{<:AbstractMutable})
    # `LinearAlgebra.diagm` doesn't work for non-one-indexed arrays in general.
    @assert _one_indexed(x)
    ZeroType = promote_operation(zero, eltype(x))
    return LinearAlgebra.diagm(0 => copyto!(similar(x, ZeroType), x))
end

################################################################################
# Interception of Base's matrix/vector arithmetic machinery
#
# Redirect calls with `eltype(ret) <: AbstractMutable` to `_mul!` to replace it
# with an implementation more efficient than `generic_matmatmul!` and
# `generic_matvecmul!` since it takes into account the mutability of the
# arithmetic. We need `args...` because SparseArrays` also gives `α` and `β`
# arguments.

function _mul!(output, A, B, α, β)
    # See SparseArrays/src/linalg.jl
    if !isone(β)
        if iszero(β)
            operate!(zero, output)
        else
            rmul!(output, scaling(β))
        end
    end
    return operate!(add_mul, output, A, B, scaling(α))
end

function _mul!(output, A, B, α)
    operate!(zero, output)
    return operate!(add_mul, output, A, B, scaling(α))
end

# LinearAlgebra uses `Base.promote_op(LinearAlgebra.matprod, ...)` to try to
# infere the return type. If the operation is not supported, it returns
# `Union{}`.
function _mul!(output::AbstractArray{Union{}}, A, B)
    # Normally, if the product is not supported, this should redirect to
    # `MA.promote_operation(*, ...)` which redirects to
    # `zero(...) * zero(...)` which should throw an appropriate error.
    # For example, in JuMP, it would say that you cannot multiply quadratic
    # expressions with an affine expression for instance.
    ProdType = promote_array_mul(typeof(A), typeof(B))
    # If we arrived here, it means that we have found a type for `output`, even
    # if LinearAlgebra couldn't. This is most probably a but so let's provide
    # extensive information to help debugging.
    return error(
        "Cannot multiply a `$(typeof(A))` with a `$(typeof(B))` because the " *
        "sum of the product of a `$(eltype(A))` and a `$(eltype(B))` could " *
        "not be inferred so a `$(typeof(output))` allocated to store the " *
        "output of the multiplication instead of a `$ProdType`.",
    )
end

_mul!(output, A, B) = operate_to!(output, *, A, B)

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::AbstractVecOrMat,
    B::AbstractVecOrMat,
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractVector{<:AbstractMutable},
    A::AbstractVecOrMat,
    B::AbstractVector,
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractVector{<:AbstractMutable},
    A::LinearAlgebra.Transpose{<:Any,<:AbstractVecOrMat},
    B::AbstractVector,
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractVector{<:AbstractMutable},
    A::LinearAlgebra.Adjoint{<:Any,<:AbstractVecOrMat},
    B::AbstractVector,
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractVector{<:AbstractMutable},
    A::LinearAlgebra.AbstractTriangular,
    B::AbstractVector,
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::LinearAlgebra.Transpose{<:Any,<:AbstractVecOrMat},
    B::AbstractMatrix,
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::LinearAlgebra.Adjoint{<:Any,<:AbstractVecOrMat},
    B::AbstractMatrix,
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::LinearAlgebra.AbstractTriangular,
    B::AbstractMatrix,
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::AbstractMatrix,
    B::LinearAlgebra.Transpose{<:Any,<:AbstractVecOrMat},
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::AbstractMatrix,
    B::LinearAlgebra.Adjoint{<:Any,<:AbstractVecOrMat},
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::LinearAlgebra.Adjoint{<:Any,<:AbstractVecOrMat},
    B::LinearAlgebra.Adjoint{<:Any,<:AbstractVecOrMat},
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::LinearAlgebra.Transpose{<:Any,<:AbstractVecOrMat},
    B::LinearAlgebra.Transpose{<:Any,<:AbstractVecOrMat},
)
    return _mul!(ret, A, B)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::AbstractVecOrMat,
    B::AbstractVecOrMat,
    α::Number,
    β::Number,
)
    return _mul!(ret, A, B, α, β)
end

function LinearAlgebra.mul!(
    ret::AbstractVector{<:AbstractMutable},
    A::AbstractVecOrMat,
    B::AbstractVector,
    α::Number,
    β::Number,
)
    return _mul!(ret, A, B, α, β)
end

function LinearAlgebra.mul!(
    ret::AbstractVector{<:AbstractMutable},
    A::LinearAlgebra.Transpose{<:Any,<:AbstractVecOrMat},
    B::AbstractVector,
    α::Number,
    β::Number,
)
    return _mul!(ret, A, B, α, β)
end

function LinearAlgebra.mul!(
    ret::AbstractVector{<:AbstractMutable},
    A::LinearAlgebra.Adjoint{<:Any,<:AbstractVecOrMat},
    B::AbstractVector,
    α::Number,
    β::Number,
)
    return _mul!(ret, A, B, α, β)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::LinearAlgebra.Transpose{<:Any,<:AbstractVecOrMat},
    B::AbstractMatrix,
    α::Number,
    β::Number,
)
    return _mul!(ret, A, B, α, β)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::LinearAlgebra.Adjoint{<:Any,<:AbstractVecOrMat},
    B::AbstractMatrix,
    α::Number,
    β::Number,
)
    return _mul!(ret, A, B, α, β)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::AbstractMatrix,
    B::LinearAlgebra.Transpose{<:Any,<:AbstractVecOrMat},
    α::Number,
    β::Number,
)
    return _mul!(ret, A, B, α, β)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::AbstractMatrix,
    B::LinearAlgebra.Adjoint{<:Any,<:AbstractVecOrMat},
    α::Number,
    β::Number,
)
    return _mul!(ret, A, B, α, β)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::LinearAlgebra.Adjoint{<:Any,<:AbstractVecOrMat},
    B::LinearAlgebra.Adjoint{<:Any,<:AbstractVecOrMat},
    α::Number,
    β::Number,
)
    return _mul!(ret, A, B, α, β)
end

function LinearAlgebra.mul!(
    ret::AbstractMatrix{<:AbstractMutable},
    A::LinearAlgebra.Transpose{<:Any,<:AbstractVecOrMat},
    B::LinearAlgebra.Transpose{<:Any,<:AbstractVecOrMat},
    α::Number,
    β::Number,
)
    return _mul!(ret, A, B, α, β)
end

# SparseArrays promotes the element types of `A` and `B` to the same type which,
# always produce quadratic expressions for JuMP even if only one of them was
# affine and the other one constant. Moreover, it does not always go through
# `LinearAlgebra.mul!` which prevents us from using mutability of the
# arithmetic. For this reason we intercept the calls and redirect them to `mul`.

const _LinearAlgebraWrappers = (
    LinearAlgebra.Adjoint,
    LinearAlgebra.Transpose,
    # TODO(odow): we could expand these overloads to other LinearAlgebra types.
    # LinearAlgebra.Symmetric,
    # LinearAlgebra.Hermitian,
    # LinearAlgebra.Diagonal,
    # LinearAlgebra.LowerTriangular,
    # LinearAlgebra.UpperTriangular,
    # LinearAlgebra.UnitLowerTriangular,
    # LinearAlgebra.UnitUpperTriangular,
)

const _MatrixLike = vcat(
    Any[T -> LA{<:T,<:_SparseMat} for LA in _LinearAlgebraWrappers],
    Any[T->_SparseMat{<:T}, T->StridedMatrix{<:T}],
)

for f_A in _MatrixLike, f_B in vcat(_MatrixLike, T -> StridedVector{<:T})
    A, mut_A = f_A(Any), f_A(AbstractMutable)
    B, mut_B = f_B(Any), f_B(AbstractMutable)
    if A <: StridedMatrix && B <: StridedMatrix
        continue
    end
    @eval begin
        Base.:*(a::$(mut_A), b::$(B)) = mul(a, b)
        Base.:*(a::$(A), b::$(mut_B)) = mul(a, b)
        Base.:*(a::$(mut_A), b::$(mut_B)) = mul(a, b)
    end
end

# See https://github.com/JuliaLang/julia/pull/37898
# The default fallback only used `promote_type` so it may get its wrong, e.g.,
# for JuMP and MultivariatePolynomials.
if VERSION >= v"1.7.0-DEV.1284"
    using LinearAlgebra: StridedMaybeAdjOrTransMat
    _mat_mat_scalar(A, B, γ) = operate!!(*, operate(*, A, B), γ)
    function LinearAlgebra.mat_mat_scalar(
        A::StridedMaybeAdjOrTransMat{<:AbstractMutable},
        B::StridedMaybeAdjOrTransMat,
        γ,
    )
        return _mat_mat_scalar(A, B, γ)
    end
    function LinearAlgebra.mat_mat_scalar(
        A::StridedMaybeAdjOrTransMat,
        B::StridedMaybeAdjOrTransMat{<:AbstractMutable},
        γ,
    )
        return _mat_mat_scalar(A, B, γ)
    end
    function LinearAlgebra.mat_mat_scalar(
        A::StridedMaybeAdjOrTransMat{<:AbstractMutable},
        B::StridedMaybeAdjOrTransMat{<:AbstractMutable},
        γ,
    )
        return _mat_mat_scalar(A, B, γ)
    end
end

# Base doesn't define efficient fallbacks for sparse array arithmetic involving
# non-`<:Number` scalar elements, so we define some of these for
# `<:AbstractMutable` scalar elements here.

for (S, T) in [
    (LinearAlgebra.UniformScaling, AbstractMutable),
    (Number, AbstractMutable),
    (AbstractMutable, Any),
]
    @eval begin
        function Base.:*(A::$S, B::_SparseMat{<:$T})
            return _SparseMat(
                B.m,
                B.n,
                copy(B.colptr),
                copy(SparseArrays.rowvals(B)),
                A .* SparseArrays.nonzeros(B),
            )
        end
        function Base.:*(B::_SparseMat{<:$T}, A::$S)
            return _SparseMat(
                B.m,
                B.n,
                copy(B.colptr),
                copy(SparseArrays.rowvals(B)),
                SparseArrays.nonzeros(B) .* A,
            )
        end
    end
end

function Base.:/(
    A::_SparseMat{<:AbstractMutable},
    B::LinearAlgebra.UniformScaling,
)
    return _SparseMat(
        A.m,
        A.n,
        copy(A.colptr),
        copy(SparseArrays.rowvals(A)),
        SparseArrays.nonzeros(A) ./ B,
    )
end

function Base.:/(A::_SparseMat{<:AbstractMutable}, B::Number)
    return _SparseMat(
        A.m,
        A.n,
        copy(A.colptr),
        copy(SparseArrays.rowvals(A)),
        SparseArrays.nonzeros(A) ./ B,
    )
end

# Base assumes that the element type is unaffected by `-`
function Base.:-(A::_SparseMat{<:AbstractMutable})
    return _SparseMat(
        A.m,
        A.n,
        copy(A.colptr),
        copy(SparseArrays.rowvals(A)),
        -SparseArrays.nonzeros(A),
    )
end

# Matrix(::SparseMatrixCSC) assumes that `zero` does not affect the element
# type of `S`.
function Base.Matrix(S::_SparseMat{T}) where {T<:AbstractMutable}
    U = promote_operation(+, promote_operation(zero, T), T)
    A = Matrix{U}(undef, size(S)...)
    operate!(zero, A)
    return operate!(+, A, S)
end

# +(::SparseMatrixCSC) is not defined for generic types in Base.
Base.:+(A::AbstractArray{<:AbstractMutable}) = A

# `Base.*(::AbstractArray, α)` is only defined if `α isa Number`
# Currently, mutable types are scalar elements (e.g. JuMP expression,
# MOI functions or polynomials) so broadcasting is the right dispatch.
# If this causes issues in the future, e.g., because a user define a non-scalar
# subtype of `AbstractMutable`, we might want to check that
# `ndims` is zero and error otherwise.

Base.:*(α::AbstractMutable, A::AbstractArray) = α .* A

Base.:*(A::AbstractArray, α::AbstractMutable) = A .* α

function operate_to!(
    output::AbstractArray,
    ::typeof(*),
    v::AbstractArray,
    α::Union{Number,AbstractMutable},
)
    return Base.broadcast!(*, output, v, α)
end

function operate_to!(
    output::AbstractArray,
    ::typeof(*),
    α::Union{Number,AbstractMutable},
    v::AbstractArray,
)
    return Base.broadcast!(*, output, α, v)
end

# Needed for Julia v1.0, otherwise, `broadcast(*, α, A)` gives a `Array` and
# not a `Symmetric`.

function _mult_triangle(
    ::Type{T},
    x,
    A::T,
) where {T<:Union{LinearAlgebra.Symmetric,LinearAlgebra.Hermitian}}
    c = LinearAlgebra.sym_uplo(A.uplo)
    B = if c == :U
        parent(x * LinearAlgebra.UpperTriangular(parent(A)))
    else
        parent(x * LinearAlgebra.LowerTriangular(parent(A)))
    end
    # Intermediate conversion to `Matrix` is needed to work around
    # https://github.com/JuliaLang/julia/issues/52895
    return T(Matrix(T(B, c)), c)
end

function Base.:*(α::Number, A::LinearAlgebra.Symmetric{<:AbstractMutable})
    return _mult_triangle(LinearAlgebra.Symmetric, α, A)
end

Base.:*(A::LinearAlgebra.Symmetric{<:AbstractMutable}, α::Number) = α * A

function Base.:*(
    α::AbstractMutable,
    A::LinearAlgebra.Symmetric{<:AbstractMutable},
)
    return _mult_triangle(LinearAlgebra.Symmetric, α, A)
end

function Base.:*(
    A::LinearAlgebra.Symmetric{<:AbstractMutable},
    α::AbstractMutable,
)
    return α * A
end

function Base.:*(α::AbstractMutable, A::LinearAlgebra.Symmetric)
    return _mult_triangle(LinearAlgebra.Symmetric, α, A)
end

Base.:*(A::LinearAlgebra.Symmetric, α::AbstractMutable) = α * A

function Base.:*(α::Real, A::LinearAlgebra.Hermitian{<:AbstractMutable})
    return _mult_triangle(LinearAlgebra.Hermitian, α, A)
end

function Base.:*(A::LinearAlgebra.Hermitian{<:AbstractMutable}, α::Real)
    return α * A
end

# These three have specific methods that just redirect to `Matrix{T}` which
# does not work, e.g. if `zero(T)` has a different type than `T`.

function Base.Matrix(x::LinearAlgebra.Tridiagonal{T}) where {T<:AbstractMutable}
    return Matrix{promote_type(promote_operation(zero, T), T)}(x)
end

function Base.Matrix(
    x::LinearAlgebra.UpperTriangular{T},
) where {T<:AbstractMutable}
    return Matrix{promote_type(promote_operation(zero, T), T)}(x)
end

function Base.Matrix(
    x::LinearAlgebra.LowerTriangular{T},
) where {T<:AbstractMutable}
    return Matrix{promote_type(promote_operation(zero, T), T)}(x)
end

# Needed for Julia v1.1 only. If `parent(A)` is for instance `Diagonal`, the
# `eltype` of `B` might be different form the `eltype` of `A`.
function Base.Matrix(A::LinearAlgebra.Symmetric{<:AbstractMutable})
    B = LinearAlgebra.copytri!(convert(Matrix, copy(A.data)), A.uplo)
    for i in axes(A, 1)
        # `B[i, i]` is used instead of `A[i, i]` on Julia v1.1 hence the need
        # to overwrite it for `AbstractMutable`.
        B[i, i] = LinearAlgebra.symmetric(
            A[i, i],
            LinearAlgebra.sym_uplo(A.uplo),
        )::LinearAlgebra.symmetric_type(eltype(A.data))
    end
    return B
end

function Base.Matrix(A::LinearAlgebra.Hermitian{<:AbstractMutable})
    B = LinearAlgebra.copytri!(convert(Matrix, copy(A.data)), A.uplo, true)
    for i in axes(A, 1)
        # `B[i, i]` is used instead of `A[i, i]` on Julia v1.1 hence the need
        # to overwrite it for `AbstractMutable`.
        B[i, i] = LinearAlgebra.hermitian(
            A[i, i],
            LinearAlgebra.sym_uplo(A.uplo),
        )::LinearAlgebra.hermitian_type(eltype(A.data))
    end
    return B
end

# Called in `getindex` of `LinearAlgebra.LowerTriangular` and
# `LinearAlgebra.UpperTriangular` as the elements may be `Array` for which
# `zero` is only defined for instances but not for the type. For
# `AbstractMutable` we assume that `zero` for the instance is the same than for
# the type by default.
Base.zero(x::AbstractMutable) = zero(typeof(x))

# This was fixed in https://github.com/JuliaLang/julia/pull/36194 but then
# reverted. Fixed again in https://github.com/JuliaLang/julia/pull/38789/.
if VERSION >= v"1.7.0-DEV.872"
    # `AbstractMutable` objects are more likely to implement `iszero` than `==`
    # with `Int`.
    LinearAlgebra.iszerodefined(::Type{<:AbstractMutable}) = true
else
    # To determine whether the funtion is zero preserving, `LinearAlgebra` calls
    # `zero` on the `eltype` of the broadcasted object and then check `_iszero`.
    # `_iszero(x)` redirects to `iszero(x)` for numbers and to `x == 0`
    # otherwise.
    # `x == 0` returns false for types that implement `iszero` but not `==` such
    # as `DummyBigInt` and MOI functions.
    LinearAlgebra._iszero(x::AbstractMutable) = iszero(x)
end
