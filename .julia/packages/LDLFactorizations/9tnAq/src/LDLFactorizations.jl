module LDLFactorizations

export ldl, ldl_analyze, ldl_factorize!, factorized, \, ldiv!, lmul!, mul!, nnz, -

using AMD, LinearAlgebra, SparseArrays

mutable struct SQDException <: Exception
  msg::String
end

const error_string = "LDL' factorization was not computed or failed"

"""
    col_symb!(n, Ap, Ai, Cp, w, Pinv)
Compute the sparse structure of missing elements of the upper triangle of PAPt. Nonzero elements have to verify Pinv[i] < Pinv[j] where i is the
row index and j the column index. Those elements are the nonzeros of the lower triangle of A that will be in the upper triangle of PAPt (after permutation)
# Arguments
- `n::Ti`: number of columns of the matrix
- `Ap::Vector{Ti}`: colptr of the matrix to factorize (CSC format)
- `Ai::Vector{Ti}`: rowval of the matrix to factorize (CSC format)
- `Cp::Vector{Ti}`: colptr of the lower triangle (to be modified)
- `w::Vector{Ti}`: work array
- `Pinv::Vector{Ti}`: inverse permutation of P. PAPt is the matrix to factorize (CSC format)
"""
function col_symb!(n, Ap, Ai, Cp, w, Pinv)
  fill!(w, 0)
  @inbounds for j = 1:n
    @inbounds for p = Ap[j]:(Ap[j + 1] - 1)
      i = Ai[p]
      i >= j && break  # only upper part
      Pinv[i] < Pinv[j] && continue  # store only what will be used during the factorization
      w[i] += 1  # count entries
    end
  end

  Cp[1] = 1
  @inbounds for i = 1:n  # cumulative sum
    Cp[i + 1] = w[i] + Cp[i]
    w[i] = Cp[i]
  end
end

"""
    col_num!(n, Ap, Ai, Ci, w, Pinv)
    
Compute the rowval and values of missing elements of the upper triangle of PAPt. Nonzero elements have to verify Pinv[i] ≥ Pinv[j] where i is the
row index and j the column index. Those elements are the nonzeros of the lower triangle of A that will be in the upper triangle of PAPt (after permutation)
# Arguments
- `n::Ti`: number of columns of the matrix
- `Ap::Vector{Ti}`: colptr of the matrix to factorize (CSC format)
- `Ai::Vector{Ti}`: rowval of the matrix to factorize (CSC format)
- `Ci::Vector{Ti}`: rowval of the lower triangle
- `w::Vector{Ti}`: work array
- `Pinv::Vector{Ti}`: inverse permutation of P. PAPt is the matrix to factorize (CSC format)
"""
function col_num!(n, Ap, Ai, Ci, w, Pinv)
  @inbounds for j = 1:n
    @inbounds for p = Ap[j]:(Ap[j + 1] - 1)
      i = Ai[p]
      i >= j && break  # only upper part
      Pinv[i] < Pinv[j] && continue  # store only what will be used during the factorization
      Ci[w[i]] = j
      w[i] += 1
    end
  end
end

function ldl_symbolic_upper!(n, Ap, Ai, Cp, Ci, Lp, parent, Lnz, flag, P, Pinv)
  @inbounds for k = 1:n
    parent[k] = -1
    flag[k] = k
    Lnz[k] = 0
    pk = P[k]
    @inbounds for p = Ap[pk]:(Ap[pk + 1] - 1)
      i = Pinv[Ai[p]]
      i ≥ k && continue
      @inbounds while flag[i] != k
        if parent[i] == -1
          parent[i] = k
        end
        Lnz[i] += 1
        flag[i] = k
        i = parent[i]
      end
    end

    # Missing nonzero elements of the upper triangle
    @inbounds for ind = Cp[pk]:(Cp[pk + 1] - 1)
      i = Pinv[Ci[ind]]
      i > k && continue
      @inbounds while flag[i] != k
        if parent[i] == -1
          parent[i] = k
        end
        Lnz[i] += 1
        flag[i] = k
        i = parent[i]
      end
    end
  end
  Lp[1] = 1
  @inbounds for k = 1:n
    Lp[k + 1] = Lp[k] + Lnz[k]
  end
end

function ldl_symbolic!(n, Ap, Ai, Lp, parent, Lnz, flag, P, Pinv)
  @inbounds for k = 1:n
    parent[k] = -1
    flag[k] = k
    Lnz[k] = 0
    pk = P[k]
    @inbounds for p = Ap[pk]:(Ap[pk + 1] - 1)
      i = Pinv[Ai[p]]
      i ≥ k && continue
      @inbounds while flag[i] != k
        if parent[i] == -1
          parent[i] = k
        end
        Lnz[i] += 1
        flag[i] = k
        i = parent[i]
      end
    end
  end
  Lp[1] = 1
  @inbounds for k = 1:n
    Lp[k + 1] = Lp[k] + Lnz[k]
  end
end

function ldl_numeric_upper!(
  n,
  Ap,
  Ai,
  Ax,
  Cp,
  Ci,
  Lp,
  parent,
  Lnz,
  Li,
  Lx,
  D,
  Y::T,
  pattern,
  flag,
  P,
  Pinv,
  r1,
  r2,
  tol,
  n_d,
) where {T}
  dynamic_reg = r1 != 0 || r2 != 0
  Y_zero = zero(eltype(T))
  @inbounds for k = 1:n
    Y[k] = Y_zero
    top = n + 1
    flag[k] = k
    Lnz[k] = 0
    pk = P[k]
    @inbounds for p = Ap[pk]:(Ap[pk + 1] - 1)
      i = Pinv[Ai[p]]
      i > k && continue
      Y[i] += Ax[p]
      len = 1
      @inbounds while flag[i] != k
        pattern[len] = i
        len += 1
        flag[i] = k
        i = parent[i]
      end
      @inbounds while len > 1
        top -= 1
        len -= 1
        pattern[top] = pattern[len]
      end
    end
    # missing non zero elements of the upper triangle
    @inbounds for ind = Cp[pk]:(Cp[pk + 1] - 1)
      i2 = Ci[ind]
      i = Pinv[i2]
      i > k && continue
      @inbounds for p = Ap[i2]:(Ap[i2 + 1] - 1)
        Ai[p] < pk && continue
        Y[i] += conj(Ax[p])
        len = 1
        @inbounds while flag[i] != k
          pattern[len] = i
          len += 1
          flag[i] = k
          i = parent[i]
        end
        @inbounds while len > 1
          top -= 1
          len -= 1
          pattern[top] = pattern[len]
        end
        break
      end
    end
    D[k] = Y[k]
    Y[k] = Y_zero
    @inbounds while top ≤ n
      i = pattern[top]
      yi = Y[i]
      Y[i] = Y_zero
      @inbounds for p = Lp[i]:(Lp[i] + Lnz[i] - 1)
        Y[Li[p]] -= Lx[p] * yi
      end
      p = Lp[i] + Lnz[i]
      l_ki = yi / D[i]
      D[k] -= conj(l_ki) * yi
      Li[p] = k
      Lx[p] = conj(l_ki)
      Lnz[i] += 1
      top += 1
    end
    if dynamic_reg && abs(D[k]) < tol
      r = P[k] <= n_d ? r1 : r2
      D[k] = sign(r) * max(abs(D[k] + r), abs(r))
    end
    D[k] == 0 && return false
  end
  return true
end

function ldl_numeric!(n, Ap, Ai, Ax, Lp, parent, Lnz, Li, Lx, D, Y, pattern, flag, P, Pinv)
  @inbounds for k = 1:n
    Y[k] = 0
    top = n + 1
    flag[k] = k
    Lnz[k] = 0
    pk = P[k]
    @inbounds for p = Ap[pk]:(Ap[pk + 1] - 1)
      i = Pinv[Ai[p]]
      i > k && continue
      Y[i] += Ax[p]
      len = 1
      @inbounds while flag[i] != k
        pattern[len] = i
        len += 1
        flag[i] = k
        i = parent[i]
      end
      @inbounds while len > 1
        top -= 1
        len -= 1
        pattern[top] = pattern[len]
      end
    end
    D[k] = Y[k]
    Y[k] = 0
    @inbounds while top ≤ n
      i = pattern[top]
      yi = Y[i]
      Y[i] = 0
      @inbounds for p = Lp[i]:(Lp[i] + Lnz[i] - 1)
        Y[Li[p]] -= Lx[p] * yi
      end
      p = Lp[i] + Lnz[i]
      l_ki = yi / D[i]
      D[k] -= conj(l_ki) * yi
      Li[p] = k
      Lx[p] = conj(l_ki)
      Lnz[i] += 1
      top += 1
    end
    D[k] == 0 && return false
  end
  return true
end

# solve functions for a single rhs
function ldl_lsolve!(n, x::AbstractVector, Lp, Li, Lx)
  @inbounds for j = 1:n
    xj = x[j]
    @simd for p = Lp[j]:(Lp[j + 1] - 1)
      x[Li[p]] -= Lx[p] * xj
    end
  end
  return x
end

function ldl_dsolve!(n, x::AbstractVector, D)
  x ./= D
  return x
end

function ldl_ltsolve!(n, x::AbstractVector, Lp, Li, Lx)
  @inbounds for j = n:-1:1
    xj = x[j]
    @simd for p = Lp[j]:(Lp[j + 1] - 1)
      xj -= conj(Lx[p]) * x[Li[p]]
    end
    x[j] = xj
  end
  return x
end

function ldl_solve!(n, b::AbstractVector, Lp, Li, Lx, D, P)
  @views y = b[P]
  ldl_lsolve!(n, y, Lp, Li, Lx)
  ldl_dsolve!(n, y, D)
  ldl_ltsolve!(n, y, Lp, Li, Lx)
  return b
end

# solve functions for multiple rhs
function ldl_lsolve!(n, X::AbstractMatrix, Lp, Li, Lx)
  @inbounds for j = 1:n
    @inbounds for p = Lp[j]:(Lp[j + 1] - 1)
      for k ∈ axes(X, 2)
        X[Li[p], k] -= Lx[p] * X[j, k]
      end
    end
  end
  return X
end

function ldl_dsolve!(n, X::AbstractMatrix, D)
  @inbounds for j = 1:n
    for k ∈ axes(X, 2)
      X[j, k] /= D[j]
    end
  end
  return X
end

function ldl_ltsolve!(n, X::AbstractMatrix, Lp, Li, Lx)
  @inbounds for j = n:-1:1
    @inbounds for p = Lp[j]:(Lp[j + 1] - 1)
      for k ∈ axes(X, 2)
        X[j, k] -= conj(Lx[p]) * X[Li[p], k]
      end
    end
  end
  return X
end

function ldl_solve!(n, B::AbstractMatrix, Lp, Li, Lx, D, P)
  @views Y = B[P, :]
  ldl_lsolve!(n, Y, Lp, Li, Lx)
  ldl_dsolve!(n, Y, D)
  ldl_ltsolve!(n, Y, Lp, Li, Lx)
  return B
end

# compute L*D*L'*x where x is a vector
function ldl_ltmul!(n, x::AbstractVector, Lp, Li, Lx)
  @inbounds for j = 1:n
    xj = x[j]
    @inbounds for p = Lp[j]:(Lp[j + 1] - 1)
      xj += conj(Lx[p]) * x[Li[p]]
    end
    x[j] = xj
  end
  return x
end

function ldl_dmul!(n, x::AbstractVector, D)
  @inbounds for j = 1:n
    x[j] *= D[j]
  end
  return x
end

function ldl_lmul!(n, x::AbstractVector, Lp, Li, Lx)
  @inbounds for j = n:-1:1
    xj = x[j]
    @inbounds for p = Lp[j]:(Lp[j + 1] - 1)
      x[Li[p]] += Lx[p] * xj
    end
  end
  return x
end

function ldl_mul!(n, x::AbstractVector, Lp, Li, Lx, D, P)
  @views y = x[P]
  ldl_ltmul!(n, y, Lp, Li, Lx)
  ldl_dmul!(n, y, D)
  ldl_lmul!(n, y, Lp, Li, Lx)
  return x
end

# compute L*D*L'*X where X is a matrix
function ldl_ltmul!(n, X::AbstractMatrix, Lp, Li, Lx)
  @inbounds for j = 1:n
    @inbounds for p = Lp[j]:(Lp[j + 1] - 1)
      for k ∈ axes(X, 2)
        X[j, k] += conj(Lx[p]) * X[Li[p], k]
      end
    end
  end
  return X
end

function ldl_dmul!(n, X::AbstractMatrix, D)
  @inbounds for j = 1:n
    @inbounds for k ∈ axes(X, 2)
      X[j, k] *= D[j]
    end
  end
  return X
end

function ldl_lmul!(n, X::AbstractMatrix, Lp, Li, Lx)
  @inbounds for j = n:-1:1
    @inbounds for p = Lp[j]:(Lp[j + 1] - 1)
      for k ∈ axes(X, 2)
        X[Li[p], k] += Lx[p] * X[j, k]
      end
    end
  end
  return X
end

function ldl_mul!(n, X::AbstractMatrix, Lp, Li, Lx, D, P)
  @views Y = X[P, :]
  ldl_ltmul!(n, Y, Lp, Li, Lx)
  ldl_dmul!(n, Y, D)
  ldl_lmul!(n, Y, Lp, Li, Lx)
  return X
end

"""
Type that contains the LDLᵀ factorization of a matrix.

The components of the factorization can be accessed via `getproperty`:

- `LDL.L`: `L` sparse lower triangular factor of the factorization without the diagonal 
    of ones that is removed to save space
- `LDL.D`: `D` diagonal matrix of the factorization.

In order to avoid zero pivots during the factorization, the user can regularize the matrix by modifying 
`LDL.r1` for the `LDL.n_d` first pivots and `LDL.r2` for the other pivots with tolerance `LDL.tol`.
"""
mutable struct LDLFactorization{T <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
  __analyzed::Bool
  __factorized::Bool
  __upper::Bool
  n::Tn
  # fields related to symbolic analysis
  parent::Vector{Ti}
  Lnz::Vector{Ti}
  flag::Vector{Ti}
  P::Vector{Tp}
  pinv::Vector{Tp}
  Lp::Vector{Ti}
  Cp::Vector{Ti}
  Ci::Vector{Ti}
  # fields related to numerical factorization
  Li::Vector{Ti}
  Lx::Vector{T}
  d::Vector{T}
  Y::Vector{T}
  pattern::Vector{Ti}
  # fields related to dynamic regularization
  r1::T
  r2::T
  tol::T
  n_d::Tn
end

"""
    isfact = factorized(LDL)

Returns true if the most recent factorization stored in `LDL` [`LDLFactorization`](@ref) succeeded.
"""
factorized(
  LDL::LDLFactorization{T, Ti, Tn, Tp},
) where {T <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer} = LDL.__factorized

"""
    LDL = ldl_analyze(A, Tf; P = amd(A))
    LDL = ldl_analyze(A; P = amd(A))
    LDL = ldl_analyze(A)

Perform symbolic analysis of the matrix `A` with permutation vector `P` (uses an
AMD permutation by default) so it can be reused.
`Tf` should be the element type of the factors, and is set to `eltype(A)` if not provided.
`A` should be a upper triangular matrix wrapped with LinearAlgebra's Symmetric / Hermitian type.

# Example
    A = sprand(Float64, 10, 10, 0.2)
    As = Symmetric(triu(A * A' + I), :U)
    LDL = ldl_analyze(As) # LDL in Float64
    LDL = ldl_analyze(As, Float32) # LDL in Float64
"""
function ldl_analyze end

"""
    ldl_factorize!(A, S)

Factorize A into the S [`LDLFactorization`](@ref) struct.
"""
function ldl_factorize! end

"""
    S = ldl(A, Tf; P = amd(A))
    S = ldl(A; P = amd(A))
    S = ldl(A)

Compute the LDLᵀ factorization of the matrix A with permutation vector P (uses an
AMD permutation by default).
`Tf` should be the element type of the factors, and is set to `eltype(A)` if not provided.
This function is equivalent to:

    S = ldl_analyze(A)
    ldl_factorize!(A, S)

A should either be a upper triangular matrix wrapped with LinearAlgebra's Symmetric / Hermitian
type, or a symmetric / hermitian matrix (not wrapped with Symmetric / Hermitian).

Using a non upper triangular matrix wrapped with Symmetric or Hermitian will not give the LDLᵀ factorization
of A.

# Example
    A = sprand(Float64, 10, 10, 0.2)
    As = Symmetric(triu(A * A' + I), :U)
    LDL = ldl(As) # LDL in Float64
    LDL = ldl(As, Float32) # LDL in Float64
"""
function ldl end

for (wrapper) in (:Symmetric, :Hermitian)
  @eval begin
    function ldl_analyze(
      A::$wrapper{T, SparseMatrixCSC{T, Ti}},
      ::Type{Tf};
      P::Vector{Tp} = amd(A),
    ) where {T <: Number, Ti <: Integer, Tp <: Integer, Tf <: Number}
      A.uplo == 'U' || error("upper triangle must be supplied")
      n = size(A, 1)
      n == size(A, 2) || throw(DimensionMismatch("matrix must be square"))
      n == length(P) || throw(DimensionMismatch("permutation size mismatch"))

      # allocate space for symbolic analysis
      parent = Vector{Ti}(undef, n)
      Lnz = Vector{Ti}(undef, n)
      flag = Vector{Ti}(undef, n)
      pinv = Vector{Tp}(undef, n)
      Lp = Vector{Ti}(undef, n + 1)

      # Compute inverse permutation
      @inbounds for k = 1:n
        pinv[P[k]] = k
      end

      Cp = Vector{Ti}(undef, n + 1)
      col_symb!(n, A.data.colptr, A.data.rowval, Cp, Lp, pinv)
      Ci = Vector{Ti}(undef, Cp[end] - 1)
      col_num!(n, A.data.colptr, A.data.rowval, Ci, Lp, pinv)

      # perform symbolic analysis
      ldl_symbolic_upper!(n, A.data.colptr, A.data.rowval, Cp, Ci, Lp, parent, Lnz, flag, P, pinv)

      # space for numerical factorization will be allocated later
      Li = Vector{Ti}(undef, Lp[n] - 1)
      Lx = Vector{Tf}(undef, Lp[n] - 1)
      d = Vector{Tf}(undef, n)
      Y = Vector{Tf}(undef, n)
      pattern = Vector{Ti}(undef, n)
      return LDLFactorization(
        true,
        false,
        true,
        n,
        parent,
        Lnz,
        flag,
        P,
        pinv,
        Lp,
        Cp,
        Ci,
        Li,
        Lx,
        d,
        Y,
        pattern,
        zero(Tf),
        zero(Tf),
        zero(Tf),
        n,
      )
    end

    ldl_analyze(
      A::$wrapper{T, SparseMatrixCSC{T, Ti}};
      kwargs...,
    ) where {T <: Number, Ti <: Integer} = ldl_analyze(A, T, kwargs...)

    # convert dense to sparse
    ldl_analyze(
      A::$wrapper{T, Matrix{T}},
      ::Type{Tf};
      kwargs...,
    ) where {T <: Number, Tf <: Number} = ldl_analyze($wrapper(sparse(A.data)), Tf; kwargs...)
    ldl_analyze(A::$wrapper{T, Matrix{T}}; kwargs...) where {T <: Number} =
      ldl_analyze($wrapper(sparse(A.data)); kwargs...)

    function ldl_factorize!(
      A::$wrapper{T, SparseMatrixCSC{T, Ti}},
      S::LDLFactorization{Tf, Ti, Tn, Tp},
    ) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
      S.__analyzed || error("perform symbolic analysis prior to numerical factorization")
      n = size(A, 1)
      n == S.n ||
        throw(DimensionMismatch("matrix size is inconsistent with symbolic analysis object"))

      # perform numerical factorization
      S.__factorized = ldl_numeric_upper!(
        S.n,
        A.data.colptr,
        A.data.rowval,
        A.data.nzval,
        S.Cp,
        S.Ci,
        S.Lp,
        S.parent,
        S.Lnz,
        S.Li,
        S.Lx,
        S.d,
        S.Y,
        S.pattern,
        S.flag,
        S.P,
        S.pinv,
        S.r1,
        S.r2,
        S.tol,
        S.n_d,
      )
      return S
    end

    # convert dense to sparse
    ldl_factorize!(A::$wrapper{T, Matrix{T}}, S::LDLFactorization) where {T <: Number} =
      ldl_factorize!($wrapper(sparse(A.data)), S)

    # symmetric or hermitian matrix input
    function ldl(
      sA::$wrapper{T, SparseMatrixCSC{T, Ti}},
      ::Type{Tf};
      P::Vector{Tp} = amd(sA),
    ) where {T <: Number, Ti <: Integer, Tp <: Integer, Tf <: Number}
      sA.uplo == 'L' && error("matrix must contain the upper triangle")
      # ldl(sA.data, args...; upper = true )
      S = ldl_analyze(sA, Tf; P = P)
      ldl_factorize!(sA, S)
    end

    ldl(sA::$wrapper{T, SparseMatrixCSC{T, Ti}}; kwargs...) where {T <: Number, Ti <: Integer} =
      ldl(sA, T; kwargs...)

    ldl(sA::$wrapper{T, Matrix{T}}, ::Type{Tf}; kwargs...) where {T <: Number, Tf <: Number} =
      ldl($wrapper(sparse(sA.data)), Tf; kwargs...)
    ldl(sA::$wrapper{T, Matrix{T}}; kwargs...) where {T <: Number} =
      ldl($wrapper(sparse(sA.data)), T; kwargs...)
  end
end

# use ldl(A, P = collect(1:n)) to suppress permutation
function ldl_analyze(
  A::SparseMatrixCSC{T, Ti},
  ::Type{Tf};
  P::Vector{Tp} = amd(A),
) where {T <: Number, Ti <: Integer, Tp <: Integer, Tf <: Number}
  n = size(A, 1)
  n == size(A, 2) || throw(DimensionMismatch("matrix must be square"))
  n == length(P) || throw(DimensionMismatch("permutation size mismatch"))

  # allocate space for symbolic analysis
  parent = Vector{Ti}(undef, n)
  Lnz = Vector{Ti}(undef, n)
  flag = Vector{Ti}(undef, n)
  pinv = Vector{Tp}(undef, n)
  Lp = Vector{Ti}(undef, n + 1)
  Cp = Ti[]
  Ci = Ti[]

  # Compute inverse permutation
  @inbounds for k = 1:n
    pinv[P[k]] = k
  end

  ldl_symbolic!(n, A.colptr, A.rowval, Lp, parent, Lnz, flag, P, pinv)

  # space for numerical factorization will be allocated later
  Li = Vector{Ti}(undef, Lp[n] - 1)
  Lx = Vector{Tf}(undef, Lp[n] - 1)
  d = Vector{Tf}(undef, n)
  Y = Vector{Tf}(undef, n)
  pattern = Vector{Ti}(undef, n)
  return LDLFactorization(
    true,
    false,
    false,
    n,
    parent,
    Lnz,
    flag,
    P,
    pinv,
    Lp,
    Cp,
    Ci,
    Li,
    Lx,
    d,
    Y,
    pattern,
    zero(Tf),
    zero(Tf),
    zero(Tf),
    n,
  )
end

ldl_analyze(A::SparseMatrixCSC{T}; kwargs...) where {T <: Number} = ldl_analyze(A, T; kwargs...)

# convert dense to sparse
ldl_analyze(A::Matrix{T}, ::Type{Tf}; kwargs...) where {T <: Number, Tf <: Number} =
  ldl_analyze(sparse(A), Tf; kwargs...)
ldl_analyze(A::Matrix{T}, kwargs...) where {T <: Number} = ldl_analyze(sparse(A), T, kwargs...)

function ldl_factorize!(
  A::SparseMatrixCSC{T, Ti},
  S::LDLFactorization{Tf, Ti, Tn, Tp},
) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
  S.__upper && error("symbolic analysis was performed for a Symmetric{} / Hermitian{} matrix")
  S.__analyzed || error("perform symbolic analysis prior to numerical factorization")
  n = size(A, 1)
  n == S.n || throw(DimensionMismatch("matrix size is inconsistent with symbolic analysis object"))

  S.__factorized = ldl_numeric!(
    S.n,
    A.colptr,
    A.rowval,
    A.nzval,
    S.Lp,
    S.parent,
    S.Lnz,
    S.Li,
    S.Lx,
    S.d,
    S.Y,
    S.pattern,
    S.flag,
    S.P,
    S.pinv,
  )
  return S
end

# convert dense to sparse
ldl_factorize!(A::Matrix{T}, S::LDLFactorization) where {T <: Number} = ldl_factorize!(sparse(A), S)

function ldl(A::SparseMatrixCSC, ::Type{Tf}; kwargs...) where {Tf <: Number}
  S = ldl_analyze(A, Tf; kwargs...)
  ldl_factorize!(A, S)
end

ldl(A::SparseMatrixCSC{T}; kwargs...) where {T <: Number} = ldl(sparse(A), T; kwargs...)

ldl(A::Matrix{T}, ::Type{Tf}; kwargs...) where {T <: Number, Tf <: Number} =
  ldl(sparse(A), Tf; kwargs...)
ldl(A::Matrix{T}; kwargs...) where {T <: Number} = ldl(sparse(A), T; kwargs...)

import Base.(\)

"""
    x = LDL \\ b 

If LDL is the LDLᵀ factorization of A, solves ``A x = b``.
"""
function (\)(
  LDL::LDLFactorization{Tf, Ti, Tn, Tp},
  b::AbstractVector{T},
) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
  y = copy(b)
  LDL.__factorized || throw(SQDException(error_string))
  ldl_solve!(LDL.n, y, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

function (\)(
  LDL::LDLFactorization{Tf, Ti, Tn, Tp},
  B::AbstractMatrix{T},
) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
  Y = copy(B)
  LDL.__factorized || throw(SQDException(error_string))
  ldl_solve!(LDL.n, Y, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

import LinearAlgebra.ldiv!

"""
    ldiv!(LDL, b) 

If LDL is the LDLᵀ factorization of A, solves `A x = b` and overwrites b with x.
"""
@inline ldiv!(
  LDL::LDLFactorization{Tf, Ti, Tn, Tp},
  b::AbstractVector{T},
) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer} =
  LDL.__factorized ? ldl_solve!(LDL.n, b, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P) :
  throw(SQDException(error_string))

@inline ldiv!(
  LDL::LDLFactorization{Tf, Ti, Tn, Tp},
  B::AbstractMatrix{T},
) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer} =
  LDL.__factorized ? ldl_solve!(LDL.n, B, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P) :
  throw(SQDException(error_string))

"""
    ldiv!(y, LDL, b) 

If LDL is the LDLᵀ factorization of A, solves `A x = b` In place.
"""
function ldiv!(
  y::AbstractVector{T},
  LDL::LDLFactorization{Tf, Ti, Tn, Tp},
  b::AbstractVector{T},
) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
  y .= b
  LDL.__factorized || throw(SQDException(error_string))
  ldl_solve!(LDL.n, y, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

function ldiv!(
  Y::AbstractMatrix{T},
  LDL::LDLFactorization{Tf, Ti, Tn, Tp},
  B::AbstractMatrix{T},
) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
  Y .= B
  LDL.__factorized || throw(SQDException(error_string))
  ldl_solve!(LDL.n, Y, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

import LinearAlgebra.lmul!, LinearAlgebra.mul!
function lmul!(
  LDL::LDLFactorization{Tf, Ti, Tn, Tp},
  x::AbstractVector{T},
) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
  LDL.__factorized || throw(SQDException(error_string))
  ldl_mul!(LDL.n, x, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

function lmul!(
  LDL::LDLFactorization{Tf, Ti, Tn, Tp},
  X::AbstractMatrix{T},
) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
  LDL.__factorized || throw(SQDException(error_string))
  ldl_mul!(LDL.n, X, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

function mul!(
  y::AbstractVector{T},
  LDL::LDLFactorization{Tf, Ti, Tn, Tp},
  x::AbstractVector{T},
) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
  y .= x
  LDL.__factorized || throw(SQDException(error_string))
  ldl_mul!(LDL.n, y, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

function mul!(
  Y::AbstractMatrix{T},
  LDL::LDLFactorization{Tf, Ti, Tn, Tp},
  X::AbstractMatrix{T},
) where {T <: Number, Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
  Y .= X
  LDL.__factorized || throw(SQDException(error_string))
  ldl_mul!(LDL.n, Y, LDL.Lp, LDL.Li, LDL.Lx, LDL.d, LDL.P)
end

Base.eltype(LDL::LDLFactorization) = eltype(LDL.d)
Base.size(LDL::LDLFactorization) = (LDL.n, LDL.n)
SparseArrays.nnz(LDL::LDLFactorization) = length(LDL.Lx) + length(LDL.d)

# user-friendly factors
@inline function Base.getproperty(LDL::LDLFactorization, prop::Symbol)
  if prop == :L
    # TODO: permute and return UnitLowerTriangular()
    return SparseMatrixCSC(LDL.n, LDL.n, LDL.Lp, LDL.Li, LDL.Lx)
  elseif prop == :D
    return Diagonal(LDL.d)
  else
    getfield(LDL, prop)
  end
end

Base.propertynames(LDL::LDLFactorization) = (:L, :D, :P)

import Base.(-)

"""
    -(LDL)

Unary minus operator returns an `LDLFactorization` with `-LDL.d`.
"""
function (-)(
  LDL::LDLFactorization{Tf, Ti, Tn, Tp},
) where {Tf <: Number, Ti <: Integer, Tn <: Integer, Tp <: Integer}
  LDL.d .= -LDL.d
  return LDL
end

end  # module
