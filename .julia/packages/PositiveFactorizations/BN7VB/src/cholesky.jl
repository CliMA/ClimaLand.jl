import Base: *, \, unsafe_getindex
using LinearAlgebra.BLAS: syr!, ger!, syrk!, syr2k!
using LinearAlgebra: BlasInt, BlasFloat, Cholesky, CholeskyPivoted
import LinearAlgebra: cholesky, cholesky!, ldlt, ldlt!

cholesky(::Type{Positive{T}}, A::AbstractMatrix, pivot=Val{false}; tol=default_tol(A), blocksize=default_blocksize(T)) where {T} = ldlt(Positive{T}, A, pivot; tol=tol, blocksize=blocksize)[1]
cholesky(::Type{Positive}, A::AbstractMatrix, pivot=Val{false}; tol=default_tol(A), blocksize=default_blocksize(floattype(eltype(A)))) = cholesky(Positive{floattype(eltype(A))}, A, pivot; tol=tol, blocksize=blocksize)

function ldlt(::Type{Positive{T}}, A::AbstractMatrix, pivot=Val{false}; tol=default_tol(A), blocksize=default_blocksize(T)) where {T}
    size(A, 1) == size(A, 2) || throw(DimensionMismatch("A must be square"))
    A0 = Array{floattype(T)}(undef, size(A))
    copyto!(A0, A)
    ldlt!(Positive{T}, A0, pivot; tol=tol, blocksize=blocksize)
end
ldlt(::Type{Positive}, A::AbstractMatrix, pivot=Val{false}; tol=default_tol(A), blocksize=default_blocksize(floattype(eltype(A)))) = ldlt(Positive{floattype(eltype(A))}, A, pivot; tol=tol, blocksize=blocksize)

cholesky!(::Type{Positive{T}}, A::AbstractMatrix{T}, pivot=Val{false}; tol=default_tol(A), blocksize=default_blocksize(T)) where {T<:AbstractFloat} = ldlt!(Positive{T}, A, pivot; tol=tol, blocksize=blocksize)[1]
cholesky!(::Type{Positive}, A::AbstractMatrix{T}, pivot=Val{false}; tol=default_tol(A), blocksize=default_blocksize(T)) where {T<:AbstractFloat} = cholesky!(Positive{T}, A; tol=tol, blocksize=blocksize)

# Blocked, cache-friendly algorithm (unpivoted)
function ldlt!(::Type{Positive{T}}, A::AbstractMatrix{T}, pivot::Type{Val{false}}=Val{false}; tol=default_tol(A), blocksize=default_blocksize(T)) where {T}
    size(A,1) == size(A,2) || error("A must be square")
    eltype(A)<:Real || error("element type $(eltype(A)) not yet supported")
    K = size(A, 1)
    d = Array{signtype(T)}(undef, K)
    for j = 1:blocksize:K
        # Split A into
        #            |
        #       B11  |
        #            |
        # A = ----------------
        #            |
        #       B21  |   B22
        #            |
        jend = min(K, j+blocksize-1)
        B11 = view(A, j:jend, j:jend)
        d1 = view(d, j:jend)
        solve_diagonal!(B11, d1, tol)
        if jend < K
            B21 = view(A, jend+1:K, j:jend)
            solve_columns!(B21, d1, B11)
            B22 = view(A, jend+1:K, jend+1:K)
            update_columns!(B22, d1, B21)
        end
    end
    @static if VERSION >= v"0.7.0-DEV.393"
        return Cholesky(A, :L, BLAS.BlasInt(0)), d
    else
        return Cholesky(A, :L), d
    end
end

# Version with pivoting
function ldlt!(::Type{Positive{T}}, A::AbstractMatrix{T}, pivot::Type{Val{true}}; tol=default_tol(A), blocksize=default_blocksize(T)) where {T}
    size(A,1) == size(A,2) || error("A must be square")
    eltype(A)<:Real || error("element type $(eltype(A)) not yet supported")
    K = size(A, 1)
    d = Array{signtype(T)}(undef, K)
    piv = convert(Vector{BlasInt}, 1:K)
    Ad = diag(A)
    for j = 1:blocksize:K
        jend = min(K, j+blocksize-1)
        solve_columns_pivot!(A, d, piv, Ad, tol, j:jend)
        if jend < K
            d1 = view(d, j:jend)
            B21 = view(A, jend+1:K, j:jend)
            B22 = view(A, jend+1:K, jend+1:K)
            update_columns!(B22, d1, B21)
        end
    end
    CholeskyPivoted(A, 'L', piv, BLAS.BlasInt(K), tol, BLAS.BlasInt(0)), d
end

ldlt!(::Type{Positive}, A::AbstractMatrix{T}, pivot=Val{false}; tol=default_tol(A), blocksize=default_blocksize(T)) where {T} = ldlt!(Positive{T}, A; tol=tol, blocksize=blocksize)


function solve_diagonal!(B, d, tol)
    K = size(B, 1)
    for j = 1:K
        Bjj = B[j,j]
        if abs(Bjj) > tol
            # compute ℓ (as the jth column of B)
            d[j] = sign(Bjj)
            s = sqrt(abs(Bjj))
            B[j,j] = s
            f = d[j]/s
            for i = j+1:K
                B[i,j] *= f
            end
            # subtract ℓ[j+1:end]⊗ℓ[j+1:end] from the lower right quadrant
            update_columns!(view(B, j+1:K, j+1:K), d[j], view(B, j+1:K, j))
        else
            # For the zero diagonals, replace them with 1. In a Newton step,
            # this corresponds to following the gradient (i.e., H = eye).
            d[j] = 0
            B[j,j] = 1
            for i = j+1:K
                B[i,j] = 0
            end
        end
    end
    B
end

function solve_columns!(B21, d, B11)
    I, J = size(B21)
    for j = 1:J
        dj = d[j]
        dj == 0 && continue
        s = B11[j,j]
        f = dj/s
        for i = 1:I
            B21[i,j] *= f
        end
        update_columns!(view(B21, :, j+1:J), dj, view(B21, :, j), view(B11, j+1:J, j))
    end
    B21
end

# Here, pivoting applies to the whole matrix, so we don't pass in a view.
# The jrange input describes the columns we're supposed to handle now.
function solve_columns_pivot!(A, d, piv, Ad, tol, jrange)
    K, KA = last(jrange), size(A, 1)
    jmin = first(jrange)
    for j in jrange
        # Find the remaining diagonal with largest magnitude
        Amax = zero(eltype(A))
        jmax = j-1
        for jj = j:KA
            tmp = abs(Ad[jj])
            if tmp > Amax
                Amax = tmp
                jmax = jj
            end
        end
        if jmax > j
            pivot!(A, j, jmax)
            Ad[j], Ad[jmax] = Ad[jmax], Ad[j]
            piv[j], piv[jmax] = piv[jmax], piv[j]
        end
        Ajj = A[j,j]
        for k = jmin:j-1
            tmp = A[j,k]
            Ajj -= d[k]*tmp*tmp
        end
        if abs(Ajj) > tol
            # compute ℓ (as the jth column of A)
            d[j] = sign(Ajj)
            s = sqrt(abs(Ajj))
            A[j,j] = s
            f = d[j]/s
        else
            d[j] = 0
            A[j,j] = 1
            f = zero(eltype(A))
        end
        for k = jmin:j-1
            @inbounds ck = d[k]*A[j,k]
            @simd for i = j+1:KA
                @inbounds A[i,j] -= ck*A[i,k]
            end
        end
        dj = d[j]
        @simd for i = j+1:KA
            @inbounds tmp = A[i,j]
            tmp *= f
            @inbounds A[i,j] = tmp
            @inbounds Ad[i] -= dj*tmp*tmp
        end
    end
    A
end

# Computes dest -= d*c*c', in the lower diagonal
@inline function update_columns!(dest::StridedMatrix{T}, d::Number, c::StridedVector{T}) where {T<:BlasFloat}
    syr!('L', convert(T, -d), c, dest)
    dest
end

# Computes dest -= d*x*y'
@inline function update_columns!(dest::StridedMatrix{T}, d::Number, x::StridedVector{T}, y::StridedVector{T}) where {T<:BlasFloat}
    ger!(convert(T, -d), x, y, dest)
    dest
end

# Computes dest -= C*diagm(d)*C', in the lower diagonal
function update_columns!(dest::StridedMatrix{T}, d::AbstractVector, C::StridedMatrix{T}) where {T<:BlasFloat}
    isempty(d) && return dest
    # If d is homogeneous, we can use syr rather than syr2
    allsame = true
    d1 = d[1]
    for i = 2:length(d)
        allsame &= (d[i] == d1)
    end
    allsame && d1 == 0 && return dest
    if allsame
        syrk!('L', 'N', convert(T, -d1), C, one(T), dest)
    else
        Cd = C*Diagonal(d)
        syr2k!('L', 'N', -one(T)/2, C, Cd, one(T), dest)
    end
    dest
end

# Pure-julia fallbacks for the above routines
# Computes dest -= d*c*c', in the lower diagonal
function update_columns!(dest, d::Number, c::AbstractVector)
    K = length(c)
    for j = 1:K
        dcj = d*c[j]
        @simd for i = j:K
            @inbounds dest[i,j] -= dcj*c[i]
        end
    end
    dest
end

# Computes dest -= d*x*y'
function update_columns!(dest, d::Number, x::AbstractVector, y::AbstractVector)
    I, J = size(dest)
    for j = 1:J
        dyj = d*y[j]
        @simd for i = 1:I
            @inbounds dest[i,j] -= dyj*x[i]
        end
    end
    dest
end

# Computes dest -= C*diagm(d)*C', in the lower diagonal
function update_columns!(dest, d::AbstractVector, C::AbstractMatrix)
    Ct = C'
    Cdt = Diagonal(d) * Ct
    K = size(dest, 1)
    nc = size(C, 2)
    for j = 1:K
        for i = j:K
            tmp = zero(eltype(dest))
            @simd for k = 1:nc
                @inbounds tmp += Ct[k,i]*Cdt[k,j]
            end
            @inbounds dest[i,j] -= tmp
        end
    end
    dest
end

# Diagonal pivoting (row&column swap) for a lower triangular matrix
function pivot!(A, i::Integer, j::Integer)
    i, j = min(i,j), max(i,j)
    for k = 1:i-1
        A[i,k], A[j,k] = A[j,k], A[i,k]
    end
    A[i,i], A[j,j] = A[j,j], A[i,i]
    for k = i+1:j-1
        A[k,i], A[j,k] = A[j,k], A[k,i]
    end
    for k = j+1:size(A,1)
        A[k,i], A[k,j] = A[k,j], A[k,i]
    end
    A
end

floattype(::Type{T}) where {T<:AbstractFloat} = T
floattype(::Type{T}) where {T<:Integer} = Float64
floattype(::Type{T}) where {T<:Real} = T

signtype(::Type{T}) where {T<:AbstractFloat} = Int8
signtype(::Type{T}) where {T<:Integer} = Int8
signtype(::Type{T}) where {T<:Real} = T

const cachesize = 2^15

default_δ(A) = 10 * size(A, 1) * eps(floattype(real(eltype(A))))
default_tol(A) = default_δ(A) * maximum(abs,A)
default_blocksize(::Type{T}) where {T} = max(4, floor(Int, sqrt(cachesize/sizeof(T)/4)))
