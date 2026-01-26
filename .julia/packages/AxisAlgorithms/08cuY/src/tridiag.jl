_axes(F::Factorization, dim::Int) = hasmethod(axes, Tuple{typeof(F), Int}) ? axes(F, dim) : Base.OneTo(size(F, dim))
_axes(A::AbstractArray, dim::Int) = axes(A, dim)

"""
`A_ldiv_B_md!(dest, F, src, dim)` solves a tridiagonal system along dimension `dim` of `src`,
storing the result in `dest`. Currently, `F` must be an LU-factorized tridiagonal matrix.
If desired, you may safely use the same array for both `src` and `dest`, so that this becomes an
in-place algorithm.
"""
function A_ldiv_B_md!(dest, F, src, dim::Integer)
    1 <= dim <= max(ndims(dest),ndims(src)) || throw(DimensionMismatch("The chosen dimension $dim is larger than $(ndims(src)) and $(ndims(dest))"))
    ax = _axes(F, 1)
    ax == axes(src, dim) && ax == axes(dest, dim) || throw(DimensionMismatch("Axes $ax, $(axes(src,dim)), and $(axes(dest,dim)) do not match"))
    axes(dest) == axes(src) || throw(DimensionMismatch("Axes $(axes(dest)), $(axes(src)) do not match"))
    check_matrix(F)
    R1 = CartesianIndices(axes(dest)[1:dim-1])
    R2 = CartesianIndices(axes(dest)[dim+1:end])
    _A_ldiv_B_md!(dest, F, src, R1, R2)
end
_A_ldiv_B_md(F, src, R1::CartesianIndices, R2::CartesianIndices) =
    _A_ldiv_B_md!(similar(src, promote_type(eltype(F), eltype(src))), F, src, R1, R2)

# Solving along the first dimension
function _A_ldiv_B_md!(dest, F::LU{T,<:Tridiagonal{T}}, src,  R1::CartesianIndices{0}, R2::CartesianIndices) where {T}
    ax = _axes(F, 1)
    axbegin, axend = first(ax), last(ax)
    dl = F.factors.dl
    d  = F.factors.d
    du = F.factors.du
    # Forward substitution
    @inbounds for I2 in R2
        dest[axbegin, I2] = src[axbegin, I2]
        for i = axbegin+1:axend       # note: cannot use @simd here!
            dest[i, I2] = src[i, I2] - dl[i-1]*dest[i-1, I2]
        end
    end
    # Backward substitution
    dinv = 1 ./ d
    @inbounds for I2 in R2
        dest[axend, I2] /= d[axend]
        for i = axend-1:-1:axbegin  # note: cannot use @simd here!
            dest[i, I2] = (dest[i, I2] - du[i]*dest[i+1, I2])*dinv[i]
        end
    end
    dest
end

# Solving along any other dimension
function _A_ldiv_B_md!(dest, F::LU{T,<:Tridiagonal{T}}, src, R1::CartesianIndices, R2::CartesianIndices) where {T}
    ax = _axes(F, 1)
    axbegin, axend = first(ax), last(ax)
    dl = F.factors.dl
    d  = F.factors.d
    du = F.factors.du
    # Forward substitution
    @inbounds for I2 in R2
        @simd for I1 in R1
            dest[I1, axbegin, I2] = src[I1, axbegin, I2]
        end
        for i = axbegin+1:axend
            @simd for I1 in R1
                dest[I1, i, I2] = src[I1, i, I2] - dl[i-1]*dest[I1, i-1, I2]
            end
        end
    end
    # Backward substitution
    dinv = 1 ./ d
    for I2 in R2
        @simd for I1 in R1
            dest[I1, axend, I2] *= dinv[axend]
        end
        for i = axend-1:-1:axbegin
            @simd for I1 in R1
                dest[I1, i, I2] = (dest[I1, i, I2] - du[i]*dest[I1, i+1, I2])*dinv[i]
            end
        end
    end
    dest
end

function check_matrix(F::LU{T,<:Tridiagonal{T}}) where {T}
    ax = _axes(F, 1)
    for i in ax
        F.ipiv[i] == i || error("For efficiency, pivoting is not supported")
    end
    nothing
end
