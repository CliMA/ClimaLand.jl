# Consider permutedims as an alternative to direct multiplication.
# Multiplication is an O(m*N) cost compared to an O(N) cost for tridiagonal algorithms.
# However, when the multiplication is only a small fraction of the total time
# (for example, when m is small), then these can be convenient and avoid the need for calls to permutedims.

# Multiplication using permutedims
"""
`A_mul_B_perm!(dest, M, src, dim)` computes `M*x` for slices `x` of `src` along dimension `dim`,
storing the result in `dest`. `M` must be an `AbstractMatrix`. This uses `permutedims` to make dimension
`dim` into the first dimension, performs a standard matrix multiplication, and restores the original
dimension ordering. In many cases, this algorithm exhibits the best cache behavior.
"""
function A_mul_B_perm!(dest, M::AbstractMatrix, src, dim::Integer)
    check_matmul_sizes(dest, M, src, dim)
    order = [dim; setdiff(1:ndims(src), dim)]
    srcp = permutedims(src, order)
    tmp = Array{eltype(dest), 2}(undef, size(dest, dim), div(length(dest), size(dest, dim)))
    mul!(tmp, M, reshape(srcp, (size(src,dim), div(length(srcp), size(src,dim)))))
    iorder = [2:dim; 1; dim+1:ndims(src)]
    permutedims!(dest, reshape(tmp, size(dest)[order]), iorder)
    dest
end

# Direct (temporary-free) multiplication
"""
`A_mul_B_md!(dest, M, src, dim)` computes `M*x` for slices `x` of `src` along dimension `dim`,
storing the result in `dest`. `M` must be an `AbstractMatrix`. This uses an in-place naive algorithm.
"""
function A_mul_B_md!(dest, M::AbstractMatrix, src, dim::Integer)
    check_matmul_sizes(dest, M, src, dim)
    if size(M,1) == size(M,2) == 1
        return mul!(dest, src, M[1,1])
    end
    R2 = CartesianIndices(size(dest)[dim+1:end])
    if dim > 1
        R1 = CartesianIndices(size(dest)[1:dim-1])
        _A_mul_B_md!(dest, M, src, R1, R2)
    else
        _A_mul_B_md!(dest, M, src, R2)
    end
end

# Multiplication along the first dimension
# Here we expect that M will typically be small and fit in cache, whereas src and dest do not
function _A_mul_B_md!(dest, M::AbstractMatrix, src, R2::CartesianIndices)
    m, n = size(M, 1), size(M, 2)
    if m == n == 2
        return _A_mul_B_md_2x2!(dest, M, src, R2)
    end
    for I2 in R2
        @inbounds for i = 1:m
            dest[i,I2] = zero(eltype(dest))
        end
        @inbounds for j = 1:n
            b = src[j,I2]
            @simd for i = 1:m
                dest[i,I2] += M[i,j]*b
            end
        end
    end
    dest
end
_A_mul_B_md(M::AbstractMatrix, src, R2::CartesianIndices) = _A_mul_B_md!(alloc_matmul(M, src, 1), M, src, R2)

function _A_mul_B_md_2x2!(dest, M::AbstractMatrix, src, R2::CartesianIndices)
    a, b, c, d = M[1,1], M[1,2], M[2,1], M[2,2]
    @simd for I2 in R2
        @inbounds begin
            s1, s2 = src[1,I2], src[2,I2]
            dest[1,I2] = a*s1 + b*s2
            dest[2,I2] = c*s1 + d*s2
        end
    end
    dest
end

function _A_mul_B_md!(dest, M::SparseMatrixCSC, src, R2::CartesianIndices)
    m, n = size(M,1), size(M,2)
    nzv = M.nzval
    rv = M.rowval
    cp = M.colptr
    for I2 in R2
        @inbounds for i = 1:m
            dest[i,I2] = zero(eltype(dest))
        end
        for j = 1:n
            b = src[j,I2]
            @inbounds for k = cp[j]:(cp[j+1]-1)
                dest[rv[k],I2] += nzv[k]*b
            end
        end
    end
    dest
end

# Multiplication along any other dimension
function _A_mul_B_md!(dest, M::AbstractMatrix, src,  R1::CartesianIndices, R2::CartesianIndices)
    m, n = size(M, 1), size(M, 2)
    if m == n == 2
        return _A_mul_B_md_2x2!(dest, M, src, R1, R2)
    end
    fill!(dest, zero(eltype(dest)))
    for I2 in R2
        for j = 1:n
            @inbounds for i = 1:m
                Mij = M[i,j]
                @simd for I1 in R1
                    dest[I1,i,I2] += Mij*src[I1,j,I2]
                end
            end
        end
    end
    dest
end
_A_mul_B_md(M::AbstractMatrix, src,  R1::CartesianIndices, R2::CartesianIndices) = _A_mul_B_md!(alloc_matmul(M, src, ndims(R1)+1), M, src, R1, R2)

function _A_mul_B_md_2x2!(dest, M::AbstractMatrix, src,  R1::CartesianIndices, R2::CartesianIndices)
    a, b, c, d = M[1,1], M[1,2], M[2,1], M[2,2]
    for I2 in R2
        @simd for I1 in R1
            @inbounds begin
                s1, s2 = src[I1,1,I2], src[I1,2,I2]
                dest[I1,1,I2] = a*s1 + b*s2
                dest[I1,2,I2] = c*s1 + d*s2
            end
        end
    end
    dest
end

function _A_mul_B_md!(dest, M::SparseMatrixCSC, src, R1::CartesianIndices, R2::CartesianIndices)
    m, n = size(M,1), size(M,2)
    nzv = M.nzval
    rv = M.rowval
    cp = M.colptr
    fill!(dest, zero(eltype(dest)))
    for I2 in R2
        for j = 1:n
            @inbounds for k = cp[j]:(cp[j+1]-1)
                i, Mij = rv[k], nzv[k]
                @simd for I1 in R1
                    dest[I1,i,I2] += Mij*src[I1,j,I2]
                end
            end
        end
    end
    dest
end

function check_matmul_sizes(dest, M::AbstractMatrix, src, dim)
    1 <= dim <= max(ndims(dest),ndims(src)) || throw(DimensionMismatch("The chosen dimension $dim is larger than $(ndims(src)) and $(ndims(dest))"))
    m, n = size(M, 1), size(M, 2)
    n == size(src, dim) && m == size(dest, dim) || throw(DimensionMismatch("Sizes $m, $n, $(size(src,dim)), and $(size(dest,dim)) do not match"))
    for i = 1:max(ndims(src), ndims(dest))
        i == dim && continue
        if size(src,i) != size(dest,i)
            throw(DimensionMismatch("Sizes $(size(dest)), $(size(src)) do not match"))
        end
    end
    nothing
end
