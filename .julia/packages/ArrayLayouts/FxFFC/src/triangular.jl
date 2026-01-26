colsupport(::TriangularLayout{'L'}, A, j) = isempty(j) ? (1:0) : colsupport(triangulardata(A), j) ∩ (minimum(j):size(A,1))
colsupport(::TriangularLayout{'U'}, A, j) = isempty(j) ? (1:0) : colsupport(triangulardata(A), j) ∩ oneto(maximum(j))
rowsupport(::TriangularLayout{'U'}, A, j) = isempty(j) ? (1:0) : rowsupport(triangulardata(A), j) ∩ (minimum(j):size(A,2))
rowsupport(::TriangularLayout{'L'}, A, j) = isempty(j) ? (1:0) : rowsupport(triangulardata(A), j) ∩ oneto(maximum(j))



###
# Lmul
###

mulreduce(M::Mul{<:TriangularLayout}) = Lmul(M)
mulreduce(M::Mul{<:TriangularLayout,<:TriangularLayout}) = Lmul(M)
mulreduce(M::Mul{<:Any,<:TriangularLayout}) = Rmul(M)
mulreduce(M::Mul{<:DiagonalLayout,<:TriangularLayout}) = Lmul(M)
mulreduce(M::Mul{<:TriangularLayout,<:DiagonalLayout}) = Rmul(M)

similar(M::Lmul{<:TriangularLayout{'U'},<:TriangularLayout{'U'}}) = UpperTriangular(Matrix{eltype(M)}(undef, size(M)))
similar(M::Lmul{<:TriangularLayout{'L'},<:TriangularLayout{'L'}}) = LowerTriangular(Matrix{eltype(M)}(undef, size(M)))
similar(M::Lmul{<:TriangularLayout{'U','U'},<:TriangularLayout{'U','U'}}) = UnitUpperTriangular(Matrix{eltype(M)}(undef, size(M)))
similar(M::Lmul{<:TriangularLayout{'L','U'},<:TriangularLayout{'L','U'}}) = UnitLowerTriangular(Matrix{eltype(M)}(undef, size(M)))


## Generic triangular multiplication
function materialize!(M::Lmul{<:TriangularLayout{'U','N'}})
    A,B = M.A,M.B
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch(LazyString("right hand side B needs first dimension of size ", size(A,1), ", has size ", m)))
    end
    Adata = triangulardata(A)
    for j = rowsupport(B)
        cs = colsupport(B,j)
        for i = cs
            Bij = Adata[i,i]*B[i,j]
            for k = (i + 1:m) ∩ cs ∩ rowsupport(Adata,i)
                Bij += Adata[i,k]*B[k,j]
            end
            B[i,j] = Bij
        end
    end
    B
end

function materialize!(M::Lmul{<:TriangularLayout{'U','U'}})
    A,B = M.A,M.B
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch(LazyString("right hand side B needs first dimension of size ", size(A,1), ", has size ", m)))
    end
    Adata = triangulardata(A)
    for j = rowsupport(B)
        cs = colsupport(B,j)
        for i = cs
            Bij = B[i,j]
            for k = (i + 1:m) ∩ cs ∩ rowsupport(Adata,i)
                Bij += Adata[i,k]*B[k,j]
            end
            B[i,j] = Bij
        end
    end
    B
end

function materialize!(M::Lmul{<:TriangularLayout{'L','N'}})
    A,B = M.A,M.B
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch(LazyString("right hand side B needs first dimension of size ", size(A,1), ", has size ", m)))
    end
    Adata = triangulardata(A)
    for j = 1:n
        for i = reverse(colsupport(A,colsupport(B,j)))
            Bij = Adata[i,i]*B[i,j]
            for k = (1:i - 1) ∩ rowsupport(Adata,i)
                Bij += Adata[i,k]*B[k,j]
            end
            B[i,j] = Bij
        end
    end
    B
end
function materialize!(M::Lmul{<:TriangularLayout{'L','U'}})
    A,B = M.A,M.B
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch(LazyString("right hand side B needs first dimension of size ", size(A,1), ", has size ", m)))
    end
    Adata = triangulardata(A)
    for j = 1:n
        for i = m:-1:1
            Bij = B[i,j]
            for k = (1:i - 1) ∩ rowsupport(Adata,i)
                Bij += Adata[i,k]*B[k,j]
            end
            B[i,j] = Bij
        end
    end
    B
end

@inline function materialize!(M::BlasMatLmulVec{<:TriangularLayout{UPLO,UNIT,<:AbstractColumnMajor},
                                         <:AbstractStridedLayout}) where {UPLO,UNIT}
    A,x = M.A,M.B
    BLAS.trmv!(UPLO, 'N', UNIT, triangulardata(A), x)
end

@inline function materialize!(M::BlasMatLmulVec{<:TriangularLayout{'U',UNIT,<:AbstractRowMajor},
                                   <:AbstractStridedLayout}) where UNIT
    A,x = M.A,M.B
    BLAS.trmv!('L', 'T', UNIT, transpose(triangulardata(A)), x)
end

@inline function materialize!(M::BlasMatLmulVec{<:TriangularLayout{'L',UNIT,<:AbstractRowMajor},
                                   <:AbstractStridedLayout}) where UNIT
    A,x = M.A,M.B
    BLAS.trmv!('U', 'T', UNIT, transpose(triangulardata(A)), x)
end

@inline function materialize!(M::BlasMatLmulVec{<:TriangularLayout{'U',UNIT,<:ConjLayout{<:AbstractRowMajor}},
                                   <:AbstractStridedLayout,<:BlasComplex}) where UNIT
    A,x = M.A,M.B
    BLAS.trmv!('L', 'C', UNIT, triangulardata(A)', x)
end

@inline function materialize!(M::BlasMatLmulVec{<:TriangularLayout{'L',UNIT,<:ConjLayout{<:AbstractRowMajor}},
                                   <:AbstractStridedLayout,<:BlasComplex}) where UNIT
    A,x = M.A,M.B
    BLAS.trmv!('U', 'C', UNIT, triangulardata(A)', x)
end
# Triangular * Matrix

@inline function materialize!(M::BlasMatLmulMat{<:TriangularLayout{UPLO,UNIT,<:AbstractColumnMajor},
                                         <:AbstractStridedLayout, T}) where {UPLO,UNIT,T<:BlasFloat}
    A,x = M.A,M.B
    BLAS.trmm!('L', UPLO, 'N', UNIT, one(T), triangulardata(A), x)
end

@inline function materialize!(M::BlasMatLmulMat{<:TriangularLayout{'L',UNIT,<:AbstractRowMajor},
                                                <:AbstractStridedLayout, T}) where {UNIT,T<:BlasFloat}
    A,x = M.A,M.B
    BLAS.trmm!('L', 'U', 'T', UNIT, one(T), transpose(triangulardata(A)), x)
end

@inline function materialize!(M::BlasMatLmulMat{<:TriangularLayout{'U',UNIT,<:AbstractRowMajor},
                                                <:AbstractStridedLayout, T}) where {UNIT,T<:BlasFloat}
    A,x = M.A,M.B
    BLAS.trmm!('L', 'L', 'T', UNIT, one(T), transpose(triangulardata(A)), x)
end

@inline function materialize!(M::BlasMatLmulMat{<:TriangularLayout{'L',UNIT,<:ConjLayout{<:AbstractRowMajor}},
                                   <:AbstractStridedLayout, T}) where {UNIT,T<:BlasComplex}
    A,x = M.A,M.B
    BLAS.trmm!('L', 'U', 'C', UNIT, one(T), triangulardata(A)', x)
end

@inline function materialize!(M::BlasMatLmulMat{<:TriangularLayout{'U',UNIT,<:ConjLayout{<:AbstractRowMajor}},
                                                <:AbstractStridedLayout, T}) where {UNIT,T<:BlasComplex}
    A,x = M.A,M.B
    BLAS.trmm!('L', 'L', 'C', UNIT, one(T), triangulardata(A)', x)
end


###
# Rmul
###

@inline function materialize!(M::BlasMatRmulMat{<:AbstractStridedLayout,
                                                <:TriangularLayout{UPLO,UNIT,<:AbstractColumnMajor},T}) where {UPLO,UNIT,T<:BlasFloat}
    x,A = M.A,M.B
    BLAS.trmm!('R', UPLO, 'N', UNIT, one(T), triangulardata(A), x)
end

@inline function materialize!(M::BlasMatRmulMat{<:AbstractStridedLayout,
                                                <:TriangularLayout{'L',UNIT,<:AbstractRowMajor},T}) where {UNIT,T<:BlasFloat}
    x,A = M.A,M.B
    BLAS.trmm!('R', 'U', 'T', UNIT, one(T), transpose(triangulardata(A)), x)
end

@inline function materialize!(M::BlasMatRmulMat{<:AbstractStridedLayout,
                                                <:TriangularLayout{'U',UNIT,<:AbstractRowMajor},T}) where {UNIT,T<:BlasFloat}
x,A = M.A,M.B
BLAS.trmm!('R', 'L', 'T', UNIT, one(T), transpose(triangulardata(A)), x)
end

@inline function materialize!(M::BlasMatRmulMat{<:AbstractStridedLayout,
                                                <:TriangularLayout{'L',UNIT,<:ConjLayout{<:AbstractRowMajor}},T}) where {UNIT,T<:BlasComplex}
    x,A = M.A,M.B
    BLAS.trmm!('R', 'U', 'C', UNIT, one(T), triangulardata(A)', x)
end

@inline function materialize!(M::BlasMatRmulMat{<:AbstractStridedLayout,
                                                <:TriangularLayout{'U',UNIT,<:ConjLayout{<:AbstractRowMajor}},T}) where {UNIT,T<:BlasComplex}
x,A = M.A,M.B
BLAS.trmm!('R', 'L', 'C', UNIT, one(T), triangulardata(A)', x)
end


materialize!(M::MatRmulMat{<:AbstractStridedLayout,<:TriangularLayout}) = lmul!(M.B', M.A')'


########
# Ldiv
########


@inline function copyto!(dest::AbstractArray, M::Ldiv{<:Union{TriangularLayout,BidiagonalLayout,DiagonalLayout}})
    A, B = M.A, M.B
    dest ≡ B || copyto!(dest, B)
    ldiv!(A, dest)
end

for UNIT in ('U', 'N')
    for UPLO in ('L', 'U')
        @eval @inline materialize!(M::BlasMatLdivVec{<:TriangularLayout{$UPLO,$UNIT,<:AbstractColumnMajor},
                                            <:AbstractStridedLayout}) =
            BLAS.trsv!($UPLO, 'N', $UNIT, triangulardata(M.A), M.B)
        @eval @inline materialize!(M::BlasMatLdivMat{<:TriangularLayout{$UPLO,$UNIT,<:AbstractColumnMajor},
                                            <:AbstractStridedLayout}) =
            LAPACK.trtrs!($UPLO, 'N', $UNIT, triangulardata(M.A), M.B)
    end

    @eval begin
        @inline materialize!(M::BlasMatLdivVec{<:TriangularLayout{'U',$UNIT,<:AbstractRowMajor},
                                                        <:AbstractStridedLayout}) =
            BLAS.trsv!('L', 'T', $UNIT, transpose(triangulardata(M.A)), M.B)
        @inline materialize!(M::BlasMatLdivMat{<:TriangularLayout{'U',$UNIT,<:AbstractRowMajor},
                                                        <:AbstractColumnMajor}) =
            LAPACK.trtrs!('L', 'T', $UNIT, transpose(triangulardata(M.A)), M.B)

        @inline materialize!(M::BlasMatLdivVec{<:TriangularLayout{'L',$UNIT,<:AbstractRowMajor},
                                                        <:AbstractStridedLayout}) =
            BLAS.trsv!('U', 'T', $UNIT, transpose(triangulardata(M.A)), M.B)
        @inline materialize!(M::BlasMatLdivMat{<:TriangularLayout{'L',$UNIT,<:AbstractRowMajor},
                                                        <:AbstractColumnMajor}) =
            LAPACK.trtrs!('U', 'T', $UNIT, transpose(triangulardata(M.A)), M.B)


        @inline materialize!(M::BlasMatLdivVec{<:TriangularLayout{'U',$UNIT,<:ConjLayout{<:AbstractRowMajor}},
                                                        <:AbstractStridedLayout}) =
            BLAS.trsv!('L', 'C', $UNIT, triangulardata(M.A)', M.B)
        @inline materialize!(M::BlasMatLdivMat{<:TriangularLayout{'U',$UNIT,<:ConjLayout{<:AbstractRowMajor}},
                                                        <:AbstractColumnMajor}) =
            LAPACK.trtrs!('L', 'C', $UNIT, triangulardata(M.A)', M.B)

        @inline materialize!(M::BlasMatLdivVec{<:TriangularLayout{'L',$UNIT,<:ConjLayout{<:AbstractRowMajor}},
                                                        <:AbstractStridedLayout}) =
            BLAS.trsv!('U', 'C', $UNIT, triangulardata(M.A)', M.B)
        @inline materialize!(M::BlasMatLdivMat{<:TriangularLayout{'L',$UNIT,<:ConjLayout{<:AbstractRowMajor}},
                                                        <:AbstractColumnMajor}) =
            LAPACK.trtrs!('U', 'C', $UNIT, triangulardata(M.A)', M.B)
    end
end

function materialize!(M::MatLdivMat{<:Union{TriangularLayout,BidiagonalLayout}})
    A,X = M.A,M.B
    size(A,2) == size(X,1) || throw(DimensionMismatch("Dimensions must match"))
    @views for j in axes(X,2)
        ldiv!(A, X[:,j])
    end
    X
end

function materialize!(M::MatLdivVec{<:TriangularLayout{'U','N'}})
    A,b = M.A,M.B
    require_one_based_indexing(A, b)
    check_mul_axes(A, b)
    data = triangulardata(A)
    @inbounds for j in reverse(colsupport(b,1))
        iszero(data[j,j]) && throw(SingularException(j))
        bj = b[j] = data[j,j] \ b[j]
        for i in (1:j-1) ∩ colsupport(data,j)
            b[i] -= data[i,j] * bj
        end
    end
    b
end

function materialize!(M::MatLdivVec{<:TriangularLayout{'U','U'}})
    A,b = M.A,M.B
    require_one_based_indexing(A, b)
    check_mul_axes(A, b)
    data = triangulardata(A)
    @inbounds for j in reverse(colsupport(b,1))
        iszero(data[j,j]) && throw(SingularException(j))
        bj = b[j]
        for i in (1:j-1) ∩ colsupport(data,j)
            b[i] -= data[i,j] * bj
        end
    end
    b
end

function materialize!(M::MatLdivVec{<:TriangularLayout{'L','N'}})
    A,b = M.A,M.B
    require_one_based_indexing(A, b)
    check_mul_axes(A, b)
    data = triangulardata(A)
    n = size(A, 2)
    @inbounds for j in 1:n
        iszero(data[j,j]) && throw(SingularException(j))
        bj = b[j] = data[j,j] \ b[j]
        for i in (j+1:n) ∩ colsupport(data,j)
            b[i] -= data[i,j] * bj
        end
    end
    b
end

function materialize!(M::MatLdivVec{<:TriangularLayout{'L','U'}})
    A,b = M.A,M.B
    require_one_based_indexing(A, b)
    check_mul_axes(A, b)
    data = triangulardata(A)
    n = size(A, 2)
    @inbounds for j in 1:n
        iszero(data[j,j]) && throw(SingularException(j))
        bj = b[j]
        for i in (j+1:n) ∩ colsupport(data,j)
            b[i] -= data[i,j] * bj
        end
    end
    b
end

function _bidiag_backsub!(M)
    A,b = M.A, M.B
    N = last(colsupport(b,1))
    dv = diagonaldata(A)
    ev = supdiagonaldata(A)
    b[N] = bj1 = dv[N]\b[N]

    @inbounds for j = (N - 1):-1:1
        bj  = b[j]
        bj -= ev[j] * bj1
        dvj = dv[j]
        if iszero(dvj)
            throw(SingularException(j))
        end
        bj   = dvj\bj
        b[j] = bj1 = bj
    end

    b
end

function _bidiag_forwardsub!(M)
    A, b = M.A, M.B
    dv = diagonaldata(A)
    ev = subdiagonaldata(A)
    N = length(b)
    b[1] = bj1 = dv[1]\b[1]
    @inbounds for j = 2:N
        bj  = b[j]
        bj -= ev[j - 1] * bj1
        dvj = dv[j]
        if iszero(dvj)
            throw(SingularException(j))
        end
        bj   = dvj\bj
        b[j] = bj1 = bj
    end
    b
end

#Generic solver using naive substitution, based on LinearAlgebra/src/bidiag.jl
function materialize!(M::MatLdivVec{<:BidiagonalLayout})
    A,b = M.A,M.B
    require_one_based_indexing(A, b)
    N = size(A, 2)
    if N != length(b)
        throw(DimensionMismatch(LazyString("second dimension of A, ", N, ", does not match one of the length of b, ", length(b))))
    end

    if N == 0
        return b
    end

    if bidiagonaluplo(A) == 'L' #do forward substitution
        _bidiag_forwardsub!(M)
    else #do backward substitution
        _bidiag_backsub!(M)
    end
end
