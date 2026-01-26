###
# Lmul
####

mulreduce(M::Mul{<:DiagonalLayout,<:DiagonalLayout}) = Lmul(M)
mulreduce(M::Mul{<:DiagonalLayout}) = Lmul(M)
mulreduce(M::Mul{<:Any,<:DiagonalLayout}) = Rmul(M)

# Diagonal multiplication never changes structure
similar(M::Lmul{<:DiagonalLayout}, ::Type{T}, axes) where T = similar(M.B, T, axes)
# equivalent to rescaling
function materialize!(M::Lmul{<:DiagonalLayout{<:AbstractFillLayout}})
    M.B .= getindex_value(M.A.diag) .* M.B
    M.B
end


copy(M::Lmul{<:DiagonalLayout,<:DiagonalLayout}) = diagonal(diagonaldata(M.A) .* diagonaldata(M.B))
copy(M::Lmul{<:DiagonalLayout}) = diagonaldata(M.A) .* M.B
copy(M::Rmul{<:Any,<:DiagonalLayout}) = M.A .* permutedims(diagonaldata(M.B))


dualadjoint(_) = adjoint
dualadjoint(::Transpose) = transpose
dualadjoint(V::SubArray) = dualadjoint(parent(V))
function copy(M::Rmul{<:DualLayout,<:DiagonalLayout})
    adj = dualadjoint(M.A)
    adj(adj(M.B) * adj(M.A))
end



# Diagonal multiplication never changes structure
similar(M::Rmul{<:Any,<:DiagonalLayout}, ::Type{T}, axes) where T = similar(M.A, T, axes)
# equivalent to rescaling
function materialize!(M::Rmul{<:Any,<:DiagonalLayout{<:AbstractFillLayout}})
    M.A .= M.A .* getindex_value(M.B.diag)
    M.A
end


function materialize!(M::Ldiv{<:DiagonalLayout})
    M.B .= M.A.diag .\ M.B
    M.B
end

copy(M::Ldiv{<:DiagonalLayout,<:DiagonalLayout}) = diagonal(M.A.diag .\ M.B.diag)
copy(M::Ldiv{<:DiagonalLayout}) = M.A.diag .\ M.B
copy(M::Ldiv{<:DiagonalLayout{<:AbstractFillLayout}}) = inv(getindex_value(M.A.diag)) .* M.B
copy(M::Ldiv{<:DiagonalLayout{<:AbstractFillLayout},<:DiagonalLayout}) = diagonal(inv(getindex_value(M.A.diag)) .* M.B.diag)

copy(M::Rdiv{<:DiagonalLayout,<:DiagonalLayout}) = diagonal(M.A.diag .* inv.(M.B.diag))
copy(M::Rdiv{<:Any,<:DiagonalLayout}) = M.A .* inv.(permutedims(M.B.diag))
copy(M::Rdiv{<:Any,<:DiagonalLayout{<:AbstractFillLayout}}) = M.A .* inv(getindex_value(M.B.diag))
copy(M::Rdiv{<:DiagonalLayout,<:DiagonalLayout{<:AbstractFillLayout}}) = diagonal(M.A.diag .* inv(getindex_value(M.B.diag)))



## bi/tridiagonal copy
# hack around the fact that a SymTridiagonal isn't fully mutable
_similar(A) = similar(A)
_similar(A::SymTridiagonal) = similar(Tridiagonal(A.ev, A.dv, A.ev))
_copy_diag(M::T, ::T) where {T<:Rmul} = copyto!(_similar(M.A), M)
_copy_diag(M::T, ::T) where {T<:Lmul} = copyto!(_similar(M.B), M)
_copy_diag(M, _) = copy(M)
_bidiagonal(A::Bidiagonal) = A
function _bidiagonal(A)
    # we assume that the matrix is indeed bidiagonal,
    # so that the conversion is lossless
    if iszero(view(A, diagind(A, -1)))
        uplo = :U
    else
        uplo = :L
    end
    Bidiagonal(A, uplo)
end
function copy(M::Rmul{<:BidiagonalLayout,<:DiagonalLayout})
    A = _bidiagonal(M.A)
    _copy_diag(Rmul(A, M.B), M)
end
function copy(M::Lmul{<:DiagonalLayout,<:BidiagonalLayout})
    B = _bidiagonal(M.B)
    _copy_diag(Lmul(M.A, B), M)
end
# we assume that the matrix is indeed tridiagonal,
# so that the conversion is lossless
_tridiagonal(A::Tridiagonal) = A
_tridiagonal(A) = Tridiagonal(A)
function copy(M::Rmul{<:TridiagonalLayout,<:DiagonalLayout})
    A = _tridiagonal(M.A)
    _copy_diag(Rmul(A, M.B), M)
end
function copy(M::Lmul{<:DiagonalLayout,<:TridiagonalLayout})
    B = _tridiagonal(M.B)
    _copy_diag(Lmul(M.A, B), M)
end
# we assume that the matrix is indeed symmetric tridiagonal,
# so that the conversion is lossless
_symtridiagonal(A::SymTridiagonal) = A
_symtridiagonal(A) = SymTridiagonal(A)
function copy(M::Rmul{<:SymTridiagonalLayout,<:DiagonalLayout})
    A = _symtridiagonal(M.A)
    _copy_diag(Rmul(A, M.B), M)
end
function copy(M::Lmul{<:DiagonalLayout,<:SymTridiagonalLayout})
    B = _symtridiagonal(M.B)
    _copy_diag(Lmul(M.A, B), M)
end

copy(M::Lmul{DiagonalLayout{OnesLayout}}) = _copy_oftype(M.B, eltype(M))
copy(M::Lmul{DiagonalLayout{OnesLayout},<:DiagonalLayout}) = Diagonal(_copy_oftype(diagonaldata(M.B), eltype(M)))
copy(M::Lmul{<:DiagonalLayout,DiagonalLayout{OnesLayout}}) = Diagonal(_copy_oftype(diagonaldata(M.A), eltype(M)))
copy(M::Lmul{DiagonalLayout{OnesLayout},DiagonalLayout{OnesLayout}}) = _copy_oftype(M.B, eltype(M))
copy(M::Rmul{<:Any,DiagonalLayout{OnesLayout}}) = _copy_oftype(M.A, eltype(M))
copy(M::Rmul{<:DualLayout,DiagonalLayout{OnesLayout}}) = _copy_oftype(M.A, eltype(M))

copy(M::Lmul{<:DiagonalLayout{<:AbstractFillLayout}}) = getindex_value(diagonaldata(M.A)) * M.B
copy(M::Lmul{<:DiagonalLayout{<:AbstractFillLayout},<:DiagonalLayout}) = getindex_value(diagonaldata(M.A)) * M.B
copy(M::Rmul{<:Any,<:DiagonalLayout{<:AbstractFillLayout}}) = M.A * getindex_value(diagonaldata(M.B))
copy(M::Rmul{<:DualLayout,<:DiagonalLayout{<:AbstractFillLayout}}) = M.A * getindex_value(diagonaldata(M.B))

copy(M::Rmul{<:BidiagonalLayout,<:DiagonalLayout{<:AbstractFillLayout}}) = M.A * getindex_value(diagonaldata(M.B))
copy(M::Lmul{<:DiagonalLayout{<:AbstractFillLayout},<:BidiagonalLayout}) = getindex_value(diagonaldata(M.A)) * M.B
copy(M::Rmul{<:TridiagonalLayout,<:DiagonalLayout{<:AbstractFillLayout}}) = M.A * getindex_value(diagonaldata(M.B))
copy(M::Lmul{<:DiagonalLayout{<:AbstractFillLayout},<:TridiagonalLayout}) = getindex_value(diagonaldata(M.A)) * M.B
copy(M::Rmul{<:SymTridiagonalLayout,<:DiagonalLayout{<:AbstractFillLayout}}) = M.A * getindex_value(diagonaldata(M.B))
copy(M::Lmul{<:DiagonalLayout{<:AbstractFillLayout},<:SymTridiagonalLayout}) = getindex_value(diagonaldata(M.A)) * M.B


copy(M::Rmul{<:BidiagonalLayout,DiagonalLayout{OnesLayout}}) = _copy_oftype(M.A, eltype(M))
copy(M::Lmul{DiagonalLayout{OnesLayout},<:BidiagonalLayout}) =  _copy_oftype(M.B, eltype(M))
copy(M::Rmul{<:TridiagonalLayout,DiagonalLayout{OnesLayout}}) = _copy_oftype(M.A, eltype(M))
copy(M::Lmul{DiagonalLayout{OnesLayout},<:TridiagonalLayout}) =  _copy_oftype(M.B, eltype(M))
copy(M::Rmul{<:SymTridiagonalLayout,DiagonalLayout{OnesLayout}}) = _copy_oftype(M.A, eltype(M))
copy(M::Lmul{DiagonalLayout{OnesLayout},<:SymTridiagonalLayout}) =  _copy_oftype(M.B, eltype(M))
