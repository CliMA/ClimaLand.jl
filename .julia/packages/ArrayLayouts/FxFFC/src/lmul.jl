for Typ in (:Lmul, :Rmul)
    @eval begin
        struct $Typ{StyleA, StyleB, TypeA, TypeB}
            A::TypeA
            B::TypeB
        end

        $Typ(A::TypeA, B::TypeB) where {TypeA,TypeB} = $Typ{typeof(MemoryLayout(TypeA)),typeof(MemoryLayout(TypeB)),TypeA,TypeB}(A,B)

        $Typ(M::Mul) = $Typ(M.A, M.B)

        BroadcastStyle(::Type{<:$Typ}) = ApplyBroadcastStyle()
        broadcastable(M::$Typ) = M

        eltype(::$Typ{<:Any,<:Any,A,B}) where {A,B} = promote_type(eltype(A), eltype(B))
        size(M::$Typ, p::Int) = size(M)[p]
        axes(M::$Typ, p::Int) = axes(M)[p]
        length(M::$Typ) = prod(size(M))
        size(M::$Typ) = map(length,axes(M))
        axes(M::$Typ) = (axes(M.A,1),axes(M.B,2))

        similar(M::$Typ, ::Type{T}, axes) where {T} = similar(Array{T}, axes)
        similar(M::$Typ, ::Type{T}) where T = similar(M, T, axes(M))
        similar(M::$Typ) = similar(M, eltype(M))
    end
end



const MatLmulVec{StyleA,StyleB} = Lmul{StyleA,StyleB,<:Union{AbstractMatrix,AbstractQ},<:AbstractVector}
const MatLmulMat{StyleA,StyleB} = Lmul{StyleA,StyleB,<:Union{AbstractMatrix,AbstractQ},<:Union{AbstractMatrix,AbstractQ}}

const BlasMatLmulVec{StyleA,StyleB,T<:BlasFloat} = Lmul{StyleA,StyleB,<:Union{AbstractMatrix{T},AbstractQ{T}},<:AbstractVector{T}}
const BlasMatLmulMat{StyleA,StyleB,T<:BlasFloat} = Lmul{StyleA,StyleB,<:Union{AbstractMatrix{T},AbstractQ{T}},<:Union{AbstractMatrix{T},AbstractQ{T}}}

const MatRmulMat{StyleA,StyleB} = Rmul{StyleA,StyleB,<:Union{AbstractMatrix,AbstractQ},<:Union{AbstractMatrix,AbstractQ}}
const BlasMatRmulMat{StyleA,StyleB,T<:BlasFloat} = Rmul{StyleA,StyleB,<:Union{AbstractMatrix{T},AbstractQ{T}},<:Union{AbstractMatrix{T},AbstractQ{T}}}


####
# LMul materialize
####

axes(M::MatLmulVec) = (axes(M.A,1),)

lmul(A, B; kwds...) = materialize(Lmul(A, B); kwds...)
lmul!(A, B; kwds...) = materialize!(Lmul(A, B); kwds...)
rmul!(A, B; kwds...) = materialize!(Rmul(A, B); kwds...)

materialize(L::Lmul) = copy(instantiate(L))

# needed since in orthogonal case dimensions might mismatch and
# so need to make sure extra entries are zero
function _zero_copyto!(dest, A)
    if axes(dest) == axes(A)
        copyto!(dest, A)
    else
        copyto!(zero!(dest), A)
    end
end

copy(M::Lmul) = lmul!(M.A, _zero_copyto!(similar(M), M.B))
copy(M::Rmul) = rmul!(_zero_copyto!(similar(M), M.A), M.B)

@inline function _lmul_copyto!(dest, M)
    M.B ≡ dest || copyto!(dest, M.B)
    lmul!(M.A,dest)
end

@inline function _rmul_copyto!(dest, M::Rmul)
    M.A ≡ dest || copyto!(dest, M.A)
    rmul!(dest,M.B)
end

copyto!(dest, M::Lmul) = _lmul_copyto!(dest, M)
copyto!(dest::AbstractArray, M::Lmul) = _lmul_copyto!(dest, M)
copyto!(dest, M::Rmul) = _rmul_copyto!(dest, M)
copyto!(dest::AbstractArray, M::Rmul) = _rmul_copyto!(dest, M)

materialize!(M::Lmul) = LinearAlgebra.lmul!(M.A,M.B)
materialize!(M::Rmul) = LinearAlgebra.rmul!(M.A,M.B)

materialize!(M::Lmul{ScalarLayout}) = invoke(LinearAlgebra.lmul!, Tuple{Number,AbstractArray}, M.A, M.B)
materialize!(M::Rmul{<:Any,ScalarLayout}) = invoke(LinearAlgebra.rmul!, Tuple{AbstractArray,Number}, M.A, M.B)

function materialize!(M::Lmul{ScalarLayout,<:SymmetricLayout})
    lmul!(M.A, symmetricdata(M.B))
    M.B
end
function materialize!(M::Lmul{ScalarLayout,<:HermitianLayout})
    lmul!(M.A, hermitiandata(M.B))
    M.B
end
function materialize!(M::Rmul{<:SymmetricLayout,ScalarLayout})
    rmul!(symmetricdata(M.A), M.B)
    M.A
end
function materialize!(M::Rmul{<:HermitianLayout,ScalarLayout})
    rmul!(hermitiandata(M.A), M.B)
    M.A
end

macro _layoutlmul(Typ)
    esc(quote
        LinearAlgebra.lmul!(A::$Typ, x::AbstractVector) = ArrayLayouts.lmul!(A,x)
        LinearAlgebra.lmul!(A::$Typ, x::AbstractMatrix) = ArrayLayouts.lmul!(A,x)
        LinearAlgebra.lmul!(A::$Typ, x::StridedVector) = ArrayLayouts.lmul!(A,x)
        LinearAlgebra.lmul!(A::$Typ, x::StridedMatrix) = ArrayLayouts.lmul!(A,x)
    end)
end

macro layoutlmul(Typ)
    esc(quote
        ArrayLayouts.@_layoutlmul UpperTriangular{T, <:$Typ{T}} where T
        ArrayLayouts.@_layoutlmul UnitUpperTriangular{T, <:$Typ{T}} where T
        ArrayLayouts.@_layoutlmul LowerTriangular{T, <:$Typ{T}} where T
        ArrayLayouts.@_layoutlmul UnitLowerTriangular{T, <:$Typ{T}} where T

        ArrayLayouts.@_layoutlmul UpperTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutlmul UnitUpperTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutlmul LowerTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutlmul UnitLowerTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T

        ArrayLayouts.@_layoutlmul UpperTriangular{T, <:Transpose{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutlmul UnitUpperTriangular{T, <:Transpose{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutlmul LowerTriangular{T, <:Transpose{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutlmul UnitLowerTriangular{T, <:Transpose{T,<:$Typ{T}}} where T         
    end)
end

macro _layoutrmul(Typ)
    esc(quote
        LinearAlgebra.rmul!(A::AbstractMatrix, B::$Typ) = ArrayLayouts.rmul!(A, B)
        LinearAlgebra.rmul!(A::StridedMatrix, B::$Typ) = ArrayLayouts.rmul!(A, B)
    end)
end

macro layoutrmul(Typ)
    esc(quote
        ArrayLayouts.@_layoutrmul UpperTriangular{T, <:$Typ{T}} where T
        ArrayLayouts.@_layoutrmul UnitUpperTriangular{T, <:$Typ{T}} where T
        ArrayLayouts.@_layoutrmul LowerTriangular{T, <:$Typ{T}} where T
        ArrayLayouts.@_layoutrmul UnitLowerTriangular{T, <:$Typ{T}} where T

        ArrayLayouts.@_layoutrmul UpperTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutrmul UnitUpperTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutrmul LowerTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutrmul UnitLowerTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T

        ArrayLayouts.@_layoutrmul UpperTriangular{T, <:Transpose{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutrmul UnitUpperTriangular{T, <:Transpose{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutrmul LowerTriangular{T, <:Transpose{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutrmul UnitLowerTriangular{T, <:Transpose{T,<:$Typ{T}}} where T        
    end)
end


# TODO: uncomment out. commented for now due to ambiguities
# LinearAlgebra.lmul!(Q::AbstractQ, v::LayoutVecOrMats) = lmul!(Q, v)
# LinearAlgebra.rmul!(A::LayoutMatrices, Q::AbstractQ) = rmul!(A, Q)

for Typ in (:QRCompactWYQ, :QRPackedQ)
    @eval begin
        LinearAlgebra.lmul!(Q::$Typ{<:Any,<:LayoutMatrix}, v::LayoutVecOrMats) = lmul!(Q, v)
        LinearAlgebra.lmul!(Q::$Typ{<:Any,<:LayoutMatrix}, v::AbstractVecOrMat) = lmul!(Q, v)
        LinearAlgebra.lmul!(Q::AdjointQ{<:Any,<:$Typ{<:Any,<:LayoutMatrix}}, v::LayoutVecOrMats) = lmul!(Q, v)
        LinearAlgebra.lmul!(Q::AdjointQ{<:Any,<:$Typ{<:Any,<:LayoutMatrix}}, v::AbstractVecOrMat) = lmul!(Q, v)
        LinearAlgebra.rmul!(V::LayoutMatrices, Q::$Typ{<:Any,<:LayoutMatrix}) = rmul!(V, Q)
        LinearAlgebra.rmul!(V::AbstractMatrix, Q::$Typ{<:Any,<:LayoutMatrix}) = rmul!(V, Q)
        LinearAlgebra.rmul!(V::LayoutMatrices, Q::AdjointQ{<:Any,<:$Typ{<:Any,<:LayoutMatrix}}) = rmul!(V, Q)
        LinearAlgebra.rmul!(V::AbstractMatrix, Q::AdjointQ{<:Any,<:$Typ{<:Any,<:LayoutMatrix}}) = rmul!(V, Q)
        LinearAlgebra.rmul!(V::StridedVecOrMat{T}, Q::AdjointQ{<:Any,<:$Typ{T,<:LayoutMatrix}}) where T<:BlasReal = rmul!(V, Q)
    end
end
