for Typ in (:Ldiv, :Rdiv)
    @eval begin
        struct $Typ{StyleA, StyleB, AType, BType}
            A::AType
            B::BType
        end

        @inline $Typ{StyleA, StyleB}(A::AType, B::BType) where {StyleA,StyleB,AType,BType} =
            $Typ{StyleA,StyleB,AType,BType}(A,B)

        @inline $Typ(A::AType, B::BType) where {AType,BType} =
            $Typ{typeof(MemoryLayout(AType)),typeof(MemoryLayout(BType)),AType,BType}(A, B)

        @inline BroadcastStyle(::Type{<:$Typ}) = ApplyBroadcastStyle()
        @inline broadcastable(M::$Typ) = M

        similar(A::$Typ, ::Type{T}, axes) where T = similar(Array{T}, axes)
        similar(A::$Typ, ::Type{T}) where T = similar(A, T, axes(A))
        similar(A::$Typ) = similar(A, eltype(A))

        @inline copy(M::$Typ; kwds...) = copyto!(similar(M), M; kwds...)
        @inline materialize(M::$Typ; kwds...) = copy(instantiate(M); kwds...)
    end
end

similar(A::Rdiv{<:DualLayout}, ::Type{T}, (ax1,ax2)) where T = dualadjoint(A.A)(similar(Array{T}, (ax2,)))

@inline _ldivaxes(::Tuple{}, ::Tuple{}) = ()
@inline _ldivaxes(::Tuple{}, Bax::Tuple) = Bax
@inline _ldivaxes(::Tuple{<:Any}, ::Tuple{<:Any}) = ()
@inline _ldivaxes(::Tuple{<:Any}, Bax::Tuple{<:Any,<:Any}) = (OneTo(1),last(Bax))
@inline _ldivaxes(Aax::Tuple{<:Any,<:Any}, ::Tuple{<:Any}) = (last(Aax),)
@inline _ldivaxes(Aax::Tuple{<:Any,<:Any}, Bax::Tuple{<:Any,<:Any}) = (last(Aax),last(Bax))

@inline ldivaxes(A, B) = _ldivaxes(axes(A), axes(B))

@inline axes(L::Ldiv) = ldivaxes(L.A, L.B)
@inline size(L::Ldiv) = map(length, axes(L))
@inline length(L::Ldiv{<:Any,<:Any,<:Any,<:AbstractVector}) =size(L.A, 2)

@inline size(L::Rdiv) = (size(L.A, 1),size(L.B,1))
@inline axes(L::Rdiv) = (axes(L.A, 1),axes(L.B,1))

@inline ndims(L::Ldiv) = ndims(L.B)
@inline ndims(L::Rdiv) = 2
@inline eltype(M::Ldiv) = promote_type(Base.promote_op(inv, eltype(M.A)), eltype(M.B))
@inline eltype(M::Rdiv) = promote_type(eltype(M.A), Base.promote_op(inv, eltype(M.B)))

# Lazy getindex

getindex(L::Ldiv, k...) = _getindex(indextype(L), L, k)
concretize(L::AbstractArray) = convert(Array,L)
concretize(L::Ldiv) = ldiv(concretize(L.A), concretize(L.B))
_getindex(::Type{Tuple{I}}, L::Ldiv, (k,)::Tuple{I}) where I = concretize(L)[k]
_getindex(::Type{Tuple{I,J}}, L::Ldiv, (k,j)::Tuple{Colon,J}) where {I,J} = Ldiv(L.A, L.B[:,j])
_getindex(::Type{Tuple{I,J}}, L::Ldiv, (k,j)::Tuple{I,J}) where {I,J} = L[:,j][k]

check_ldiv_axes(A, B) =
    axes(A,1) == axes(B,1) || throw(DimensionMismatch(LazyString("First axis of A, ", axes(A,1), ", and first axis of B, ", axes(B,1), " must match")))

check_rdiv_axes(A, B) =
    axes(A,2) == axes(B,2) || throw(DimensionMismatch(LazyString("Second axis of A, ", axes(A,2), ", and second axis of B, ", axes(B,2), " must match")))



@inline function instantiate(L::Ldiv)
    check_ldiv_axes(L.A, L.B)
    Ldiv(instantiate(L.A), instantiate(L.B))
end

@inline function instantiate(L::Rdiv)
    check_rdiv_axes(L.A, L.B)
    Rdiv(instantiate(L.A), instantiate(L.B))
end

__ldiv!(::Mat, ::Mat, B) where Mat = error(LazyString("Overload materialize!(::Ldiv{", typeof(MemoryLayout(Mat)), ",", typeof(MemoryLayout(B)), "})"))
__ldiv!(::Mat, ::Mat, B::LayoutArray) where Mat = error(LazyString("Overload materialize!(::Ldiv{", typeof(MemoryLayout(Mat)), ",", typeof(MemoryLayout(B)), "})"))
__ldiv!(_, F, B) = LinearAlgebra.ldiv!(F, B)
@inline _ldiv!(A, B) = __ldiv!(A, factorize(A), B)
@inline _ldiv!(A::Factorization, B) = LinearAlgebra.ldiv!(A, B)
@inline _ldiv!(A::Factorization, B::LayoutArray) = error(LazyString("Overload materialize!(::Ldiv{", typeof(MemoryLayout(A)), ",", typeof(MemoryLayout(B)), "})"))

@inline _ldiv!(dest, A, B; kwds...) = ldiv!(dest, factorize(A), B; kwds...)
@inline _ldiv!(dest, A::Factorization, B; kwds...) = LinearAlgebra.ldiv!(dest, A, B; kwds...)

if VERSION ≥ v"1.10-"
    using LinearAlgebra: TransposeFactorization, AdjointFactorization
else
    const TransposeFactorization = Transpose
    const AdjointFactorization = Adjoint

end
@inline _ldiv!(dest, A::TransposeFactorization{<:Any,<:Factorization}, B; kwds...) = LinearAlgebra.ldiv!(dest, A, B; kwds...)
@inline _ldiv!(dest, A::AdjointFactorization{<:Any,<:Factorization}, B; kwds...) = LinearAlgebra.ldiv!(dest, A, B; kwds...)



@inline ldiv(A, B; kwds...) = materialize(Ldiv(A,B); kwds...)
@inline rdiv(A, B; kwds...) = materialize(Rdiv(A,B); kwds...)

@inline ldiv!(A, B; kwds...) = materialize!(Ldiv(A,B); kwds...)
@inline rdiv!(A, B; kwds...) = materialize!(Rdiv(A,B); kwds...)

@inline ldiv!(C, A, B; kwds...) = copyto!(C, Ldiv(A,B); kwds...)
@inline rdiv!(C, A, B; kwds...) = copyto!(C, Rdiv(A,B); kwds...)

@inline materialize!(M::Ldiv) = _ldiv!(M.A, M.B)
@inline materialize!(M::Rdiv) = ldiv!(M.B', M.A')'
@inline function copyto!(dest::AbstractArray, M::Rdiv; kwds...)
    adj = dualadjoint(dest)
    adj(copyto!(adj(dest), Ldiv(adj(M.B), adj(M.A)); kwds...))
end
@inline copyto!(dest::AbstractArray, M::Ldiv; kwds...) = _ldiv!(dest, M.A, copy(M.B); kwds...)

const MatLdivVec{styleA, styleB, T, V} = Ldiv{styleA, styleB, <:AbstractMatrix{T}, <:AbstractVector{V}}
const MatLdivMat{styleA, styleB, T, V} = Ldiv{styleA, styleB, <:AbstractMatrix{T}, <:AbstractMatrix{V}}
const BlasMatLdivVec{styleA, styleB, T<:BlasFloat} = MatLdivVec{styleA, styleB, T, T}
const BlasMatLdivMat{styleA, styleB, T<:BlasFloat} = MatLdivMat{styleA, styleB, T, T}

const MatRdivMat{styleA, styleB, T, V} = Rdiv{styleA, styleB, <:AbstractMatrix{T}, <:AbstractMatrix{V}}
const BlasMatRdivMat{styleA, styleB, T<:BlasFloat} = MatRdivMat{styleA, styleB, T, T}

materialize!(M::Ldiv{ScalarLayout}) = invoke(LinearAlgebra.ldiv!, Tuple{Number,AbstractArray}, M.A, M.B)
materialize!(M::Rdiv{<:Any,ScalarLayout}) = invoke(LinearAlgebra.rdiv!, Tuple{AbstractArray,Number}, M.A, M.B)

function materialize!(M::Ldiv{ScalarLayout,<:SymmetricLayout})
    ldiv!(M.A, symmetricdata(M.B))
    M.B
end
function materialize!(M::Ldiv{ScalarLayout,<:HermitianLayout})
    ldiv!(M.A, hermitiandata(M.B))
    M.B
end
function materialize!(M::Rdiv{<:SymmetricLayout,ScalarLayout})
    rdiv!(symmetricdata(M.A), M.B)
    M.A
end
function materialize!(M::Rdiv{<:HermitianLayout,ScalarLayout})
    rdiv!(hermitiandata(M.A), M.B)
    M.A
end

macro _layoutldiv(Typ)
    ret = quote
        LinearAlgebra.ldiv!(A::$Typ, x::AbstractVector; kwds...) = ArrayLayouts.ldiv!(A,x; kwds...)
        LinearAlgebra.ldiv!(A::$Typ, x::AbstractMatrix; kwds...) = ArrayLayouts.ldiv!(A,x; kwds...)
        LinearAlgebra.ldiv!(A::$Typ, x::StridedVector; kwds...) = ArrayLayouts.ldiv!(A,x; kwds...)
        LinearAlgebra.ldiv!(A::$Typ, x::StridedMatrix; kwds...) = ArrayLayouts.ldiv!(A,x; kwds...)

        LinearAlgebra.ldiv!(A::Factorization, x::$Typ; kwds...) = ArrayLayouts.ldiv!(A,x; kwds...)
        LinearAlgebra.ldiv!(A::LU, x::$Typ; kwds...) = ArrayLayouts.ldiv!(A,x; kwds...)
        LinearAlgebra.ldiv!(A::Cholesky, x::$Typ; kwds...) = ArrayLayouts.ldiv!(A,x; kwds...)
        LinearAlgebra.ldiv!(A::LinearAlgebra.QRCompactWY, x::$Typ; kwds...) = ArrayLayouts.ldiv!(A,x; kwds...)

        LinearAlgebra.ldiv!(A::Bidiagonal, B::$Typ; kwds...) = ArrayLayouts.ldiv!(A,B; kwds...)

        LinearAlgebra.rdiv!(A::AbstractMatrix, B::$Typ; kwds...) = ArrayLayouts.rdiv!(A,B; kwds...)

        # Fix ambiguity issue
        LinearAlgebra.rdiv!(A::StridedMatrix, B::$Typ; kwds...) = ArrayLayouts.rdiv!(A,B; kwds...)

        (\)(A::$Typ, x::AbstractVector; kwds...) = ArrayLayouts.ldiv(A,x; kwds...)
        (\)(A::$Typ, x::AbstractMatrix; kwds...) = ArrayLayouts.ldiv(A,x; kwds...)

        (\)(x::AbstractMatrix, A::$Typ; kwds...) = ArrayLayouts.ldiv(x,A; kwds...)
        (\)(x::LinearAlgebra.HermOrSym, A::$Typ; kwds...) = ArrayLayouts.ldiv(x,A; kwds...)
        if VERSION < v"1.9-" # disambiguation
            \(x::LinearAlgebra.HermOrSym{<:Any,<:StridedMatrix}, A::$Typ; kwds...) = ArrayLayouts.ldiv(x,A; kwds...)
        end
        (\)(x::UpperTriangular, A::$Typ; kwds...) = ArrayLayouts.ldiv(x,A; kwds...)
        (\)(x::UnitUpperTriangular, A::$Typ; kwds...) = ArrayLayouts.ldiv(x,A; kwds...)
        (\)(x::LowerTriangular, A::$Typ; kwds...) = ArrayLayouts.ldiv(x,A; kwds...)
        (\)(x::UnitLowerTriangular, A::$Typ; kwds...) = ArrayLayouts.ldiv(x,A; kwds...)
        (\)(x::Diagonal, A::$Typ; kwds...) = ArrayLayouts.ldiv(x,A; kwds...)

        (\)(A::Bidiagonal{<:Number}, B::$Typ{<:Number}; kwds...) = ArrayLayouts.ldiv(A,B; kwds...)
        (\)(A::Bidiagonal, B::$Typ; kwds...) = ArrayLayouts.ldiv(A,B; kwds...)
        (\)(transA::Transpose{<:Number,<:Bidiagonal{<:Number}}, B::$Typ{<:Number}; kwds...) = ArrayLayouts.ldiv(transA,B; kwds...)
        (\)(transA::Transpose{<:Any,<:Bidiagonal}, B::$Typ; kwds...) = ArrayLayouts.ldiv(transA,B; kwds...)
        (\)(adjA::Adjoint{<:Number,<:Bidiagonal{<:Number}}, B::$Typ{<:Number}; kwds...) = ArrayLayouts.ldiv(adjA,B; kwds...)
        (\)(adjA::Adjoint{<:Any,<:Bidiagonal}, B::$Typ; kwds...) = ArrayLayouts.ldiv(adjA,B; kwds...)

        (\)(x::$Typ, A::$Typ; kwds...) = ArrayLayouts.ldiv(x,A; kwds...)

        (/)(A::$Typ, x::AbstractVector; kwds...) = ArrayLayouts.rdiv(A,x; kwds...)
        (/)(A::$Typ, x::AbstractMatrix; kwds...) = ArrayLayouts.rdiv(A,x; kwds...)

        (/)(x::AbstractMatrix, A::$Typ; kwds...) = ArrayLayouts.rdiv(x,A; kwds...)
        (/)(D::Diagonal, A::$Typ; kwds...) = ArrayLayouts.rdiv(D,A; kwds...)
        (/)(A::$Typ, D::Diagonal; kwds...) = ArrayLayouts.rdiv(A,D; kwds...)

        (/)(x::$Typ, A::$Typ; kwds...) = ArrayLayouts.rdiv(x,A; kwds...)
        (/)(D::Adjoint{<:Any,<:AbstractVector}, A::$Typ; kwds...) = ArrayLayouts.rdiv(D,A; kwds...)
        (/)(D::Transpose{<:Any,<:AbstractVector}, A::$Typ; kwds...) = ArrayLayouts.rdiv(D,A; kwds...)
    end
    if Typ ≠ :LayoutVector
        ret = quote
            $ret
            (\)(A::$Typ, x::LayoutVector; kwds...) = ArrayLayouts.ldiv(A,x; kwds...)
            (\)(x::Diagonal{<:Any,<:LayoutVector}, A::$Typ; kwds...) = ArrayLayouts.ldiv(x,A; kwds...)
            \(A::$Typ, B::Diagonal{<:Any,<:LayoutVector}; kwds...) = ArrayLayouts.ldiv(A, B; kwds...)
        end
    end
    if Typ ≠ :LayoutMatrix
        ret = quote
            $ret
            (\)(A::$Typ, x::LayoutMatrix; kwds...) = ArrayLayouts.ldiv(A,x; kwds...)
            (/)(x::LayoutMatrix, A::$Typ; kwds...) = ArrayLayouts.rdiv(x,A; kwds...)
        end
    end
    esc(ret)
end

macro layoutldiv(Typ)
    esc(quote
        ArrayLayouts.@_layoutldiv $Typ
        ArrayLayouts.@_layoutldiv UpperTriangular{T, <:$Typ{T}} where T
        ArrayLayouts.@_layoutldiv UnitUpperTriangular{T, <:$Typ{T}} where T
        ArrayLayouts.@_layoutldiv LowerTriangular{T, <:$Typ{T}} where T
        ArrayLayouts.@_layoutldiv UnitLowerTriangular{T, <:$Typ{T}} where T

        ArrayLayouts.@_layoutldiv UpperTriangular{T, <:SubArray{T,2,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutldiv UnitUpperTriangular{T, <:SubArray{T,2,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutldiv LowerTriangular{T, <:SubArray{T,2,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutldiv UnitLowerTriangular{T, <:SubArray{T,2,<:$Typ{T}}} where T

        ArrayLayouts.@_layoutldiv UpperTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutldiv UnitUpperTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutldiv LowerTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutldiv UnitLowerTriangular{T, <:Adjoint{T,<:$Typ{T}}} where T

        ArrayLayouts.@_layoutldiv UpperTriangular{T, <:Transpose{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutldiv UnitUpperTriangular{T, <:Transpose{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutldiv LowerTriangular{T, <:Transpose{T,<:$Typ{T}}} where T
        ArrayLayouts.@_layoutldiv UnitLowerTriangular{T, <:Transpose{T,<:$Typ{T}}} where T
    end)
end

@_layoutldiv LayoutVector
