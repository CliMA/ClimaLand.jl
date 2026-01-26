module ArrayLayouts
using Base: _typed_hcat
using Base.Broadcast, LinearAlgebra, FillArrays
using LinearAlgebra.BLAS

using Base: AbstractCartesianIndex, OneTo, oneto, RangeIndex, ReinterpretArray, ReshapedArray,
            Slice, tuple_type_tail, unalias,
            @propagate_inbounds

import Base: axes, size, length, eltype, ndims, first, last, diff, isempty, union, sort!, sort,
                ==, *, +, -, /, \, copy, copyto!, similar, getproperty, getindex, strides,
                reverse, unsafe_convert, convert, view

using Base.Broadcast: Broadcasted

import Base.Broadcast: BroadcastStyle, broadcastable, instantiate, materialize, materialize!

using LinearAlgebra: AbstractQ, AbstractTriangular, AdjOrTrans, AdjointAbsVec, HermOrSym, HessenbergQ, QRCompactWYQ,
                     QRPackedQ, RealHermSymComplexHerm, TransposeAbsVec, _apply_ipiv_rows!, checknonsingular,
                     checksquare, chkfullrank, cholcopy, ipiv2perm, QRCompactWYQ, AdjointQ, QRPackedQ

using LinearAlgebra.BLAS: BlasComplex, BlasFloat, BlasReal

const AdjointQtype{T} = isdefined(LinearAlgebra, :AdjointQ) ? LinearAlgebra.AdjointQ{T} : Adjoint{T,<:AbstractQ}

using FillArrays: AbstractFill, axes_print_matrix_row, getindex_value
using StaticArrays

using Base: require_one_based_indexing

export materialize, materialize!, MulAdd, muladd!, Ldiv, Rdiv, Lmul, Rmul, Dot,
        lmul, mul, ldiv, rdiv, mul, MemoryLayout, AbstractStridedLayout,
        DenseColumnMajor, ColumnMajor, ZerosLayout, FillLayout, AbstractColumnMajor, RowMajor, AbstractRowMajor, UnitStride,
        DiagonalLayout, ScalarLayout, SymTridiagonalLayout, TridiagonalLayout, BidiagonalLayout,
        HermitianLayout, SymmetricLayout, TriangularLayout,
        UnknownLayout, AbstractBandedLayout, ApplyBroadcastStyle, ConjLayout, AbstractFillLayout, DualLayout,
        colsupport, rowsupport, layout_getindex, AbstractQLayout, LayoutArray, LayoutMatrix, LayoutVector,
        RangeCumsum

if VERSION < v"1.7-"
    const ColumnNorm = Val{true}
    const RowMaximum = Val{true}
    const NoPivot = Val{false}
end

if VERSION < v"1.8-"
    const CRowMaximum = Val{true}
    const CNoPivot = Val{false}
else
    const CRowMaximum = RowMaximum
    const CNoPivot = NoPivot
end

if VERSION ≥ v"1.11.0-DEV.21"
    using LinearAlgebra: UpperOrLowerTriangular
else
    const UpperOrLowerTriangular{T,S} = Union{LinearAlgebra.UpperTriangular{T,S},
                                              LinearAlgebra.UnitUpperTriangular{T,S},
                                              LinearAlgebra.LowerTriangular{T,S},
                                              LinearAlgebra.UnitLowerTriangular{T,S}}
end

@static if VERSION ≥ v"1.8.0"
    import Base: LazyString
else
    const LazyString = string
end

# Originally defined in FillArrays
_copy_oftype(A::AbstractArray, ::Type{S}) where {S} = eltype(A) == S ? copy(A) : AbstractArray{S}(A)
_copy_oftype(A::AbstractRange, ::Type{S}) where {S} = eltype(A) == S ? copy(A) : map(S, A)

struct ApplyBroadcastStyle <: BroadcastStyle end
@inline function copyto!(dest::AbstractArray, bc::Broadcasted{ApplyBroadcastStyle})
    @assert length(bc.args) == 1
    copyto!(dest, first(bc.args))
end

# Subtypes of LayoutArray default to
# ArrayLayouts routines
abstract type LayoutArray{T,N} <: AbstractArray{T,N} end
const LayoutMatrix{T} = LayoutArray{T,2}
const LayoutVector{T} = LayoutArray{T,1}
const LayoutVecOrMat{T} = Union{LayoutVector{T},LayoutMatrix{T}}
const LayoutMatrices{T} = Union{LayoutMatrix{T}, SubArray{T,2,<:LayoutMatrix}}
const LayoutVecOrMats{T} = Union{LayoutVecOrMat{T}, SubArray{T,1,<:LayoutVecOrMat}, SubArray{T,2,<:LayoutMatrix}}


## TODO: Following are type piracy which may be removed in Julia v1.5
## No, it can't, because strides(::AdjointAbsMat) is defined only for real eltype!
_transpose_strides(a) = (a,1)
_transpose_strides(a,b) = (b,a)
strides(A::Adjoint) = _transpose_strides(strides(parent(A))...)
strides(A::Transpose) = _transpose_strides(strides(parent(A))...)

"""
    ConjPtr{T}

A memory address referring to complex conjugated data of type T. However, there is no guarantee
that the memory is actually valid, or that it actually represents the complex conjugate of data of
the specified type.
"""
struct ConjPtr{T}
    ptr::Ptr{T}
end


# work-around issue with complex conjugation of pointer
unsafe_convert(::Type{Ptr{T}}, Ac::Adjoint{<:Complex}) where T<:Complex = unsafe_convert(ConjPtr{T}, parent(Ac))
unsafe_convert(::Type{ConjPtr{T}}, Ac::Adjoint{<:Complex}) where T<:Complex = unsafe_convert(Ptr{T}, parent(Ac))
function unsafe_convert(::Type{ConjPtr{T}}, V::SubArray{T,2}) where {T}
    kr, jr = parentindices(V)
    unsafe_convert(Ptr{T}, view(parent(V)', jr, kr))
end

Base.elsize(::Type{<:Adjoint{<:Complex,P}}) where P<:AbstractVecOrMat = conjelsize(P)
conjelsize(::Type{<:Adjoint{<:Complex,P}}) where P<:AbstractVecOrMat = Base.elsize(P)
conjelsize(::Type{<:Transpose{<:Any, P}}) where {P<:AbstractVecOrMat} = conjelsize(P)
conjelsize(::Type{<:PermutedDimsArray{<:Any, <:Any, <:Any, <:Any, P}}) where {P} = conjelsize(P)
conjelsize(::Type{<:ReshapedArray{<:Any,<:Any,P}}) where {P} = conjelsize(P)
conjelsize(::Type{<:SubArray{<:Any,<:Any,P}}) where {P} = conjelsize(P)
conjelsize(A::AbstractArray) = conjelsize(typeof(A))

include("memorylayout.jl")
include("mul.jl")
include("muladd.jl")
include("lmul.jl")
include("ldiv.jl")
include("diagonal.jl")
include("triangular.jl")
include("factorizations.jl")

# Extend this function if you're only looking to dispatch on the axes
@inline sub_materialize_axes(V, _) = copyto!(similar(V, axes(V)), V)
@inline sub_materialize(_, V, ax) = sub_materialize_axes(V, ax)
@inline sub_materialize(L, V) = sub_materialize(L, V, axes(V))
@inline sub_materialize(V::SubArray) = sub_materialize(MemoryLayout(V), V)
@inline sub_materialize(V) = copy(V) # Anything not a SubArray is already materialized.

copy(A::SubArray{<:Any,N,<:LayoutArray}) where N = sub_materialize(A)
copy(A::SubArray{<:Any,N,<:AdjOrTrans{<:Any,<:LayoutArray}}) where N = sub_materialize(A)

@inline layout_getindex(A, I...) = sub_materialize(view(A, I...))
@propagate_inbounds function layout_getindex(A::AbstractArray, k::Int...)
    Base.error_if_canonical_getindex(IndexStyle(A), A, k...)
    Base._getindex(IndexStyle(A), A, k...)
end


# avoid 0-dimensional views
Base.@propagate_inbounds layout_getindex(A::AbstractArray, I::CartesianIndex) = A[to_indices(A, (I,))...]

macro _layoutgetindex(Typ)
    esc(quote
        @inline getindex(A::$Typ, kr::Colon, jr::Colon) = ArrayLayouts.layout_getindex(A, kr, jr)
        @inline getindex(A::$Typ, kr::Colon, jr::AbstractUnitRange) = ArrayLayouts.layout_getindex(A, kr, jr)
        @inline getindex(A::$Typ, kr::AbstractUnitRange, jr::Colon) = ArrayLayouts.layout_getindex(A, kr, jr)
        @inline getindex(A::$Typ, kr::AbstractUnitRange, jr::AbstractUnitRange) = ArrayLayouts.layout_getindex(A, kr, jr)
        @inline getindex(A::$Typ, kr::AbstractVector, jr::AbstractVector) = ArrayLayouts.layout_getindex(A, kr, jr)
        @inline getindex(A::$Typ, kr::Colon, jr::AbstractVector) = ArrayLayouts.layout_getindex(A, kr, jr)
        @inline getindex(A::$Typ, kr::Colon, jr::Integer) = ArrayLayouts.layout_getindex(A, kr, jr)
        @inline getindex(A::$Typ, kr::AbstractVector, jr::Colon) = ArrayLayouts.layout_getindex(A, kr, jr)
        @inline getindex(A::$Typ, kr::Integer, jr::Colon) = ArrayLayouts.layout_getindex(A, kr, jr)
        @inline getindex(A::$Typ, kr::Integer, jr::AbstractVector) = ArrayLayouts.layout_getindex(A, kr, jr)
    end)
end

macro layoutgetindex(Typ)
    esc(quote
        ArrayLayouts.@_layoutgetindex $Typ
        ArrayLayouts.@_layoutgetindex ArrayLayouts.UpperOrLowerTriangular{<:Any,<:$Typ}
        ArrayLayouts.@_layoutgetindex LinearAlgebra.Symmetric{<:Any,<:$Typ}
        ArrayLayouts.@_layoutgetindex LinearAlgebra.Hermitian{<:Any,<:$Typ}
        ArrayLayouts.@_layoutgetindex LinearAlgebra.Adjoint{<:Any,<:$Typ}
        ArrayLayouts.@_layoutgetindex LinearAlgebra.Transpose{<:Any,<:$Typ}
        ArrayLayouts.@_layoutgetindex LinearAlgebra.SubArray{<:Any,2,<:$Typ}
    end)
end


macro layoutmatrix(Typ)
    esc(quote
        ArrayLayouts.@layoutldiv $Typ
        ArrayLayouts.@layoutmul $Typ
        ArrayLayouts.@layoutlmul $Typ
        ArrayLayouts.@layoutrmul $Typ
        ArrayLayouts.@layoutfactorizations $Typ
        ArrayLayouts.@layoutgetindex $Typ
    end)
end

@layoutmatrix LayoutMatrix

for Typ in (:LayoutArray, :(Transpose{<:Any,<:LayoutMatrix}), :(Adjoint{<:Any,<:LayoutMatrix}), :(Symmetric{<:Any,<:LayoutMatrix}), :(Hermitian{<:Any,<:LayoutMatrix}))
    @eval begin
        LinearAlgebra.lmul!(α::Number, A::$Typ) = lmul!(α, A)
        LinearAlgebra.rmul!(A::$Typ, α::Number) = rmul!(A, α)
        LinearAlgebra.ldiv!(α::Number, A::$Typ) = ldiv!(α, A)
        LinearAlgebra.rdiv!(A::$Typ, α::Number) = rdiv!(A, α)
    end
end

getindex(A::LayoutVector, kr::AbstractVector) = layout_getindex(A, kr)
getindex(A::LayoutVector, kr::Colon) = layout_getindex(A, kr)
getindex(A::AdjOrTrans{<:Any,<:LayoutVector}, kr::Integer, jr::Colon) = layout_getindex(A, kr, jr)
getindex(A::AdjOrTrans{<:Any,<:LayoutVector}, kr::Integer, jr::AbstractVector) = layout_getindex(A, kr, jr)

*(a::Zeros{<:Any,2}, b::LayoutMatrix) = FillArrays.mult_zeros(a, b)
*(a::LayoutMatrix, b::Zeros{<:Any,2}) = FillArrays.mult_zeros(a, b)
*(a::LayoutMatrix, b::Zeros{<:Any,1}) = FillArrays.mult_zeros(a, b)
*(a::Transpose{T, <:LayoutMatrix{T}} where T, b::Zeros{<:Any, 2}) = FillArrays.mult_zeros(a, b)
*(a::Adjoint{T, <:LayoutMatrix{T}} where T, b::Zeros{<:Any, 2}) = FillArrays.mult_zeros(a, b)
*(A::Adjoint{<:Any, <:Zeros{<:Any,1}}, B::Diagonal{<:Any,<:LayoutVector}) = (B' * A')'
*(A::Transpose{<:Any, <:Zeros{<:Any,1}}, B::Diagonal{<:Any,<:LayoutVector}) = transpose(transpose(B) * transpose(A))
*(a::Adjoint{<:Number,<:LayoutVector}, b::Zeros{<:Number,1})= FillArrays._adjvec_mul_zeros(a, b)
function *(a::Transpose{T, <:LayoutVector{T}}, b::Zeros{T, 1}) where T<:Real
    la, lb = length(a), length(b)
    if la ≠ lb
        throw(DimensionMismatch(LazyString("dot product arguments have lengths ", la, " and ", lb)))
    end
    return zero(T)
end

*(A::Diagonal{<:Any,<:LayoutVector}, B::Diagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::Diagonal{<:Any,<:LayoutVector}, B::AbstractMatrix) = mul(A, B)
*(A::AbstractMatrix, B::Diagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::Union{Bidiagonal, Tridiagonal}, B::Diagonal{<:Any, <:LayoutVector}) = mul(A, B) # ambiguity
*(A::Diagonal{<:Any, <:LayoutVector}, B::Union{Bidiagonal, Tridiagonal}) = mul(A, B) # ambiguity
*(A::Diagonal{<:Any,<:LayoutVector}, B::LayoutMatrix) = mul(A, B)
*(A::LayoutMatrix, B::Diagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::UpperTriangular, B::Diagonal{<:Any,<:LayoutVector}) = mul(A, B)
*(A::Diagonal{<:Any,<:LayoutVector}, B::UpperTriangular) = mul(A, B)
*(A::Diagonal{<:Any,<:LayoutVector}, B::Diagonal) = mul(A, B)
*(A::Diagonal, B::Diagonal{<:Any,<:LayoutVector}) = mul(A, B)

for Mod in (:Adjoint, :Transpose, :Symmetric, :Hermitian)
    @eval begin
        *(A::Diagonal{<:Any,<:LayoutVector}, B::$Mod{<:Any,<:LayoutMatrix}) = mul(A,B)
        *(A::$Mod{<:Any,<:LayoutMatrix}, B::Diagonal{<:Any,<:LayoutVector}) = mul(A,B)
        *(A::Diagonal{<:Any,<:LayoutVector}, B::$Mod{<:Any,<:AbstractMatrix}) = mul(A,B)
        *(A::$Mod{<:Any,<:AbstractMatrix}, B::Diagonal{<:Any,<:LayoutVector}) = mul(A,B)
    end
end

*(A::LayoutMatrix, B::Adjoint{<:Any,<:AbstractTriangular}) = mul(A, B)
*(A::LayoutMatrix, B::Transpose{<:Any,<:AbstractTriangular}) = mul(A, B)
*(A::Adjoint{<:Any,<:AbstractTriangular}, B::LayoutMatrix) = mul(A, B)
*(A::Transpose{<:Any,<:AbstractTriangular}, B::LayoutMatrix) = mul(A, B)

*(A::Transpose{<:Any,<:AbstractTriangular}, B::Adjoint{<:Any,<:LayoutMatrix}) = mul(A, B)
*(A::Adjoint{<:Any,<:AbstractTriangular}, B::Transpose{<:Any,<:LayoutMatrix}) = mul(A, B)
*(A::Adjoint{<:Any,<:AbstractTriangular}, B::Adjoint{<:Any,<:LayoutMatrix}) = mul(A, B)
*(A::Transpose{<:Any,<:AbstractTriangular}, B::Transpose{<:Any,<:LayoutMatrix}) = mul(A, B)
*(A::Adjoint{<:Any,<:LayoutMatrix}, B::Transpose{<:Any,<:AbstractTriangular}) = mul(A, B)
*(A::Transpose{<:Any,<:LayoutMatrix}, B::Adjoint{<:Any,<:AbstractTriangular}) = mul(A, B)
*(A::Adjoint{<:Any,<:LayoutMatrix}, B::Adjoint{<:Any,<:AbstractTriangular}) = mul(A, B)
*(A::Transpose{<:Any,<:LayoutMatrix}, B::Transpose{<:Any,<:AbstractTriangular}) = mul(A, B)

\(A::Diagonal{<:Any,<:LayoutVector}, B::Diagonal{<:Any,<:LayoutVector}) = ldiv(A, B)
\(A::Diagonal{<:Any,<:LayoutVector}, B::AbstractMatrix) = ldiv(A, B)
\(A::AbstractMatrix, B::Diagonal{<:Any,<:LayoutVector}) = ldiv(A, B)
\(A::Diagonal{<:Any,<:LayoutVector}, B::Diagonal) = ldiv(A, B)
\(A::Diagonal, B::Diagonal{<:Any,<:LayoutVector}) = ldiv(A, B)


copyto!_layout(_, _, dest::AbstractArray{T,N}, src::AbstractArray{V,N}) where {T,V,N} =
    invoke(copyto!, Tuple{AbstractArray{T,N},AbstractArray{V,N}}, dest, src)

const _copyto! = copyto!_layout # TODO: deprecate


copyto!_layout(dest, src) = copyto!_layout(MemoryLayout(dest), MemoryLayout(src), dest, src)
copyto!(dest::LayoutArray{<:Any,N}, src::LayoutArray{<:Any,N}) where N = copyto!_layout(dest, src)
copyto!(dest::AbstractArray{<:Any,N}, src::LayoutArray{<:Any,N}) where N = copyto!_layout(dest, src)
copyto!(dest::LayoutArray{<:Any,N}, src::AbstractArray{<:Any,N}) where N = copyto!_layout(dest, src)

copyto!(dest::SubArray{<:Any,N,<:LayoutArray}, src::SubArray{<:Any,N,<:LayoutArray}) where N = copyto!_layout(dest, src)
copyto!(dest::SubArray{<:Any,N,<:LayoutArray}, src::LayoutArray{<:Any,N}) where N = copyto!_layout(dest, src)
copyto!(dest::LayoutArray{<:Any,N}, src::SubArray{<:Any,N,<:LayoutArray}) where N = copyto!_layout(dest, src)
copyto!(dest::SubArray{<:Any,N,<:LayoutArray}, src::AbstractArray{<:Any,N}) where N = copyto!_layout(dest, src)
copyto!(dest::AbstractArray{<:Any,N}, src::SubArray{<:Any,N,<:LayoutArray}) where N = copyto!_layout(dest, src)

copyto!(dest::LayoutMatrix, src::AdjOrTrans{<:Any,<:LayoutArray}) = copyto!_layout(dest, src)
copyto!(dest::LayoutMatrix, src::SubArray{<:Any,2,<:AdjOrTrans{<:Any,<:LayoutArray}}) = copyto!_layout(dest, src)
copyto!(dest::AbstractMatrix, src::AdjOrTrans{<:Any,<:LayoutArray}) = copyto!_layout(dest, src)
copyto!(dest::SubArray{<:Any,2,<:LayoutArray}, src::AdjOrTrans{<:Any,<:LayoutArray}) = copyto!_layout(dest, src)
copyto!(dest::SubArray{<:Any,2,<:LayoutMatrix}, src::SubArray{<:Any,2,<:AdjOrTrans{<:Any,<:LayoutArray}}) = copyto!_layout(dest, src)
copyto!(dest::AbstractMatrix, src::SubArray{<:Any,2,<:AdjOrTrans{<:Any,<:LayoutArray}}) = copyto!_layout(dest, src)
if isdefined(LinearAlgebra, :copymutable_oftype)
    LinearAlgebra.copymutable_oftype(A::Union{LayoutArray,Symmetric{<:Any,<:LayoutMatrix},Hermitian{<:Any,<:LayoutMatrix},
                                                                UpperOrLowerTriangular{<:Any,<:LayoutMatrix},
                                                                AdjOrTrans{<:Any,<:LayoutMatrix}}, ::Type{T}) where T = copymutable_oftype_layout(MemoryLayout(A), A, T)
end
copymutable_oftype_layout(_, A, ::Type{S}) where S = copyto!(similar(A, S), A)

# avoid bad copy in Base
Base.map(::typeof(copy), D::Diagonal{<:LayoutArray}) = Diagonal(map(copy, D.diag))
Base.permutedims(D::Diagonal{<:Any,<:LayoutVector}) = D


zero!(A) = zero!(MemoryLayout(A), A)
zero!(_, A) = fill!(A,zero(eltype(A)))
function zero!(::DualLayout, A)
    zero!(parent(A))
    A
end
function zero!(_, A::AbstractArray{<:AbstractArray})
    for a in A
        zero!(a)
    end
    A
end

zero!(_, A::AbstractArray{<:SArray}) = fill!(A,zero(eltype(A)))

_norm(_, A, p) = invoke(norm, Tuple{Any,Real}, A, p)
LinearAlgebra.norm(A::LayoutArray, p::Real=2) = _norm(MemoryLayout(A), A, p)
LinearAlgebra.norm(A::SubArray{<:Any,N,<:LayoutArray}, p::Real=2) where N = _norm(MemoryLayout(A), A, p)


_fill_lmul!(β, A::AbstractArray) = iszero(β) ? zero!(A) : lmul!(β, A)


# Elementary reflection similar to LAPACK. The reflector is not Hermitian but
# ensures that tridiagonalization of Hermitian matrices become real. See lawn72
@inline function reflector!(x::AbstractVector)
    require_one_based_indexing(x)
    n = length(x)
    n == 0 && return zero(eltype(x))
    @inbounds begin
        ξ1 = x[1]
        normu = abs2(ξ1)
        for i = 2:n
            normu += abs2(x[i])
        end
        if iszero(normu)
            return zero(ξ1/normu)
        end
        normu = sqrt(normu)
        ν = copysign(normu, real(ξ1))
        ξ1 += ν
        x[1] = -ν
        for i = 2:n
            x[i] /= ξ1
        end
    end
    ξ1/ν
end

# apply reflector from left
@inline function reflectorApply!(x::AbstractVector, τ::Number, A::AbstractVecOrMat)
    m,n = size(A,1),size(A,2)
    if length(x) != m
        throw(DimensionMismatch(LazyString("reflector has length ", length(x), ", which must match the first dimension of matrix A, ", m)))
    end
    m == 0 && return A
    @inbounds begin
        for j = 1:n
            # dot
            vAj = A[1, j]
            for i = 2:m
                vAj += x[i]'*A[i, j]
            end

            vAj = conj(τ)*vAj

            # ger
            A[1, j] -= vAj
            for i = 2:m
                A[i, j] -= x[i]*vAj
            end
        end
    end
    return A
end

###
# printing
###

layout_replace_in_print_matrix(_, A, i, j, s) =
    i in colsupport(A,j) ? s : Base.replace_with_centered_mark(s)

Base.replace_in_print_matrix(A::Union{LayoutVector,
                                      LayoutMatrix,
                                      UpperTriangular{<:Any,<:LayoutMatrix},
                                      UnitUpperTriangular{<:Any,<:LayoutMatrix},
                                      LowerTriangular{<:Any,<:LayoutMatrix},
                                      UnitLowerTriangular{<:Any,<:LayoutMatrix},
                                      AdjOrTrans{<:Any,<:LayoutVecOrMat},
                                      HermOrSym{<:Any,<:LayoutMatrix},
                                      SubArray{<:Any,2,<:LayoutMatrix},
                                      UpperTriangular{<:Any,<:AdjOrTrans{<:Any,<:LayoutMatrix}},
                                      UnitUpperTriangular{<:Any,<:AdjOrTrans{<:Any,<:LayoutMatrix}},
                                      LowerTriangular{<:Any,<:AdjOrTrans{<:Any,<:LayoutMatrix}},
                                      UnitLowerTriangular{<:Any,<:AdjOrTrans{<:Any,<:LayoutMatrix}}}, i::Integer, j::Integer, s::AbstractString) =
    layout_replace_in_print_matrix(MemoryLayout(A), A, i, j, s)

Base.print_matrix_row(io::IO,
        X::Union{LayoutMatrix,
        LayoutVector,
        UpperOrLowerTriangular{<:Any,<:LayoutMatrix},
        AdjOrTrans{<:Any,<:LayoutMatrix},
        AdjOrTrans{<:Any,<:LayoutVector},
        HermOrSym{<:Any,<:LayoutMatrix},
        SubArray{<:Any,2,<:LayoutMatrix},
        Diagonal{<:Any,<:LayoutVector}}, A::Vector,
        i::Integer, cols::AbstractVector, sep::AbstractString, idxlast::Integer=last(axes(X, 2))) =
        axes_print_matrix_row(axes(X), io, X, A, i, cols, sep)


include("cumsum.jl")

###
# support overloading hcat/vcat for ∞-arrays
###



typed_hcat(::Type{T}, ::Tuple, A...) where T = Base._typed_hcat(T, A)
typed_vcat(::Type{T}, ::Tuple, A...) where T = Base._typed_vcat(T, A)
typed_hcat(::Type{T}, A...) where T = typed_hcat(T, (size(A[1],1), sum(size.(A,2))), A...)
typed_vcat(::Type{T}, A...) where T = typed_vcat(T, (sum(size.(A,1)),size(A[1],2)), A...)

Base.typed_vcat(::Type{T}, A::LayoutVecOrMats, B::AbstractVecOrMat...) where T = typed_vcat(T, A, B...)
Base.typed_hcat(::Type{T}, A::LayoutVecOrMats, B::AbstractVecOrMat...) where T = typed_hcat(T, A, B...)
Base.typed_vcat(::Type{T}, A::LayoutVecOrMats, B::LayoutVecOrMats, C::AbstractVecOrMat...) where T = typed_vcat(T, A, B, C...)
Base.typed_hcat(::Type{T}, A::LayoutVecOrMats, B::LayoutVecOrMats, C::AbstractVecOrMat...) where T = typed_hcat(T, A, B, C...)
Base.typed_vcat(::Type{T}, A::AbstractVecOrMat, B::LayoutVecOrMats, C::AbstractVecOrMat...) where T = typed_vcat(T, A, B, C...)
Base.typed_hcat(::Type{T}, A::AbstractVecOrMat, B::LayoutVecOrMats, C::AbstractVecOrMat...) where T = typed_hcat(T, A, B, C...)

end # module
