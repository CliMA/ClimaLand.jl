
tuple_type_reverse(::Type{T}) where T<:Tuple = Tuple{reverse(tuple(T.parameters...))...}
tuple_type_reverse(::Type{Tuple{}}) = Tuple{}
tuple_type_reverse(::Type{Tuple{A}}) where A = Tuple{A}
tuple_type_reverse(::Type{Tuple{A,B}}) where {A,B} = Tuple{B,A}


abstract type MemoryLayout end
struct UnknownLayout <: MemoryLayout end
abstract type AbstractStridedLayout <: MemoryLayout end
abstract type AbstractIncreasingStrides <: AbstractStridedLayout end
abstract type AbstractColumnMajor <: AbstractIncreasingStrides end
struct DenseColumnMajor <: AbstractColumnMajor end
struct ColumnMajor <: AbstractColumnMajor end
struct IncreasingStrides <: AbstractIncreasingStrides end
abstract type AbstractDecreasingStrides <: AbstractStridedLayout end
abstract type AbstractRowMajor <: AbstractDecreasingStrides end
struct DenseRowMajor <: AbstractRowMajor end
struct RowMajor <: AbstractRowMajor end
struct DecreasingStrides <: AbstractIncreasingStrides end
struct UnitStride{D} <: AbstractStridedLayout end
struct StridedLayout <: AbstractStridedLayout end
struct ScalarLayout <: MemoryLayout end

"""
    UnknownLayout()

is returned by `MemoryLayout(A)` if it is unknown how the entries of an array `A`
are stored in memory.
"""
UnknownLayout

"""
    AbstractStridedLayout

is an abstract type whose subtypes are returned by `MemoryLayout(A)`
if an array `A` has storage laid out at regular offsets in memory,
and which can therefore be passed to external C and Fortran functions expecting
this memory layout.

Julia's internal linear algebra machinery will automatically (and invisibly)
dispatch to BLAS and LAPACK routines if the memory layout is BLAS compatible and
the element type is a `Float32`, `Float64`, `ComplexF32`, or `ComplexF64`.
In this case, one must implement the strided array interface, which requires
overrides of `strides(A::MyMatrix)` and `unknown_convert(::Type{Ptr{T}}, A::MyMatrix)`.

The complete list of more specialised types is as follows:
```
julia> using ArrayLayouts, AbstractTrees

julia> AbstractTrees.children(x::Type) = subtypes(x)

julia> print_tree(AbstractStridedLayout)
AbstractStridedLayout
├─ AbstractDecreasingStrides
│  └─ AbstractRowMajor
│     ├─ DenseRowMajor
│     └─ RowMajor
├─ AbstractIncreasingStrides
│  ├─ AbstractColumnMajor
│  │  ├─ ColumnMajor
│  │  └─ DenseColumnMajor
│  ├─ DecreasingStrides
│  └─ IncreasingStrides
├─ StridedLayout
└─ UnitStride

julia> Base.show_supertypes(AbstractStridedLayout)
AbstractStridedLayout <: MemoryLayout <: Any
```
"""
AbstractStridedLayout

"""
    DenseColumnMajor()

is returned by `MemoryLayout(A)` if an array `A` has storage in memory
equivalent to an `Array`, so that `stride(A,1) == 1` and
`stride(A,i) ≡ size(A,i-1) * stride(A,i-1)` for `2 ≤ i ≤ ndims(A)`. In particular,
if `A` is a matrix then `strides(A) == `(1, size(A,1))`.

Arrays with `DenseColumnMajor` memory layout must conform to the `DenseArray` interface.
"""
DenseColumnMajor

"""
    ColumnMajor()

is returned by `MemoryLayout(A)` if an array `A` has storage in memory
as a column major array, so that `stride(A,1) == 1` and
`stride(A,i) ≥ size(A,i-1) * stride(A,i-1)` for `2 ≤ i ≤ ndims(A)`.

Arrays with `ColumnMajor` memory layout must conform to the `DenseArray` interface.
"""
ColumnMajor

"""
    IncreasingStrides()

is returned by `MemoryLayout(A)` if an array `A` has storage in memory
as a strided array with  increasing strides, so that `stride(A,1) ≥ 1` and
`stride(A,i) ≥ size(A,i-1) * stride(A,i-1)` for `2 ≤ i ≤ ndims(A)`.
"""
IncreasingStrides

"""
    DenseRowMajor()

is returned by `MemoryLayout(A)` if an array `A` has storage in memory
as a row major array with dense entries, so that `stride(A,ndims(A)) == 1` and
`stride(A,i) ≡ size(A,i+1) * stride(A,i+1)` for `1 ≤ i ≤ ndims(A)-1`. In particular,
if `A` is a matrix then `strides(A) == `(size(A,2), 1)`.
"""
DenseRowMajor

"""
    RowMajor()

is returned by `MemoryLayout(A)` if an array `A` has storage in memory
as a row major array, so that `stride(A,ndims(A)) == 1` and
stride(A,i) ≥ size(A,i+1) * stride(A,i+1)` for `1 ≤ i ≤ ndims(A)-1`.

If `A` is a matrix  with `RowMajor` memory layout, then
`transpose(A)` should return a matrix whose layout is `ColumnMajor`.
"""
RowMajor

"""
    DecreasingStrides()

is returned by `MemoryLayout(A)` if an array `A` has storage in memory
as a strided array with decreasing strides, so that `stride(A,ndims(A)) ≥ 1` and
stride(A,i) ≥ size(A,i+1) * stride(A,i+1)` for `1 ≤ i ≤ ndims(A)-1`.
"""
DecreasingStrides

"""
    StridedLayout()

is returned by `MemoryLayout(A)` if an array `A` has storage laid out at regular
offsets in memory. `Array`s with `StridedLayout` must conform to the `DenseArray` interface.
"""
StridedLayout

"""
    ScalarLayout()

is returned by `MemoryLayout(A)` if A is a scalar, which does not live in memory
"""
ScalarLayout

"""
    MemoryLayout(A)

specifies the layout in memory for an array `A`. When
you define a new `AbstractArray` type, you can choose to override
`MemoryLayout` to indicate how an array is stored in memory.
For example, if your matrix is column major with `stride(A,2) == size(A,1)`,
then override as follows:

    MemoryLayout(::MyMatrix) = DenseColumnMajor()

The default is `UnknownLayout()` to indicate that the layout
in memory is unknown.

Julia's internal linear algebra machinery will automatically (and invisibly)
dispatch to BLAS and LAPACK routines if the memory layout is compatible.
"""
@inline MemoryLayout(A) = MemoryLayout(typeof(A))

@inline MemoryLayout(::Type) = UnknownLayout()

@inline MemoryLayout(::Type{<:Number}) = ScalarLayout()
@inline MemoryLayout(::Type{<:DenseArray}) = DenseColumnMajor()

@inline MemoryLayout(::Type{<:ReinterpretArray{T,N,S,P}}) where {T,N,S,P} = reinterpretedlayout(MemoryLayout(P))
@inline reinterpretedlayout(::MemoryLayout) = UnknownLayout()
@inline reinterpretedlayout(::DenseColumnMajor) = DenseColumnMajor()


MemoryLayout(::Type{<:ReshapedArray{T,N,A,DIMS}}) where {T,N,A,DIMS} = reshapedlayout(MemoryLayout(A), DIMS)
@inline reshapedlayout(_, _) = UnknownLayout()
@inline reshapedlayout(::DenseColumnMajor, _) = DenseColumnMajor()


@inline MemoryLayout(A::Type{<:SubArray{T,N,P,I}}) where {T,N,P,I} =
    sublayout(MemoryLayout(P), I)
sublayout(_1, _2) = UnknownLayout()
sublayout(_1, _2, _3)= UnknownLayout()
sublayout(::DenseColumnMajor, ::Type{<:Tuple{<:Union{AbstractUnitRange{Int},Int,AbstractCartesianIndex}}}) =
    DenseColumnMajor()  # A[:] is DenseColumnMajor if A is DenseColumnMajor
sublayout(ml::AbstractColumnMajor, inds) = _column_sublayout1(ml, inds)
sublayout(::AbstractRowMajor, ::Type{<:Tuple{<:Any}}) =
    UnknownLayout()  # A[:] does not have any structure if A is AbstractRowMajor
sublayout(ml::AbstractRowMajor, inds) = _row_sublayout1(ml, tuple_type_reverse(inds))
sublayout(ml::AbstractStridedLayout, inds) = _strided_sublayout(ml, inds)

_column_sublayout1(::DenseColumnMajor, inds::Type{<:Tuple{I,Vararg{Int}}}) where I<:Union{Int,AbstractCartesianIndex} =
    DenseColumnMajor() # view(A,1,1,2) is a scalar, which we include in DenseColumnMajor
_column_sublayout1(::DenseColumnMajor, inds::Type{<:Tuple{I,Vararg{Int}}}) where I<:Slice =
    DenseColumnMajor() # view(A,:,1,2) is a DenseColumnMajor vector
_column_sublayout1(::DenseColumnMajor, inds::Type{<:Tuple{I,Vararg{Int}}}) where I<:AbstractUnitRange{Int} =
    DenseColumnMajor() # view(A,1:3,1,2) is a DenseColumnMajor vector
_column_sublayout1(par, inds::Type{<:Tuple{I,Vararg{Int}}}) where I<:Union{Int,AbstractCartesianIndex} =
    DenseColumnMajor() # view(A,1,1,2) is a scalar, which we include in DenseColumnMajor
_column_sublayout1(par, inds::Type{<:Tuple{I,Vararg{Int}}}) where I<:AbstractUnitRange{Int} =
    DenseColumnMajor() # view(A,1:3,1,2) is a DenseColumnMajor vector
_column_sublayout1(::DenseColumnMajor, inds::Type{<:Tuple{I,Vararg{Any}}}) where I<:Slice =
    _column_sublayout(DenseColumnMajor(), DenseColumnMajor(), tuple_type_tail(inds))
_column_sublayout1(par::DenseColumnMajor, inds::Type{<:Tuple{I,Vararg{Any}}}) where I<:AbstractUnitRange{Int} =
    _column_sublayout(par, ColumnMajor(), tuple_type_tail(inds))
_column_sublayout1(par, inds::Type{<:Tuple{I,Vararg{Any}}}) where I<:AbstractUnitRange{Int} =
    _column_sublayout(par, ColumnMajor(), tuple_type_tail(inds))
_column_sublayout1(par::DenseColumnMajor, inds::Type{<:Tuple{I,Vararg{Any}}}) where I<:Union{RangeIndex,AbstractCartesianIndex} =
    _column_sublayout(par, StridedLayout(), tuple_type_tail(inds))
_column_sublayout1(par, inds::Type{<:Tuple{I,Vararg{Any}}}) where I<:Union{RangeIndex,AbstractCartesianIndex} =
    _column_sublayout(par, StridedLayout(), tuple_type_tail(inds))
_column_sublayout1(par, inds) = UnknownLayout()
_column_sublayout(par, ret, ::Type{<:Tuple{}}) = ret
_column_sublayout(par, ret, ::Type{<:Tuple{I}}) where I = UnknownLayout()
_column_sublayout(::DenseColumnMajor, ::DenseColumnMajor, inds::Type{<:Tuple{I,Vararg{Int}}}) where I<:Union{AbstractUnitRange{Int},Int,AbstractCartesianIndex} =
    DenseColumnMajor() # A[:,1:3,1,2] is DenseColumnMajor if A is DenseColumnMajor
_column_sublayout(par::DenseColumnMajor, ::DenseColumnMajor, inds::Type{<:Tuple{I, Vararg{Int}}}) where I<:Slice =
    DenseColumnMajor()
_column_sublayout(par::DenseColumnMajor, ::DenseColumnMajor, inds::Type{<:Tuple{I, Vararg{Any}}}) where I<:Slice =
    _column_sublayout(par, DenseColumnMajor(), tuple_type_tail(inds))
_column_sublayout(par, ::AbstractColumnMajor, inds::Type{<:Tuple{I, Vararg{Any}}}) where I<:Union{AbstractUnitRange{Int},Int,AbstractCartesianIndex} =
    _column_sublayout(par, ColumnMajor(), tuple_type_tail(inds))
_column_sublayout(par, ::AbstractStridedLayout, inds::Type{<:Tuple{I, Vararg{Any}}}) where I<:Union{RangeIndex,AbstractCartesianIndex} =
    _column_sublayout(par, StridedLayout(), tuple_type_tail(inds))

_row_sublayout1(par, inds::Type{<:Tuple{I,Vararg{Int}}}) where I<:Union{Int,AbstractCartesianIndex} =
    DenseColumnMajor() # view(A,1,1,2) is a scalar, which we include in DenseColumnMajor
_row_sublayout1(::DenseRowMajor, inds::Type{<:Tuple{I,Vararg{Int}}}) where I<:Slice =
    DenseColumnMajor() # view(A,1,2,:) is a DenseColumnMajor vector
_row_sublayout1(par, inds::Type{<:Tuple{I,Vararg{Int}}}) where I<:AbstractUnitRange{Int} =
    DenseColumnMajor() # view(A,1,2,1:3) is a DenseColumnMajor vector
_row_sublayout1(::DenseRowMajor, inds::Type{<:Tuple{I,Vararg{Any}}}) where I<:Slice =
    _row_sublayout(DenseRowMajor(), DenseRowMajor(), tuple_type_tail(inds))
_row_sublayout1(par, inds::Type{<:Tuple{I,Vararg{Any}}}) where I<:AbstractUnitRange{Int} =
    _row_sublayout(par, RowMajor(), tuple_type_tail(inds))
_row_sublayout1(par, inds::Type{<:Tuple{I,Vararg{Any}}}) where I<:Union{RangeIndex,AbstractCartesianIndex} =
    _row_sublayout(par, StridedLayout(), tuple_type_tail(inds))
_row_sublayout1(par, inds) = UnknownLayout()
_row_sublayout(par, ret, ::Type{<:Tuple{}}) = ret
_row_sublayout(par, ret, ::Type{<:Tuple{I}}) where I = UnknownLayout()
_row_sublayout(::DenseRowMajor, ::DenseRowMajor, inds::Type{<:Tuple{I,Vararg{Int}}}) where I<:Union{AbstractUnitRange{Int},Int,AbstractCartesianIndex} =
    DenseRowMajor() # A[1,2,1:3,:] is DenseRowMajor if A is DenseRowMajor
_row_sublayout(par::DenseRowMajor, ::DenseRowMajor, inds::Type{<:Tuple{I, Vararg{Int}}}) where I<:Slice =
    DenseRowMajor()
_row_sublayout(par::DenseRowMajor, ::DenseRowMajor, inds::Type{<:Tuple{I, Vararg{Any}}}) where I<:Slice =
    _row_sublayout(par, DenseRowMajor(), tuple_type_tail(inds))
_row_sublayout(par::AbstractRowMajor, ::AbstractRowMajor, inds::Type{<:Tuple{I, Vararg{Any}}}) where I<:Union{AbstractUnitRange{Int},Int,AbstractCartesianIndex} =
    _row_sublayout(par, RowMajor(), tuple_type_tail(inds))
_row_sublayout(par::AbstractRowMajor, ::AbstractStridedLayout, inds::Type{<:Tuple{I, Vararg{Any}}}) where I<:Union{RangeIndex,AbstractCartesianIndex} =
    _row_sublayout(par, StridedLayout(), tuple_type_tail(inds))

_strided_sublayout(par, inds) = UnknownLayout()
_strided_sublayout(par, ::Type{<:Tuple{}}) = StridedLayout()
_strided_sublayout(par, inds::Type{<:Tuple{I, Vararg{Any}}}) where I<:Union{RangeIndex,AbstractCartesianIndex} =
    _strided_sublayout(par, tuple_type_tail(inds))

# MemoryLayout of transposed and adjoint matrices
struct ConjLayout{ML<:MemoryLayout} <: MemoryLayout end

conjlayout(_1, _2) = UnknownLayout()
conjlayout(::Type{<:Complex}, ::ConjLayout{ML}) where ML = ML()
conjlayout(::Type{<:Complex}, ::ML) where ML<:AbstractStridedLayout = ConjLayout{ML}()
conjlayout(::Type{<:Real}, M::MemoryLayout) = M


sublayout(::ConjLayout{ML}, t::Type{<:Tuple}) where ML = ConjLayout{typeof(sublayout(ML(), t))}()

"""
   DualLayout{ML<:MemoryLayout}()

represents a row-vector that should behave like a dual-vector, that is
multiplication times a column-vector returns a scalar.
"""
struct DualLayout{ML<:MemoryLayout} <: MemoryLayout end

MemoryLayout(::Type{Transpose{T,P}}) where {T,P} = transposelayout(MemoryLayout(P))
MemoryLayout(::Type{Adjoint{T,P}}) where {T,P} = adjointlayout(T, MemoryLayout(P))
MemoryLayout(::Type{AdjointAbsVec{T,P}}) where {T,P<:AbstractVector} = DualLayout{typeof(adjointlayout(T,MemoryLayout(P)))}()
MemoryLayout(::Type{TransposeAbsVec{T,P}}) where {T,P<:AbstractVector} = DualLayout{typeof(transposelayout(MemoryLayout(P)))}()
if isdefined(LinearAlgebra, :AdjointQ)
    MemoryLayout(::Type{LinearAlgebra.AdjointQ{T,P}}) where {T,P} = adjointlayout(T, MemoryLayout(P))
end

transposelayout(_) = UnknownLayout()
transposelayout(::StridedLayout) = StridedLayout()
transposelayout(::ColumnMajor) = RowMajor()
transposelayout(::RowMajor) = ColumnMajor()
transposelayout(::DenseColumnMajor) = DenseRowMajor()
transposelayout(::DenseRowMajor) = DenseColumnMajor()
transposelayout(::ConjLayout{ML}) where ML = ConjLayout{typeof(transposelayout(ML()))}()
adjointlayout(::Type{T}, ::DualLayout{ML}) where {T,ML} = adjointlayout(T, ML())
adjointlayout(::Type{T}, M::MemoryLayout) where T = transposelayout(conjlayout(T, M))
sublayout(::DualLayout{ML}, ::Type{<:Tuple{KR,JR}}) where {ML,KR<:Slice,JR} = DualLayout{typeof(sublayout(ML(),Tuple{KR,JR}))}()
sublayout(::DualLayout{ML}, INDS::Type) where ML = sublayout(ML(), INDS)

# try to maintain "type" of parent ie adjoint/transpose when materialising, even though layouts are equivalent
# to preserve special overloads
_dual_adjoint(a::SubArray{<:Any, 2, <:Any, <:Tuple{Slice,Any}}) = view(parent(a)', parentindices(a)[2])
_dual_transpose(a::SubArray{<:Any, 2, <:Any, <:Tuple{Slice,Any}}) = view(transpose(parent(a)), parentindices(a)[2])
sub_materialize(::DualLayout{ML}, A::AbstractMatrix{<:Real}) where ML = sub_materialize(adjointlayout(eltype(A), ML()), _dual_adjoint(A))'
sub_materialize(::DualLayout{ML}, A::AbstractMatrix) where ML<:ConjLayout = sub_materialize(adjointlayout(eltype(A), ML()), _dual_adjoint(A))'
sub_materialize(::DualLayout{ML}, A::AbstractMatrix) where ML = transpose(sub_materialize(transposelayout(ML()), _dual_transpose(A)))

copyto!_layout(dlay, ::DualLayout{ML}, dest::AbstractArray{T,N}, src::AbstractArray{V,N}) where {T,V,N,ML} =
    copyto!_layout(dlay, ML(), dest, src)

# Layouts of PermutedDimsArrays
"""
    UnitStride{D}()

is returned by `MemoryLayout(A)` for arrays of `ndims(A) >= 3` which have `stride(A,D) == 1`.

`UnitStride{1}` is weaker than `ColumnMajor` in that it does not demand that the other
strides are increasing, hence it is not a subtype of `AbstractIncreasingStrides`.
To ensure that `stride(A,1) == 1`, you may dispatch on `Union{UnitStride{1}, AbstractColumnMajor}`
to allow for both options. (With complex numbers, you may also need their `ConjLayout` versions.)

Likewise, both `UnitStride{ndims(A)}` and `AbstractRowMajor` have `stride(A, ndims(A)) == 1`.
"""
UnitStride

MemoryLayout(::Type{PermutedDimsArray{T,N,P,Q,S}}) where {T,N,P,Q,S} = permutelayout(MemoryLayout(S), Val(P))

permutelayout(::Any, perm) = UnknownLayout()
permutelayout(::StridedLayout, perm) = StridedLayout()
permutelayout(::ConjLayout{ML}, perm) where ML = ConjLayout{typeof(permutelayout(ML(), perm))}()

function permutelayout(layout::AbstractColumnMajor, ::Val{perm}) where {perm}
    issorted(perm) && return layout
    issorted(reverse(perm)) && return reverse(layout)
    D = sum(ntuple(dim -> perm[dim] == 1 ? dim : 0, length(perm)))
    return UnitStride{D}()
end
function permutelayout(layout::AbstractRowMajor, ::Val{perm}) where {perm}
    issorted(perm) && return layout
    issorted(reverse(perm)) && return reverse(layout)
    N = length(perm) # == ndims(A)
    D = sum(ntuple(dim -> perm[dim] == N ? dim : 0, N))
    return UnitStride{D}()
end
function permutelayout(layout::UnitStride{D0}, ::Val{perm}) where {D0, perm}
    D = sum(ntuple(dim -> perm[dim] == D0 ? dim : 0, length(perm)))
    return UnitStride{D}()
end
function permutelayout(layout::Union{IncreasingStrides,DecreasingStrides}, ::Val{perm}) where {perm}
    issorted(perm) && return layout
    issorted(reverse(perm)) && return reverse(layout)
    return StridedLayout()
end

reverse(::DenseRowMajor) = DenseColumnMajor()
reverse(::RowMajor) = ColumnMajor()
reverse(::DenseColumnMajor) = DenseRowMajor()
reverse(::ColumnMajor) = RowMajor()
reverse(::IncreasingStrides) = DecreasingStrides()
reverse(::DecreasingStrides) = IncreasingStrides()
reverse(::AbstractStridedLayout) = StridedLayout()


# MemoryLayout of Symmetric/Hermitian
"""
    SymmetricLayout{layout}()


is returned by `MemoryLayout(A)` if a matrix `A` has storage in memory
as a symmetrized version of `layout`, where the entries used are dictated by the
`uplo`, which can be `'U'` or `L'`.

A matrix that has memory layout `SymmetricLayout(layout, uplo)` must overrided
`symmetricdata(A)` to return a matrix `B` such that `MemoryLayout(B) == layout` and
`A[k,j] == B[k,j]` for `j ≥ k` if `uplo == 'U'` (`j ≤ k` if `uplo == 'L'`) and
`A[k,j] == B[j,k]` for `j < k` if `uplo == 'U'` (`j > k` if `uplo == 'L'`).
"""
struct SymmetricLayout{ML<:MemoryLayout} <: MemoryLayout end

"""
    HermitianLayout(layout, uplo)


is returned by `MemoryLayout(A)` if a matrix `A` has storage in memory
as a hermitianized version of `layout`, where the entries used are dictated by the
`uplo`, which can be `'U'` or `L'`.

A matrix that has memory layout `HermitianLayout(layout, uplo)` must overrided
`hermitiandata(A)` to return a matrix `B` such that `MemoryLayout(B) == layout` and
`A[k,j] == B[k,j]` for `j ≥ k` if `uplo == 'U'` (`j ≤ k` if `uplo == 'L'`) and
`A[k,j] == conj(B[j,k])` for `j < k` if `uplo == 'U'` (`j > k` if `uplo == 'L'`).
"""
struct HermitianLayout{ML<:MemoryLayout} <: MemoryLayout end

MemoryLayout(::Type{Hermitian{T,P}}) where {T,P} = hermitianlayout(T, MemoryLayout(P))
MemoryLayout(::Type{Symmetric{T,P}}) where {T,P} = symmetriclayout(MemoryLayout(P))
hermitianlayout(_1, _2) = UnknownLayout()
hermitianlayout(::Type{<:Complex}, ::ML) where ML<:AbstractColumnMajor = HermitianLayout{ML}()
hermitianlayout(::Type{<:Real}, ::ML) where ML<:AbstractColumnMajor = SymmetricLayout{ML}()
hermitianlayout(::Type{<:Complex}, ::ML) where ML<:AbstractRowMajor = HermitianLayout{ML}()
hermitianlayout(::Type{<:Real}, ::ML) where ML<:AbstractRowMajor = SymmetricLayout{ML}()
symmetriclayout(_1) = UnknownLayout()
symmetriclayout(::ML) where ML<:AbstractColumnMajor = SymmetricLayout{ML}()
symmetriclayout(::ML) where ML<:AbstractRowMajor = SymmetricLayout{ML}()
transposelayout(S::SymmetricLayout) = S
adjointlayout(::Type{T}, S::SymmetricLayout) where T<:Real = S
adjointlayout(::Type{T}, S::HermitianLayout) where T = S
sublayout(S::SymmetricLayout, ::Type{<:Tuple{<:Slice,<:Slice}}) = S
sublayout(S::HermitianLayout, ::Type{<:Tuple{<:Slice,<:Slice}}) = S

symmetricdata(A::Symmetric) = A.data
symmetricdata(A::Hermitian{<:Real}) = A.data
symmetricdata(V::SubArray{<:Any, 2, <:Any, <:Tuple{<:Slice,<:Slice}}) = symmetricdata(parent(V))
symmetricdata(V::Adjoint{<:Real}) = symmetricdata(parent(V))
symmetricdata(V::Transpose) = symmetricdata(parent(V))
hermitiandata(A::Hermitian) = A.data
hermitiandata(V::SubArray{<:Any, 2, <:Any, <:Tuple{<:Slice,<:Slice}}) = hermitiandata(parent(V))
hermitiandata(V::Adjoint) = hermitiandata(parent(V))
hermitiandata(V::Transpose{<:Real}) = hermitiandata(parent(V))

symmetricuplo(A::Symmetric) = A.uplo
symmetricuplo(A::Hermitian) = A.uplo
symmetricuplo(A::Adjoint) = symmetricuplo(parent(A))
symmetricuplo(A::Transpose) = A.uplo
symmetricuplo(A::SubArray{<:Any, 2, <:Any, <:Tuple{<:Slice,<:Slice}}) = symmetricuplo(parent(A))

# MemoryLayout of triangular matrices
struct TriangularLayout{UPLO,UNIT,ML} <: MemoryLayout end


"""
    LowerTriangularLayout(layout)

is returned by `MemoryLayout(A)` if a matrix `A` has storage in memory
equivalent to a `LowerTriangular(B)` where `B` satisfies `MemoryLayout(B) == layout`.

A matrix that has memory layout `LowerTriangularLayout(layout)` must overrided
`triangulardata(A)` to return a matrix `B` such that `MemoryLayout(B) == layout` and
`A[k,j] ≡ zero(eltype(A))` for `j > k` and
`A[k,j] ≡ B[k,j]` for `j ≤ k`.

Moreover, `transpose(A)` and `adjoint(A)` must return a matrix that has memory
layout `UpperTriangularLayout`.
"""
const LowerTriangularLayout{ML} = TriangularLayout{'L','N',ML}

"""
    UnitLowerTriangularLayout(ML::MemoryLayout)

is returned by `MemoryLayout(A)` if a matrix `A` has storage in memory
equivalent to a `UnitLowerTriangular(B)` where `B` satisfies `MemoryLayout(B) == layout`.

A matrix that has memory layout `UnitLowerTriangularLayout(layout)` must overrided
`triangulardata(A)` to return a matrix `B` such that `MemoryLayout(B) == layout` and
`A[k,j] ≡ zero(eltype(A))` for `j > k`,
`A[k,j] ≡ one(eltype(A))` for `j == k`,
`A[k,j] ≡ B[k,j]` for `j < k`.

Moreover, `transpose(A)` and `adjoint(A)` must return a matrix that has memory
layout `UnitUpperTriangularLayout`.
"""
const UnitLowerTriangularLayout{ML} = TriangularLayout{'L','U',ML}

"""
    UpperTriangularLayout(ML::MemoryLayout)

is returned by `MemoryLayout(A)` if a matrix `A` has storage in memory
equivalent to a `UpperTriangularLayout(B)` where `B` satisfies `MemoryLayout(B) == ML`.

A matrix that has memory layout `UpperTriangularLayout(layout)` must overrided
`triangulardata(A)` to return a matrix `B` such that `MemoryLayout(B) == layout` and
`A[k,j] ≡ B[k,j]` for `j ≥ k` and
`A[k,j] ≡ zero(eltype(A))` for `j < k`.

Moreover, `transpose(A)` and `adjoint(A)` must return a matrix that has memory
layout `LowerTriangularLayout`.
"""
const UpperTriangularLayout{ML} = TriangularLayout{'U','N',ML}

"""
    UnitUpperTriangularLayout(ML::MemoryLayout)

is returned by `MemoryLayout(A)` if a matrix `A` has storage in memory
equivalent to a `UpperTriangularLayout(B)` where `B` satisfies `MemoryLayout(B) == ML`.

A matrix that has memory layout `UnitUpperTriangularLayout(layout)` must overrided
`triangulardata(A)` to return a matrix `B` such that `MemoryLayout(B) == layout` and
`A[k,j] ≡ B[k,j]` for `j > k`,
`A[k,j] ≡ one(eltype(A))` for `j == k`,
`A[k,j] ≡ zero(eltype(A))` for `j < k`.

Moreover, `transpose(A)` and `adjoint(A)` must return a matrix that has memory
layout `UnitLowerTriangularLayout`.
"""
const UnitUpperTriangularLayout{ML} = TriangularLayout{'U','U',ML}


MemoryLayout(A::Type{UpperTriangular{T,P}}) where {T,P} = triangularlayout(UpperTriangularLayout, MemoryLayout(P))
MemoryLayout(A::Type{UnitUpperTriangular{T,P}}) where {T,P} = triangularlayout(UnitUpperTriangularLayout, MemoryLayout(P))
MemoryLayout(A::Type{LowerTriangular{T,P}}) where {T,P} = triangularlayout(LowerTriangularLayout, MemoryLayout(P))
MemoryLayout(A::Type{UnitLowerTriangular{T,P}}) where {T,P} = triangularlayout(UnitLowerTriangularLayout, MemoryLayout(P))
triangularlayout(::Type{Tri}, ::MemoryLayout) where Tri = Tri{UnknownLayout}()
triangularlayout(::Type{Tri}, ::ML) where {Tri, ML<:AbstractColumnMajor} = Tri{ML}()
triangularlayout(::Type{Tri}, ::ML) where {Tri, ML<:AbstractRowMajor} = Tri{ML}()
triangularlayout(::Type{Tri}, ::ML) where {Tri, ML<:ConjLayout{<:AbstractRowMajor}} = Tri{ML}()
sublayout(layout::TriangularLayout, ::Type{<:Tuple{<:Union{Slice,OneTo},<:Union{Slice,OneTo}}}) = layout
conjlayout(::Type{<:Complex}, ::TriangularLayout{UPLO,UNIT,ML}) where {UPLO,UNIT,ML} =
    TriangularLayout{UPLO,UNIT,ConjLayout{ML}}()

for (TriLayout, TriLayoutTrans) in ((UpperTriangularLayout,     LowerTriangularLayout),
                                    (UnitUpperTriangularLayout, UnitLowerTriangularLayout),
                                    (LowerTriangularLayout,     UpperTriangularLayout),
                                    (UnitLowerTriangularLayout, UnitUpperTriangularLayout))
    @eval transposelayout(::$TriLayout{ML}) where ML = $TriLayoutTrans{typeof(transposelayout(ML()))}()
end

triangulardata(A::AbstractTriangular) = parent(A)
triangulardata(A::Adjoint) = Adjoint(triangulardata(parent(A)))
triangulardata(A::Transpose) = Transpose(triangulardata(parent(A)))
triangulardata(A::SubArray{<:Any,2,<:Any,<:Tuple{<:Union{Slice,OneTo},<:Union{Slice,OneTo}}}) =
    view(triangulardata(parent(A)), parentindices(A)...)

###
# Fill
####
abstract type AbstractFillLayout <: MemoryLayout end
struct FillLayout <: AbstractFillLayout end
struct ZerosLayout <: AbstractFillLayout end
struct OnesLayout <: AbstractFillLayout end
struct EyeLayout <: MemoryLayout end

MemoryLayout(::Type{<:AbstractFill}) = FillLayout()
MemoryLayout(::Type{<:Zeros}) = ZerosLayout()
MemoryLayout(::Type{<:Ones}) = OnesLayout()

# all sub arrays are same
sublayout(L::AbstractFillLayout, inds::Type) = L
reshapedlayout(L::AbstractFillLayout, _) = L
adjointlayout(::Type, L::AbstractFillLayout) = L
transposelayout(L::AbstractFillLayout) = L

copyto!_layout(_, ::AbstractFillLayout, dest::AbstractArray{<:Any,N}, src::AbstractArray{<:Any,N}) where N =
    fill!(dest, getindex_value(src))

sub_materialize(::AbstractFillLayout, V, ax) = Fill(getindex_value(V), ax)
sub_materialize(::ZerosLayout, V, ax) = Zeros{eltype(V)}(ax)
sub_materialize(::OnesLayout, V, ax) = Ones{eltype(V)}(ax)

abstract type AbstractBandedLayout <: MemoryLayout end
abstract type AbstractTridiagonalLayout <: AbstractBandedLayout end

struct DiagonalLayout{ML} <: AbstractBandedLayout end
struct BidiagonalLayout{DV,EV} <: AbstractBandedLayout end
struct SymTridiagonalLayout{DV,EV} <: AbstractTridiagonalLayout end
struct TridiagonalLayout{DL,D,DU} <: AbstractTridiagonalLayout end

bidiagonallayout(dv, ev) = BidiagonalLayout{UnknownLayout,UnknownLayout}()
tridiagonallayout(dl, d, du) = TridiagonalLayout{UnknownLayout,UnknownLayout,UnknownLayout}()

symtridiagonallayout(d, ev) = SymTridiagonalLayout{UnknownLayout,UnknownLayout}()
bidiagonallayout(d) = bidiagonallayout(d, d)
tridiagonallayout(d) = tridiagonallayout(d,d,d)
symtridiagonallayout(d) = symtridiagonallayout(d,d)
diagonallayout(_) = DiagonalLayout{UnknownLayout}()
diagonal(d::AbstractVector) = Diagonal(d) # support non-array diagonal objects like QuasiDiagonal

diagonallayout(::Lay) where Lay<:Union{AbstractStridedLayout, AbstractFillLayout} = DiagonalLayout{Lay}()
bidiagonallayout(::Lay, ::Lay) where Lay<:Union{AbstractStridedLayout, AbstractFillLayout} = BidiagonalLayout{Lay,Lay}()
tridiagonallayout(::Lay, ::Lay, ::Lay) where Lay<:Union{AbstractStridedLayout, AbstractFillLayout} = TridiagonalLayout{Lay,Lay,Lay}()
symtridiagonallayout(::Lay, ::Lay) where Lay<:Union{AbstractStridedLayout, AbstractFillLayout} = SymTridiagonalLayout{Lay,Lay}()


MemoryLayout(D::Type{Diagonal{T,P}}) where {T,P} = diagonallayout(MemoryLayout(P))
MemoryLayout(::Type{Bidiagonal{T,V}}) where {T,V} = bidiagonallayout(MemoryLayout(V))
MemoryLayout(::Type{SymTridiagonal{T,P}}) where {T,P} = symtridiagonallayout(MemoryLayout(P))
MemoryLayout(::Type{Tridiagonal{T,P}}) where {T,P} = tridiagonallayout(MemoryLayout(P))

bidiagonaluplo(A::Bidiagonal) = A.uplo
bidiagonaluplo(A::AdjOrTrans) = bidiagonaluplo(parent(A)) == 'L' ? 'U' : 'L'

diagonaldata(D::Diagonal) = parent(D)
diagonaldata(D::Bidiagonal) = D.dv
diagonaldata(D::SymTridiagonal) = D.dv
diagonaldata(D::Tridiagonal) = D.d

supdiagonaldata(D::Bidiagonal) = D.uplo == 'U' ? D.ev : throw(ArgumentError(LazyString(D, " is lower-bidiagonal")))
subdiagonaldata(D::Bidiagonal) = D.uplo == 'L' ? D.ev : throw(ArgumentError(LazyString(D, " is upper-bidiagonal")))

supdiagonaldata(D::SymTridiagonal) = D.ev
subdiagonaldata(D::SymTridiagonal) = D.ev

subdiagonaldata(D::Tridiagonal) = D.dl
supdiagonaldata(D::Tridiagonal) = D.du

transposelayout(ml::DiagonalLayout) = ml
transposelayout(ml::BidiagonalLayout) = ml
transposelayout(ml::SymTridiagonalLayout) = ml
transposelayout(ml::TridiagonalLayout) = ml
transposelayout(ml::ConjLayout{DiagonalLayout}) = ml

triangularlayout(::Type{<:TriangularLayout{'L','N'}}, ::TridiagonalLayout{DL,D,DU}) where {DL,D,DU} = BidiagonalLayout{D,DL}()
triangularlayout(::Type{<:TriangularLayout{'U','N'}}, ::TridiagonalLayout{DL,D,DU}) where {DL,D,DU} = BidiagonalLayout{D,DU}()
triangularlayout(::Type{<:TriangularLayout{'L','N'}}, ::TridiagonalLayout{FillLayout,FillLayout,FillLayout}) = BidiagonalLayout{FillLayout,FillLayout}()
triangularlayout(::Type{<:TriangularLayout{'U','N'}}, ::TridiagonalLayout{FillLayout,FillLayout,FillLayout}) = BidiagonalLayout{FillLayout,FillLayout}()
triangularlayout(::Type{<:TriangularLayout{UPLO,'U'}}, ::TridiagonalLayout{FillLayout,FillLayout,FillLayout}) where UPLO = BidiagonalLayout{FillLayout,FillLayout}()

bidiagonaluplo(::Union{UpperTriangular,UnitUpperTriangular}) = 'U'
bidiagonaluplo(::Union{LowerTriangular,UnitLowerTriangular}) = 'L'
diagonaldata(U::Union{UnitUpperTriangular{T},UnitLowerTriangular{T}}) where T = Ones{T}(size(U,1))
diagonaldata(U::Union{UpperTriangular{T},LowerTriangular{T}}) where T = diagonaldata(triangulardata(U))
supdiagonaldata(U::Union{UnitUpperTriangular,UpperTriangular}) = supdiagonaldata(triangulardata(U))
subdiagonaldata(U::Union{UnitLowerTriangular,LowerTriangular}) = subdiagonaldata(triangulardata(U))

adjointlayout(::Type{<:Real}, ml::SymTridiagonalLayout) = ml
adjointlayout(::Type{<:Real}, ::TridiagonalLayout{DL,D,DU}) where {DL,D,DU} = TridiagonalLayout{DU,D,DL}()
adjointlayout(::Type{<:Real}, ml::BidiagonalLayout) = ml

symmetriclayout(B::BidiagonalLayout{DV,EV}) where {DV,EV} = SymTridiagonalLayout{DV,EV}()
hermitianlayout(::Type{<:Real}, B::BidiagonalLayout{DV,EV}) where {DV,EV} = SymTridiagonalLayout{DV,EV}()
hermitianlayout(_, B::BidiagonalLayout) = HermitianLayout{typeof(B)}()

diagonaldata(D::Transpose) = diagonaldata(parent(D))
subdiagonaldata(D::Transpose) = supdiagonaldata(parent(D))
supdiagonaldata(D::Transpose) = subdiagonaldata(parent(D))

diagonaldata(D::Adjoint{<:Real}) = diagonaldata(parent(D))
subdiagonaldata(D::Adjoint{<:Real}) = supdiagonaldata(parent(D))
supdiagonaldata(D::Adjoint{<:Real}) = subdiagonaldata(parent(D))

diagonaldata(S::HermOrSym) = diagonaldata(parent(S))
subdiagonaldata(S::HermOrSym) = symmetricuplo(S) == 'L' ? subdiagonaldata(parent(S)) : supdiagonaldata(parent(S))
supdiagonaldata(S::HermOrSym) = symmetricuplo(S) == 'L' ? subdiagonaldata(parent(S)) : supdiagonaldata(parent(S))





rowsupport(_, A, k) = axes(A,2)
"""
    rowsupport(A, k)

Return an iterator containing the column indices of the possible non-zero entries in the `k`-th row of `A`.
"""
rowsupport(A, k) = rowsupport(MemoryLayout(A), A, k)
"""
    rowsupport(A)

Return an iterator containing the column indices of the possible non-zero entries in `A`.
"""
rowsupport(A) = rowsupport(A, axes(A,1))

colsupport(_, A, j) = axes(A,1)


"""
    colsupport(A, j)

Return an iterator containing the row indices of the possible non-zero entries in the `j`-th column of `A`.
"""
colsupport(A, j) = colsupport(MemoryLayout(A), A, j)
"""
    colsupport(A)

Return an iterator containing the row indices of the possible non-zero entries in `A`.
"""
colsupport(A) = colsupport(A, axes(A,2))

# TODO: generalise to other subarrays
function colsupport(A::SubArray{<:Any,N,<:Any,<:Tuple{Slice,AbstractVector}}, j) where N
    _, jr = parentindices(A)
    colsupport(parent(A), jr[j])
end

rowsupport(::ZerosLayout, A, _) = 1:0
colsupport(::ZerosLayout, A, _) = 1:0

colsupport(::UnknownLayout, A::OneElement{<:Any,1}, _) =
    intersect(axes(A,1), A.ind[1]:A.ind[1])
function colsupport(::UnknownLayout, A::OneElement{<:Any,2}, j)
    intersect(axes(A,1), range(A.ind[1], length = Int(A.ind[2] ∈ j)))
end
function rowsupport(::UnknownLayout, A::OneElement{<:Any,2}, k)
    intersect(axes(A,2), range(A.ind[2], length = Int(A.ind[1] ∈ k)))
end

rowsupport(::DiagonalLayout, _, k) = k
colsupport(::DiagonalLayout, _, j) = j

function colsupport(::BidiagonalLayout, A, j)
    isempty(j) && return 1:0
    bidiagonaluplo(A) == 'L' ? (minimum(j):min(size(A,1),maximum(j)+1)) : (max(minimum(j)-1,1):maximum(j))
end
function rowsupport(::BidiagonalLayout, A, j)
    isempty(j) && return 1:0
    bidiagonaluplo(A) == 'U' ? (minimum(j):min(size(A,2),maximum(j)+1)) : (max(minimum(j)-1,1):maximum(j))
end

function colsupport(::AbstractTridiagonalLayout, A, j)
	isempty(j) && return 1:0
	max(minimum(j)-1,1):min(size(A,1),maximum(j)+1)
end
function rowsupport(::AbstractTridiagonalLayout, A, j)
	isempty(j) && return 1:0
	max(minimum(j)-1,1):min(size(A,2),maximum(j)+1)
end

function colsupport(::SymmetricLayout, A, j)
    if symmetricuplo(A) == 'U'
        first(colsupport(symmetricdata(A),j)):last(rowsupport(symmetricdata(A),j))
    else
        first(rowsupport(symmetricdata(A),j)):last(colsupport(symmetricdata(A),j))
    end
end
function colsupport(::HermitianLayout, A, j)
    if symmetricuplo(A) == 'U'
        first(colsupport(hermitiandata(A),j)):last(rowsupport(hermitiandata(A),j))
    else
        first(rowsupport(hermitiandata(A),j)):last(colsupport(hermitiandata(A),j))
    end
end
rowsupport(::Union{SymmetricLayout,HermitianLayout}, A, j) = colsupport(A, j)



function _sym_axes(A)
    ax = axes(parent(A),2)
    (ax, ax)
end

###
# axes overloads to support block indexing
###

axes(A::HermOrSym{<:Any,<:LayoutMatrix}) = _sym_axes(A)
axes(A::HermOrSym{<:Any,<:SubArray{<:Any,2,<:LayoutMatrix}}) = _sym_axes(A)
axes(A::UpperOrLowerTriangular{<:Any,<:LayoutMatrix}) = axes(parent(A))
axes(A::UpperOrLowerTriangular{<:Any,<:SubArray{<:Any,2,<:LayoutMatrix}}) = axes(parent(A))

function axes(D::Diagonal{<:Any,<:LayoutVector})
    a = axes(parent(D),1)
    (a,a)
end
