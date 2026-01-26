# Infinite Arrays implementation from
# https://github.com/JuliaLang/julia/blob/master/test/testhelpers/InfiniteArrays.jl
module InfiniteArrays
    using Infinities, LinearAlgebra, Random
    using ArrayLayouts: ArrayLayouts, LayoutVector, LayoutMatrix, Mul, DenseColumnMajor
    export OneToInf,
        InfSymTridiagonal,
        InfTridiagonal,
        InfBidiagonal,
        InfUnitUpperTriangular,
        InfUnitLowerTriangular,
        InfUpperTriangular,
        InfLowerTriangular,
        InfDiagonal 

    abstract type AbstractInfUnitRange{T<:Real} <: AbstractUnitRange{T} end
    Base.length(r::AbstractInfUnitRange) = ℵ₀
    Base.size(r::AbstractInfUnitRange) = (ℵ₀,)
    Base.last(r::AbstractInfUnitRange) = ℵ₀
    Base.axes(r::AbstractInfUnitRange) = (OneToInf(),)

    Base.IteratorSize(::Type{<:AbstractInfUnitRange}) = Base.IsInfinite()
    Base.sum(r::AbstractInfUnitRange) = last(r)

    """
        OneToInf(n)
    Define an `AbstractInfUnitRange` that behaves like `1:∞`, with the added
    distinction that the limits are guaranteed (by the type system) to
    be 1 and ∞.
    """
    struct OneToInf{T<:Integer} <: AbstractInfUnitRange{T} end

    OneToInf() = OneToInf{Int}()

    Base.axes(r::OneToInf) = (r,)
    Base.first(r::OneToInf{T}) where {T} = oneunit(T)
    Base.oneto(::InfiniteCardinal{0}) = OneToInf()


    struct InfUnitRange{T<:Real} <: AbstractInfUnitRange{T}
        start::T
    end
    Base.first(r::InfUnitRange) = r.start
    InfUnitRange(a::InfUnitRange) = a
    InfUnitRange{T}(a::AbstractInfUnitRange) where T<:Real = InfUnitRange{T}(first(a))
    InfUnitRange(a::AbstractInfUnitRange{T}) where T<:Real = InfUnitRange{T}(first(a))
    Base.:(:)(start::T, stop::InfiniteCardinal{0}) where {T<:Integer} = InfUnitRange{T}(start)
    function Base.getindex(v::InfUnitRange{T}, i::Integer) where T
        @boundscheck i > 0 || throw(BoundsError(v, i))
        convert(T, first(v) + i - 1)
    end
    function Base.getindex(v::InfUnitRange{T}, i::AbstractUnitRange{<:Integer}) where T
        @boundscheck checkbounds(v, i)
        v[first(i)]:v[last(i)]
    end
    function Base.getindex(v::InfUnitRange{T}, i::AbstractInfUnitRange{<:Integer}) where T
        @boundscheck checkbounds(v, first(i))
        v[first(i)]:ℵ₀
    end

    ## Methods for testing infinite arrays 
    struct InfVec{RNG} <: LayoutVector{Float64} # show is broken for InfVec
        rng::RNG
        data::Vector{Float64}
    end
    InfVec() = InfVec(copy(Random.seed!(Random.default_rng(), rand(UInt64))), Float64[])
    function resizedata!(v::InfVec, i)
        n = length(v.data)
        i ≤ n && return v
        resize!(v.data, i)
        for j in (n+1):i
            v[j] = rand(v.rng)
        end
        return v
    end
    Base.getindex(v::InfVec, i::Int) = (resizedata!(v, i); v.data[i])
    Base.setindex!(v::InfVec, r, i::Int) = setindex!(v.data, r, i)
    Base.size(v::InfVec) = (ℵ₀,)
    Base.axes(v::InfVec) = (OneToInf(),)
    ArrayLayouts.MemoryLayout(::Type{<:InfVec}) = DenseColumnMajor()
    Base.similar(v::InfVec, ::Type{T}, ::Tuple{OneToInf{Int}}) where {T} = InfVec()
    Base.copy(v::InfVec) = InfVec(copy(v.rng), copy(v.data))

    struct InfMat{RNG} <: LayoutMatrix{Float64} # show is broken for InfMat
        vec::InfVec{RNG}
    end
    InfMat() = InfMat(InfVec())
    function diagtrav_idx(i, j)
        band = i + j - 1
        nelm = (band * (band - 1)) ÷ 2
        return nelm + i
    end
    Base.getindex(A::InfMat, i::Int, j::Int) = A.vec[diagtrav_idx(i, j)]
    Base.setindex!(A::InfMat, r, i::Int, j::Int) = setindex!(A.vec, r, diagtrav_idx(i, j))
    Base.size(A::InfMat) = (ℵ₀, ℵ₀)
    Base.axes(v::InfMat) = (OneToInf(), OneToInf())
    ArrayLayouts.MemoryLayout(::Type{<:InfMat}) = DenseColumnMajor()
    Base.copy(A::InfMat) = InfMat(copy(A.vec))

    const InfSymTridiagonal = SymTridiagonal{Float64,<:InfVec}
    const InfTridiagonal = Tridiagonal{Float64,<:InfVec}
    const InfBidiagonal = Bidiagonal{Float64,<:InfVec}
    const InfUnitUpperTriangular = UnitUpperTriangular{Float64,<:InfMat}
    const InfUnitLowerTriangular = UnitLowerTriangular{Float64,<:InfMat}
    const InfUpperTriangular = UpperTriangular{Float64,<:InfMat}
    const InfLowerTriangular = LowerTriangular{Float64,<:InfMat}
    const InfDiagonal = Diagonal{Float64,<:InfVec}
    InfSymTridiagonal() = SymTridiagonal(InfVec(), InfVec())
    InfTridiagonal() = Tridiagonal(InfVec(), InfVec(), InfVec())
    InfBidiagonal(uplo) = Bidiagonal(InfVec(), InfVec(), uplo)
    InfUnitUpperTriangular() = UnitUpperTriangular(InfMat())
    InfUnitLowerTriangular() = UnitLowerTriangular(InfMat())
    InfUpperTriangular() = UpperTriangular(InfMat())
    InfLowerTriangular() = LowerTriangular(InfMat())
    InfDiagonal() = Diagonal(InfVec())
    Base.copy(D::InfDiagonal) = Diagonal(copy(D.diag))

    # Without LazyArrays we have no access to the lazy machinery, so we must define copy(::Mul) to leave mul(A, B) as a lazy Mul(A, B)
    const InfNamedMatrix = Union{InfSymTridiagonal,InfTridiagonal,InfBidiagonal,
        InfUnitUpperTriangular,InfUnitLowerTriangular,
        InfUpperTriangular,InfLowerTriangular,
        InfDiagonal}
    const InfMul{L1,L2} = Mul{L1,L2,<:InfNamedMatrix,<:InfNamedMatrix}
    Base.copy(M::InfMul{L1,L2}) where {L1,L2} = Mul{L1,L2}(copy(M.A), copy(M.B))
end
