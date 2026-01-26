module StaticArrayInterfaceStaticArraysExt

using StaticArrayInterface
using StaticArrayInterface.Static
using Static: StaticInt

if isdefined(Base, :get_extension) 
    using StaticArrays
    using LinearAlgebra
else 
    using ..StaticArrays
    using ..LinearAlgebra
end

const CanonicalInt = Union{Int,StaticInt}

function Static.OptionallyStaticUnitRange(::StaticArrays.SOneTo{N}) where {N}
    Static.OptionallyStaticUnitRange(StaticInt(1), StaticInt(N))
end
StaticArrayInterface.known_first(::Type{<:StaticArrays.SOneTo}) = 1
StaticArrayInterface.known_last(::Type{StaticArrays.SOneTo{N}}) where {N} = N
StaticArrayInterface.known_length(::Type{StaticArrays.SOneTo{N}}) where {N} = N
StaticArrayInterface.known_length(::Type{StaticArrays.Length{L}}) where {L} = L
function StaticArrayInterface.known_length(::Type{A}) where {A<:StaticArrays.StaticArray}
    StaticArrayInterface.known_length(StaticArrays.Length(A))
end

@inline StaticArrayInterface.static_length(x::StaticArrays.StaticArray) = Static.maybe_static(StaticArrayInterface.known_length, Base.length, x)
StaticArrayInterface.device(::Type{<:StaticArrays.MArray}) = StaticArrayInterface.CPUPointer()
StaticArrayInterface.device(::Type{<:StaticArrays.SArray}) = StaticArrayInterface.CPUTuple()
StaticArrayInterface.contiguous_axis(::Type{<:StaticArrays.StaticArray}) = StaticInt{1}()
StaticArrayInterface.contiguous_batch_size(::Type{<:StaticArrays.StaticArray}) = StaticInt{0}()
function StaticArrayInterface.stride_rank(::Type{T}) where {N,T<:StaticArray{<:Any,<:Any,N}}
    ntuple(static, StaticInt(N))
end
function StaticArrayInterface.dense_dims(::Type{<:StaticArray{S,T,N}}) where {S,T,N}
    StaticArrayInterface._all_dense(Val(N))
end
StaticArrayInterface.defines_strides(::Type{<:StaticArrays.SArray}) = true
StaticArrayInterface.defines_strides(::Type{<:StaticArrays.MArray}) = true

@generated function StaticArrayInterface.axes_types(::Type{<:StaticArrays.StaticArray{S}}) where {S}
    Tuple{[StaticArrays.SOneTo{s} for s in S.parameters]...}
end
@generated function StaticArrayInterface.static_size(A::StaticArrays.StaticArray{S}) where {S}
    t = Expr(:tuple)
    Sp = S.parameters
    for n = 1:length(Sp)
        push!(t.args, Expr(:call, Expr(:curly, :StaticInt, Sp[n])))
    end
    return t
end
@generated function StaticArrayInterface.static_strides(A::StaticArrays.StaticArray{S}) where {S}
    t = Expr(:tuple, Expr(:call, Expr(:curly, :StaticInt, 1)))
    Sp = S.parameters
    x = 1
    for n = 1:(length(Sp)-1)
        push!(t.args, Expr(:call, Expr(:curly, :StaticInt, (x *= Sp[n]))))
    end
    return t
end
if StaticArrays.SizedArray{Tuple{8,8},Float64,2,2} isa UnionAll
    @inline StaticArrayInterface.static_strides(B::StaticArrays.SizedArray{S,T,M,N,A}) where {S,T,M,N,A<:SubArray} = StaticArrayInterface.static_strides(B.data)
    StaticArrayInterface.parent_type(::Type{<:StaticArrays.SizedArray{S,T,M,N,A}}) where {S,T,M,N,A} = A
else
    StaticArrayInterface.parent_type(::Type{<:StaticArrays.SizedArray{S,T,M,N}}) where {S,T,M,N} = Array{T,N}
end

end # module
