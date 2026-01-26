module SIMDTypes

# using Static: StaticInt

struct Bit; data::Bool; end # Dummy for Ptr
# @inline Base.convert(::Type{Bool}, b::Bit) = getfield(b, :data)

const FloatingTypes = Union{Float16,Float32,Float64}
const SignedHW = Union{Int8,Int16,Int32,Int64}
const UnsignedHW = Union{UInt8,UInt16,UInt32,UInt64}
const IntegerTypesHW = Union{SignedHW,UnsignedHW}
# const IntegerTypes = Union{StaticInt,IntegerTypesHW}

const NativeTypesExceptBitandFloat16 = Union{Bool,Base.HWReal}
const NativeTypesExceptBit = Union{Bool,Base.HWReal,Float16}
const NativeTypesExceptFloat16 = Union{Bool,Base.HWReal,Bit}
const NativeTypes = Union{NativeTypesExceptBit, Bit}

const _Vec{W,T<:Number} = NTuple{W,Core.VecElement{T}}

end
