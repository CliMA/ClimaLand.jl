module Bits

export AbstractBits, asbits

abstract type AbstractBits end

function asuint end
function asbits end

for T in [UInt8, UInt16, UInt32, UInt64, UInt128]
    C = :(Base.$(nameof(T)))
    nbits = sizeof(T) * 8
    B = Symbol("Bits$(nbits)")
    @eval begin
        export $B
        primitive type $B <: AbstractBits $nbits end
        asbits(::Type{$T}) = $B
        asuint(::Type{$B}) = $T
        $B(x::$T) = reinterpret($B, x)
        $C(x::$B) = reinterpret($T, x)
    end
end

Base.rand(::Type{B}) where {B<:AbstractBits} = B(rand(asuint(B)))

end  # module
