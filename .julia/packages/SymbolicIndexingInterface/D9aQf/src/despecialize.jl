struct __InternalInvalidator1 <: Real end
struct __InternalInvalidator2 <: Real end
struct __InternalInvalidator3 <: AbstractVector{Real} end
symbolic_type(::Type{__InternalInvalidator1}) = ScalarSymbolic()
symbolic_type(::Type{__InternalInvalidator2}) = ArraySymbolic()
symbolic_type(::Type{__InternalInvalidator3}) = ScalarSymbolic()
symbolic_type(::Type{Int}) = NotSymbolic()
symbolic_type(::Type{UInt}) = NotSymbolic()
symbolic_type(::Type{Float64}) = NotSymbolic()
symbolic_type(::Type{Vector{Int}}) = NotSymbolic()
symbolic_type(::Type{Matrix{Float64}}) = NotSymbolic()
symbolic_type(::Type{Array{UInt, 3}}) = NotSymbolic()
