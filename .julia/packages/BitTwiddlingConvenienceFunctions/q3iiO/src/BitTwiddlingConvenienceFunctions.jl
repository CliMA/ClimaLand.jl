module BitTwiddlingConvenienceFunctions

using Static

intlog2(N::I) where {I <: Integer} = (8sizeof(I) - one(I) - leading_zeros(N)) % I
intlog2(::Type{T}) where {T} = intlog2(sizeof(T))
nextpow2(W::T) where {T<:Base.BitInteger} = (one(T) << ((T(8sizeof(W)) - leading_zeros((W - one(T))))%Unsigned))
prevpow2(W::T) where {T<:Base.BitInteger} = (one(T) << ((((T(8sizeof(W))) - one(T)) - leading_zeros(W))%Unsigned))
@generated nextpow2(::StaticInt{N}) where {N} = Expr(:call, Expr(:curly, :StaticInt, nextpow2(N)))
@generated prevpow2(::StaticInt{N}) where {N} = Expr(:call, Expr(:curly, :StaticInt, prevpow2(N)))
@generated intlog2(::StaticInt{N}) where {N} = Expr(:call, Expr(:curly, :StaticInt, intlog2(N)))


end
