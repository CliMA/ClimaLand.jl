function Base.:*(A::Matrix{S}, b::Vector{T}) where {S,T}
    c = Vector{U}(undef, size(A, 1)) # What is U ?
    return mul_to!(c, A, b)
end
