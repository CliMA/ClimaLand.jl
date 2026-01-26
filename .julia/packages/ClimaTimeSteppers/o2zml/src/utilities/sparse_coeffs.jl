"""
    SparseCoeffs(coefficients)

A mask for coefficients. Supports `getindex(::SparseCoeffs, ijk...)`
that forwards to coefficients, and `getindex(::Type{SparseCoeffs}, ijk...)`
that forwards (at compile time) to the mask, which behaves as a
BitArray of the coefficients.
"""
struct SparseCoeffs{S, m, C}
    coeffs::C
    function SparseCoeffs(coeffs::C) where {C}
        m = BitArray(iszero.(coeffs))
        return new{size(m), Tuple(m), C}(coeffs)
    end
end

# Forward array behavior:
Base.@propagate_inbounds Base.getindex(sc::SparseCoeffs, inds...) = @inbounds sc.coeffs[inds...]
Base.length(sc::SparseCoeffs) = length(sc.coeffs)
import LinearAlgebra
LinearAlgebra.diag(sc::SparseCoeffs, args...) = LinearAlgebra.diag(sc.coeffs, args...)
LinearAlgebra.adjoint(sc::SparseCoeffs) = LinearAlgebra.adjoint(sc.coeffs)

get_S(::SparseCoeffs{S}) where {S} = S

# Special behavior of SparseCoeffs:
Base.@propagate_inbounds zero_coeff(::Type{SparseCoeffs{S, m, C}}, i::Int, j::Int) where {S, m, C} =
    @inbounds m[i + S[1] * (j - 1)]
Base.@propagate_inbounds zero_coeff(::Type{SparseCoeffs{S, m, C}}, j::Int) where {S, m, C} = @inbounds m[j]

Base.convert(::Type{T}, x::SArray) where {T <: SparseCoeffs} = SparseCoeffs(x)
