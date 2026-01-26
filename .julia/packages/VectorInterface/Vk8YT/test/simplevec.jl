using LinearAlgebra: LinearAlgebra, rmul!, mul!, axpy!, axpby!, dot, norm

struct SimpleVec{V<:AbstractArray}
    vec::V
end

Base.getindex(v::SimpleVec) = v.vec # for convience, should not interfere
Base.eltype(::Type{SimpleVec{V}}) where {V} = V
Base.iterate(v::SimpleVec, args...) = iterate((v[],), args...)
Base.length(v::SimpleVec) = 1

Base.:+(v::SimpleVec, w::SimpleVec) = SimpleVec(v[] + w[])
Base.:*(v::SimpleVec, a::Number) = SimpleVec(a * v[])

Base.similar(v::SimpleVec, S) = SimpleVec(similar(v[], S))

LinearAlgebra.rmul!(v::SimpleVec, α) = (rmul!(v[], α); return v)
LinearAlgebra.mul!(w::SimpleVec, v::SimpleVec, α) = (mul!(w[], v[], α); return w)

LinearAlgebra.axpy!(α, v::SimpleVec, w::SimpleVec) = (axpy!(α, v[], w[]); return w)
LinearAlgebra.axpby!(α, v::SimpleVec, β, w::SimpleVec) = (axpby!(α, v[], β, w[]); return w)

LinearAlgebra.dot(v::SimpleVec, w::SimpleVec) = dot(v[], w[])
LinearAlgebra.norm(v::SimpleVec) = norm(v[])
