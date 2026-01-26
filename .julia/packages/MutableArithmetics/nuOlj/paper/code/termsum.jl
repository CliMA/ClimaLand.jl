function sum(x)
    acc = zero(eltype(x))
    for el in x
        acc = add!!(acc, el)
    end
    return acc
end
add!!(a, b) = a + b # default fallback
add!!(a::BigInt, b::BigInt) = Base.GMP.MPZ.add!(a, b)
function add!!(s::Sum, t::Term)
    push!(s.terms, t)
    return s
end
