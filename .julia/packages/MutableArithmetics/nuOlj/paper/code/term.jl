struct Term{T}
    coef::T
    sym::SymbolicVariable
end
struct Sum{T}
    terms::Vector{Term{T}}
end
Base.:+(s::Sum, t::Term) = Sum(push!(copy(s.terms), t))
Base.zero(::Type{Term{T}}) where {T} = Sum(Term{T}[])
