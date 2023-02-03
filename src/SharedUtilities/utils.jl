"""
     heaviside(x::FT)::FT where {FT}

Computes the heaviside function.
"""
function heaviside(x::FT)::FT where {FT}
    if x >= FT(0.0)
        return FT(1.0)
    else
        return FT(0.0)
    end
end
