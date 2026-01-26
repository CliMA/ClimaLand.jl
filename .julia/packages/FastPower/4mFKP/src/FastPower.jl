module FastPower

# From David Goldberg's blog post with minor modifications
# https://tech.ebayinc.com/engineering/fast-approximate-logarithms-part-iii-the-formulas/
@inline function fastlog2(x::Float32)::Float32
    # (x-1)*(a*(x-1) + b)/((x-1) + c) (line 8 of table 2)
    # b = b-a and c = c-1 so we can delay subtracting 1 from signif
    # in the computation of quot to get better ipc
    a = 0.338953f0
    b = 1.8596461f0
    c = 0.523692f0
    # IEEE is sgn(1):exp(8):frac(23) representing
    # (1+frac)*2^(exp-127). 1+frac is called the significand
    
    # get exponent
    ux1i = reinterpret(UInt32, x)
    exp = Int32((ux1i & 0x7F800000) >> 23)
    # actual exponent is exp-127, will subtract later

    if iszero(ux1i & 0x00400000) # if signif > 1.5
        ux2i = (ux1i & 0x007FFFFF) | 0x3f800000
        exp -= 0x7f # 127
    else
        ux2i = (ux1i & 0x007FFFFF) | 0x3f000000
        exp -= 0x7e # 126 instead of 127 compensates for division by 2
    end
    signif = reinterpret(Float32, ux2i) 
    quot = muladd(signif, a, b) / (signif + c)
    return muladd(signif - 1.0f0, quot, exp)
end

"""
    fastpower(x::T, y::T) where {T} -> float(T)
    Trips through Float32 for performance.
"""
@inline function fastpower(x::T, y::T) where {T <: Real}
    outT = float(T)
    if iszero(x)
        return zero(outT)
    elseif isinf(x) && isinf(y)
        return convert(outT, Inf)
    else
        return convert(
            outT,
            @fastmath exp2(convert(Float32, y) * fastlog2(convert(Float32, x)))
        )
    end
end

fastpower(x, y) = x^y

end
