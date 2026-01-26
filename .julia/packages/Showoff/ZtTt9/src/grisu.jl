# This is used on Julia version that have the Base.Grisu module

import Base.Grisu

function grisu(v::AbstractFloat, mode, requested_digits)
    return tuple(Grisu.grisu(v, mode, requested_digits)..., Grisu.DIGITS)
end

function plain_precision_heuristic(xs::AbstractArray{<:AbstractFloat})
    ys = filter(isfinite, xs)
    precision = 0
    for y in ys
        len, point, neg, digits = grisu(convert(Float32, y), Grisu.SHORTEST, 0)
        precision = max(precision, len - point)
    end
    return max(precision, 0)
end

# Print a floating point number at fixed precision. Pretty much equivalent to
# @sprintf("%0.$(precision)f", x), without the macro issues.
function format_fixed(x::AbstractFloat, precision::Integer)
    @assert precision >= 0

    if x == Inf
        return "∞"
    elseif x == -Inf
        return "-∞"
    elseif isnan(x)
        return "NaN"
    end

    len, point, neg, digits = grisu(x, Grisu.FIXED, precision)

    buf = IOBuffer()
    if x < 0
        print(buf, '-')
    end

    for c in digits[1:min(point, len)]
        print(buf, convert(Char, c))
    end

    if point > len
        for _ in len:point-1
            print(buf, '0')
        end
    elseif point < len
        if point <= 0
            print(buf, '0')
        end
        print(buf, '.')
        if point < 0
            for _ in 1:-point
                print(buf, '0')
            end
            for c in digits[1:len]
                print(buf, convert(Char, c))
            end
        else
            for c in digits[point+1:len]
                print(buf, convert(Char, c))
            end
        end
    end

    trailing_zeros = precision - max(0, len - point)
    if trailing_zeros > 0 && point >= len
        print(buf, '.')
    end

    for _ in 1:trailing_zeros
        print(buf, '0')
    end

    String(take!(buf))
end

# Print a floating point number in scientific notation at fixed precision. Sort of equivalent
# to @sprintf("%0.$(precision)e", x), but prettier printing.
function format_fixed_scientific(x::AbstractFloat, precision::Integer,
                                 engineering::Bool)
    if x == 0.0
        return "0"
    elseif x == Inf
        return "∞"
    elseif x == -Inf
        return "-∞"
    elseif isnan(x)
        return "NaN"
    end

    mag = floor(Int, log10(abs(x)))
    grisu_precision = precision

    len, point, neg, digits = grisu((x / 10.0^mag), Grisu.FIXED, grisu_precision)
    point += mag

    @assert len > 0

    buf = IOBuffer()
    if x < 0
        print(buf, '-')
    end

    print(buf, convert(Char, digits[1]))
    nextdigit = 2
    if engineering
        while (point - 1) % 3 != 0
            if nextdigit <= len
                print(buf, convert(Char, digits[nextdigit]))
            else
                print(buf, '0')
            end
            nextdigit += 1
            point -= 1
        end
    end

    if precision > 1
        print(buf, '.')
    end

    for i in nextdigit:len
        print(buf, convert(Char, digits[i]))
    end

    for i in (len+1):precision
        print(buf, '0')
    end

    print(buf, "×10")
    for c in string(point - 1)
        if '0' <= c <= '9'
            print(buf, superscript_numerals[c - '0' + 1])
        elseif c == '-'
            print(buf, '⁻')
        end
    end

    return String(take!(buf))
end

