# This is used on Julia version that have the Base.Ryu module.

using Base.Ryu

function plain_precision_heuristic(xs::AbstractArray{<:AbstractFloat})
    ys = filter(isfinite, xs)
    e10max = -(e10min = typemax(Int))
    for y in ys
        if isapprox(y, 0, atol=1e-16)
            e10 = min(e10min, 0)
        else
            _, e10 = Ryu.reduce_shortest(convert(Float32, y))
        end
        e10min = min(e10min, e10)
        e10max = max(e10max, e10)
    end
    return precision = min(-e10min, -e10max+16)
end

# Print a floating point number at fixed precision. Pretty much equivalent to
# @sprintf("%0.$(precision)f", x), without the macro issues.
function format_fixed(x::AbstractFloat, precision::Integer)
    if x == Inf
        return "∞"
    elseif x == -Inf
        return "-∞"
    elseif isnan(x)
        return "NaN"
    end

    return Ryu.writefixed(x, precision)
end

# Print a floating point number in scientific notation at fixed precision. Sort of equivalent
# to @sprintf("%0.$(precision)e", x), but prettier printing.
function format_fixed_scientific(x::AbstractFloat, precision::Integer,
                                 engineering::Bool)
    if iszero(x)
        return "0"
    elseif isinf(x)
        return signbit(x) ? "-∞" : "∞"
    elseif isnan(x)
        return "NaN"
    end

    if engineering
        # precision applies before convertic to engineering format
        # that is, number of digits would be the same as the with non-engineering format
        base_digits, power = get_engineering_string(x, precision)
    else
        e_format_number = Ryu.writeexp(x, precision)
        base_digits, power = split(e_format_number, 'e')
    end


    buf = IOBuffer()

    print(buf, base_digits)
    print(buf, "×10")

    if power[1] == '-'
        print(buf, '⁻')
    end
    leading_index = findfirst(c -> '1' <= c <= '9', power)

    if leading_index === nothing
        print(buf, superscript_numerals[1])
        return String(take!(buf))
    end

    for digit in power[leading_index:end]
        if digit == '-'
            print(buf, '⁻')
        elseif '0' <= digit <= '9'
            print(buf, superscript_numerals[digit - '0' + 1])
        end

    end

    return String(take!(buf))
end


function get_engineering_string(x::AbstractFloat, precision::Integer)
    e_format_number = Ryu.writeexp(x, precision)
    base_digits, power = split(e_format_number, 'e')

    int_power = parse(Int, power)
    positive = int_power >= 0

    # round the power to the nearest multiple of 3
    # positive power -> move the "." to the right by mode, round the power to the higher power
    # negative power -> move the "." to the right by mode, round the power to the lower power
    # ex:
    # 1.2334e5 = 123.334e3
    # 1.2334-5 = 12.3334e-6

    if positive
        indices_to_move = int_power % 3
    else
        indices_to_move = 3 - abs(int_power) % 3
    end

    buf = IOBuffer()
    negative_base_compensation = base_digits[1] == '-' ? 1 : 0
    for i in eachindex(base_digits)
        if base_digits[i] != '.'
            print(buf, base_digits[i])
        end
        if i == 2 + indices_to_move + negative_base_compensation
            print(buf, '.')
        end
    end

    return String(take!(buf)), string(int_power - indices_to_move)
end
