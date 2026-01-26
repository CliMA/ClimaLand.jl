using Base.Rounding: setrounding_raw, to_fenv
using Printf

function compare_calc_raw(op, updown, calc, raw, args...)
    if isequal(calc, raw)
        true
    else
        @info("Erorr", op, updown)
        for (i, v) in enumerate(args...)
            @info(@sprintf("a%d = %0.18e, bit rep : %s", i, v, bitstring(v)))
        end

        @info(@sprintf("calc = %0.18e, bit rep : %s", calc, bitstring(calc)))
        @info(@sprintf("raw = %0.18e, bit rep : %s", raw, bitstring(raw)))
        false
    end
end

function rounding_check(op, base_op, arrays...)
    elt = eltype(first(arrays))
    setrounding_raw(elt, to_fenv(RoundNearest))

    @eval begin
        up_calc = broadcast($(Symbol(op, "_up")), $(arrays...))
        down_calc = broadcast($(Symbol(op, "_down")), $(arrays...))

        setrounding_raw($elt, to_fenv(RoundUp))
        up_raw = broadcast($base_op, $(arrays...))
    
        setrounding_raw($elt, to_fenv(RoundDown))
        down_raw = broadcast($base_op, $(arrays...))
    end

    # Compare
    for (calc, raw, args) in zip(up_calc, up_raw, zip(arrays...))
        @test compare_calc_raw(op, "up", calc, raw, args)
    end

    for (calc, raw, args) in zip(down_calc, down_raw, zip(arrays...))
        @test compare_calc_raw(op, "down", calc, raw, args)
    end

    setrounding_raw(elt, to_fenv(RoundNearest)) 
end

rounding_check_unary(a::AbstractVector) = rounding_check(:sqrt, :sqrt, a)
rounding_check_unary(a) = rounding_check_unary([a])

function rounding_check_binary(a::T, b::T) where {T<:AbstractVector}
    for (op, base_op) in zip((:add, :sub, :mul, :div), (:+, :-, :*, :/))
        rounding_check(op, base_op, a, b)
    end
end
rounding_check_binary(a, b) = rounding_check_binary([a], [b])

special_value_list(T::Type) = [
    zero(T), -zero(T),                                  # 0.0, -0.0
    one(T), -one(T),                                    # 1.0, -1.0
    nextfloat(zero(T)), prevfloat(zero(T)),             # N_min^s, -N_min^s
    prevfloat(floatmin(T)), nextfloat(-floatmin(T)),    # N_max^s, -N_max^s
    floatmin(T), -floatmin(T),                          # N_min^n, -N_min^n
    floatmax(T), -floatmax(T),                          # N_max^n, -N_max^n
    eps(T), -eps(T),                                    # machine epsilon
    typemax(T), typemin(T),                             # Inf, -Inf
    T(NaN)                                              # NaN
]

for T in (Float64, Float32)
    @testset "$(T), Special Cases" begin
        special_values = special_value_list(T)
        len = Base.length(special_values)
        a = repeat(special_values, len)
        b = sort(a)
        rounding_check_unary(filter(x->x ≥ zero(x), special_values)) # sqrt
        rounding_check_binary(a, b)
    end
end

for n in 3:6
    N = 10^n
    for T in (Float64, Float32)
        @testset "$(T), Random Sampling, 10^$(n)" begin
            rand_a = reinterpret.(T, rand(Base.uinttype(T), N))
            rand_b = reinterpret.(T, rand(Base.uinttype(T), N))
            rounding_check_unary(abs.(rand_a))
            rounding_check_unary(abs.(rand_b))
            rounding_check_binary(rand_a, rand_b)
            rounding_check_binary(rand_b, rand_a)
        end
    end
end
