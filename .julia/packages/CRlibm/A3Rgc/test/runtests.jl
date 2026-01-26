using Test

using CRlibm

# Float64
@test CRlibm.cos(0.5, RoundDown) == 0.8775825618903726
@test CRlibm.cos(0.5, RoundUp) == 0.8775825618903728
@test CRlibm.cos(0.5, RoundNearest) == cos(0.5) == 0.8775825618903728
@test CRlibm.cos(1.6, RoundToZero) == -0.029199522301288812
@test CRlibm.cos(1.6, RoundDown) == -0.029199522301288815
@test CRlibm.cos(0.5) == CRlibm.cos(0.5, RoundNearest)

# Float32
@test CRlibm.cos(0.5f0, RoundDown) == 0.87758255f0
@test CRlibm.cos(0.5f0, RoundUp) == 0.8775826f0
@test CRlibm.cos(0.5f0, RoundNearest) == cos(0.5f0) == 0.87758255f0
@test CRlibm.cos(1.6f0, RoundToZero) == -0.029199544f0
@test CRlibm.cos(1.6f0, RoundDown) == -0.029199546f0
@test CRlibm.cos(0.5f0) == CRlibm.cos(0.5f0, RoundNearest)

# Float16
@test CRlibm.cos(Float16(0.5), RoundDown) == Float16(0.8774)
@test CRlibm.cos(Float16(0.5), RoundUp) == Float16(0.878)
@test CRlibm.cos(Float16(0.5), RoundNearest) == cos(Float16(0.5)) == Float16(0.8774)
@test CRlibm.cos(Float16(1.6), RoundToZero) == Float16(-0.02881)
@test CRlibm.cos(Float16(1.6), RoundDown) == Float16(-0.02882)
@test CRlibm.cos(Float16(0.5)) == CRlibm.cos(Float16(0.5), RoundNearest)

# One failing example from #45
@test CRlibm.log10(10.0, RoundUp) == 1.0

function my_eps(prec::Int)
    ldexp(eps(Float64), 53-prec)
end

function my_eps(x::BigFloat)  # only works for precision >= 10?
    my_eps(precision(x) - exponent(x))
end


is_log(f) = string(f)[1:3] == "log"

function do_test(f, val)
    a = f(val, RoundDown)
    b = f(val, RoundUp)
    @test b - a == eps(a) || b - a == eps(b) || b - a == 0
end

function test_CRlibm(functions)
    @show functions
    for f in functions
        #println("Testing CRlibm.$f")

        ff = getfield(CRlibm, f)  # the actual Julia function

        for val in (0.51, 103.2, -17.1, -0.00005, 1)

            #@show f, val

            val <= 0.0 && is_log(f) && continue
            abs(val) > 1 && f ∈ (:asin, :acos) && continue

            do_test(ff, val)
        end
    end
end

function test_MPFR()
    for f in CRlibm.MPFR_functions
        #println("Testing CRlibm.$f")

        ff = getfield(CRlibm, f)  # the actual Julia function

        for val in (0.51, 103.2, -17.1, -0.00005)
            #print(val, " ")

            #@show f, val

            val <= 0.0 && is_log(f) && continue
            abs(val) > 1 && f ∈ (:asin, :acos) && continue

            for prec in (20, 100, 1000)

                setprecision(BigFloat, prec) do
                    val = BigFloat(val)
                    do_test(ff, val)
                end

            end
        end
    end
end


if !(CRlibm.is_32_bit)
    println("Testing CRlibm")
    test_CRlibm(CRlibm.functions)
    # This will currently fail on the :sinpi etc. functions (that are not defined in MPFR) if MPFR is already enabled because the CRlibm library could not be found


    println("Testing shadowing MPFR")
    CRlibm.setup(true)
    test_CRlibm(CRlibm.MPFR_functions)
end

println("Testing MPFR")
test_MPFR()
