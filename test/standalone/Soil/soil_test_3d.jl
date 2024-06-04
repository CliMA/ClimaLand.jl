FT = Float64
using Test
x = (10 ^ 11) * eps(FT)
# x - 2eps() ≤ x ≤ x + 2eps()
# @test 2.220446049205904e-5 ≤ x ≤ 2.220446049294722e-5

@show eps(FT)
@show (10 ^ 1) * eps(FT)
@show (10 ^ 2) * eps(FT)
@show (10 ^ 3) * eps(FT)
@show (10 ^ 4) * eps(FT)
@show (10 ^ 5) * eps(FT)
@show (10 ^ 6) * eps(FT)
@show (10 ^ 7) * eps(FT)
@show (10 ^ 8) * eps(FT)
@show (10 ^ 9) * eps(FT)
@show (10 ^ 10) * eps(FT)
@show (10 ^ 11) * eps(FT)

x = (10 ^ 11) * eps(Float64)
# x - 2eps() ≤ x ≤ x + 2eps()
# @test 2.220446049205904e-5 ≤ x ≤ 2.220446049294722e-5

@show eps(Float64)
@show (10 ^ 1) * eps(Float64)
@show (10 ^ 2) * eps(Float64)
@show (10 ^ 3) * eps(Float64)
@show (10 ^ 4) * eps(Float64)
@show (10 ^ 5) * eps(Float64)
@show (10 ^ 6) * eps(Float64)
@show (10 ^ 7) * eps(Float64)
@show (10 ^ 8) * eps(Float64)
@show (10 ^ 9) * eps(Float64)
@show (10 ^ 10) * eps(Float64)
@show (10 ^ 11) * eps(Float64)
