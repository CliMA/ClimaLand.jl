using Showoff
using Test
using Dates
using Printf

const drops0s = !isdefined(Base, :Ryu)

@testset "Internals" begin
    @test Showoff.@grisu_ccall(1, 2, 3) === nothing
    if isdefined(Showoff, :Grisu)
        @test Showoff.grisu(1.0, Showoff.Grisu.SHORTEST, 2) == (1, 1, false, Showoff.Grisu.DIGITS)
    end

    let x = [1.0, Inf, 2.0, NaN]
        @test Showoff.concrete_minimum(x) == 1.0
        @test Showoff.concrete_maximum(x) == 2.0
    end

    @test_throws ArgumentError Showoff.concrete_minimum([])
    @test_throws ArgumentError Showoff.concrete_maximum([])

    let x = [1.12345, 4.5678]
        @test Showoff.plain_precision_heuristic(x) == 5
        @test Showoff.scientific_precision_heuristic(x) == 6
    end
end

@testset "Formatting" begin
    @test Showoff.format_fixed(-10.0, 0) == "-10"
    @test Showoff.format_fixed(0.012345, 3) == "0.012"
    @test Showoff.format_fixed(Inf, 1) == "∞"
    @test Showoff.format_fixed(-Inf, 1) == "-∞"
    @test Showoff.format_fixed(NaN, 1) == "NaN"
    @test Showoff.format_fixed_scientific(0.0, 1, false) == "0"
    @test Showoff.format_fixed_scientific(Inf, 1, false) == "∞"
    @test Showoff.format_fixed_scientific(-Inf, 1, false) == "-∞"
    @test Showoff.format_fixed_scientific(NaN, 1, false) == "NaN"
    @test Showoff.format_fixed_scientific(0.012345678, 4, true) == "12.346×10⁻³"
    @test Showoff.format_fixed_scientific(0.012345678, 4, false) == "1.2346×10⁻²"
    @test Showoff.format_fixed_scientific(-10.0, 4, false) == (drops0s ? "-1.000×10¹" : "-1.0000×10¹")
    @test Showoff.format_fixed_scientific(-10.0, 4, true) == "-10.000×10⁰"
    @test Showoff.format_fixed_scientific(-10.0, 4, false)[1:end-5] == (drops0s ? @sprintf("%0.4e", -10.0)[1:end-5] : @sprintf("%0.4e", -10.0)[1:end-4])
    @test Showoff.format_fixed_scientific(1.23456e7, 3, false)[1:end-5] == @sprintf("%0.3e", 1.23456e7)[1:end-4]
    @test Showoff.format_fixed_scientific(2.99999999999999956E-16, 2, false) == (drops0s ? "3.0×10⁻¹⁶" : "3.00×10⁻¹⁶")
end

@testset "Showoff" begin
    x = [1.12345, 4.5678]
    @test showoff(x) == ["1.12345", "4.56780"]
    @test showoff([0.0, 50000.0]) == (drops0s ? ["0", "5×10⁴"] : ["0", "5.0×10⁴"])
    @test showoff(x, :plain) == ["1.12345", "4.56780"]
    @test showoff([0.0], :plain) == ["0"]
    @test showoff(x, :scientific) == (drops0s ? ["1.12345×10⁰", "4.56780×10⁰"] : ["1.123450×10⁰", "4.567800×10⁰"])
    @test showoff(x, :engineering) == showoff(x, :scientific)
    @test showoff([DateTime("2017-04-11", "yyyy-mm-dd")]) == ["Apr 11, 2017"]
    @test showoff(["a", "b"]) == ["\"a\"", "\"b\""]
    @test showoff([1, 1e39]) == (drops0s ? ["1×10⁰", "1×10³⁹"] : ["1.0×10⁰", "1.0×10³⁹"])
    @test_throws ArgumentError showoff(x, :nevergonnagiveyouup)
    @test showoff([Inf, Inf, NaN]) == ["Inf", "Inf", "NaN"]
end
