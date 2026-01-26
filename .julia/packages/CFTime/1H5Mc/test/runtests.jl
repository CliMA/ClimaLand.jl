using CFTime
using Test

@testset verbose = true "All tests" begin
    @testset verbose = true "Basic functionality" begin
        include("test_time.jl")
    end

    @testset verbose = true "Comparision with reference algorithm" begin
        include("test_validity.jl")
    end

    @testset "Operators (+,-,<,>...)" begin
        include("test_operators.jl")
    end

    @testset verbose = true "CF conventions" begin
        include("test_cf.jl")
    end

    @testset "Different time resolutions" begin
        include("test_resolution.jl")
    end

    @testset "Rounding" begin
        include("test_rounding.jl")
    end

    @testset "Year 0" begin
        include("test_year0.jl")
    end

    @testset "Reported issues" begin
        include("test_issues.jl")
    end

    @testset "Aqua" begin
        include("test_aqua.jl")
    end
end
