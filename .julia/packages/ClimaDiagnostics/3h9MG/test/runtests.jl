using SafeTestsets
using Test

#! format: off
@safetestset "Aqua" begin @time include("aqua.jl") end
@safetestset "Format" begin @time include("format.jl") end
@safetestset "Doctest" begin @time include("doctest.jl") end

@safetestset "Writers" begin @time include("writers.jl") end

@safetestset "Schedules" begin @time include("schedules.jl") end
@safetestset "DiagnosticVariable" begin @time include("diagnostic_variable.jl") end
@safetestset "SchduledDiagnostics" begin @time include("diagnostics.jl") end

@safetestset "Integration test" begin @time include("integration_test.jl") end
#! format: on

return nothing
