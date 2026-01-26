using SciMLStructures, Test, SafeTestsets

@testset "SciMLStructures" begin
    @safetestset "Quality Assurance" include("qa.jl")
end
