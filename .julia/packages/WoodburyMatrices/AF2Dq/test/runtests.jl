using WoodburyMatrices
using LinearAlgebra, SparseArrays, Test
using Aqua
using Random: seed!

@testset "Aqua" begin
    Aqua.test_all(WoodburyMatrices)
end

# helper function for testing logdet
function randpsd(sidelength)
    Q = randn(sidelength, sidelength)
    return Q * Q'
end

include("woodbury.jl")
include("symwoodbury.jl")
