using NLSolversBase, Test
using Random
using LinearAlgebra: Diagonal, I
using ComponentArrays
using SparseArrays
using OptimTestProblems
using RecursiveArrayTools
using ADTypes
MVP = OptimTestProblems.MultivariateProblems

# TODO: Use OptimTestProblems (but it does not have exponential_gradient_hession etc.)
# TODO: MultivariateProblems.UnconstrainedProblems.examples["Exponential"]

# Test example
function exponential(x)
    return exp((2.0 - x[1])^2) + exp((3.0 - x[2])^2)
end

function exponential_gradient!(storage, x)
    storage[1] = -2.0 * (2.0 - x[1]) * exp((2.0 - x[1])^2)
    storage[2] = -2.0 * (3.0 - x[2]) * exp((3.0 - x[2])^2)
end


function exponential_gradient(x)
    storage = similar(x)
    storage[1] = -2.0 * (2.0 - x[1]) * exp((2.0 - x[1])^2)
    storage[2] = -2.0 * (3.0 - x[2]) * exp((3.0 - x[2])^2)
    storage
end


function exponential_value_gradient!(storage, x)
    storage[1] = -2.0 * (2.0 - x[1]) * exp((2.0 - x[1])^2)
    storage[2] = -2.0 * (3.0 - x[2]) * exp((3.0 - x[2])^2)
    return exp((2.0 - x[1])^2) + exp((3.0 - x[2])^2)
end

function exponential_value_gradient(x)
    return exponential(x), exponential_gradient(x)
end

function exponential_gradient_hessian(x)
    F = similar(x)
    F[1] = -2.0 * (2.0 - x[1]) * exp((2.0 - x[1])^2)
    F[2] = -2.0 * (3.0 - x[2]) * exp((3.0 - x[2])^2)

    nx = length(x)
    J = zeros(nx, nx)
    J[1, 1] = 2.0 * exp((2.0 - x[1])^2) * (2.0 * x[1]^2 - 8.0 * x[1] + 9)
    J[1, 2] = 0.0
    J[2, 1] = 0.0
    J[2, 2] = 2.0 * exp((3.0 - x[1])^2) * (2.0 * x[2]^2 - 12.0 * x[2] + 19)
    F, J
end

function exponential_hessian(x)
    nx = length(x)
    storage = zeros(nx, nx)
    storage[1, 1] = 2.0 * exp((2.0 - x[1])^2) * (2.0 * x[1]^2 - 8.0 * x[1] + 9)
    storage[1, 2] = 0.0
    storage[2, 1] = 0.0
    storage[2, 2] = 2.0 * exp((3.0 - x[1])^2) * (2.0 * x[2]^2 - 12.0 * x[2] + 19)
    storage
end

function exponential_hessian!(storage, x)
    storage[1, 1] = 2.0 * exp((2.0 - x[1])^2) * (2.0 * x[1]^2 - 8.0 * x[1] + 9)
    storage[1, 2] = 0.0
    storage[2, 1] = 0.0
    storage[2, 2] = 2.0 * exp((3.0 - x[1])^2) * (2.0 * x[2]^2 - 12.0 * x[2] + 19)
end

function exponential_hessian_product!(storage, x)
    storage[1, 1] = 2.0 * exp((2.0 - x[1])^2) * (2.0 * x[1]^2 - 8.0 * x[1] + 9)
    storage[1, 2] = 0.0
    storage[2, 1] = 0.0
    storage[2, 2] = 2.0 * exp((3.0 - x[1])^2) * (2.0 * x[2]^2 - 12.0 * x[2] + 19)
end

@testset verbose=true "NLSolversBase.jl" begin
    include("objective_types.jl")
    include("interface.jl")
    include("incomplete.jl")
    include("constraints.jl")
    include("abstractarrays.jl")
    include("autodiff.jl")
    include("sparse.jl")
    include("kwargs.jl")
end
