#=
julia --project=test
using Revise; include("test/tabulate_convergence_orders_rosenbrock.jl")
=#
using ClimaTimeSteppers, LinearAlgebra, Test
import PrettyTables
import SciMLBase

if !@isdefined(IntegratorTestCase)
    include(joinpath(@__DIR__, "convergence_orders.jl"))
    include(joinpath(@__DIR__, "convergence_utils.jl"))
    include(joinpath(@__DIR__, "utils.jl"))
    include(joinpath(@__DIR__, "problems.jl"))
end

function tabulate_convergence_orders_rosenbrock()
    tabs = [SSPKnoth]
    tabs = map(t -> t(), tabs)
    test_cases = all_test_cases(Float64)
    results = convergence_order_results(tabs, test_cases)

    prob_names = map(t -> t.test_name, test_cases)
    expected_orders = SciMLBase.alg_order.(tabs)
    algs = algorithm.(tabs)

    tabulate_convergence_orders(prob_names, algs, results, expected_orders; tabs)
    return results
end

tabulate_convergence_orders_rosenbrock()
