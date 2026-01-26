using ArgParse, Test, BenchmarkTools, DiffEqBase, ClimaTimeSteppers
import ClimaTimeSteppers as CTS
function parse_commandline()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        "--job_id"
        help = "Job ID"
        arg_type = String
    end
    parsed_args = ArgParse.parse_args(ARGS, s)
    return (s, parsed_args)
end
(s, parsed_args) = parse_commandline()
cts = joinpath(dirname(@__DIR__));
include(joinpath(cts, "test", "problems.jl"))
function do_work!(integrator, cache)
    for _ in 1:100_000
        CTS.step_u!(integrator, cache)
    end
end
test_case = climacore_1Dheat_test_cts(Float64)
prob = test_case.split_prob
algorithm = CTS.IMEXAlgorithm(ARS343(), NewtonsMethod(; max_iters = 2))
dt = 0.01
integrator = DiffEqBase.init(prob, algorithm; dt)
cache = CTS.init_cache(prob, algorithm)
do_work!(integrator, cache) # compile first
import Profile
Profile.clear_malloc_data()
Profile.clear()
prof = Profile.@profile begin
    do_work!(integrator, cache)
end

p = @allocated do_work!(integrator, cache)
using Test
@testset "Allocations" begin
    @test p == 0
end

import ProfileCanvas

if haskey(ENV, "BUILDKITE_COMMIT") || haskey(ENV, "BUILDKITE_BRANCH")
    output_dir = joinpath(cts, parsed_args["job_id"])
    mkpath(output_dir)
    ProfileCanvas.html_file(joinpath(output_dir, "flame.html"))
else
    ProfileCanvas.view(Profile.fetch())
end
