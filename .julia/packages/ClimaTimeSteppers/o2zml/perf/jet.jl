# using Revise; include("perf/jet.jl")
using ArgParse, JET, Test, BenchmarkTools, SciMLBase, ClimaTimeSteppers
import ClimaTimeSteppers as CTS
function parse_commandline()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        "--problem"
        help = "Problem type [`ode_fun`, `fe`]"
        arg_type = String
        default = "diffusion2d"
    end
    parsed_args = ArgParse.parse_args(ARGS, s)
    return (s, parsed_args)
end
(s, parsed_args) = parse_commandline()
cts = joinpath(dirname(@__DIR__));
include(joinpath(cts, "test", "problems.jl"))
config_integrators(itc::IntegratorTestCase) = config_integrators(itc.prob)

struct Foo end
foo!(integrator) = nothing
(::Foo)(integrator) = foo!(integrator)
struct Bar end
bar!(integrator) = nothing
(::Bar)(integrator) = bar!(integrator)

function discrete_cb(cb!, n)
    cond = if n == 1
        (u, t, integrator) -> isnothing(cb!(integrator))
    else
        (u, t, integrator) -> isnothing(cb!(integrator)) || rand() ≤ 0.5
    end
    SciMLBase.DiscreteCallback(cond, cb!;)
end
function config_integrators(problem)
    algorithm = CTS.IMEXAlgorithm(ARS343(), NewtonsMethod(; max_iters = 2))
    dt = 0.01
    discrete_callbacks = (discrete_cb(Foo(), 0), discrete_cb(Bar(), 0), discrete_cb(Foo(), 1), discrete_cb(Bar(), 1))
    callback = SciMLBase.CallbackSet((), discrete_callbacks)

    integrator = SciMLBase.init(problem, algorithm; dt, callback)
    integrator.cache = CTS.init_cache(problem, algorithm)
    return (; integrator)
end
prob = if parsed_args["problem"] == "diffusion2d"
    climacore_2Dheat_test_cts(Float64)
elseif parsed_args["problem"] == "ode_fun"
    split_linear_prob_wfact_split()
elseif parsed_args["problem"] == "fe"
    split_linear_prob_wfact_split_fe()
else
    error("Bad option")
end
(; integrator) = config_integrators(prob)

@testset "JET / allocations" begin
    CTS.step_u!(integrator, integrator.cache) # compile first, and make sure it runs
    step_allocs = @allocated CTS.step_u!(integrator, integrator.cache)
    @show step_allocs
    JET.@test_opt CTS.step_u!(integrator, integrator.cache)

    CTS.__step!(integrator) # compile first, and make sure it runs
    JET.@test_opt CTS.__step!(integrator)
end
