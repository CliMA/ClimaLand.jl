using ArgParse, JET, Test, DiffEqBase, ClimaTimeSteppers
using ClimaComms
using CUDA, BenchmarkTools, OrderedCollections, StatsBase, PrettyTables # needed for CTS.benchmark_step
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
function config_integrators(problem)
    algorithm = CTS.IMEXAlgorithm(ARS343(), NewtonsMethod(; max_iters = 2))
    dt = 0.01
    integrator = DiffEqBase.init(problem, algorithm; dt)
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

device = ClimaComms.device()
CTS.benchmark_step(integrator, device)
