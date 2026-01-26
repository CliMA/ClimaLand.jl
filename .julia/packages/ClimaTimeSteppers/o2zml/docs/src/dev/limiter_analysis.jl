using ClimaTimeSteppers
import ArgParse
import JLD2

include(joinpath(pkgdir(ClimaTimeSteppers), "test", "problems.jl"))

function parse_commandline()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        "--alg"
        help = "Algorithm to test convergence"
        arg_type = String
        default = "ARS343"
        "--use_limiter"
        help = "Bool indicating to use limiter"
        arg_type = Bool
        default = true
        "--use_hyperdiffusion"
        help = "Bool indicating to use hyperdiffusion"
        arg_type = Bool
        default = true
    end
    parsed_args = ArgParse.parse_args(ARGS, s)
    return parsed_args
end
parsed_args = parse_commandline()
alg_str = parsed_args["alg"]
use_limiter = parsed_args["use_limiter"]
use_hyperdiffusion = parsed_args["use_hyperdiffusion"]
alg_name = getproperty(ClimaTimeSteppers, Symbol(alg_str))()

# This also runs with num_steps = 1000, but with larger under/overshoots; 4800
# is the value used in the paper.
num_steps = 4800

test_case = horizontal_deformational_flow_test(Float64; use_limiter, use_hyperdiffusion)
prob = test_case.split_prob
dt = test_case.t_end / num_steps
algorithm =
    alg_name isa ClimaTimeSteppers.IMEXARKAlgorithmName ? IMEXAlgorithm(alg_name, NewtonsMethod()) :
    ExplicitAlgorithm(alg_name)
solution = solve(deepcopy(prob), algorithm; dt).u

limiter_results = Dict()
limiter_results[alg_str, use_hyperdiffusion, use_limiter] = solution
mkpath("output")
filename = "limiter_$(alg_str)_hyperdiff_$(use_hyperdiffusion)_lim_$(use_limiter).jld2"
JLD2.save_object(joinpath("output", filename), limiter_results)
