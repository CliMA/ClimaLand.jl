import ArgParse
import JLD2

include(joinpath(@__DIR__, "compute_convergence.jl"))
include(joinpath(pkgdir(ClimaTimeSteppers), "test", "problems.jl"))

function parse_commandline()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        "--alg"
        help = "Algorithm to test convergence"
        arg_type = String
        default = "ARS343"
    end
    parsed_args = ArgParse.parse_args(ARGS, s)
    return parsed_args
end
parsed_args = parse_commandline()
alg_str = parsed_args["alg"]
alg_name = getproperty(ClimaTimeSteppers, Symbol(alg_str))()

convergence_results = Dict()

test_convergence!(convergence_results, alg_name, ark_analytic_nonlin_test_cts(Float64), 450)
test_convergence!(convergence_results, alg_name, ark_analytic_sys_test_cts(Float64), 200)
test_convergence!(
    convergence_results,
    alg_name,
    ark_analytic_test_cts(Float64),
    25000;
    super_convergence_algorithm_names = (ARS121(),),
)
test_convergence!(convergence_results, alg_name, onewaycouple_mri_test_cts(Float64), 9000; high_order_sample_shifts = 3)

test_convergence!(
    convergence_results,
    alg_name,
    climacore_1Dheat_test_cts(Float64),
    200;
    numerical_reference_algorithm_name = ARK548L2SA2(),
)
test_convergence!(
    convergence_results,
    alg_name,
    climacore_1Dheat_test_implicit_cts(Float64),
    200;
    high_order_sample_shifts = 2,
    numerical_reference_algorithm_name = ARK548L2SA2(),
    broken_tests = (SSP22Heuns(), SSP33ShuOsher(), RK4(), IMKG253b(), IMKG254a(), IMKG254b(), IMKG343a()),
) # This problem exceeds the CFL bounds of ERK methods, and it is unstable for several IMKG methods.
test_convergence!(
    convergence_results,
    alg_name,
    climacore_2Dheat_test_cts(Float64),
    200;
    numerical_reference_algorithm_name = ARK548L2SA2(),
)

mkpath("output")
JLD2.save_object(joinpath("output", "convergence_$alg_str.jld2"), convergence_results)

# Unconstrained vs SSP tests
if alg_name isa ClimaTimeSteppers.IMEXSSPRKAlgorithmName
    test_unconstrained_vs_ssp_without_limiters(alg_name, ark_analytic_nonlin_test_cts(Float64), 450)
    test_unconstrained_vs_ssp_without_limiters(alg_name, ark_analytic_sys_test_cts(Float64), 200)
    test_unconstrained_vs_ssp_without_limiters(alg_name, ark_analytic_test_cts(Float64), 25000)
    test_unconstrained_vs_ssp_without_limiters(alg_name, onewaycouple_mri_test_cts(Float64), 9000)
    test_unconstrained_vs_ssp_without_limiters(alg_name, climacore_1Dheat_test_cts(Float64), 200)
    test_unconstrained_vs_ssp_without_limiters(alg_name, climacore_1Dheat_test_implicit_cts(Float64), 200)
    test_unconstrained_vs_ssp_without_limiters(alg_name, climacore_2Dheat_test_cts(Float64), 200)
end
