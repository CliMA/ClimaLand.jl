using ClimaTimeSteppers
using InteractiveUtils: subtypes
using Distributions: quantile, TDist
using LinearAlgebra: norm

all_subtypes(::Type{T}) where {T} = isabstracttype(T) ? vcat(all_subtypes.(subtypes(T))...) : [T]

"""
    imex_convergence_orders(algorithm_name)

Return a tuple containing the expected convergence order of the tableau when
using only an implicit tendency, the order when using only an explicit tendency,
and the order when using both tendencies.
"""
function imex_convergence_orders end
imex_convergence_orders(::ARS111) = (1, 1, 1)
imex_convergence_orders(::ARS121) = (1, 1, 1)
imex_convergence_orders(::ARS122) = (2, 2, 2)
imex_convergence_orders(::ARS222) = (2, 2, 2)
imex_convergence_orders(::ARS232) = (2, 3, 2)
imex_convergence_orders(::ARS233) = (3, 3, 3)
imex_convergence_orders(::ARS343) = (3, 4, 3)
imex_convergence_orders(::ARS443) = (3, 3, 3)
imex_convergence_orders(::Union{IMKG232a, IMKG232b}) = (2, 2, 2)
imex_convergence_orders(::Union{IMKG242a, IMKG242b}) = (2, 4, 2)
imex_convergence_orders(::IMKG243a) = (2, 4, 2)
imex_convergence_orders(::Union{IMKG252a, IMKG252b}) = (2, 2, 2)
imex_convergence_orders(::Union{IMKG253a, IMKG253b}) = (2, 2, 2)
imex_convergence_orders(::Union{IMKG254a, IMKG254b, IMKG254c}) = (2, 2, 2)
imex_convergence_orders(::IMKG342a) = (3, 4, 3)
imex_convergence_orders(::IMKG343a) = (3, 4, 3)
imex_convergence_orders(::SSP222) = (2, 2, 2)
imex_convergence_orders(::SSP322) = (2, 2, 2)
imex_convergence_orders(::SSP332) = (2, 3, 2)
imex_convergence_orders(::SSP333) = (3, 3, 3)
imex_convergence_orders(::SSP433) = (3, 3, 3)
imex_convergence_orders(::DBM453) = (3, 3, 3)
imex_convergence_orders(::HOMMEM1) = (2, 3, 2)
imex_convergence_orders(::ARK2GKC) = (2, 2, 2)
imex_convergence_orders(::ARK437L2SA1) = (4, 4, 4)
imex_convergence_orders(::ARK548L2SA2) = (5, 5, 5)
imex_convergence_orders(::SSP22Heuns) = (2, 2, 2)
imex_convergence_orders(::SSP33ShuOsher) = (3, 3, 3)
imex_convergence_orders(::RK4) = (4, 4, 4)
imex_convergence_orders(::SSPKnoth) = (2, 3, 2)

# Compute a confidence interval for the convergence order, returning the
# estimated convergence order and its uncertainty.
function convergence_order(dts, errs, confidence)
    log_dts = log10.(dts)
    log_errs = log10.(errs)
    n_dts = length(dts)

    # slope and vertical intercept of linear regression (for each log_dt and
    # log_err, log_err ≈ log_dt * order + log_err_of_dt1)
    order, log_err_of_dt1 = hcat(log_dts, ones(n_dts)) \ log_errs

    # number of degrees of freedom of linear regression (number of data
    # points minus number of fitted parameters)
    n_dof = n_dts - 2

    # critical value of Student's t-distribution for two-sided confidence
    # interval
    confidence_interval_t_value = -quantile(TDist(n_dof), (1 - confidence) / 2)

    # standard deviation of linear regression
    regression_standard_deviation = sqrt(sum((log_errs .- (log_dts .* order .+ log_err_of_dt1)) .^ 2) / n_dof)

    # standard deviation of slope
    order_standard_deviation = regression_standard_deviation / sqrt(sum((log_dts .- sum(log_dts) / n_dts) .^ 2))

    # "uncertainty" in slope (half of width of confidence interval)
    order_uncertainty = confidence_interval_t_value * order_standard_deviation

    return order, order_uncertainty
end

"""
    predicted_convergence_order(algorithm_name, ode_function)

Return the predicted convergence order of the algorithm for the given ODE
function (assuming that the algorithm converges).
"""
function predicted_convergence_order(algorithm_name::AbstractAlgorithmName, ode_function::AbstractClimaODEFunction)
    (imp_order, exp_order, combined_order) = imex_convergence_orders(algorithm_name)
    has_imp = !isnothing(ode_function.T_imp!)
    has_exp = ClimaTimeSteppers.has_T_exp(ode_function)
    has_imp && !has_exp && return imp_order
    !has_imp && has_exp && return exp_order
    has_imp && has_exp && return combined_order
    return 0
end

"""
    algorithm(algorithm_name, [linear_implicit])

Generates an appropriate `DistributedODEAlgorithm` from an `AbstractAlgorithmName`.
For `IMEXAlgorithmNames`, `linear_implicit` must also be specified. One Newton
iteration is used for linear implicit problems, and two iterations are used for
nonlinear implicit problems.
"""
algorithm(algorithm_name, _) = algorithm(algorithm_name)
algorithm(algorithm_name::ClimaTimeSteppers.ERKAlgorithmName) = ExplicitAlgorithm(algorithm_name)
algorithm(algorithm_name::ClimaTimeSteppers.SSPKnoth) =
    ClimaTimeSteppers.RosenbrockAlgorithm(ClimaTimeSteppers.tableau(ClimaTimeSteppers.SSPKnoth()))
algorithm(algorithm_name::ClimaTimeSteppers.IMEXARKAlgorithmName, linear_implicit) =
    IMEXAlgorithm(algorithm_name, NewtonsMethod(; max_iters = linear_implicit ? 1 : 2))

rms(array) = norm(array) / sqrt(length(array))
rms_error(u, t, sol) = rms(abs.(u .- sol(t)))

function test_convergence!(
    convergence_results,
    algorithm_name,
    test_case,
    default_num_steps;
    super_convergence_algorithm_names = (),
    num_steps_scaling_factor = 4,
    high_order_sample_shifts = 1,
    numerical_reference_algorithm_name = nothing,
    numerical_reference_num_steps = num_steps_scaling_factor^3 * default_num_steps,
    broken_tests = (),
    error_on_failure = true,
    verbose = false,
)
    (; test_name, t_end, linear_implicit, analytic_sol) = test_case
    prob = test_case.split_prob

    default_dt = t_end / default_num_steps
    ref_sol = if isnothing(numerical_reference_algorithm_name)
        analytic_sol
    else
        # TODO: Do not regenerate the reference solution for every algorithm!!
        ref_alg_str = string(nameof(typeof(numerical_reference_algorithm_name)))
        ref_alg = algorithm(numerical_reference_algorithm_name, linear_implicit)
        ref_dt = t_end / numerical_reference_num_steps
        verbose && @info "Generating reference solution for $test_name with $ref_alg_str and dt = $ref_dt"
        solve(deepcopy(prob), ref_alg; dt = ref_dt, save_everystep = true)
    end
    numerical_reference_info = if isnothing(numerical_reference_algorithm_name)
        nothing
    else
        ref_average_rms_error = rms(rms_error.(ref_sol.u, ref_sol.t, (analytic_sol,)))
        (ref_alg_str, ref_dt, ref_average_rms_error)
    end

    alg_str = string(nameof(typeof(algorithm_name)))
    alg = algorithm(algorithm_name, linear_implicit)
    verbose && @info "Testing convergence of $alg_str for $test_name"

    predicted_order = predicted_convergence_order(algorithm_name, prob.f)
    predicted_super_convergence = algorithm_name in super_convergence_algorithm_names
    num_steps_powers = (-1:0.5:1) .- high_order_sample_shifts * max(0, predicted_order - 3) / 2
    sampled_num_steps = default_num_steps .* num_steps_scaling_factor .^ num_steps_powers
    sampled_dts = t_end ./ round.(Int, sampled_num_steps)
    average_rms_errors = map(sampled_dts) do dt
        sol = solve(deepcopy(prob), alg; dt = dt, save_everystep = true)
        rms(rms_error.(sol.u, sol.t, (ref_sol,)))
    end
    verbose && @info "Sampled timesteps = $sampled_dts"
    verbose && @info "Average RMS errors = $average_rms_errors"

    # Compute a 99% confidence interval for the convergence order
    order, order_uncertainty = convergence_order(sampled_dts, average_rms_errors, 0.99)
    verbose && @info "Convergence order = $order ± $order_uncertainty"
    actual_predicted_order = predicted_order + Bool(predicted_super_convergence)
    convergence_test_error = if isnan(order)
        "Timestepper does not converge for $alg_str ($test_name)"
    elseif abs(order - actual_predicted_order) > order_uncertainty
        "Predicted order outside error bars for $alg_str ($test_name)"
    elseif order_uncertainty > actual_predicted_order / 10
        "Order uncertainty too large for $alg_str ($test_name)"
    else
        nothing
    end
    if isnothing(convergence_test_error)
        @assert !(algorithm_name in broken_tests)
    elseif error_on_failure && !(algorithm_name in broken_tests)
        error(convergence_test_error)
    else
        @warn convergence_test_error
    end

    default_dt_sol = solve(deepcopy(prob), alg; dt = default_dt, save_everystep = true)
    default_dt_times = default_dt_sol.t
    default_dt_solutions = rms.(default_dt_sol.u)
    default_dt_errors = rms_error.(default_dt_sol.u, default_dt_sol.t, (ref_sol,))

    convergence_results[test_name] = Dict()
    convergence_results[test_name]["default_dt"] = default_dt
    convergence_results[test_name]["numerical_reference_info"] = numerical_reference_info
    convergence_results[test_name]["all_alg_results"] = Dict()
    convergence_results[test_name]["all_alg_results"][alg_str] = Dict()
    alg_results = convergence_results[test_name]["all_alg_results"][alg_str]
    alg_results["predicted_order"] = predicted_order
    alg_results["predicted_super_convergence"] = predicted_super_convergence
    alg_results["sampled_dts"] = sampled_dts
    alg_results["average_rms_errors"] = average_rms_errors
    alg_results["default_dt_times"] = default_dt_times
    alg_results["default_dt_solutions"] = default_dt_solutions
    alg_results["default_dt_errors"] = default_dt_errors
    return convergence_results
end

function test_unconstrained_vs_ssp_without_limiters(algorithm_name, test_case, num_steps)
    prob = test_case.split_prob
    dt = test_case.t_end / num_steps
    newtons_method = NewtonsMethod(; max_iters = test_case.linear_implicit ? 1 : 2)
    algorithm = IMEXAlgorithm(algorithm_name, newtons_method)
    reference_algorithm = IMEXAlgorithm(algorithm_name, newtons_method, Unconstrained())
    solution = solve(deepcopy(prob), algorithm; dt).u[end]
    reference_solution = solve(deepcopy(prob), reference_algorithm; dt).u[end]
    relative_error = norm(solution .- reference_solution) / norm(reference_solution)
    if relative_error > 100 * eps(Float64)
        error("Unconstrained and SSP versions of $algorithm_name give \
               different results for $(test_case.test_name): relative \
               error = $(round(Int, relative_error / eps(Float64))) * eps")
    end
end
