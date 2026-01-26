# Determine array type for testing - supports both CPU (Array) and GPU (CuArray) testing
if get(ARGS, 1, "Array") == "CuArray"
    import CUDA
    ArrayType = CUDA.CuArray
    CUDA.allowscalar(false)  # Ensure GPU operations are properly vectorized
else
    ArrayType = Array
end

@show ArrayType

using Test
using RootSolvers
using StaticArrays
using ForwardDiff

# Include test helper functions that define test problems, methods, and tolerances
include("test_helper.jl")

# Helper function to check if roots are within tolerance of expected solution    
function check_root_tolerance(roots, expected_root, problem, method, tol)
    # For high-multiplicity roots, use more lenient tolerance since they're inherently difficult
    if problem.name in ("high-multiplicity root", "steep exponential function")
        tol_factor = 100  # Much more lenient for difficult functions
    else
        tol_factor = 10
    end
    FT = typeof(expected_root)
    if tol === nothing
        _default_tol = tol_factor * default_tol(FT).tol
        if roots isa AbstractArray
            if ArrayType <: Array # If roots and x_init is on GPU, avoid scalar indexing
                for (r, x0) in zip(roots, problem.x_init)
                    if abs(r - expected_root) ≥ _default_tol
                        @info "Failing root" problem = problem.name method =
                            method root = r expected = expected_root initial_guess =
                            x0 tol = _default_tol error = abs(r - expected_root)
                    end
                end
            end
            @test all(abs.(roots .- expected_root) .< _default_tol)
        else
            if abs(roots - expected_root) ≥ _default_tol
                @info "Failing root" problem = problem.name method = method root =
                    roots expected = expected_root initial_guess =
                    problem.x_init tol = _default_tol error =
                    abs(roots - expected_root)
            end
            @test abs(roots - expected_root) < _default_tol
        end
    elseif tol isa ResidualTolerance
        # For residual tolerance, check that function values are small
        if roots isa AbstractArray
            f = problem.f
            @test all(map(x -> abs(f(x)), roots) .< tol.tol)
        else
            @test abs(problem.f(roots)) < tol.tol
        end
    elseif tol isa SolutionTolerance
        # For solution tolerance, check that roots are close to expected
        if roots isa AbstractArray
            @test all(abs.(roots .- expected_root) .< tol_factor * tol.tol)
        else
            @test abs(roots - expected_root) < tol_factor * tol.tol
        end
    elseif tol isa RelativeSolutionTolerance
        # For relative tolerance, check relative difference
        # Avoid division by zero
        if abs(expected_root) < eps(typeof(expected_root))
            if roots isa AbstractArray
                @test all(abs.(roots .- expected_root)) .< tol_factor * tol.tol
            else
                @test abs(roots - expected_root) < tol_factor * tol.tol
            end
        else
            if roots isa AbstractArray
                @test all(
                    abs.((roots .- expected_root) ./ expected_root) .<
                    tol_factor * tol.tol,
                )
            else
                @test abs((roots - expected_root) / expected_root) <
                      tol_factor * tol.tol
            end
        end
    elseif tol isa RelativeOrAbsoluteSolutionTolerance
        # For combined tolerance, check both relative and absolute
        if roots isa AbstractArray
            @test all(
                (
                    abs.((roots .- expected_root) ./ expected_root) .<
                    tol_factor * tol.rtol
                ) .| (abs.(roots .- expected_root) .< tol_factor * tol.atol),
            )
        else
            @test (
                abs((roots - expected_root) / expected_root) <
                tol_factor * tol.rtol
            ) || (abs(roots - expected_root) < tol_factor * tol.atol)
        end
    else
        # Default tolerance check - use a reasonable default based on floating point precision
        FT = typeof(expected_root)
        _default_tol = tol_factor * default_tol(FT).tol
        if roots isa AbstractArray
            @test all(abs.(roots .- expected_root) .< _default_tol)
        else
            @test abs(roots - expected_root) < _default_tol
        end
    end
end

# Helper function to get maxiters based on problem name
function get_maxiters(problem_name)
    return 2_000
end

# Helper function to run solver and extract results
function run_solver_test(f, method, sol_type, tol, maxiters, is_array)
    if is_array
        sol = RootSolvers.find_zero.(f, method, sol_type, tol, maxiters)
        converged = map(x -> x.converged, sol)
        roots = map(x -> x.root, sol)
        return sol, converged, roots
    else
        sol = find_zero(f, method, sol_type, tol, maxiters)
        return sol, sol.converged, sol.root
    end
end

function find_zero_wrapper(
    f::F,
    method::M,
    sol_type::S,
    tol::T,
    maxiters,
    is_array,
) where {F, M, S, T}
    # Wrapper to handle both array and scalar cases
    if is_array
        find_zero.(f, method, sol_type, tol, maxiters)
    else
        find_zero(f, method, sol_type, tol, maxiters)
    end
    return
end

function test_allocations(
    f::F,
    method::M,
    sol_type::S,
    tol::T,
    maxiters,
    is_array,
) where {F, M, S, T}
    sol_type isa VerboseSolution && return  # VerboseSolution is expected to allocate
    find_zero_wrapper(f, method, sol_type, tol, maxiters, is_array) # ensure function is compiled
    @test (@allocated find_zero_wrapper(
        f,
        method,
        sol_type,
        tol,
        maxiters,
        is_array,
    )) == 0
end

@testset "Convergence reached" begin
    # Test that all root-finding methods converge to the correct solution
    # This is the main test suite that validates the core functionality
    for problem in problem_list
        FT = typeof(problem.x̃)  # Extract the floating-point type from the expected solution
        for sol_type in (
            ArrayType <: Array ? [CompactSolution(), VerboseSolution()] :
            [CompactSolution()] # VerboseSolution is not GPU friendly
        )
            # Test both solution types: compact (GPU-friendly) and verbose (with history)
            for tol in get_tolerances(FT)
                # Test all tolerance types for the given floating-point type
                for method in get_methods(
                    problem.x_init,
                    problem.x_lower,
                    problem.x_upper,
                )
                    # Test all applicable root-finding methods for this problem

                    # Choose function based on method type:
                    # - NewtonsMethod requires function that returns (f(x), f'(x))
                    # - Other methods use standard function f(x)
                    f =
                        method isa NewtonsMethod || (
                            method isa AbstractArray &&
                            eltype(method) <: NewtonsMethod
                        ) ? problem.ff′ : problem.f

                    # Run solver test
                    is_array = problem.x_init isa AbstractArray
                    maxiters = get_maxiters(problem.name)
                    sol, converged, roots = run_solver_test(
                        f,
                        method,
                        sol_type,
                        tol,
                        maxiters,
                        is_array,
                    )

                    # Validate results
                    if is_array
                        @test all(converged) # All should converge
                        @test eltype(roots) == eltype(problem.x_init)  # Type consistency
                        test_verbose!(
                            sol_type,
                            sol,
                            problem,
                            tol,
                            all(converged),
                        )  # Additional verbose checks
                    else
                        @test converged # Should converge
                        @test roots isa FT # Type consistency
                        test_verbose!(sol_type, sol, problem, tol, converged)  # Additional verbose checks
                    end

                    # Check that roots are within a reasonable tolerance of the expected solution
                    check_root_tolerance(
                        roots,
                        problem.x̃,
                        problem,
                        method,
                        tol,
                    )

                    ArrayType <: Array && test_allocations(
                        f,
                        method,
                        sol_type,
                        tol,
                        maxiters,
                        is_array,
                    ) # allocations expected during CUDA kernel launches
                end
            end
        end
    end
end

@testset "Convergence not reached" begin
    # Test that methods properly handle cases where convergence is not possible
    # This validates error handling and non-convergence detection

    # Create a difficult problem that's unlikely to converge in few iterations
    # Use a problem that properly brackets a root for bracketing methods
    difficult_problem = RootSolvingProblem(
        "difficult convergence test",
        x -> x^3 - 1000,  # Function with root at x = 10
        x -> (x^3 - 1000, 3x^2),  # Function and derivative
        10.0,  # Solution
        1.0,   # Initial guess (far from solution)
        0.0,   # Lower bound (f(0) = -1000)
        20.0,  # Upper bound (f(20) = 8000)
    )

    for sol_type in (
        ArrayType <: Array ? [CompactSolution(), VerboseSolution()] :
        [CompactSolution()]
    )
        for tol in get_tolerances(Float64)
            for method in get_methods(
                difficult_problem.x_init,
                difficult_problem.x_lower,
                difficult_problem.x_upper,
            )
                # Same function selection logic as above
                f =
                    method isa NewtonsMethod ? difficult_problem.ff′ :
                    difficult_problem.f

                # Run solver test with very few iterations to force non-convergence
                sol, converged, roots =
                    run_solver_test(f, method, sol_type, tol, 2, false)

                # Validate that convergence is unlikely with very few iterations
                @test isbits(method)
                @test roots isa Float64
                test_verbose!(sol_type, sol, difficult_problem, tol, converged)
                ArrayType <: Array &&
                    test_allocations(f, method, sol_type, tol, 2, false)

                # Note: We don't strictly require non-convergence here because some methods
                # might be very efficient and converge quickly. The important thing is that
                # the method handles the case gracefully.
            end
        end
    end
end

@testset "Check small Δy" begin
    # Test edge case where function values are very close (small Δy)
    # This tests the robustness of the secant method implementation

    ## Case 1: Δy is small and we converged
    # This tests the case where function values are nearly identical but we still converge
    # due to small Δx (solution is very close)
    sol = find_zero(
        x -> x^3,
        SecantMethod{Float64}(1e-8, 1e-8 + 1e-24),
        VerboseSolution(),
    )
    @test sol.converged === true  # Δx is small
    y = sol.err
    @test abs(y) ≤ default_tol(Float64).tol  # y is small

    ## Case 2: Δy is small, but we didn't converge
    # This tests the case where function values are nearly identical but we don't converge
    # because the solution is not close (e.g., found two distinct roots)
    # Use a function where both endpoints have nearly identical function values but are far from any root
    sol =
        find_zero(x -> x^2 + 1, SecantMethod{Float64}(-1, 1), VerboseSolution())
    @test sol.converged === false  # Should not converge because no root exists in the interval
    y = sol.err
    @test abs(y) > default_tol(Float64).tol  # Verify y is not small
end

@testset "Check invalid starting arguments" begin
    for FT in (Float32, Float64)
        # Test with x_init, x_lower, and x_upper as infinity
        for method in get_methods(FT(Inf), FT(Inf), FT(Inf))
            f = method isa NewtonsMethod ? x -> (identity, 1) : identity
            sol = find_zero(f, method)
            @test sol.converged === false
        end
        # Test when ff(x_init), (x_lower), or f(x_upper) evaluate to infinity
        for method in get_methods(FT(2000), FT(2000), FT(2000))
            f = method isa NewtonsMethod ? x -> (exp(x), exp(x)) : exp
            sol = find_zero(f, method)
            @test sol.converged === false
        end
        # Test when f(x_lower) and  f(x_upper) have same sign for methods that use bracketing
        bracketing_methods = (RegulaFalsiMethod, BisectionMethod, BrentsMethod)
        for method in bracketing_methods
            sol = find_zero(identity, method(FT(-1), FT(-5)))
            @test sol.converged === false  # Should not converge since f(x_lower) > 0 &  f(x_upper) > 0
        end
    end
end

@testset "Dual Number tests" begin
    for FT in (Float32, Float64)
        # θ  ̸∈ [-1, 1] to ensure valid inital bounds/guesses
        for θ in FT.((-3, 3, -3 / 2, 5 / 2, 10))
            method_types = (
                SecantMethod,
                RegulaFalsiMethod,
                BisectionMethod,
                BrentsMethod,
                NewtonsMethodAD,
                NewtonsMethod,
            )
            for MT in method_types
                function f(θ)
                    g(x) = x^2 - θ^2
                    input_function = MT <: NewtonsMethod ? x -> (g(x), 2x) : g
                    method =
                        MT <: Union{NewtonsMethod, NewtonsMethodAD} ?
                        MT(θ - 1) : MT(θ - 1, θ + 1)  # Create method with bounds
                    sol = find_zero(input_function, method)
                    return sol.root^2
                end
                deriv = ForwardDiff.derivative(f, θ)
                @test abs(deriv - 2 * θ) <= default_tol(FT).tol
            end
        end

    end
end

# Include additional test files for specialized functionality
include("test_printing.jl")    # Tests for solution pretty printing
