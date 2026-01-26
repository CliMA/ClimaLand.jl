using ClimaTimeSteppers, LinearAlgebra, Random, Test
import ClimaTimeSteppers as CTS

function linear_equation(FT, n)
    rng = MersenneTwister(1)
    A = rand(rng, FT, n, n)
    b = rand(rng, FT, n)
    f!(f, x) = f .= A * x .- b
    j!(j, x) = j .= A
    x_exact = A \ b
    x_init = zeros(FT, n)
    return (f!, j!, x_exact, x_init)
end

function nonlinear_equation(FT, n)
    rng = MersenneTwister(1)
    A = rand(rng, FT, n, n)
    b = rand(rng, FT, n)
    f!(f, x) = f .= A * sinh.(A * x .- b) .- b
    j!(j, x) = j .= A * diagm(cosh.(A * x .- b)) * A
    x_exact = A \ (asinh.(A \ b) .+ b)
    x_init = zeros(FT, n)
    return (f!, j!, x_exact, x_init)
end

@testset "Newton's Method" begin
    # For each test case, step_adjustment is set to the value that minimizes the
    # number of Newton iterations and the normed error of the last Newton
    # iteration (for the specified value of rtol). In general, it is large when
    # the corresponding function f is very smooth or has a large roundoff error.
    for (is_linear, FT, n, step_adjustment) in (
        (true, Float32, 1, 100000),
        (true, Float32, 10, 1000000),
        (true, Float64, 1, 1000000000),
        (true, Float64, 10, 10000000000),
        (false, Float32, 1, 1),
        (false, Float32, 10, 1),
        (false, Float64, 1, 0.1),
        (false, Float64, 10, 10),
    )
        equation = is_linear ? linear_equation : nonlinear_equation
        max_iters = is_linear ? 1 : 10
        rtol = 100 * eps(FT)
        convergence_checker = ConvergenceChecker(; norm_condition = MaximumRelativeError(rtol))
        alg1 = NewtonsMethod(; max_iters, convergence_checker)
        alg2 = NewtonsMethod(;
            max_iters,
            krylov_method = KrylovMethod(; forcing_term = ConstantForcing(rtol)),
            convergence_checker,
        )
        alg3 = NewtonsMethod(;
            max_iters,
            krylov_method = KrylovMethod(;
                jacobian_free_jvp = ForwardDiffJVP(; step_adjustment),
                forcing_term = ConstantForcing(rtol),
            ),
            convergence_checker,
        )
        f!, j!, x_exact, x_init = equation(FT, n)
        for (alg, use_j) in ((alg1, true), (alg2, true), (alg3, false))
            x = copy(x_init)
            j_prototype = similar(x, length(x), length(x))
            cache = CTS.allocate_cache(alg, x, use_j ? j_prototype : nothing)
            CTS.solve_newton!(alg, cache, x, f!, use_j ? j! : nothing)
            @test norm(x .- x_exact) / norm(x_exact) < rtol
        end
    end
end

@testset "Dense Matrix From Operator" begin
    n = 10
    rng = MersenneTwister(1)
    A = rand(rng, n, n)
    vector = similar(A, n)
    matrix = similar(A, n, n)
    ClimaTimeSteppers.dense_matrix_from_operator!(matrix, vector, A)
    @test matrix == A
    ClimaTimeSteppers.dense_inverse_matrix_from_operator!(matrix, vector, lu(A))
    @test matrix ≈ inv(A)
end
