export NewtonsMethod, KrylovMethod
export JacobianFreeJVP, ForwardDiffJVP, ForwardDiffStepSize
export ForwardDiffStepSize1, ForwardDiffStepSize2, ForwardDiffStepSize3
export ForcingTerm, ConstantForcing, EisenstatWalkerForcing

# TODO: Implement AutoDiffJVP after ClimaAtmos's cache is moved from f! to x (so
#       that we only need to define Dual.(x), and not also make_dual(f!)).
# TODO: Define ktypeof(::FieldVector) so that it returns CuVector for a
#       FieldVector that contains CuArrays.
# TODO: Consider implementing line search backtracking to adjust Δx in
#       NewtonsMethod (scale Δx[n] by some λ[n] that reduces ‖f(x[n + 1])‖).
# TODO: Consider implementing BroydensMethod, BadBroydensMethod, and other
#       automatic Jacobian approximations.

abstract type AbstractVerbosity end
struct Verbose <: AbstractVerbosity end
struct Silent <: AbstractVerbosity end
is_verbose(v::AbstractVerbosity) = v isa Verbose

const KrylovWorkspace = @static pkgversion(Krylov) < v"0.10" ? Krylov.KrylovSolver : Krylov.KrylovWorkspace
const krylov_solve! = @static pkgversion(Krylov) < v"0.10" ? Krylov.solve! : Krylov.krylov_solve!
const GmresWorkspace = @static pkgversion(Krylov) < v"0.10" ? Krylov.GmresSolver : Krylov.GmresWorkspace

"""
    ForwardDiffStepSize

Computes a default step size for the forward difference approximation of the
Jacobian-vector product. This is done by calling `default_step(Δx, x)`, where
`default_step` is a `ForwardDiffStepSize`.
"""
abstract type ForwardDiffStepSize end

"""
    ForwardDiffStepSize1()

A `ForwardDiffStepSize` that is derived based on the notes here:
https://web.engr.oregonstate.edu/~webbky/MAE4020_5020_files/Section%204%20Roundoff%20and%20Truncation%20Error.pdf. Although it is not often
used with Newton-Krylov methods in practice, it can provide some intuition for
for how to set the value of `step_adjustment` in a `ForwardDiffJVP`.

The first-order Taylor series expansion of `f(x + ε * Δx)` around `x` is
    `f(x + ε * Δx) = f(x) + j(x) * (ε * Δx) + e_trunc(x, ε * Δx)`,
where `j(x) = f'(x)` and `e_trunc` is the expansion's truncation error.
Due to roundoff error, we are unable to directly compute the value of `f(x)`;
instead, we can only determine `f̂(x)`, where
    `f(x) = f̂(x) + e_round(x)`.
Substituting this into the expansion tells us that
    `f̂(x + ε * Δx) + e_round(x + ε * Δx) =
        f̂(x) + e_round(x) + j(x) * (ε * Δx) + e_trunc(x, ε * Δx)`.
Rearranging this gives us the Jacobian-vector product
    `j(x) * Δx = (f̂(x + ε * Δx) - f̂(x)) / ε - e_trunc(x, ε * Δx) / ε +
                 (e_round(x + ε * Δx) - e_round(x)) / ε`.
So, the normed error of the forward difference approximation of this product is
    `‖error‖ = ‖(f̂(x + ε * Δx) - f̂(x)) / ε - j(x) * Δx‖ =
             = ‖e_trunc(x, ε * Δx) - e_round(x + ε * Δx) + e_round(x)‖ / ε`.
We can use the triangle inequality to get the upper bound
    `‖error‖ ≤
        (‖e_trunc(x, ε * Δx)‖ + ‖e_round(x + ε * Δx)‖ + ‖e_round(x)‖) / ε`.
If `ε` is sufficiently small, we can approximate
    `‖e_round(x + ε * Δx)‖ ≈ ‖e_round(x)‖`.
This simplifies the upper bound to
    `‖error‖ ≤ (‖e_trunc(x, ε * Δx)‖ + 2 * ‖e_round(x)‖) / ε`.

From Taylor's theorem (for multivariate vector-valued functions), the truncation
error of the first-order expansion is bounded by
    `‖e_trunc(x, ε * Δx)‖ ≤ (sup_{x̂ ∈ X} ‖f''(x̂)‖) / 2 * ‖ε * Δx‖^2`,
where `X` is a closed ball around `x` that contains `x + ε * Δx` (see
https://math.stackexchange.com/questions/3478229 for a proof of this).
Let us define the value
    `S = ‖f(x)‖ / sup_{x̂ ∈ X} ‖f''(x̂)‖`.
By default, we will assume that `S ≈ 1`, but we will let users pass other values
to indicate the "smoothness" of `f(x)` (a large value of `S` should indicate
that the Hessian tensor of `f(x)` has a small norm compared to `f(x)` itself).
We then have that
    `‖e_trunc(x, ε * Δx)‖| ≤ ε^2 / (2 * S) * ‖Δx‖^2 * ‖f(x)‖`.

If only the last bit in each component of `f(x)` can be altered by roundoff
error, then the `i`-th component of `e_round(x)` is bounded by
    `|e_round(x)[i]| ≤ eps(f(x)[i])`.
More generally, we can assume that there is some constant `R` (by default, we
will assume that `R ≈ 1`) such that
    `|e_round(x)[i]| ≤ R * eps(f(x)[i])`.
We can also make the approximation (which is accurate to within `eps(FT)`)
    `eps(f(x)[i]) ≈ eps(FT) * |f(x)[i]|`.
This implies that
    `|e_round(x)[i]| ≤ R * eps(FT) * |f(x)[i]|`.
Since this is true for every component of `e_round(x)` and `f(x)`, we find that
    `‖e_round(x)‖ ≤ R * eps(FT) * ‖f(x)‖`.

Substituting the bounds on the truncation and roundoff errors into the bound on
the overall error gives us
    `‖error‖ ≤ ε / (2 * S) * ‖Δx‖^2 * ‖f(x)‖ + 2 / ε * R * eps(FT) * ‖f(x)‖`.
Differentiating the right-hand side with respect to `ε` and setting the result
equal to 0 (and noting that the second derivative is always positive) tells us
that this upper bound is minimized when
    `ε = step_adjustment * sqrt(eps(FT)) / ‖Δx‖`,
where `step_adjustment = 2 * sqrt(S * R)`.
By default, we will assume that `step_adjustment = 1`, but it should be made
larger when `f` is very smooth or has a large roundoff error.

Note that, if we were to replace the forward difference approximation in the
derivation above with a central difference approximation, the square root would
end up being replaced with a cube root (or, more generally, with an `n`-th root
for a finite difference approximation of order `n - 1`).
"""
struct ForwardDiffStepSize1 <: ForwardDiffStepSize end
(::ForwardDiffStepSize1)(Δx, x) = sqrt(eps(eltype(Δx))) / norm(Δx)

"""
    ForwardDiffStepSize2()

A `ForwardDiffStepSize` that is described in the paper "Jacobian-free
Newton–Krylov methods: a survey of approaches and applications" by D.A. Knoll
and D.E. Keyes. According to the paper, this is the step size used by the
Fortran package NITSOL.
"""
struct ForwardDiffStepSize2 <: ForwardDiffStepSize end
(::ForwardDiffStepSize2)(Δx, x) = sqrt(eps(eltype(Δx)) * (1 + norm(x))) / norm(Δx)

"""
    ForwardDiffStepSize3()

A `ForwardDiffStepSize` that is described in the paper "Jacobian-free
Newton–Krylov methods: a survey of approaches and applications" by D.A. Knoll
and D.E. Keyes. According to the paper, this is the average step size one gets
when using a certain forward difference approximation for each Jacobian element.
"""
struct ForwardDiffStepSize3 <: ForwardDiffStepSize end
(::ForwardDiffStepSize3)(Δx, x) = sqrt(eps(eltype(Δx))) * sum(x_i -> 1 + abs(x_i), x) / (length(x) * norm(Δx))

"""
    JacobianFreeJVP

Computes the Jacobian-vector product `j(x[n]) * Δx[n]` for a Newton-Krylov
method without directly using the Jacobian `j(x[n])`, and instead only using
`x[n]`, `f(x[n])`, and other function evaluations `f(x′)`. This is done by
calling `jvp!(::JacobianFreeJVP, cache, jΔx, Δx, x, f!, f, prepare_for_f!)`.
The `jΔx` passed to a Jacobian-free JVP is modified in-place. The `cache` can
be obtained with `allocate_cache(::JacobianFreeJVP, x_prototype)`, where
`x_prototype` is `similar` to `x` (and also to `Δx` and `f`).
"""
abstract type JacobianFreeJVP end

"""
    ForwardDiffJVP(; default_step = ForwardDiffStepSize3(), step_adjustment = 1)

Computes the Jacobian-vector product using the forward difference approximation
`j(x) * Δx = (f(x + ε * Δx) - f(x)) / ε`, where
`ε = step_adjustment * default_step(Δx, x)`.
"""
Base.@kwdef struct ForwardDiffJVP{S <: ForwardDiffStepSize, T} <: JacobianFreeJVP
    default_step::S = ForwardDiffStepSize3()
    step_adjustment::T = 1
end

allocate_cache(::ForwardDiffJVP, x_prototype) = (; x2 = zero(x_prototype), f2 = zero(x_prototype))

function jvp!(alg::ForwardDiffJVP, cache, jΔx, Δx, x, f!, f, prepare_for_f!)
    (; default_step, step_adjustment) = alg
    (; x2, f2) = cache
    FT = eltype(x)
    ε = FT(step_adjustment) * default_step(Δx, x)
    @. x2 = x + ε * Δx
    isnothing(prepare_for_f!) || prepare_for_f!(x2)
    f!(f2, x2)
    @. jΔx = (f2 - f) / ε
end

"""
    ForcingTerm

Computes the value of `rtol[n]` for a Newton-Krylov method. This is done by
calling `get_rtol!(::ForcingTerm, cache, f, n)`, which returns `rtol[n]`. The `cache`
can be obtained with `allocate_cache(::ForcingTerm, x_prototype)`, where
`x_prototype` is `similar` to `f`.

For a detailed discussion of forcing terms and their convergence guarantees, see
"Choosing the Forcing Terms in an Inexact Newton Method" by S.C. Eisenstat and
H.F. Walker (http://softlib.rice.edu/pub/CRPC-TRs/reports/CRPC-TR94463.pdf).
"""
abstract type ForcingTerm end

"""
    ConstantForcing(rtol)

A `ForcingTerm` that always returns the value `rtol`, which must be in the
interval `[0, 1)`. If `x` and `f!` satisfy certain assumptions, this forcing
term guarantees that the Newton-Krylov method will converge linearly with an
asymptotic rate of at most `rtol`. If `rtol` is 0 (or `eps(FT)`), this forces
the approximation of `Δx[n]` to be exact (or exact to within machine precision)
and guarantees that the Newton-Krylov method will converge quadratically. Note
that, although a smaller value of `rtol` guarantees faster asymptotic
convergence, it also leads to a higher probability of oversolving.
"""
struct ConstantForcing{T} <: ForcingTerm
    rtol::T
end

allocate_cache(::ConstantForcing, x_prototype) = (;)

function get_rtol!(alg::ConstantForcing, cache, f, n)
    FT = eltype(f)
    return FT(alg.rtol)
end

"""
    EisenstatWalkerForcing(;
        initial_rtol = 0.5,
        γ = 1,
        α = 2,
        min_rtol_threshold = 0.1,
        max_rtol = 0.9,
    )

The `ForcingTerm` called "Choice 2" in the paper "Choosing the Forcing Terms in
an Inexact Newton Method" by S.C. Eisenstat and H.F. Walker. The values of
`initial_rtol`, `min_rtol_threshold`, and `max_rtol` must be in the interval
`[0, 1)`, the value of `γ` must be in the interval `[0, 1]`, and the value of
`α` must be in the interval `(1, 2]`. These values can all be tuned to prevent
the Newton-Krylov method from oversolving. If `x` and `f!` satisfy certain
assumptions, this forcing term guarantees that the Newton-Krylov method will
converge with order `α`. Note that, although a larger value of `α` guarantees a
higher convergence order, it also leads to a higher probability of oversolving.

This forcing term was implemented instead of the one called "Choice 1" because
it has a significantly simpler implementation---it only requires the value of
`‖f(x[n])‖` to compute `rtol[n]`, whereas "Choice 1" also requires the norm of
the final residual from the Krylov solver.
"""
Base.@kwdef struct EisenstatWalkerForcing{T1, T2, T3, T4, T5} <: ForcingTerm
    initial_rtol::T1 = 0.5
    γ::T2 = 1
    α::T3 = 2
    min_rtol_threshold::T4 = 0.1
    max_rtol::T5 = 0.9
end

function allocate_cache(::EisenstatWalkerForcing, x_prototype)
    FT = eltype(x_prototype)
    return (; prev_norm_f = Ref{FT}(), prev_rtol = Ref{FT}())
end

function get_rtol!(alg::EisenstatWalkerForcing, cache, f, n)
    (; initial_rtol, γ, α, min_rtol_threshold, max_rtol) = alg
    (; prev_norm_f, prev_rtol) = cache
    FT = eltype(f)
    norm_f = norm(f)
    if n == 0
        rtol = FT(initial_rtol)
    else
        α = α isa Integer ? α : FT(α)
        rtol = FT(γ) * (norm_f / prev_norm_f[])^α
        min_rtol = FT(γ) * prev_rtol[]^α
        if min_rtol > FT(min_rtol_threshold)
            rtol = max(rtol, min_rtol)
        end
    end
    rtol = min(rtol, FT(max_rtol))
    prev_norm_f[] = norm_f
    prev_rtol[] = rtol
    return rtol
end

"""
    KrylovMethodDebugger

Prints information about the Jacobian matrix `j` and the preconditioner `M` (if
it is available) that are passed to a Krylov method. This is done by calling
`print_debug!(::KrylovMethodDebugger, cache, j, M)`. The `cache` can be obtained with
`allocate_cache(::KrylovMethodDebugger, x_prototype)`, where `x_prototype` is
`similar` to `x`.
"""
abstract type KrylovMethodDebugger end

"""
    PrintConditionNumber()

Prints the condition number of the Jacobian matrix `j`, and, if a preconditioner
`M` is available, also prints the condition number of `inv(M) * j` (i.e., the
matrix that actually gets "inverted" by the Krylov method). This requires
computing dense representations of `j` and `inv(M) * j`, which is likely to be
significantly slower than the Krylov method itself.
"""
struct PrintConditionNumber <: KrylovMethodDebugger end

function allocate_cache(::PrintConditionNumber, x_prototype)
    l = length(x_prototype)
    FT = eltype(x_prototype)
    return (;
        dense_vector = Array{FT}(undef, l),
        dense_j = Array{FT}(undef, l, l),
        dense_inv_M = Array{FT}(undef, l, l),
        dense_inv_M_j = Array{FT}(undef, l, l),
    )
end

print_debug!(::Nothing, cache, j, M) = nothing

function print_debug!(::PrintConditionNumber, cache, j, M)
    (; dense_vector, dense_j, dense_inv_M, dense_inv_M_j) = cache
    dense_matrix_from_operator!(dense_j, dense_vector, j)
    if M === I
        @info "Condition number = $(cond(dense_j))"
    else
        dense_inverse_matrix_from_operator!(dense_inv_M, dense_vector, M)
        mul!(dense_inv_M_j, dense_inv_M, dense_j)
        @info "Condition number = $(cond(dense_inv_M_j)) ($(cond(dense_j)) \
               without the preconditioner)"
    end
end

# Like Matrix(op::AbstractLinearOperator) from LinearOperators.jl, but in-place.
function dense_matrix_from_operator!(dense_matrix, dense_vector, op)
    n_columns = size(dense_matrix)[2]
    dense_vector .= 0
    for column in 1:n_columns
        dense_vector[column] = 1
        mul!(view(dense_matrix, :, column), op, dense_vector)
        dense_vector[column] = 0
    end
end

# Same as dense_matrix_from_operator!, but with ldiv! instead of mul!.
function dense_inverse_matrix_from_operator!(dense_inv_matrix, dense_vector, op)
    n_columns = size(dense_inv_matrix)[2]
    dense_vector .= 0
    for column in 1:n_columns
        dense_vector[column] = 1
        ldiv!(view(dense_inv_matrix, :, column), op, dense_vector)
        dense_vector[column] = 0
    end
end

"""
    KrylovMethod(;
        type = Val(Krylov.GmresSolver),
        jacobian_free_jvp = nothing,
        forcing_term = ConstantForcing(0),
        args = (20,),
        kwargs = (;),
        solve_kwargs = (;),
        disable_preconditioner = false,
        verbose = Silent(),
        debugger = nothing,
    )

Finds an approximation `Δx[n] ≈ j(x[n]) \\ f(x[n])` for Newton's method such
that `‖f(x[n]) - j(x[n]) * Δx[n]‖ ≤ rtol[n] * ‖f(x[n])‖`, where `rtol[n]` is the
value of the forcing term on iteration `n`. This is done by calling
`solve_krylov!(::KrylovMethod, cache, Δx, x, f!, f, n, prepare_for_f!, j = nothing)`,
where `f` is `f(x[n])` and, if it is specified, `j` is either `j(x[n])` or an
approximation of `j(x[n])`. The `Δx` passed to a Krylov method is modified in-place.
The `cache` can be obtained with `allocate_cache(::KrylovMethod, x_prototype)`,
where `x_prototype` is `similar` to `x` (and also to `Δx` and `f`).

This is primarily a wrapper for a `Krylov.KrylovSolver` from `Krylov.jl`. In
`allocate_cache`, the solver is constructed with
`solver = type(l, l, args..., Krylov.ktypeof(x_prototype); kwargs...)`, where
`l = length(x_prototype)` and `Krylov.ktypeof(x_prototype)` is a subtype of
`DenseVector` that can be used to store `x_prototype`. By default, the solver
is a `Krylov.GmresSolver` with a Krylov subspace size of 20 (the default Krylov
subspace size for this solver in `Krylov.jl`). In `solve_krylov!`, the solver is run with
`Krylov.solve!(solver, opj, f; M, ldiv, atol, rtol, verbose, solve_kwargs...)`.
The solver's type can be changed by specifying a different value for `type`,
though this value has to be wrapped in a `Val` to avoid runtime compilation.

In the call to `Krylov.solve!`, `opj` is a `LinearOperator` that represents
`j(x[n])`, which the solver uses by evaluating `mul!(jΔx, opj, Δx)`. If a
Jacobian-free JVP (Jacobian-vector product) is specified, it gets used to
construct `opj` and to evaluate the calls to `mul!`; otherwise, `j` itself gets
used to construct `opj`, and the calls to `mul!` simplify to `mul!(jΔx, j, Δx)`.

If a Jacobian-free JVP and `j` are both specified, and if
`disable_preconditioner` is set to `false`, `j` is treated as an approximation
of `j(x[n])` and is used as the (left) preconditioner `M` in order to speed up
the solver; otherwise, the preconditioner is simply set to the identity matrix
`I`. The keyword argument `ldiv` is set to `true` so that the solver calls
`ldiv!(Δx′, M, f′)` instead of `mul!(Δx′, M, f′)`, where `Δx′` and `f′` denote
internal variables of the solver that roughly correspond to `Δx` and `f`. In
other words, setting `ldiv` to `true` makes the solver treat `M` as an
approximation of `j` instead of as the inverse of an approximation of `j`.

The keyword argument `atol` is set to 0 because, if it is set to some other
value, the inequality `‖f(x[n]) - j(x[n]) * Δx[n]‖ ≤ rtol[n] * ‖f(x[n])‖`
changes to `‖f(x[n]) - j(x[n]) * Δx[n]‖ ≤ rtol[n] * ‖f(x[n])‖ + atol`, which
eliminates any convergence guarantees provided by the forcing term (in order for
the Newton-Krylov method to converge, the right-hand side of this inequality
must approach 0 as `n` increases, which cannot happen if `atol` is not 0).

All of the arguments and keyword arguments used to construct and run the solver
can be modified using `args`, `kwargs`, and `solve_kwargs`. So, the default
behavior of this wrapper can be easily overwritten, and any features of
`Krylov.jl` that are not explicitly covered by this wrapper can still be used.

If `verbose` is `true`, the residual `‖f(x[n]) - j(x[n]) * Δx[n]‖` is printed on
each iteration of the Krylov method. If a debugger is specified, it is run
before the call to `Kyrlov.solve!`.
"""
Base.@kwdef struct KrylovMethod{
    T <: Val{<:KrylovWorkspace},
    J <: Union{Nothing, JacobianFreeJVP},
    F <: ForcingTerm,
    A <: Tuple,
    K <: NamedTuple,
    S <: NamedTuple,
    V <: AbstractVerbosity,
    D <: Union{Nothing, KrylovMethodDebugger},
}
    type::T = Val(GmresWorkspace)
    jacobian_free_jvp::J = nothing
    forcing_term::F = ConstantForcing(0)
    args::A = ()
    kwargs::K = (; memory = 20)
    solve_kwargs::S = (;)
    disable_preconditioner::Bool = false
    verbose::V = Silent()
    debugger::D = nothing
end

solver_type(::KrylovMethod{Val{T}}) where {T} = T

function allocate_cache(alg::KrylovMethod, x_prototype)
    (; jacobian_free_jvp, forcing_term, args, kwargs, debugger) = alg
    type = solver_type(alg)
    l = length(x_prototype)

    # Version 0.10 changed how the memory is set
    if pkgversion(Krylov) < v"0.10"
        args = isempty(args) ? (20,) : ()
        kwargs = haskey(kwargs, :memory) ? Base.structdiff(kwargs, NamedTuple{(:memory,)}((kwargs[:memory],))) : kwargs
    end

    return (;
        jacobian_free_jvp_cache = isnothing(jacobian_free_jvp) ? nothing :
                                  allocate_cache(jacobian_free_jvp, x_prototype),
        forcing_term_cache = allocate_cache(forcing_term, x_prototype),
        solver = type(l, l, args..., Krylov.ktypeof(x_prototype); kwargs...),
        debugger_cache = isnothing(debugger) ? nothing : allocate_cache(debugger, x_prototype),
    )
end

NVTX.@annotate function solve_krylov!(alg::KrylovMethod, cache, Δx, x, f!, f, n, prepare_for_f!, j = nothing)
    (; jacobian_free_jvp, forcing_term, solve_kwargs) = alg
    (; disable_preconditioner, debugger) = alg
    type = solver_type(alg)
    (; jacobian_free_jvp_cache, forcing_term_cache, solver, debugger_cache) = cache
    jΔx!(jΔx, Δx) =
        isnothing(jacobian_free_jvp) ? mul!(jΔx, j, Δx) :
        jvp!(jacobian_free_jvp, jacobian_free_jvp_cache, jΔx, Δx, x, f!, f, prepare_for_f!)
    opj = LinearOperator(eltype(x), length(x), length(x), false, false, jΔx!)
    M = disable_preconditioner || isnothing(j) || isnothing(jacobian_free_jvp) ? I : j
    print_debug!(debugger, debugger_cache, opj, M)
    ldiv = true
    atol = zero(eltype(Δx))
    rtol = get_rtol!(forcing_term, forcing_term_cache, f, n)
    verbose = Int(is_verbose(alg.verbose))
    krylov_solve!(solver, opj, f; M, ldiv, atol, rtol, verbose, solve_kwargs...)
    iter = solver.stats.niter
    if !solver.stats.solved
        str1 = isnothing(j) ? () : ("the Jacobian",)
        str2 = isnothing(jacobian_free_jvp) ? () : ("the Jacobian-vector product",)
        str = join((str1..., str2...), " and/or ")
        if solver.stats.inconsistent
            @warn "$type detected that the Jacobian is singular on iteration \
                   $iter; if possible, try improving the approximation of $str"
        else
            @warn "$type did not converge within $iter iterations; if \
                   possible, try improving the approximation of $str, or try \
                   increasing the forcing term"
        end
    elseif iter == 0 && solver.stats.status != "x = 0 is a zero-residual solution"
        @warn "$type set Δx to 0 without running any iterations; if possible, \
               try decreasing the forcing term"
    end
    Δx .= Krylov.solution(solver)
end

"""
    NewtonsMethod(;
        max_iters = 1,
        update_j = UpdateEvery(NewNewtonIteration),
        krylov_method = nothing,
        convergence_checker = nothing,
        verbose = Silent(),
    )

Solves the equation `f(x) = 0`, using the Jacobian (or an approximation of the
Jacobian) `j(x) = f'(x)` if it is available. This is done by calling
`solve_newton!(::NewtonsMethod, cache, x, f!, j! = nothing)`, where `f!(f, x)` is a
function that sets `f(x)` in-place and, if it is specified, `j!(j, x)` is a
function that sets `j(x)` in-place. The `x` passed to Newton's method is
modified in-place, and its initial value is used as a starting guess for the
root. The `cache` can be obtained with
`allocate_cache(::NewtonsMethod, x_prototype, j_prototype = nothing)`, where
`x_prototype` is `similar` to `x` and `f(x)`, and, if it is specified,
`j_prototype` is `similar` to `j(x)`. In order for `j(x)` to be invertible, it
must be a square matrix, which implies that `x` and `f(x)` must be `similar` to
to each other.

Let `x[n]` denote the value of `x` on the `n`-th Newton iteration (with `x[0]`
denoting the initial value of `x`), and suppose that `x[n]` is sufficiently
close to some root `x̂` of `f(x)` to make the first-order approximation
    `f(x̂) ≈ f(x[n]) + j(x[n]) * (x̂ - x[n])`.
Since `f(x̂) = 0`, the error on the `n`-th iteration is roughly
    `x[n] - x̂ ≈ Δx[n]`, where `Δx[n] = j(x[n]) \\ f(x[n])`.
Newton's method sets `x[n + 1]` to be the value of `x̂` given by this
approximation:
    `x[n + 1] = x[n] - Δx[n]`.

If a Krylov method is specified, it gets used to compute the error
`Δx[n] = j(x[n]) \\ f(x[n])`; otherwise, the error is directly computed by
calling `ldiv!(Δx, j, f)`. If the Krylov method uses a Jacobian-free JVP
(Jacobian-vector product), `j_prototype` and `j!` do not need to be specified.
When Newton's method uses a Krylov method, it is called a "Newton-Krylov
method"; furthermore, when the Krylov method uses a Jacobian-free JVP, it is
called a "Jacobian-free Newton-Krylov method".

If `j_prototype` is specified, it should not be an `DenseMatrix`. If it is, it
has to be factorized with `lu` before `ldiv!` is called, which requires the
allocation of additional memory. Instead, `j_prototype` should be an object that
can directly be passed to `ldiv!`. For convenience, though, the use of an
`DenseMatrix` is supported. However, `Krylov.jl` does not provide such support
for its preconditioners, so, since the value computed with `j!` is used as a
preconditioner in Krylov methods with a Jacobian-free JVP, using such a Krylov
method requires specifying a `j_prototype` that can be passed to `ldiv!`.

If `j(x)` changes sufficiently slowly, `update_j` may be changed from
`UpdateEvery(NewNewtonIteration)` to some other `UpdateSignalHandler` that
gets triggered less frequently, such as `UpdateEvery(NewNewtonSolve)`. This
can be used to make the approximation `j(x[n]) ≈ j(x₀)`, where `x₀` is a
previous value of `x[n]` (possibly even a value from a previous `solve_newton!` of
Newton's method). When Newton's method uses such an approximation, it is called
the "chord method".

In addition, `update_j` can be set to an `UpdateSignalHandler` that gets
triggered by signals that originate outside of Newton's method, such as
`UpdateEvery(NewTimeStep)`. It is possible to send any signal for updating `j`
to Newton's method while it is not running by calling
`update!(::NewtonsMethod, cache, ::UpdateSignal, j!)`, where in this case
`j!(j)` is a function that sets `j` in-place without any dependence on `x`
(since `x` is not necessarily defined while Newton's method is not running, this
version of `j!` does not take `x` as an argument). This can be used to make the
approximation `j(x[n]) ≈ j₀`, where `j₀` can have an arbitrary value.

If a convergence checker is provided, it gets used to determine whether to stop
iterating on iteration `n` based on the value `x[n]` and its error `Δx[n]`;
otherwise, Newton's method iterates from `n = 0` to `n = max_iters`. If the
convergence checker determines that `x[n]` has not converged by the time
`n = max_iters`, a warning gets printed.

If `verbose` is set to `true`, the norms of `x[n]` and `Δx[n]` get printed on
every iteration. If there is no convergence checker, `Δx[n]` is not computed on
the last iteration, so its final norm is not printed.
"""
Base.@kwdef struct NewtonsMethod{
    U <: UpdateSignalHandler,
    K <: Union{Nothing, KrylovMethod},
    C <: Union{Nothing, ConvergenceChecker},
    V <: AbstractVerbosity,
}
    max_iters::Int = 1
    update_j::U = UpdateEvery(NewNewtonIteration)
    krylov_method::K = nothing
    convergence_checker::C = nothing
    verbose::V = Silent()
end

function allocate_cache(alg::NewtonsMethod, x_prototype, j_prototype = nothing)
    (; update_j, krylov_method, convergence_checker) = alg
    @assert !(isnothing(j_prototype) && (isnothing(krylov_method) || isnothing(krylov_method.jacobian_free_jvp)))
    return (;
        krylov_method_cache = isnothing(krylov_method) ? nothing : allocate_cache(krylov_method, x_prototype),
        convergence_checker_cache = isnothing(convergence_checker) ? nothing :
                                    allocate_cache(convergence_checker, x_prototype),
        Δx = zero(x_prototype),
        f = zero(x_prototype),
        j = isnothing(j_prototype) ? nothing : zero(j_prototype),
    )
end

solve_newton!(alg::NewtonsMethod, cache::Nothing, x, f!, j! = nothing, prepare_for_f! = nothing) = nothing

NVTX.@annotate function solve_newton!(alg::NewtonsMethod, cache, x, f!, j! = nothing, prepare_for_f! = nothing)
    (; max_iters, update_j, krylov_method, convergence_checker, verbose) = alg
    (; krylov_method_cache, convergence_checker_cache) = cache
    (; Δx, f, j) = cache
    if (!isnothing(j)) && needs_update!(update_j, NewNewtonSolve())
        j!(j, x)
    end
    for n in 1:max_iters
        # Compute Δx[n].
        if (!isnothing(j)) && needs_update!(update_j, NewNewtonIteration())
            j!(j, x)
        end
        f!(f, x)
        if isnothing(krylov_method)
            if j isa DenseMatrix
                ldiv!(Δx, lu(j), f) # Highly inefficient! Only used for testing.
            else
                ldiv!(Δx, j, f)
            end
        else
            solve_krylov!(krylov_method, krylov_method_cache, Δx, x, f!, f, n, prepare_for_f!, j)
        end
        is_verbose(verbose) && @info "Newton iteration $n: ‖x‖ = $(norm(x)), ‖Δx‖ = $(norm(Δx))"

        x .-= Δx
        # Update x[n] with Δx[n - 1], and exit the loop if Δx[n] is not needed.
        # Check for convergence if necessary.
        if is_converged!(convergence_checker, convergence_checker_cache, x, Δx, n)
            break
        elseif n < max_iters
            isnothing(prepare_for_f!) || prepare_for_f!(x)
        elseif is_verbose(verbose)
            @warn "Newton's method did not converge within $n iterations: ‖x‖ = $(norm(x)), ‖Δx‖ = $(norm(Δx))"
        end
    end
end
