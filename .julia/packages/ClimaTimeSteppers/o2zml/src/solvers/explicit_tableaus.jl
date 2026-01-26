export ExplicitTableau, ExplicitAlgorithm
export SSP22Heuns, SSP33ShuOsher, RK4

abstract type ERKAlgorithmName <: AbstractAlgorithmName end

"""
    ExplicitTableau(; a, b, c)

A wrapper for an explicit Butcher tableau. Only `a` is a required argument; the
default value for `b` assumes that the algorithm is FSAL (first same as last),
and the default value for `c` assumes that it is internally consistent. The
matrix `a` must be strictly lower triangular.
"""
struct ExplicitTableau{A <: SPCO, B <: SPCO, C <: SPCO}
    a::A # matrix of size s×s
    b::B # vector of length s
    c::C # vector of length s
end
ExplicitTableau(args...) = ExplicitTableau(map(x -> SparseCoeffs(x), args)...)
function ExplicitTableau(; a, b = a[end, :], c = vec(sum(a; dims = 2)))
    @assert all(iszero, UpperTriangular(a))
    b, c = promote(b, c) # TODO: add generic promote_eltype
    return ExplicitTableau(a, b, c)
end

"""
    ExplicitAlgorithm(tableau, [constraint])
    ExplicitAlgorithm(name, [constraint])

Constructs an explicit algorithm for solving ODEs, with an optional name and
constraint. The first constructor accepts any `ExplicitTableau` and an optional
constraint, leaving the algorithm unnamed. The second constructor automatically
determines the tableau and the default constraint from the algorithm name,
which must be an `ERKAlgorithmName`.

Note that using an `ExplicitAlgorithm` is merely a shorthand for using an
`IMEXAlgorithm` with the same tableau for explicit and implicit tendencies (and
without Newton's method).
"""
ExplicitAlgorithm(tableau::ExplicitTableau, constraint = Unconstrained()) =
    IMEXAlgorithm(constraint, nothing, IMEXTableau(tableau), nothing)
ExplicitAlgorithm(name::ERKAlgorithmName, constraint = default_constraint(name)) =
    IMEXAlgorithm(constraint, name, IMEXTableau(name), nothing)

IMEXTableau(name::ERKAlgorithmName) = IMEXTableau(ExplicitTableau(name))
IMEXTableau((; a, b, c)::ExplicitTableau) = IMEXTableau(a, b, c, a, b, c)

################################################################################

abstract type SSPRKAlgorithmName <: ERKAlgorithmName end

default_constraint(::SSPRKAlgorithmName) = SSP()

"""
    SSP22Heuns

An SSPRK algorithm from [SO1988](@cite), with 2 stages and 2nd order accuracy.
Also called Heun's method ([Heun1900](@cite)).
"""
struct SSP22Heuns <: SSPRKAlgorithmName end
function ExplicitTableau(::SSP22Heuns)
    return ExplicitTableau(; a = @SArray([
        0 0
        1 0
    ]), b = @SArray([1 / 2, 1 / 2]))
end

"""
    SSP33ShuOsher

An SSPRK algorithm from [SO1988](@cite), with 3 stages and 3rd order accuracy.
"""
struct SSP33ShuOsher <: SSPRKAlgorithmName end
function ExplicitTableau(::SSP33ShuOsher)
    return ExplicitTableau(; a = @SArray([
        0 0 0
        1 0 0
        1/4 1/4 0
    ]), b = @SArray([1 / 6, 1 / 6, 2 / 3]))
end

"""
    RK4

The RK4 algorithm from [SM2003](@cite), a Runge-Kutta method with
4 stages and 4th order accuracy.
"""
struct RK4 <: ERKAlgorithmName end
function ExplicitTableau(::RK4)
    return ExplicitTableau(; a = @SArray([
        0 0 0 0
        1/2 0 0 0
        0 1/2 0 0
        0 0 1 0
    ]), b = @SArray([1 / 6, 1 / 3, 1 / 3, 1 / 6]))
end
