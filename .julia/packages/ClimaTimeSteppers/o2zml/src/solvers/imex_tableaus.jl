export IMEXTableau, IMEXAlgorithm
export ARS111, ARS121, ARS122, ARS233, ARS232, ARS222, ARS343, ARS443
export IMKG232a, IMKG232b, IMKG242a, IMKG242b, IMKG243a, IMKG252a, IMKG252b
export IMKG253a, IMKG253b, IMKG254a, IMKG254b, IMKG254c, IMKG342a, IMKG343a
export SSP222, SSP322, SSP332, SSP333, SSP433
export DBM453, HOMMEM1, ARK2GKC, ARK437L2SA1, ARK548L2SA2

abstract type IMEXARKAlgorithmName <: AbstractAlgorithmName end

"""
    IMEXTableau(; a_exp, b_exp, c_exp, a_imp, b_imp, c_imp)

A wrapper for an IMEX Butcher tableau (or, more accurately, a pair of Butcher
tableaus, one for explicit tendencies and the other for implicit tendencies).
Only `a_exp` and `a_imp` are required arguments; the default values for `b_exp`
and `b_imp` assume that the algorithm is FSAL (first same as last), and the
default values for `c_exp` and `c_imp` assume that it is internally consistent.

The explicit tableau must be strictly lower triangular, and the implicit tableau
must be lower triangular (only DIRK algorithms are currently supported).
"""
struct IMEXTableau{AE <: SPCO, BE <: SPCO, CE <: SPCO, AI <: SPCO, BI <: SPCO, CI <: SPCO}
    a_exp::AE # matrix of size s×s
    b_exp::BE # vector of length s
    c_exp::CE # vector of length s
    a_imp::AI # matrix of size s×s
    b_imp::BI # vector of length s
    c_imp::CI # vector of length s
end
IMEXTableau(args...) = IMEXTableau(map(x -> SparseCoeffs(x), args)...)

function IMEXTableau(;
    a_exp,
    b_exp = a_exp[end, :],
    c_exp = vec(sum(a_exp; dims = 2)),
    a_imp,
    b_imp = a_imp[end, :],
    c_imp = vec(sum(a_imp; dims = 2)),
)
    @assert all(iszero, UpperTriangular(a_exp))
    @assert all(iszero, UpperTriangular(a_imp) - Diagonal(a_imp))

    # TODO: add generic promote_eltype
    a_exp, a_imp = promote(a_exp, a_imp)
    b_exp, b_imp, c_exp, c_imp = promote(b_exp, b_imp, c_exp, c_imp)
    return IMEXTableau(a_exp, b_exp, c_exp, a_imp, b_imp, c_imp)
end

"""
    IMEXAlgorithm(tableau, newtons_method, [constraint])
    IMEXAlgorithm(name, newtons_method, [constraint])

Constructs an IMEX algorithm for solving ODEs, with an optional name and
constraint. The first constructor accepts any `IMEXTableau` and an optional
constraint, leaving the algorithm unnamed. The second constructor automatically
determines the tableau and the default constraint from the algorithm name,
which must be an `IMEXARKAlgorithmName`.
"""
struct IMEXAlgorithm{
    C <: AbstractAlgorithmConstraint,
    N <: Union{Nothing, AbstractAlgorithmName},
    T <: IMEXTableau,
    NM <: Union{Nothing, NewtonsMethod},
} <: DistributedODEAlgorithm
    constraint::C
    name::N
    tableau::T
    newtons_method::NM
end
IMEXAlgorithm(tableau::IMEXTableau, newtons_method, constraint = Unconstrained()) =
    IMEXAlgorithm(constraint, nothing, tableau, newtons_method)
IMEXAlgorithm(name::IMEXARKAlgorithmName, newtons_method, constraint = default_constraint(name)) =
    IMEXAlgorithm(constraint, name, IMEXTableau(name), newtons_method)

################################################################################

# ARS algorithms

# The naming convention is ARSsσp, where s is the number of implicit stages,
# σ is the number of explicit stages, and p is the order of accuracy.

"""
    ARS111

An IMEX ARK algorithm from [ARS1997](@cite), section 2, with 1 implicit stage,
1 explicit stage and 1st order accuracy. Also called *IMEX Euler* or
*forward-backward Euler*; equivalent to `OrdinaryDiffEq.IMEXEuler`.
"""
struct ARS111 <: IMEXARKAlgorithmName end
function IMEXTableau(::ARS111)
    IMEXTableau(; a_exp = @SArray([0 0; 1 0]), a_imp = @SArray([0 0; 0 1]))
end

"""
    ARS121

An IMEX ARK algorithm from [ARS1997](@cite), section 2, with 1 implicit stage, 2
explicit stages, and 1st order accuracy. Also called *IMEX Euler* or
*forward-backward Euler*; equivalent to `OrdinaryDiffEq.IMEXEulerARK`.
"""
struct ARS121 <: IMEXARKAlgorithmName end
function IMEXTableau(::ARS121)
    IMEXTableau(; a_exp = @SArray([0 0; 1 0]), b_exp = @SArray([0, 1]), a_imp = @SArray([0 0; 0 1]))
end

"""
    ARS122

An IMEX ARK algorithm from [ARS1997](@cite), section 2, with 1 implicit stage, 2
explicit stages, and 2nd order accuracy. Also called *IMEX midpoint*.
"""
struct ARS122 <: IMEXARKAlgorithmName end
function IMEXTableau(::ARS122)
    IMEXTableau(;
        a_exp = @SArray([0 0; 1/2 0]),
        b_exp = @SArray([0, 1]),
        a_imp = @SArray([0 0; 0 1/2]),
        b_imp = @SArray([0, 1])
    )
end

"""
    ARS233

An IMEX ARK algorithm from [ARS1997](@cite), section 2, with 2 implicit stages,
3 explicit stages, and 3rd order accuracy.
"""
struct ARS233 <: IMEXARKAlgorithmName end
function IMEXTableau(::ARS233)
    γ = 1 / 2 + √3 / 6
    IMEXTableau(;
        a_exp = @SArray([
            0 0 0
            γ 0 0
            (γ-1) (2-2γ) 0
        ]),
        b_exp = @SArray([0, 1 / 2, 1 / 2]),
        a_imp = @SArray([
            0 0 0
            0 γ 0
            0 (1-2γ) γ
        ]),
        b_imp = @SArray([0, 1 / 2, 1 / 2])
    )
end

"""
    ARS232

An IMEX ARK algorithm from [ARS1997](@cite), section 2, with 2 implicit stages,
3 explicit stages, and 2nd order accuracy.
"""
struct ARS232 <: IMEXARKAlgorithmName end
function IMEXTableau(::ARS232)
    γ = 1 - √2 / 2
    δ = -2√2 / 3
    IMEXTableau(;
        a_exp = @SArray([
            0 0 0
            γ 0 0
            δ (1-δ) 0
        ]),
        b_exp = @SArray([0, 1 - γ, γ]),
        a_imp = @SArray([
            0 0 0
            0 γ 0
            0 (1-γ) γ
        ])
    )
end

"""
    ARS222

An IMEX ARK algorithm from [ARS1997](@cite), section 2, with 2 implicit stages,
2 explicit stages, and 2nd order accuracy.
"""
struct ARS222 <: IMEXARKAlgorithmName end
function IMEXTableau(::ARS222)
    γ = 1 - √2 / 2
    δ = 1 - 1 / 2γ
    IMEXTableau(; a_exp = @SArray([
        0 0 0
        γ 0 0
        δ (1-δ) 0
    ]), a_imp = @SArray([
        0 0 0
        0 γ 0
        0 (1-γ) γ
    ]))
end

"""
    ARS343

An IMEX ARK algorithm from [ARS1997](@cite), section 2, with 3 implicit stages,
4 explicit stages, and 3rd order accuracy.
"""
struct ARS343 <: IMEXARKAlgorithmName end
function IMEXTableau(::ARS343)
    γ = 0.4358665215084590
    a42 = 0.5529291480359398
    a43 = a42
    b1 = -3 / 2 * γ^2 + 4 * γ - 1 / 4
    b2 = 3 / 2 * γ^2 - 5 * γ + 5 / 4
    a31 =
        (1 - 9 / 2 * γ + 3 / 2 * γ^2) * a42 + (11 / 4 - 21 / 2 * γ + 15 / 4 * γ^2) * a43 - 7 / 2 + 13 * γ - 9 / 2 * γ^2
    a32 =
        (-1 + 9 / 2 * γ - 3 / 2 * γ^2) * a42 + (-11 / 4 + 21 / 2 * γ - 15 / 4 * γ^2) * a43 + 4 - 25 / 2 * γ +
        9 / 2 * γ^2
    a41 = 1 - a42 - a43
    IMEXTableau(;
        a_exp = @SArray([
            0 0 0 0
            γ 0 0 0
            a31 a32 0 0
            a41 a42 a43 0
        ]),
        b_exp = @SArray([0, b1, b2, γ]),
        a_imp = @SArray([
            0 0 0 0
            0 γ 0 0
            0 (1 - γ)/2 γ 0
            0 b1 b2 γ
        ])
    )
end

"""
    ARS443

An IMEX ARK algorithm from [ARS1997](@cite), section 2, with 4 implicit stages,
4 explicit stages, and 3rd order accuracy.
"""
struct ARS443 <: IMEXARKAlgorithmName end
function IMEXTableau(::ARS443)
    IMEXTableau(;
        a_exp = @SArray([
            0 0 0 0 0
            1/2 0 0 0 0
            11/18 1/18 0 0 0
            5/6 -5/6 1/2 0 0
            1/4 7/4 3/4 -7/4 0
        ]),
        a_imp = @SArray([
            0 0 0 0 0
            0 1/2 0 0 0
            0 1/6 1/2 0 0
            0 -1/2 1/2 1/2 0
            0 3/2 -3/2 1/2 1/2
        ])
    )
end

################################################################################

# IMKG algorithms

# The naming convention is IMKGpfjl, where p is the order of accuracy, f is the
# number of explicit stages, j is the number of implicit stages, and l is an
# identifying letter.

# TODO: Tables 3 and 4 are riddled with typos, but most of these can be easily
# identified and corrected by referencing the implementations in HOMME:
# https://github.com/E3SM-Project/E3SM/blob/v2.0.0/components/homme/src/arkode/arkode_tables.F90
# Unfortunately, the implementations of IMKG353a and IMKG354a in HOMME also
# appear to be wrong, so they are left unimplemented. Eventually, we should get
# the official implementations from the paper's authors.

# a_exp:                                   a_imp:
# 0      0      ⋯      0      0      0     0      0      ⋯      0      0      0
# α[1]   0      ⋱      ⋮      ⋮      ⋮     α̂[1]   δ̂[1]   ⋱      ⋮      ⋮       0
# β[1]   α[2]   ⋱      0      0      0    β[1]   α̂[2]   ⋱      0      0      0
# ⋮      0      ⋱      0      0      0     ⋮      0      ⋱      δ̂[s-3] 0      0
# ⋮      ⋮       ⋱     α[s-2] 0      0     ⋮      ⋮       ⋱      α̂[s-2] δ̂[s-2] 0
# β[s-2] 0      ⋯      0      α[s-1] 0     β[s-2] 0      ⋯       0      α̂[s-1] 0
imkg_exp(i, j, α, β) = i == j + 1 ? α[j] : (i > 2 && j == 1 ? β[i - 2] : 0)
imkg_imp(i, j, α̂, β, δ̂) =
    i == j + 1 ? α̂[j] : (i > 2 && j == 1 ? β[i - 2] : (1 < i <= length(α̂) && i == j ? δ̂[i - 1] : 0))
function IMKGTableau(α, α̂, δ̂, β = ntuple(_ -> 0, length(δ̂)))
    s = length(α̂) + 1
    type = SMatrix{s, s}
    return IMEXTableau(;
        a_exp = StaticArrays.sacollect(type, imkg_exp(i, j, α, β) for i in 1:s, j in 1:s),
        a_imp = StaticArrays.sacollect(type, imkg_imp(i, j, α̂, β, δ̂) for i in 1:s, j in 1:s),
    )
end

"""
    IMKG232a

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 2 implicit stages,
3 explicit stages, and 2nd order accuracy.
"""
struct IMKG232a <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG232a)
    IMKGTableau((1 / 2, 1 / 2, 1), (0, -1 / 2 + √2 / 2, 1), (1 - √2 / 2, 1 - √2 / 2))
end

"""
    IMKG232b

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 2 implicit stages,
3 explicit stages, and 2nd order accuracy.
"""
struct IMKG232b <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG232b)
    IMKGTableau((1 / 2, 1 / 2, 1), (0, -1 / 2 - √2 / 2, 1), (1 + √2 / 2, 1 + √2 / 2))
end

"""
    IMKG242a

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 2 implicit stages,
4 explicit stages, and 2nd order accuracy.
"""
struct IMKG242a <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG242a)
    IMKGTableau((1 / 4, 1 / 3, 1 / 2, 1), (0, 0, -1 / 2 + √2 / 2, 1), (0, 1 - √2 / 2, 1 - √2 / 2))
end

"""
    IMKG242b

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 2 implicit stages,
4 explicit stages, and 2nd order accuracy.
"""
struct IMKG242b <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG242b)
    IMKGTableau((1 / 4, 1 / 3, 1 / 2, 1), (0, 0, -1 / 2 - √2 / 2, 1), (0, 1 + √2 / 2, 1 + √2 / 2))
end

"""
    IMKG243a

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 3 implicit stages,
4 explicit stages, and 2nd order accuracy.
"""
struct IMKG243a <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG243a)
    IMKGTableau((1 / 4, 1 / 3, 1 / 2, 1), (0, 1 / 6, -√3 / 6, 1), (1 / 2 + √3 / 6, 1 / 2 + √3 / 6, 1 / 2 + √3 / 6))
end
# The paper uses √3/6 for α̂[3], which also seems to work.

"""
    IMKG252a

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 2 implicit stages,
5 explicit stages, and 2nd order accuracy.
"""
struct IMKG252a <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG252a)
    IMKGTableau((1 / 4, 1 / 6, 3 / 8, 1 / 2, 1), (0, 0, 0, -1 / 2 + √2 / 2, 1), (0, 0, 1 - √2 / 2, 1 - √2 / 2))
end

"""
    IMKG252b

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 2 implicit stages,
5 explicit stages, and 2nd order accuracy.
"""
struct IMKG252b <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG252b)
    IMKGTableau((1 / 4, 1 / 6, 3 / 8, 1 / 2, 1), (0, 0, 0, -1 / 2 - √2 / 2, 1), (0, 0, 1 + √2 / 2, 1 + √2 / 2))
end

"""
    IMKG253a

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 3 implicit stages,
5 explicit stages, and 2nd order accuracy.
"""
struct IMKG253a <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG253a)
    IMKGTableau(
        (1 / 4, 1 / 6, 3 / 8, 1 / 2, 1),
        (0, 0, √3 / 4 * (1 - √3 / 3) * ((1 + √3 / 3)^2 - 2), √3 / 6, 1),
        (0, 1 / 2 - √3 / 6, 1 / 2 - √3 / 6, 1 / 2 - √3 / 6),
    )
end
# The paper uses 0.08931639747704086 for α̂[3], which also seems to work.

"""
    IMKG253b

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 3 implicit stages,
5 explicit stages, and 2nd order accuracy.
"""
struct IMKG253b <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG253b)
    IMKGTableau(
        (1 / 4, 1 / 6, 3 / 8, 1 / 2, 1),
        (0, 0, √3 / 4 * (1 + √3 / 3) * ((1 - √3 / 3)^2 - 2), -√3 / 6, 1),
        (0, 1 / 2 + √3 / 6, 1 / 2 + √3 / 6, 1 / 2 + √3 / 6),
    )
end
# The paper uses 1.2440169358562922 for α̂[3], which also seems to work.

"""
    IMKG254a

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 4 implicit stages,
5 explicit stages, and 2nd order accuracy.
"""
struct IMKG254a <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG254a)
    IMKGTableau((1 / 4, 1 / 6, 3 / 8, 1 / 2, 1), (0, -3 / 10, 5 / 6, -3 / 2, 1), (-1 / 2, 1, 1, 2))
end

"""
    IMKG254b

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 4 implicit stages,
5 explicit stages, and 2nd order accuracy.
"""
struct IMKG254b <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG254b)
    IMKGTableau((1 / 4, 1 / 6, 3 / 8, 1 / 2, 1), (0, -1 / 20, 5 / 4, -1 / 2, 1), (-1 / 2, 1, 1, 1))
end

"""
    IMKG254c

An IMEX ARK algorithm from [SVTG2019](@cite), Table 3, with 4 implicit stages,
5 explicit stages, and 2nd order accuracy.
"""
struct IMKG254c <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG254c)
    IMKGTableau((1 / 4, 1 / 6, 3 / 8, 1 / 2, 1), (0, 1 / 20, 5 / 36, 1 / 3, 1), (1 / 6, 1 / 6, 1 / 6, 1 / 6))
end

"""
    IMKG342a

An IMEX ARK algorithm from [SVTG2019](@cite), Table 4, with 2 implicit stages,
4 explicit stages, and 3rd order accuracy.
"""
struct IMKG342a <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG342a)
    IMKGTableau(
        (1 / 4, 2 / 3, 1 / 3, 3 / 4),
        (0, 1 / 6 - √3 / 6, -1 / 6 - √3 / 6, 3 / 4),
        (0, 1 / 2 + √3 / 6, 1 / 2 + √3 / 6),
        (0, 1 / 3, 1 / 4),
    )
end
# The paper and HOMME completely disagree on IMKG342a. Since the version in the
# paper is not "342" (it appears to be "332"), the version from HOMME is used
# here. The paper's version is
# IMKGTableau(
#     (0, 1/3, 1/3, 3/4),
#     (0, -1/6 - √3/6, -1/6 - √3/6, 3/4),
#     (0, 1/2 + √3/6, 1/2 + √3/6),
#     (1/3, 1/3, 1/4),
# )

"""
    IMKG343a

An IMEX ARK algorithm from [SVTG2019](@cite), Table 4, with 3 implicit stages,
4 explicit stages, and 3rd order accuracy.
"""
struct IMKG343a <: IMEXARKAlgorithmName end
function IMEXTableau(::IMKG343a)
    IMKGTableau((1 / 4, 2 / 3, 1 / 3, 3 / 4), (0, -1 / 3, -2 / 3, 3 / 4), (-1 / 3, 1, 1), (0, 1 / 3, 1 / 4))
end

# The paper and HOMME completely disagree on IMKG353a, but neither version
# is "353" (they appear to be "343" and "354", respectively). The paper's
# version is
# IMKGTableau(
#     (1/4, 2/3, 1/3, 3/4),
#     (0, -359/600, -559/600, 3/4),
#     (-1.1678009811335388, 253/200, 253/200),
#     (0, 1/3, 1/4),
# )
# HOMME's version is
# IMKGTableau(
#     (-0.017391304347826087, -23/25, 5/3, 1/3, 3/4),
#     (0.3075640504095504, -1.2990164859879263, 751/600, -49/60, 3/4),
#     (-0.2981612530370581, 83/200, 83/200, 23/20),
#     (1, -1, 1/3, 1/4),
# )

# The version of IMKG354a in the paper is not "354" (it appears to be "253"),
# and IMKG354a is missing from HOMME (or, more precisely, the tableau for
# IMKG353a is mistakenly used to define IMKG354a, and the tableau for IMKG354a
# is not specified). The paper's version is
# IMKGTableau(
#     (1/5, 1/5, 2/3, 1/3, 3/4),
#     (0, 0, 11/30, -2/3, 3/4),
#     (0, 2/4, 2/5, 1),
#     (0, 0, 1/3, 1/4),
# )

################################################################################

# IMEX SSPRK algorithms

# The naming convention is SSPsσp, where s is the number of implicit stages,
# σ is the number of explicit stages, and p is the order of accuracy.

abstract type IMEXSSPRKAlgorithmName <: IMEXARKAlgorithmName end

default_constraint(::IMEXSSPRKAlgorithmName) = SSP()

"""
    SSP222

An IMEX SSPRK algorithm from [PR2005](@cite), with 2 implicit stages, 2 explicit
stages, and 2nd order accuracy. Also called *SSP2(222)* in [GGHRUW2018](@cite).
"""
struct SSP222 <: IMEXSSPRKAlgorithmName end
function IMEXTableau(::SSP222)
    γ = 1 - √2 / 2
    return IMEXTableau(;
        a_exp = @SArray([
            0 0
            1 0
        ]),
        b_exp = @SArray([1 / 2, 1 / 2]),
        a_imp = @SArray([
            γ 0
            (1-2γ) γ
        ]),
        b_imp = @SArray([1 / 2, 1 / 2])
    )
end

"""
    SSP322

An IMEX SSPRK algorithm from [PR2005](@cite), with 3 implicit stages, 2 explicit
stages, and 2nd order accuracy.
"""
struct SSP322 <: IMEXSSPRKAlgorithmName end
function IMEXTableau(::SSP322)
    return IMEXTableau(;
        a_exp = @SArray([
            0 0 0
            0 0 0
            0 1 0
        ]),
        b_exp = @SArray([0, 1 / 2, 1 / 2]),
        a_imp = @SArray([
            1/2 0 0
            -1/2 1/2 0
            0 1/2 1/2
        ]),
        b_imp = @SArray([0, 1 / 2, 1 / 2])
    )
end

"""
    SSP332

An IMEX SSPRK algorithm from [PR2005](@cite), with 3 implicit stages, 3 explicit
stages, and 2nd order accuracy. Also called *SSP2(332)a* in [GGHRUW2018](@cite).
"""
struct SSP332 <: IMEXSSPRKAlgorithmName end
function IMEXTableau(::SSP332)
    γ = 1 - √2 / 2
    return IMEXTableau(;
        a_exp = @SArray([
            0 0 0
            1 0 0
            1/4 1/4 0
        ]),
        b_exp = @SArray([1 / 6, 1 / 6, 2 / 3]),
        a_imp = @SArray([
            γ 0 0
            (1-2γ) γ 0
            (1 / 2-γ) 0 γ
        ]),
        b_imp = @SArray([1 / 6, 1 / 6, 2 / 3])
    )
end

"""
    SSP333(; β = 1/2 + √3/6)

Family of IMEX SSPRK algorithms parametrized by the value β from
[CGGS2017](@cite), Section 3.2, with 3 implicit stages, 3 explicit stages, and
3rd order accuracy. The default value of β results in an SDIRK algorithm, which
is also called *SSP3(333)c* in [GGHRUW2018](@cite).
"""
Base.@kwdef struct SSP333{FT <: AbstractFloat} <: IMEXSSPRKAlgorithmName
    β::FT = 1 / 2 + √3 / 6
end
function IMEXTableau((; β)::SSP333)
    @assert β > 1 / 2
    γ = (2β^2 - 3β / 2 + 1 / 3) / (2 - 4β)
    return IMEXTableau(;
        a_exp = @SArray([
            0 0 0
            1 0 0
            1/4 1/4 0
        ]),
        b_exp = @SArray([1 / 6, 1 / 6, 2 / 3]),
        a_imp = @SArray([
            0 0 0
            (4γ+2β) (1 - 4γ-2β) 0
            (1 / 2 - β-γ) γ β
        ]),
        b_imp = @SArray([1 / 6, 1 / 6, 2 / 3])
    )
end

"""
    SSP433

An IMEX SSPRK algorithm from [PR2005](@cite), with 4 implicit stages, 3 explicit
stages, and 3rd order accuracy. Also called *SSP3(433)* in [GGHRUW2018](@cite).
"""
struct SSP433 <: IMEXSSPRKAlgorithmName end
function IMEXTableau(::SSP433)
    α = 0.24169426078821
    β = 0.06042356519705
    η = 0.12915286960590
    return IMEXTableau(;
        a_exp = @SArray([
            0 0 0 0
            0 0 0 0
            0 1 0 0
            0 1/4 1/4 0
        ]),
        b_exp = @SArray([0, 1 / 6, 1 / 6, 2 / 3]),
        a_imp = @SArray([
            α 0 0 0
            -α α 0 0
            0 (1-α) α 0
            β η (1 / 2 - α - β-η) α
        ]),
        b_imp = @SArray([0, 1 / 6, 1 / 6, 2 / 3])
    )
end

################################################################################

# Miscellaneous algorithms

"""
    DBM453

An IMEX ARK algorithm from [VSRUW2019](@cite), Appendix A, with 4 implicit
stages, 5 explicit stages, and 3rd order accuracy.
"""
struct DBM453 <: IMEXARKAlgorithmName end
function IMEXTableau(::DBM453)
    γ = 0.32591194130117247
    IMEXTableau(;
        a_exp = @SArray(
            [
                0 0 0 0 0
                0.10306208811591838 0 0 0 0
                -0.94124866143519894 1.6626399742527356 0 0 0
                -1.3670975201437765 1.3815852911016873 1.2673234025619065 0 0
                -0.81287582068772448 0.81223739060505738 0.90644429603699305 0.094194134045674111 0
            ]
        ),
        b_exp = @SArray([0.87795339639076675, -0.72692641526151547, 0.7520413715737272, -0.22898029400415088, γ]),
        a_imp = @SArray(
            [
                0 0 0 0 0
                -0.2228498531852541 γ 0 0 0
                -0.46801347074080545 0.86349284225716961 γ 0 0
                -0.46509906651927421 0.81063103116959553 0.61036726756832357 γ 0
                0.87795339639076675 -0.72692641526151547 0.7520413715737272 -0.22898029400415088 γ
            ]
        )
    )
end

"""
    HOMMEM1

An IMEX ARK algorithm from [GTBBS2020](@cite), section 4.1, with 5 implicit
stages, 6 explicit stages, and 2nd order accuracy.
"""
struct HOMMEM1 <: IMEXARKAlgorithmName end
function IMEXTableau(::HOMMEM1)
    IMEXTableau(;
        a_exp = @SArray([
            0 0 0 0 0 0
            1/5 0 0 0 0 0
            0 1/5 0 0 0 0
            0 0 1/3 0 0 0
            0 0 0 1/2 0 0
            0 0 0 0 1 0
        ]),
        a_imp = @SArray([
            0 0 0 0 0 0
            0 1/5 0 0 0 0
            0 0 1/5 0 0 0
            0 0 0 1/3 0 0
            0 0 0 0 1/2 0
            5/18 5/18 0 0 0 8/18
        ])
    )
end

"""
    ARK2GKC(; paper_version = false)

An IMEX ARK algorithm from [GKC2013](@cite) with 2 implicit stages, 3 explicit
stages, and 2nd order accuracy. If `paper_version = true`, the algorithm uses
coefficients from the paper. Otherwise, it uses coefficients that make it more
stable but less accurate.
"""
Base.@kwdef struct ARK2GKC <: IMEXARKAlgorithmName
    paper_version::Bool = false
end
function IMEXTableau((; paper_version)::ARK2GKC)
    a32 = paper_version ? 1 / 2 + √2 / 3 : 1 / 2
    IMEXTableau(;
        a_exp = @SArray([
            0 0 0
            (2-√2) 0 0
            (1-a32) a32 0
        ]),
        b_exp = @SArray([√2 / 4, √2 / 4, 1 - √2 / 2]),
        a_imp = @SArray([
            0 0 0
            (1-√2 / 2) (1-√2 / 2) 0
            √2/4 √2/4 (1-√2 / 2)
        ])
    )
end

"""
    ARK437L2SA1

An IMEX ARK algorithm from [KC2019](@cite), Table 8, with 6 implicit stages, 7
explicit stages, and 4th order accuracy. Written as *ARK4(3)7L[2]SA₁* in the
paper.
"""
struct ARK437L2SA1 <: IMEXARKAlgorithmName end
function IMEXTableau(::ARK437L2SA1)
    a_exp = zeros(Rational{Int64}, 7, 7)
    a_imp = zeros(Rational{Int64}, 7, 7)
    b = zeros(Rational{Int64}, 7)
    c = zeros(Rational{Int64}, 7)

    γ = 1235 // 10000
    for i in 2:7
        a_imp[i, i] = γ
    end

    a_imp[3, 2] = 624185399699 // 4186980696204
    a_imp[4, 2] = 1258591069120 // 10082082980243
    a_imp[4, 3] = -322722984531 // 8455138723562
    a_imp[5, 2] = -436103496990 // 5971407786587
    a_imp[5, 3] = -2689175662187 // 11046760208243
    a_imp[5, 4] = 4431412449334 // 12995360898505
    a_imp[6, 2] = -2207373168298 // 14430576638973
    a_imp[6, 3] = 242511121179 // 3358618340039
    a_imp[6, 4] = 3145666661981 // 7780404714551
    a_imp[6, 5] = 5882073923981 // 14490790706663
    a_imp[7, 2] = 0
    a_imp[7, 3] = 9164257142617 // 17756377923965
    a_imp[7, 4] = -10812980402763 // 74029279521829
    a_imp[7, 5] = 1335994250573 // 5691609445217
    a_imp[7, 6] = 2273837961795 // 8368240463276

    a_exp[3, 1] = 247 // 4000
    a_exp[3, 2] = 2694949928731 // 7487940209513
    a_exp[4, 1] = 464650059369 // 8764239774964
    a_exp[4, 2] = 878889893998 // 2444806327765
    a_exp[4, 3] = -952945855348 // 12294611323341
    a_exp[5, 1] = 476636172619 // 8159180917465
    a_exp[5, 2] = -1271469283451 // 7793814740893
    a_exp[5, 3] = -859560642026 // 4356155882851
    a_exp[5, 4] = 1723805262919 // 4571918432560
    a_exp[6, 1] = 6338158500785 // 11769362343261
    a_exp[6, 2] = -4970555480458 // 10924838743837
    a_exp[6, 3] = 3326578051521 // 2647936831840
    a_exp[6, 4] = -880713585975 // 1841400956686
    a_exp[6, 5] = -1428733748635 // 8843423958496
    a_exp[7, 2] = 760814592956 // 3276306540349
    a_exp[7, 3] = -47223648122716 // 6934462133451
    a_exp[7, 4] = 71187472546993 // 9669769126921
    a_exp[7, 5] = -13330509492149 // 9695768672337
    a_exp[7, 6] = 11565764226357 // 8513123442827

    b[2] = 0
    b[3] = 9164257142617 // 17756377923965
    b[4] = -10812980402763 // 74029279521829
    b[5] = 1335994250573 // 5691609445217
    b[6] = 2273837961795 // 8368240463276
    b[7] = 247 // 2000

    c[2] = 247 // 1000
    c[3] = 4276536705230 // 10142255878289
    c[4] = 67 // 200
    c[5] = 3 // 40
    c[6] = 7 // 10

    for i in 2:7
        a_imp[i, 1] = a_imp[i, 2]
    end
    b[1] = b[2]
    a_exp[2, 1] = c[2]
    a_exp[7, 1] = a_exp[7, 2]
    c[1] = 0
    c[7] = 1

    IMEXTableau(;
        a_exp = SArray{Tuple{7, 7}}(a_exp),
        b_exp = SArray{Tuple{7}}(b),
        c_exp = SArray{Tuple{7}}(c),
        a_imp = SArray{Tuple{7, 7}}(a_imp),
        b_imp = SArray{Tuple{7}}(b),
        c_imp = SArray{Tuple{7}}(c),
    )
end

"""
    ARK548L2SA2

An IMEX ARK algorithm from [KC2019](@cite), Table 8, with 7 implicit stages, 8
explicit stages, and 5th order accuracy. Written as *ARK5(4)8L[2]SA₂* in the
paper.
"""
struct ARK548L2SA2 <: IMEXARKAlgorithmName end
function IMEXTableau(::ARK548L2SA2)
    a_exp = zeros(Rational{Int64}, 8, 8)
    a_imp = zeros(Rational{Int64}, 8, 8)
    b = zeros(Rational{Int64}, 8)
    c = zeros(Rational{Int64}, 8)

    γ = 2 // 9
    for i in 2:8
        a_imp[i, i] = γ
    end

    a_imp[3, 2] = 2366667076620 // 8822750406821
    a_imp[4, 2] = -257962897183 // 4451812247028
    a_imp[4, 3] = 128530224461 // 14379561246022
    a_imp[5, 2] = -486229321650 // 11227943450093
    a_imp[5, 3] = -225633144460 // 6633558740617
    a_imp[5, 4] = 1741320951451 // 6824444397158
    a_imp[6, 2] = 621307788657 // 4714163060173
    a_imp[6, 3] = -125196015625 // 3866852212004
    a_imp[6, 4] = 940440206406 // 7593089888465
    a_imp[6, 5] = 961109811699 // 6734810228204
    a_imp[7, 2] = 2036305566805 // 6583108094622
    a_imp[7, 3] = -3039402635899 // 4450598839912
    a_imp[7, 4] = -1829510709469 // 31102090912115
    a_imp[7, 5] = -286320471013 // 6931253422520
    a_imp[7, 6] = 8651533662697 // 9642993110008

    a_exp[3, 1] = 1 // 9
    a_exp[3, 2] = 1183333538310 // 1827251437969
    a_exp[4, 1] = 895379019517 // 9750411845327
    a_exp[4, 2] = 477606656805 // 13473228687314
    a_exp[4, 3] = -112564739183 // 9373365219272
    a_exp[5, 1] = -4458043123994 // 13015289567637
    a_exp[5, 2] = -2500665203865 // 9342069639922
    a_exp[5, 3] = 983347055801 // 8893519644487
    a_exp[5, 4] = 2185051477207 // 2551468980502
    a_exp[6, 1] = -167316361917 // 17121522574472
    a_exp[6, 2] = 1605541814917 // 7619724128744
    a_exp[6, 3] = 991021770328 // 13052792161721
    a_exp[6, 4] = 2342280609577 // 11279663441611
    a_exp[6, 5] = 3012424348531 // 12792462456678
    a_exp[7, 1] = 6680998715867 // 14310383562358
    a_exp[7, 2] = 5029118570809 // 3897454228471
    a_exp[7, 3] = 2415062538259 // 6382199904604
    a_exp[7, 4] = -3924368632305 // 6964820224454
    a_exp[7, 5] = -4331110370267 // 15021686902756
    a_exp[7, 6] = -3944303808049 // 11994238218192
    a_exp[8, 1] = 2193717860234 // 3570523412979
    a_exp[8, 2] = a_exp[8, 1]
    a_exp[8, 3] = 5952760925747 // 18750164281544
    a_exp[8, 4] = -4412967128996 // 6196664114337
    a_exp[8, 5] = 4151782504231 // 36106512998704
    a_exp[8, 6] = 572599549169 // 6265429158920
    a_exp[8, 7] = -457874356192 // 11306498036315

    b[2] = 0
    b[3] = 3517720773327 // 20256071687669
    b[4] = 4569610470461 // 17934693873752
    b[5] = 2819471173109 // 11655438449929
    b[6] = 3296210113763 // 10722700128969
    b[7] = -1142099968913 // 5710983926999

    c[2] = 4 // 9
    c[3] = 6456083330201 // 8509243623797
    c[4] = 1632083962415 // 14158861528103
    c[5] = 6365430648612 // 17842476412687
    c[6] = 18 // 25
    c[7] = 191 // 200

    for i in 2:8
        a_imp[i, 1] = a_imp[i, 2]
    end
    b[1] = b[2]
    b[8] = γ
    for i in 1:8
        a_imp[8, i] = b[i]
    end
    a_exp[2, 1] = c[2]
    a_exp[8, 1] = a_exp[8, 2]
    c[1] = 0
    c[8] = 1

    IMEXTableau(;
        a_exp = SArray{Tuple{8, 8}}(a_exp),
        b_exp = SArray{Tuple{8}}(b),
        c_exp = SArray{Tuple{8}}(c),
        a_imp = SArray{Tuple{8, 8}}(a_imp),
        b_imp = SArray{Tuple{8}}(b),
        c_imp = SArray{Tuple{8}}(c),
    )
end
