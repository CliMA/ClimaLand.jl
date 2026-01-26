using DiffEqBase, ClimaTimeSteppers, StaticArrays
import LinearAlgebra
import LinearAlgebra: norm, Diagonal, mul!
using ClimaCore
using ClimaComms
import ClimaCore.MatrixFields: @name, ⋅, FieldMatrixWithSolver
import ClimaCore: Domains, Geometry, Meshes, Topologies, Spaces, Fields, Operators, Limiters

@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

import Krylov # Trigger ClimaCore/ext/KrylovExt

"""
Single variable linear ODE

```math
\\frac{du}{dt} = αu
```
with initial condition ``u_0=\\frac{1}{2}``, parameter ``α=1.01``, and solution
```math
u(t) = u_0 e^{αt}
```

This is an in-place variant of the one from DiffEqProblemLibrary.jl.
"""
function linear_prob()
    ODEProblem(
        IncrementingODEFunction{true}((du, u, p, t, α = true, β = false) -> (du .= α .* p .* u .+ β .* du)),
        [1 / 2],
        (0.0, 1.0),
        1.01,
    )
end
function linear_prob_fe()
    ODEProblem(ForwardEulerODEFunction((un, u, p, t, dt) -> (un .= u .+ dt .* p .* u)), [1.0], (0.0, 1.0), -0.2)
end

function linear_prob_wfactt()
    ODEProblem(
        ODEFunction(
            (du, u, p, t) -> (du .= p .* u);
            jac_prototype = zeros(1, 1),
            Wfact_t = (W, u, p, γ, t) -> (W[1, 1] = 1 / γ - p),
        ),
        [1 / 2],
        (0.0, 1.0),
        -0.2,
    )
end

function split_linear_prob_wfact_split()
    ODEProblem(
        ClimaODEFunction(;
            T_imp! = ODEFunction(
                (du, u, p, t) -> (du .= real(p) .* u);
                jac_prototype = zeros(ComplexF64, 1, 1),
                Wfact = (W, u, p, γ, t) -> (W[1, 1] = γ * real(p) - 1),
            ),
            T_exp! = ODEFunction((du, u, p, t) -> (du .= imag(p) * im .* u)),
        ),
        [1 / 2 + 0.0 * im],
        (0.0, 1.0),
        -0.2 + 0.1 * im,
    )
end

function split_linear_prob_wfact_split_fe()
    ODEProblem(
        ClimaODEFunction(;
            T_imp! = ODEFunction(
                (du, u, p, t) -> (du .= real(p) .* u);
                jac_prototype = zeros(ComplexF64, 1, 1),
                Wfact = (W, u, p, γ, t) -> (W[1, 1] = γ * real(p) - 1),
            ),
            T_exp! = (du, u, p, t) -> (du .= imag(p) * im .* u),
        ),
        [1 / 2 + 0.0 * im],
        (0.0, 1.0),
        -0.2 + 0.1 * im,
    )
end

# DiffEqProblemLibrary.jl uses the `analytic` argument to ODEFunction to store the exact solution
# IncrementingODEFunction doesn't have that
function linear_sol(u0, p, t)
    u0 .* exp(p * t)
end


"""
Two variable linear ODE

```math
\\frac{du_1}{dt} = α u_2 \\quad \\frac{du_2}{dt} = -α u_1
```
with initial condition ``u_0=[0,1]``, parameter ``α=2``, and solution
```math
u(t) = [cos(αt) sin(αt); -sin(αt) cos(αt) ] u_0
```
"""
function sincos_prob()
    ODEProblem(
        IncrementingODEFunction{true}((du, u, p, t, α = true, β = false) -> (du[1] = α * p * u[2] + β * du[1];
        du[2] = -α * p * u[1] + β * du[2])),
        [0.0, 1.0],
        (0.0, 1.0),
        2.0,
    )
end
function sincos_prob_fe()
    ODEProblem(
        ForwardEulerODEFunction((un, u, p, t, dt) -> (un[1] = u[1] + dt * p * u[2]; un[2] = u[2] - dt * p * u[1])),
        [0.0, 1.0],
        (0.0, 1.0),
        2.0,
    )
end

function sincos_sol(u0, p, t)
    s, c = sincos(p * t)
    [c s; -s c] * u0
end

"""
IMEX problem with autonomous linear part

```math
\\frac{du}{dt} = u + cos(t)/α
```
with initial condition ``u_0 = 1/2``, parameter ``α=4``, and solution
```math
u(t) = \frac{e^t + sin(t) - cos(t)}{2 α} + u_0 e^t
```
"""
function imex_autonomous_prob(::Type{ArrayType}) where {ArrayType}
    SplitODEProblem(
        IncrementingODEFunction{true}((du, u, p, t, α = true, β = false) -> (du .= α .* u .+ β .* du)),
        IncrementingODEFunction{true}((du, u, p, t, α = true, β = false) -> (du .= α .* cos(t) / p .+ β .* du)),
        ArrayType([0.5]),
        (0.0, 1.0),
        4.0,
    )
end

function linsolve_direct(::Type{Val{:init}}, f, u0; kwargs...)
    function _linsolve!(x, A, b, update_matrix = false; kwargs...)
        x .= A \ b
    end
end

function imex_autonomous_prob_jac(::Type{ArrayType}) where {ArrayType}
    ODEProblem(
        ODEFunction(
            IncrementingODEFunction{true}(
                (du, u, p, t, α = true, β = false) -> (du .= α .* (u .+ cos(t) / p) .+ β .* du),
            ),
            jac_prototype = zeros(1, 1),
            Wfact = (W, u, p, γ, t) -> W[1, 1] = 1 - γ,
        ),
        ArrayType([0.5]),
        (0.0, 1.0),
        4.0,
    )
end

function imex_autonomous_sol(u0, p, t)
    (exp(t) + sin(t) - cos(t)) / 2p .+ exp(t) .* u0
end

function imex_nonautonomous_prob(::Type{ArrayType}) where {ArrayType}
    SplitODEProblem(
        IncrementingODEFunction{true}((du, u, p, t, α = true, β = false) -> (du .= α .* cos(t) * u .+ β .* du)),
        IncrementingODEFunction{true}((du, u, p, t, α = true, β = false) -> (du .= α .* cos(t) / p .+ β .* du)),
        ArrayType([0.5]),
        (0.0, 2.0),
        4.0,
    )
end

function imex_nonautonomous_sol(u0, p, t)
    (exp(sin(t)) .* (1 .+ p .* u0) .- 1) ./ p
end

function kpr_rhs(Q, param, t)
    yf, ys = Q
    ω, λf, λs, ξ, α = param

    ηfs = ((1 - ξ) / α) * (λf - λs)
    ηsf = -ξ * α * (λf - λs)
    Ω = @SMatrix [
        λf ηfs
        ηsf λs
    ]

    g = @SVector [(-3 + yf^2 - cos(ω * t)) / 2yf, (-2 + ys^2 - cos(t)) / 2ys]
    h = @SVector [ω * sin(ω * t) / 2yf, sin(t) / 2ys]
    return Ω * g - h
end

function kpr_fast!(dQ, Q, param, t, α = 1, β = 0)
    dQf, _ = kpr_rhs(Q, param, t)
    dQ[1] = α * dQf + β * dQ[1]
    dQ[2] = β * dQ[2]
end

function kpr_slow!(dQ, Q, param, t, α = 1, β = 0)
    _, dQs = kpr_rhs(Q, param, t)
    dQ[1] = β * dQ[1]
    dQ[2] = α * dQs + β * dQ[2]
end

function kpr_sol(u0, param, t)
    #@assert u0 == [sqrt(4) sqrt(3)]
    ω, λf, λs, ξ, α = param
    return [
        sqrt(3 + cos(ω * t))
        sqrt(2 + cos(t))
    ]
end

"""
Test problem (4.2) from [RobertsSarsharSandu2018arxiv](@cite)

Note: The actual rates are all over the place with this test and passing largely
        depends on final dt size
"""
function kpr_multirate_prob()
    kpr_param = (ω = 100, λf = -10, λs = -1, ξ = 0.1, α = 1)
    SplitODEProblem(
        IncrementingODEFunction{true}(kpr_fast!),
        IncrementingODEFunction{true}(kpr_slow!),
        [sqrt(4), sqrt(3)],
        (0.0, 5π / 2),
        kpr_param,
    )
end

function kpr_singlerate_prob()
    kpr_param = (ω = 100, λf = -10, λs = -1, ξ = 0.1, α = 1)
    ODEProblem(
        IncrementingODEFunction{true}((dQ, Q, param, t, α = 1, β = 0) -> (dQ .= α .* kpr_rhs(Q, param, t) .+ β .* dQ)),
        [sqrt(4), sqrt(3)],
        (0.0, 5π / 2),
        kpr_param,
    )
end

reverse_problem(prob, analytic_sol) = ODEProblem(prob.f, analytic_sol(prob.tspan[2]), reverse(prob.tspan), prob.p)

struct IntegratorTestCase{FT, A, P, SP}
    test_name::String
    linear_implicit::Bool
    t_end::FT
    analytic_sol::A
    prob::P
    split_prob::SP
end

function IntegratorTestCase(;
    test_name,
    linear_implicit,
    t_end,
    Y₀,
    analytic_sol,
    tendency!,
    implicit_tendency! = nothing,
    explicit_tendency! = nothing,
    Wfact!,
    tgrad! = nothing,
)
    FT = typeof(t_end)
    jac_prototype = Matrix{FT}(undef, length(Y₀), length(Y₀))
    func_args = (; jac_prototype, Wfact = Wfact!, tgrad = tgrad!)
    tendency_func = ODEFunction(tendency!; func_args...)
    split_tendency_func = if isnothing(implicit_tendency!) # assume that related args are also nothing
        no_tendency!(Yₜ, Y, _, t) = Yₜ .= 0
        SplitFunction(tendency_func, no_tendency!)
    else
        SplitFunction(ODEFunction(implicit_tendency!; func_args...), explicit_tendency!)
    end
    make_prob(func) = ODEProblem(func, Y₀, (FT(0), t_end), nothing)
    IntegratorTestCase(
        test_name,
        linear_implicit,
        t_end,
        analytic_sol,
        make_prob(tendency_func),
        make_prob(split_tendency_func),
    )
end

function ClimaIntegratorTestCase(;
    test_name,
    linear_implicit,
    t_end,
    Y₀,
    analytic_sol,
    tendency!,
    implicit_tendency! = nothing,
    explicit_tendency! = nothing,
    Wfact!,
    tgrad! = nothing,
)
    FT = typeof(t_end)
    jac_prototype = Matrix{FT}(undef, length(Y₀), length(Y₀))
    func_args = (; jac_prototype, Wfact = Wfact!, tgrad = tgrad!)
    tendency_func = ClimaODEFunction(; T_imp! = ODEFunction(tendency!; func_args...))

    T_imp! = if isnothing(implicit_tendency!)
        # assume that related args are also nothing
        ODEFunction(tendency!; func_args...)
    else
        ODEFunction(implicit_tendency!; func_args...)
    end
    split_tendency_func = ClimaODEFunction(; T_exp! = explicit_tendency!, T_imp!)
    make_prob(func) = ODEProblem(func, Y₀, (FT(0), t_end), nothing)
    IntegratorTestCase(
        test_name,
        linear_implicit,
        t_end,
        analytic_sol,
        make_prob(tendency_func),
        make_prob(split_tendency_func),
    )
end

# A trivial test case for which any value of dt will work
function constant_tendency_test(::Type{FT}) where {FT}
    tendency = FT[1, 2, 3]
    IntegratorTestCase(;
        test_name = "constant_tendency",
        linear_implicit = true,
        t_end = FT(1),
        Y₀ = FT[0, 0, 0],
        analytic_sol = (t) -> tendency .* t,
        tendency! = (Yₜ, Y, _, t) -> Yₜ .= tendency,
        Wfact! = (W, Y, _, Δt, t) -> W .= -1,
    )
end

function clima_constant_tendency_test(::Type{FT}) where {FT}
    tendency = FT[1, 2, 3]
    ClimaIntegratorTestCase(;
        test_name = "constant_tendency",
        linear_implicit = true,
        t_end = FT(1),
        Y₀ = FT[0, 0, 0],
        analytic_sol = (t) -> tendency .* t,
        tendency! = (Yₜ, Y, _, t) -> Yₜ .= tendency,
        Wfact! = (W, Y, _, Δt, t) -> W .= -1,
    )
end

# From Section 1.1 of "Example Programs for ARKode v4.4.0" by D. R. Reynolds
function ark_analytic_test_cts(::Type{FT}) where {FT}
    λ = FT(-100) # increase magnitude for more stiffness
    source(t) = 1 / (1 + t^2) - λ * atan(t)
    ClimaIntegratorTestCase(;
        test_name = "ark_analytic",
        linear_implicit = true,
        t_end = FT(10),
        Y₀ = FT[0],
        analytic_sol = (t) -> [atan(t)],
        tendency! = (Yₜ, Y, _, t) -> Yₜ .= λ .* Y .+ source(t),
        implicit_tendency! = (Yₜ, Y, _, t) -> Yₜ .= λ .* Y,
        explicit_tendency! = (Yₜ, Y, _, t) -> Yₜ .= source(t),
        Wfact! = (W, Y, _, Δt, t) -> W .= Δt * λ - 1,
        tgrad! = (∂Y∂t, Y, _, t) -> ∂Y∂t .= -(λ + 2 * t + λ * t^2) / (1 + t^2)^2,
    )
end

function analytic_linear_no_source_test_cts(::Type{FT}) where {FT}
    λ = FT(-1)
    ClimaIntegratorTestCase(;
        test_name = "analytic_linear_no_source",
        linear_implicit = true,
        t_end = FT(10),
        Y₀ = FT[1],
        analytic_sol = (t) -> [exp(λ * t)],
        tendency! = (Yₜ, Y, _, t) -> Yₜ .= λ .* Y,
        implicit_tendency! = (Yₜ, Y, _, t) -> Yₜ .= λ .* Y,
        Wfact! = (W, Y, _, dtγ, t) -> W .= dtγ * λ - 1,
    )
end

# From Section 1.2 of "Example Programs for ARKode v4.4.0" by D. R. Reynolds
function ark_analytic_nonlin_test_cts(::Type{FT}) where {FT}
    ClimaIntegratorTestCase(;
        test_name = "ark_analytic_nonlin",
        linear_implicit = false,
        t_end = FT(10),
        Y₀ = FT[0],
        analytic_sol = (t) -> [log(t^2 / 2 + t + 1)],
        tendency! = (Yₜ, Y, _, t) -> Yₜ .= (t + 1) .* exp.(.-Y),
        Wfact! = (W, Y, _, Δt, t) -> W .= (-Δt * (t + 1) .* exp.(.-Y) .- 1),
        tgrad! = (∂Y∂t, Y, _, t) -> ∂Y∂t .= exp.(.-Y),
    )
end

# From Section 5.1 of "Example Programs for ARKode v4.4.0" by D. R. Reynolds
function ark_analytic_sys_test_cts(::Type{FT}) where {FT}
    λ = FT(-100) # increase magnitude for more stiffness
    V = FT[1 -1 1; -1 2 1; 0 -1 2]
    V⁻¹ = FT[5 1 -3; 2 2 -2; 1 1 1] / 4
    D = Diagonal(FT[-1 / 2, -1 / 10, λ])
    A = V * D * V⁻¹
    I = LinearAlgebra.I(3)
    Y₀ = FT[1, 1, 1]
    ClimaIntegratorTestCase(;
        test_name = "ark_analytic_sys",
        linear_implicit = true,
        t_end = FT(1 / 20),
        Y₀,
        analytic_sol = (t) -> V * exp(D * t) * V⁻¹ * Y₀,
        tendency! = (Yₜ, Y, _, t) -> mul!(Yₜ, A, Y),
        Wfact! = (W, Y, _, Δt, t) -> W .= Δt .* A .- I,
    )
end

function onewaycouple_mri_test_cts(::Type{FT}) where {FT}
    Y₀ = FT[1, 0, 2]
    function analytic_sol(t)
        Y = similar(Y₀)
        Y[1] = cos(50 * t)
        Y[2] = sin(50 * t)
        Y[3] = 5051 / 2501 * exp(-t) - 49 / 2501 * cos(50 * t) + 51 / 2501 * sin(50 * t)
        return Y
    end
    L = FT[0 -50 0; 50 0 0; 1 1 -1]
    I = LinearAlgebra.I(3)
    ClimaIntegratorTestCase(;
        test_name = "ark_onewaycouple_mri",
        linear_implicit = true,
        t_end = FT(1),
        Y₀ = FT[1, 0, 2],
        analytic_sol,
        tendency! = (Yₜ, Y, _, t) -> mul!(Yₜ, L, Y),
        Wfact! = (W, Y, _, Δt, t) -> W .= Δt .* L .- I,
    )
end

"""
    climacore_2Dheat_test_cts(::Type{<:AbstractFloat})

2D diffusion test problem. See [`2D diffusion problem`](@ref) for more details.
"""
function climacore_2Dheat_test_cts(::Type{FT}) where {FT}
    context = ClimaComms.context()
    dss_tendency = true

    n_elem_x = 2
    n_elem_y = 2
    n_poly = 2
    n_x = 1 # denoted by n above
    n_y = 1 # denoted by m above
    f_0 = FT(0) # denoted by f̂₀ above
    Δλ = FT(1) # denoted by Δλ̂ above
    t_end = FT(0.05) # denoted by t̂ above

    domain = Domains.RectangleDomain(
        Domains.IntervalDomain(Geometry.XPoint(FT(0)), Geometry.XPoint(FT(1)), periodic = true),
        Domains.IntervalDomain(Geometry.YPoint(FT(0)), Geometry.YPoint(FT(1)), periodic = true),
    )
    mesh = Meshes.RectilinearMesh(domain, n_elem_x, n_elem_y)
    topology = Topologies.Topology2D(context, mesh)
    quadrature = Spaces.Quadratures.GLL{n_poly + 1}()
    space = Spaces.SpectralElementSpace2D(topology, quadrature)
    (; x, y) = Fields.coordinate_field(space)

    λ = (2 * FT(π))^2 * (n_x^2 + n_y^2)
    φ_sin_sin = @. sin(2 * FT(π) * n_x * x) * sin(2 * FT(π) * n_y * y)

    init_state = Fields.FieldVector(; u = φ_sin_sin)

    wdiv = Operators.WeakDivergence()
    grad = Operators.Gradient()
    function T_exp!(tendency, state, _, t)
        @. tendency.u = wdiv(grad(state.u)) + f_0 * exp(-(λ + Δλ) * t) * φ_sin_sin
        dss_tendency && Spaces.weighted_dss!(tendency.u)
    end

    function dss!(state, _, t)
        dss_tendency || Spaces.weighted_dss!(state.u)
    end

    function analytic_sol(t)
        state = similar(init_state)
        @. state.u = (1 + f_0 / Δλ * (1 - exp(-Δλ * t))) * exp(-λ * t) * φ_sin_sin
        return state
    end

    tendency_func = ClimaODEFunction(; T_exp!, dss!)
    split_tendency_func = tendency_func
    make_prob(func) = ODEProblem(func, init_state, (FT(0), t_end), nothing)
    IntegratorTestCase(
        "2D Heat Equation",
        false,
        t_end,
        analytic_sol,
        make_prob(tendency_func),
        make_prob(split_tendency_func),
    )
end

function climacore_1Dheat_test_cts(::Type{FT}) where {FT}
    n_elem_z = 10
    n_z = 1
    f_0 = FT(0) # denoted by f̂₀ above
    Δλ = FT(1) # denoted by Δλ̂ above
    t_end = FT(0.1) # denoted by t̂ above

    domain = Domains.IntervalDomain(Geometry.ZPoint(FT(0)), Geometry.ZPoint(FT(1)), boundary_names = (:bottom, :top))
    mesh = Meshes.IntervalMesh(domain, nelems = n_elem_z)
    space = Spaces.FaceFiniteDifferenceSpace(mesh)
    (; z) = Fields.coordinate_field(space)

    λ = (2 * FT(π) * n_z)^2
    φ_sin = @. sin(2 * FT(π) * n_z * z)

    init_state = Fields.FieldVector(; u = φ_sin)

    div = Operators.DivergenceC2F(; bottom = Operators.SetDivergence(FT(0)), top = Operators.SetDivergence(FT(0)))
    grad = Operators.GradientF2C()
    function T_exp!(tendency, state, _, t)
        @. tendency.u = div(grad(state.u)) + f_0 * exp(-(λ + Δλ) * t) * φ_sin
    end

    function analytic_sol(t)
        state = similar(init_state)
        @. state.u = (1 + f_0 / Δλ * (1 - exp(-Δλ * t))) * exp(-λ * t) * φ_sin
        return state
    end

    tendency_func = ClimaODEFunction(; T_exp!)
    split_tendency_func = tendency_func
    make_prob(func) = ODEProblem(func, init_state, (FT(0), t_end), nothing)
    IntegratorTestCase(
        "1D Heat Equation",
        false,
        t_end,
        analytic_sol,
        make_prob(tendency_func),
        make_prob(split_tendency_func),
    )
end

function climacore_1Dheat_test_implicit_cts(::Type{FT}) where {FT}
    n_elem_z = 100 # this should be larger than the explicit problem's n_elem_z
    n_z = 1
    f_0 = FT(0) # denoted by f̂₀ above
    Δλ = FT(1) # denoted by Δλ̂ above
    t_end = FT(0.1) # denoted by t̂ above

    domain = Domains.IntervalDomain(Geometry.ZPoint(FT(0)), Geometry.ZPoint(FT(1)), boundary_names = (:bottom, :top))
    mesh = Meshes.IntervalMesh(domain, nelems = n_elem_z)
    space = Spaces.FaceFiniteDifferenceSpace(mesh)
    (; z) = Fields.coordinate_field(space)

    λ = (2 * FT(π) * n_z)^2
    φ_sin = @. sin(2 * FT(π) * n_z * z)

    init_state = Fields.FieldVector(; u = φ_sin)

    div = Operators.DivergenceC2F(; bottom = Operators.SetDivergence(FT(0)), top = Operators.SetDivergence(FT(0)))
    grad = Operators.GradientF2C()
    function T_imp_func!(tendency, state, _, t)
        @. tendency.u = div(grad(state.u)) + f_0 * exp(-(λ + Δλ) * t) * φ_sin
    end

    div_matrix = ClimaCore.MatrixFields.operator_matrix(div)
    grad_matrix = ClimaCore.MatrixFields.operator_matrix(grad)
    function Wfact(W, Y, p, dtγ, t)
        name = @name(u)
        # NOTE: We need MatrixFields.⋅, not LinearAlgebra.⋅
        @. W.matrix[name, name] = dtγ * div_matrix() ⋅ grad_matrix() - (LinearAlgebra.I,)
        return nothing
    end

    function tgrad(∂Y∂t, state, _, t)
        @. ∂Y∂t.u = -f_0 * (λ + Δλ) * exp(-(λ + Δλ) * t) * φ_sin
    end

    jacobian = ClimaCore.MatrixFields.FieldMatrix(
        (@name(u), @name(u)) => similar(φ_sin, ClimaCore.MatrixFields.TridiagonalMatrixRow{FT}),
    )

    jac_prototype = FieldMatrixWithSolver(jacobian, init_state)

    T_imp! = SciMLBase.ODEFunction(T_imp_func!; jac_prototype = jac_prototype, Wfact = Wfact, tgrad = tgrad)

    function analytic_sol(t)
        state = similar(init_state)
        @. state.u = (1 + f_0 / Δλ * (1 - exp(-Δλ * t))) * exp(-λ * t) * φ_sin
        return state
    end

    tendency_func = ClimaODEFunction(; T_imp!)
    split_tendency_func = tendency_func
    make_prob(func) = ODEProblem(func, init_state, (FT(0), t_end), nothing)
    IntegratorTestCase(
        "1D Heat Equation Implicit",
        false,
        t_end,
        analytic_sol,
        make_prob(tendency_func),
        make_prob(split_tendency_func),
    )
end

# "Dynamical Core Model Intercomparison Project (DCMIP) Test Case Document" by
# Ullrich et al., Section 1.1
# (http://www-personal.umich.edu/~cjablono/DCMIP-2012_TestCaseDocument_v1.7.pdf)
# Implemented in flux form, with an optional limiter and hyperdiffusion.
# TODO: Use this as an integration test.
function deformational_flow_test(::Type{FT}; use_limiter = true, use_hyperdiffusion = true) where {FT}
    context = ClimaComms.context()
    # Table III
    # Note: the paper uses "a" in place of "R"
    R = FT(6371220)          # radius of Earth [m]
    g = FT(9.80616)          # gravitational acceleration [m/s^2]
    R_d = FT(287.0)          # gas constant for dry air [J/kg/K]

    # Table IX
    # Note: the paper specifies that λ_c1 = 5π/6 ∈ [0, 2π), and that
    # λ_c2 = 7π/6 ∈ [0, 2π), whereas we use λ_c1, λ_c2 ∈ [-π, π).
    z_top = FT(12000)         # altitude at model top [m]
    p_top = FT(25494.4)       # pressure at model top [Pa]
    T_0 = FT(300)             # isothermal atmospheric temperature [K]
    p_0 = FT(100000)          # reference pressure [Pa]
    τ = 60 * 60 * 24 * FT(12) # period of motion (12 days) [s]
    ω_0 = 23000 * FT(π) / τ   # maximum vertical pressure velocity [Pa/s]
    b = FT(0.2)               # normalized pressure depth of divergent layer
    λ_c1 = -FT(π) / 6         # initial longitude of first tracer
    λ_c2 = FT(π) / 6          # initial longitude of second tracer
    φ_c = FT(0)               # initial latitude of tracers
    z_c = FT(5000)            # initial altitude of tracers [m]
    R_t = R / 2               # horizontal half-width of tracers [m]
    Z_t = FT(1000)            # vertical half-width of tracers [m]

    # scale height [m] (Equation 3)
    H = R_d * T_0 / g

    # hyperviscosity coefficient [m^4] (specified in the limiter paper)
    D₄ = FT(6.6e14)

    centers = Geometry.LatLongZPoint.(rad2deg(φ_c), rad2deg.((λ_c1, λ_c2)), FT(0))

    # custom discretization (paper's discretization results in a very slow test)
    vert_nelems = 10
    horz_nelems = 4
    horz_npoly = 3

    vert_domain =
        Domains.IntervalDomain(Geometry.ZPoint(FT(0)), Geometry.ZPoint(z_top); boundary_names = (:bottom, :top))
    vert_mesh = Meshes.IntervalMesh(vert_domain, nelems = vert_nelems)
    z_topology = Topologies.IntervalTopology(context, vert_mesh)
    vert_cent_space = Spaces.CenterFiniteDifferenceSpace(z_topology)

    horz_domain = Domains.SphereDomain(R)
    horz_mesh = Meshes.EquiangularCubedSphere(horz_domain, horz_nelems)
    horz_topology = Topologies.Topology2D(context, horz_mesh)
    horz_quad = Spaces.Quadratures.GLL{horz_npoly + 1}()
    horz_space = Spaces.SpectralElementSpace2D(horz_topology, horz_quad)

    cent_space = Spaces.ExtrudedFiniteDifferenceSpace(horz_space, vert_cent_space)
    cent_coords = Fields.coordinate_field(cent_space)
    face_space = Spaces.FaceExtrudedFiniteDifferenceSpace(cent_space)

    # initial density (Equation 8)
    cent_ρ = @. p_0 / (R_d * T_0) * exp(-cent_coords.z / H)

    # initial tracer concentrations (Equations 28--35)
    cent_q = map(cent_coords) do coord
        z = coord.z
        φ = deg2rad(coord.lat)

        ds = map(centers) do center
            r = Geometry.great_circle_distance(coord, center, Spaces.global_geometry(horz_space))
            return min(1, (r / R_t)^2 + ((z - z_c) / Z_t)^2)
        end
        in_slot = z > z_c && φ_c - FT(0.125) < φ < φ_c + FT(0.125)

        q1 = (1 + cos(FT(π) * ds[1])) / 2 + (1 + cos(FT(π) * ds[2])) / 2
        q2 = FT(0.9) - FT(0.8) * q1^2
        q3 = (ds[1] < FT(0.5) || ds[2] < FT(0.5)) && !in_slot ? FT(1) : FT(0.1)
        q4 = 1 - FT(0.3) * (q1 + q2 + q3)
        q5 = FT(1)
        return (; q1, q2, q3, q4, q5)
    end

    init_state = Fields.FieldVector(; cent_ρ, cent_ρq = cent_ρ .* cent_q)

    # current wind vector (Equations 15--26)
    current_cent_wind_vector = Fields.Field(Geometry.UVWVector{FT}, cent_space)
    current_face_wind_vector = Fields.Field(Geometry.UVWVector{FT}, face_space)
    function wind_vector(coord, ρ, t)
        z = coord.z
        φ = deg2rad(coord.lat)
        λ = deg2rad(coord.long)

        p = p_0 * exp(-g * z / (R_d * T_0)) # initial pressure (Equation 1)
        λ′ = λ - 2 * FT(π) * t / τ
        k = 10 * R / τ

        u_a = k * sin(λ′)^2 * sin(2 * φ) * cos(FT(π) * t / τ) + 2 * FT(π) * R / τ * cos(φ)
        u_d =
            ω_0 * R / (b * p_top) *
            cos(λ′) *
            cos(φ)^2 *
            cos(2 * FT(π) * t / τ) *
            (-exp((p - p_0) / (b * p_top)) + exp((p_top - p) / (b * p_top)))
        u = u_a + u_d
        v = k * sin(2 * λ′) * cos(φ) * cos(FT(π) * t / τ)
        s = 1 + exp((p_top - p_0) / (b * p_top)) - exp((p - p_0) / (b * p_top)) - exp((p_top - p) / (b * p_top))
        ω = ω_0 * sin(λ′) * cos(φ) * cos(2 * FT(π) * t / τ) * s
        w = -ω / (g * ρ)

        return Geometry.UVWVector(u, v, w)
    end

    horz_div = Operators.Divergence()
    horz_wdiv = Operators.WeakDivergence()
    horz_grad = Operators.Gradient()
    cent_χ = similar(cent_q)
    function T_lim!(tendency, state, _, t)
        @. current_cent_wind_vector = wind_vector(cent_coords, state.cent_ρ, t)
        @. tendency.cent_ρ = -horz_div(state.cent_ρ * current_cent_wind_vector)
        @. tendency.cent_ρq = -horz_div(state.cent_ρq * current_cent_wind_vector)
        use_hyperdiffusion || return nothing
        @. cent_χ = horz_wdiv(horz_grad(state.cent_ρq / state.cent_ρ))
        Spaces.weighted_dss!(cent_χ)
        @. tendency.cent_ρq += -D₄ * horz_wdiv(state.cent_ρ * horz_grad(cent_χ))
        return nothing
    end

    limiter = Limiters.QuasiMonotoneLimiter(cent_q; rtol = FT(0))
    function lim!(state, _, t, ref_state)
        use_limiter || return nothing
        Limiters.compute_bounds!(limiter, ref_state.cent_ρq, ref_state.cent_ρ)
        Limiters.apply_limiter!(state.cent_ρq, state.cent_ρ, limiter)
        return nothing
    end

    vert_div = Operators.DivergenceF2C()
    vert_interp = Operators.InterpolateC2F(top = Operators.Extrapolate(), bottom = Operators.Extrapolate())
    function T_exp!(tendency, state, _, t)
        @. current_face_wind_vector = wind_vector(face_coords, vert_interp(state.cent_ρ), t)
        @. tendency.cent_ρ = -vert_div(vert_interp(state.cent_ρ) * current_face_wind_vector)
        @. tendency.cent_ρq = -vert_div(vert_interp(state.cent_ρq) * current_face_wind_vector)
    end

    function dss!(state, _, t)
        Spaces.weighted_dss!(state.q)
    end

    function analytic_sol(t)
        t ∈ (0, τ) || error("Analytic solution only defined at start and end")
        return copy(init_state)
    end

    tendency_func = ClimaODEFunction(; T_lim!, T_exp!, lim!, dss!)
    split_tendency_func = tendency_func
    make_prob(func) = ODEProblem(func, init_state, (FT(0), τ), nothing)
    IntegratorTestCase(
        "Deformational Flow",
        false,
        τ,
        analytic_sol,
        make_prob(tendency_func),
        make_prob(split_tendency_func),
    )
end

# "A standard test case suite for two-dimensional linear transport on the
# sphere" by Lauritzen et al.
# (https://gmd.copernicus.org/articles/5/887/2012/gmd-5-887-2012.pdf)
# Implemented in flux form, with an optional limiter and hyperdiffusion.
function horizontal_deformational_flow_test(::Type{FT}; use_limiter = true, use_hyperdiffusion = true) where {FT}
    context = ClimaComms.context()
    # constants (using the same notation as deformational_flow_test)
    R = FT(6371220)           # radius of Earth [m]
    τ = 60 * 60 * 24 * FT(12) # period of motion (12 days) [s]
    λ_c1 = -FT(π) / 6         # initial longitude of first tracer
    λ_c2 = FT(π) / 6          # initial longitude of second tracer
    φ_c = FT(0)               # initial latitude of tracers
    R_t = R / 2               # horizontal half-width of tracers [m]
    D₄ = FT(6.6e14)           # hyperviscosity coefficient [m^4] (specified in the limiter paper)

    centers = Geometry.LatLongPoint.(rad2deg(φ_c), rad2deg.((λ_c1, λ_c2)))

    # 1.5° resolution on the equator: 360° / (4 * nelems * npoly) = 1.5°
    nelems = 20
    npoly = 3

    domain = Domains.SphereDomain(R)
    mesh = Meshes.EquiangularCubedSphere(domain, nelems)
    topology = Topologies.Topology2D(context, mesh)
    quad = Spaces.Quadratures.GLL{npoly + 1}()
    space = Spaces.SpectralElementSpace2D(topology, quad)
    coords = Fields.coordinate_field(space)

    # initial conditions (Section 2.2)
    ρ = ones(space)
    q = map(coords) do coord
        φ = deg2rad(coord.lat)
        λ = deg2rad(coord.long)

        hs = map(centers) do center
            center′ = Geometry.CartesianPoint(center, Spaces.global_geometry(space))
            coord′ = Geometry.CartesianPoint(coord, Spaces.global_geometry(space))
            dist_squared = (coord′.x1 - center′.x1)^2 + (coord′.x2 - center′.x2)^2 + (coord′.x3 - center′.x3)^2
            # Note: the paper doesn't divide by R^2, which only works if R = 1
            return FT(0.95) * exp(-5 * dist_squared / R^2)
        end
        gaussian_hills = hs[1] + hs[2]
        rs = map(centers) do center
            return Geometry.great_circle_distance(coord, center, Spaces.global_geometry(space))
        end
        cosine_bells = if rs[1] < R_t
            FT(0.1) + FT(0.9) * (1 + cos(FT(π) * rs[1] / R_t)) / 2
        elseif rs[2] < R_t
            FT(0.1) + FT(0.9) * (1 + cos(FT(π) * rs[2] / R_t)) / 2
        else
            FT(0.1)
        end
        slotted_cylinders =
            if (
                (rs[1] <= R_t && abs(λ - λ_c1) >= R_t / 6R) ||
                (rs[2] <= R_t && abs(λ - λ_c2) >= R_t / 6R) ||
                (rs[1] <= R_t && abs(λ - λ_c1) < R_t / 6R && φ - φ_c < -5R_t / 12R) ||
                (rs[2] <= R_t && abs(λ - λ_c2) < R_t / 6R && φ - φ_c > 5R_t / 12R)
            )
                FT(1)
            else
                FT(0.1)
            end

        return (; gaussian_hills, cosine_bells, slotted_cylinders)
    end
    init_state = Fields.FieldVector(; ρ, ρq = ρ .* q)

    # current wind vector (Section 2.3)
    current_wind_vector = Fields.Field(Geometry.UVVector{FT}, space)
    function wind_vector(coord, t)
        φ = deg2rad(coord.lat)
        λ = deg2rad(coord.long)

        λ′ = λ - 2 * FT(π) * t / τ
        k = 10 * R / τ

        u = k * sin(λ′)^2 * sin(2 * φ) * cos(FT(π) * t / τ) + 2 * FT(π) * R / τ * cos(φ)
        v = k * sin(2 * λ′) * cos(φ) * cos(FT(π) * t / τ)

        return Geometry.UVVector(u, v)
    end

    div = Operators.Divergence()
    wdiv = Operators.WeakDivergence()
    grad = Operators.Gradient()
    χ = similar(q)
    function T_lim!(tendency, state, _, t)
        @. current_wind_vector = wind_vector(coords, t)
        @. tendency.ρ = -div(state.ρ * current_wind_vector)
        @. tendency.ρq = -div(state.ρq * current_wind_vector)
        use_hyperdiffusion || return nothing
        @. χ = wdiv(grad(state.ρq / state.ρ))
        Spaces.weighted_dss!(χ)
        @. tendency.ρq += -D₄ * wdiv(state.ρ * grad(χ))
        return nothing
    end

    limiter = Limiters.QuasiMonotoneLimiter(q; rtol = FT(0))
    function lim!(state, _, t, ref_state)
        use_limiter || return nothing
        Limiters.compute_bounds!(limiter, ref_state.ρq, ref_state.ρ)
        Limiters.apply_limiter!(state.ρq, state.ρ, limiter)
        return nothing
    end

    function dss!(state, _, t)
        Spaces.weighted_dss!(state.ρ)
        Spaces.weighted_dss!(state.ρq)
    end

    function analytic_sol(t)
        t ∈ (0, τ) || error("Analytic solution only defined at start and end")
        return copy(init_state)
    end

    tendency_func = ClimaODEFunction(; T_lim!, lim!, dss!)
    split_tendency_func = tendency_func
    make_prob(func) = ODEProblem(func, init_state, (FT(0), τ), nothing)
    IntegratorTestCase(
        "Horizontal Deformational Flow",
        false,
        τ,
        analytic_sol,
        make_prob(tendency_func),
        make_prob(split_tendency_func),
    )
end

function all_test_cases(::Type{FT}) where {FT}
    return [
        ark_analytic_nonlin_test_cts(FT),
        ark_analytic_sys_test_cts(FT),
        ark_analytic_test_cts(FT),
        onewaycouple_mri_test_cts(FT),
        climacore_2Dheat_test_cts(FT),
        climacore_1Dheat_test_cts(FT),
    ]
end
