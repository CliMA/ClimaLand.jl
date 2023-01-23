using ClimaCore
using ClimaCore: Operators, Spaces, Fields
using ClimaTimeSteppers
using UnPack
using LinearAlgebra
using DocStringExtensions
using OrdinaryDiffEq

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Soil
using ClimaLSM.Domains: HybridBox, SphericalShell, Column

include("../../SharedUtilities/boundary_conditions.jl")
include("../rre.jl")

FT = Float64
divf2c_op = Operators.DivergenceF2C(
    top = Operators.SetValue(FT(0.0)),
    bottom = Operators.SetValue(FT(0.0)),
)
divf2c_stencil = Operators.Operator2Stencil(divf2c_op)
gradc2f_op = Operators.GradientC2F()
gradc2f_stencil = Operators.Operator2Stencil(gradc2f_op)
interpc2f_op = Operators.InterpolateC2F()
interpc2f_stencil = Operators.Operator2Stencil(interpc2f_op)


struct TridiagonalW{R, J1, J2}
    # reference to dtγ, which is specified by the ODE solver
    dtγ_ref::R
    ∂ϑt∂ϑ::J1 # TODO change to ∂ϑₜ
    ∂ϑt∂ϑ_column_array::J2
end

# TODO add identityW

"""
    TridiagonalW(Y)

Constructor for the TridiagonalW struct.
Uses the space information from Y to allocate space for the
necessary fields.
"""
function TridiagonalW(Y)
    FT = eltype(Y.soil.ϑ_l)
    space = axes(Y.soil.ϑ_l)
    N = Spaces.nlevels(space)
    Jtype = Operators.StencilCoefs{-1, 1, NTuple{3, FT}}
    J = Fields.Field(Jtype, space)
    J_column_array = Tridiagonal(
        Array{FT}(undef, N - 1),
        Array{FT}(undef, N),
        Array{FT}(undef, N - 1),
    )
    dtγ_ref = Ref(zero(FT))

    args = (dtγ_ref, J, J_column_array)
    return TridiagonalW{typeof.(args)...}(args...)
end


"""
    Wfact!

Compute the entries of the Jacobian and overwrite W with them.
See overleaf for Jacobian derivation: https://www.overleaf.com/project/63be02f455f84a77642ef485
"""
# TODO add Wfact! method for identity input
function Wfact!(W::TridiagonalW, Y, p, dtγ, t)
    (; dtγ_ref, ∂ϑt∂ϑ) = W
    dtγ_ref[] = dtγ

    # TODO ask Kat - better way to get face_space?
    if model.domain.space isa CenterFiniteDifferenceSpace
        face_space = Spaces.FaceFiniteDifferenceSpace(axes(Y.soil.ϑ_l))
    elseif model.domain.space isa CenterExtrudedFiniteDifferenceSpace
        face_space = Spaces.FaceExtrudedFiniteDifferenceSpace(axes(Y.soil.ϑ_l))
    else
        error("invalid model space")
    end

    @. ∂ϑt∂ϑ = compose(
        divf2c_stencil(one(face_space)),
        (
            interpc2f_op(p.soil.K) *
            to_scalar_coefs(gradc2f_stencil(dψdθ(Y.soil.ϑ_l, ν, θ_r, vg_α, vg_n, vg_m, S_s)))
        )
    )
end

function Wfact!(W::UniformScaling{T}, x...) where T
    nothing
end

# Copied from https://github.com/CliMA/ClimaLSM.jl/blob/f41c497a12f91725ff23a9cd7ba8d563285f3bd8/examples/richards_implicit.jl#L152
# Checked against soiltest https://github.com/CliMA/ClimaLSM.jl/blob/e7eaf2e6fffaf64b2d824e9c5755d2f60fa17a69/test/Soil/soiltest.jl#L459
function dψdθ(θ, ν, θ_r, vg_α, vg_n, vg_m, S_s)
    S = (θ - θ_r) / (ν - θ_r)
    if S < 1.0
        return 1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r) *
               (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
               S^(-1 / vg_m - 1)
    else
        return 1.0 / S_s
    end
end

"""
    thomas_algorithm!(A, b)

Thomas algorithm for solving a linear system A x = b,
where A is a tri-diagonal matrix.
A and b are overwritten, solution is written to b.
Pass this as linsolve to ODEFunction.

Copied directly from https://github.com/CliMA/ClimaAtmos.jl/blob/99e44f4cd97307c4e8f760a16e7958d66d67e6e8/src/tendencies/implicit/schur_complement_W.jl#L410
"""
function thomas_algorithm!(A, b)
    nrows = size(A, 1)
    # first row
    @inbounds A[1, 2] /= A[1, 1]
    @inbounds b[1] /= A[1, 1]
    # interior rows
    for row in 2:(nrows - 1)
        @inbounds fac = A[row, row] - (A[row, row - 1] * A[row - 1, row])
        @inbounds A[row, row + 1] /= fac
        @inbounds b[row] = (b[row] - A[row, row - 1] * b[row - 1]) / fac
    end
    # last row
    @inbounds fac = A[nrows, nrows] - A[nrows - 1, nrows] * A[nrows, nrows - 1]
    @inbounds b[nrows] = (b[nrows] - A[nrows, nrows - 1] * b[nrows - 1]) / fac
    # back substitution
    for row in (nrows - 1):-1:1
        @inbounds b[row] -= b[row + 1] * A[row, row + 1]
    end
    return nothing
end

"""
    make_implicit_tendency(model::Soil.RichardsModel)

Construct the tendency function for the implicit terms of the RHS.
Adapted from https://github.com/CliMA/ClimaLSM.jl/blob/f41c497a12f91725ff23a9cd7ba8d563285f3bd8/examples/richards_implicit.jl#L173
"""
# compared math to make_rhs https://github.com/CliMA/ClimaLSM.jl/blob/main/src/Soil/rre.jl#L88
function make_implicit_tendency(model::Soil.RichardsModel)
    function implicit_tendency!(dY, Y, p, t)
        @unpack ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r = model.parameters
        # TODO when does aux get updated? should we be updating p.soil.K/ψ instead of just K/ψ here?
        (; K, ψ) = p.soil

        @. K = hydraulic_conductivity(
            K_sat,
            vg_m,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        )
        @. ψ = pressure_head(vg_α, vg_n, vg_m, θ_r, Y.soil.ϑ_l, ν, S_s)

        z = model.coordinates.z
        Δz_top, Δz_bottom = get_Δz(z)

        top_flux_bc = boundary_flux(
            model.boundary_conditions.water.top,
            TopBoundary(),
            Δz_top,
            p,
            t,
            model.parameters,
        )
        bot_flux_bc = boundary_flux(
            model.boundary_conditions.water.bottom,
            BottomBoundary(),
            Δz_bottom,
            p,
            t,
            model.parameters,
        )

        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector(top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector(bot_flux_bc)),
        )

        @. dY.soil.ϑ_l = -(divf2c_water(-interpc2f(K) * gradc2f_water(ψ + z)))
    end
    return implicit_tendency!
end

"""
    make_explicit_tendency(model::Soil.RichardsModel)

Construct the tendency function for the explicit terms of the RHS.
Adapted from https://github.com/CliMA/ClimaLSM.jl/blob/f41c497a12f91725ff23a9cd7ba8d563285f3bd8/examples/richards_implicit.jl#L204
"""
# compared math to make_rhs https://github.com/CliMA/ClimaLSM.jl/blob/main/src/Soil/rre.jl#L88
function make_explicit_tendency(model::Soil.RichardsModel)
    function explicit_tendency!(dY, Y, p, t)
        @unpack ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r = model.parameters

        @. p.soil.K = hydraulic_conductivity(
            K_sat,
            vg_m,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        )
        @. p.soil.ψ = pressure_head(vg_α, vg_n, vg_m, θ_r, Y.soil.ϑ_l, ν, S_s)
        hdiv = Operators.WeakDivergence()
        hgrad = Operators.Gradient()
        @. dY.soil.ϑ_l +=
            -hdiv(-p.soil.K * hgrad(p.soil.ψ + model.coordinates.z))
        Spaces.weighted_dss!(dY.soil.ϑ_l)
    end
    return explicit_tendency!
end

# TODO Need a linsolve...this doesnt do anything yet.
function linsolve!(::Type{Val{:init}}, f, u0; kwargs...)
    function _linsolve!(x, A, b, update_matrix = false; kwargs...)
        (; dtγ_ref, ∂Tt∂T) = W
        A = Tridiagonal(∂Tt∂t.upper_diagonal,
                        ∂Tt∂t.diagonal,
                        ∂Tt∂t.lower_diagonal)
        dtγ = dtγ_ref[]
    end
end

# Setup largely taken from test/Soil/soiltest.jl
function run()
    is_picard = true

    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0)
    zmax = FT(0)
    zmin = FT(-10)
    nelems = 50

    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
    top_flux_bc = FluxBC((p, t) -> eltype(t)(0.0))
    bot_flux_bc = FluxBC((p, t) -> eltype(t)(0.0))
    sources = ()
    boundary_fluxes = (; water = (top = top_flux_bc, bottom = bot_flux_bc))
    params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)

    soil = Soil.RichardsModel{FT}(;
        parameters = params,
        domain = soil_domain,
        boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil)

    # specify ICs
    function init_soil!(Ysoil, z, params)
        function hydrostatic_profile(
            z::FT,
            params::Soil.RichardsParameters{FT},
        ) where {FT}
            @unpack ν, vg_α, vg_n, vg_m, θ_r = params
            #unsaturated zone only, assumes water table starts at z_∇
            z_∇ = FT(-10)# matches zmin
            S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
            ϑ_l = S * (ν - θ_r) + θ_r
            return FT(ϑ_l)
        end
        Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
    end

    init_soil!(Y, coords.z, soil.parameters)

    t_start = 0.0
    t_end = 100.0
    dt = 10.0

    update_aux! = make_update_aux(soil)
    update_aux!(p, Y, t_start)

    W = is_picard ? -I : TridiagonalW(Y)

    implicit_tendency! = make_implicit_tendency(soil)
    explicit_tendency! = make_explicit_tendency(soil)
    jac_kwargs = (; jac_prototype = W, Wfact = Wfact!)
    alg_kwargs = (; linsolve = linsolve!)

    problem = SplitODEProblem(
        ODEFunction(
            implicit_tendency!;
            jac_kwargs...,
            tgrad = (∂Y∂t, Y, p, t) -> (∂Y∂t .= FT(0)),
        ),
        explicit_tendency!,
        Y,
        (t_start, t_end),
        p,
    )
    integrator = init(
        problem,
        Rosenbrock23(; alg_kwargs...);
        dt = dt,
        adaptive = false,
        progress = true,
        progress_steps = 1,
    )

    sol = @timev OrdinaryDiffEq.solve!(integrator)

end

run()
