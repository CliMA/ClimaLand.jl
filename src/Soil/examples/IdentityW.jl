# This Jacobian used for true picard approximation
module IdentityJacobian
using ClimaCore: Operators, Spaces, Fields, Geometry
using LinearAlgebra
using NVTX
using Colors
using ClimaLSM
using ClimaLSM.Soil
using ClimaLSM.Domains: HybridBox, SphericalShell, Column

export IdentityW, Wfact!, make_implicit_tendency, make_explicit_tendency
to_scalar_coefs(vector_coefs) = map(vector_coef -> vector_coef.u₃, vector_coefs)

"""
    IdentityW{R}

Jacobian representation used for the true Picard method.
Here, J = 0 and W = -I, so there is no need to store ∂ϑₜ∂ϑ.
W_column_arrays would contain -I for each column, so we can
avoid storing them and instead set x .= b in ldiv for this case.
"""
struct IdentityW{R}
    # reference to dtγ, which is specified by the ODE solver
    dtγ_ref::R

    # whether this struct is used to compute Wfact_t or Wfact
    transform::Bool
end

# We only use Wfact, but the implicit/IMEX solvers require us to pass
# jac_prototype, then call similar(jac_prototype) to obtain J and Wfact. Here
# is a temporary workaround to avoid unnecessary allocations.
Base.similar(w::IdentityW) = w

"""
    Wfact!(W::IdentityW, _, _, dtγ, _)

In the case of the true Picard method, where W = -I, we don't update
W when Wfact! is called.
"""
function Wfact!(W::IdentityW, _, _, dtγ, _)
    (; dtγ_ref) = W
    dtγ_ref[] = dtγ
end

"""
    LinearAlgebra.ldiv!(x, A::IdentityW, b)

In the case of true picard, we have x = W*b where W = -I,
so we can set x = -b.
"""
function LinearAlgebra.ldiv!(x, A::IdentityW, b)
    (; dtγ_ref, transform) = A
    dtγ = dtγ_ref[]

    x .= -b

    # Apply transform (if needed) - used with ODE
    if transform
        Fields.bycolumn(axes(x.soil.ϑ_l)) do colidx
            x.soil.ϑ_l[colidx] .*= dtγ
        end
    end
end

# Keep linsolve! function for compatibility with OrdinaryDiffEq
linsolve!(::Type{Val{:init}}, f, u0; kwargs...) = _linsolve!
_linsolve!(x, A, b, update_matrix = false; kwargs...) =
    LinearAlgebra.ldiv!(x, A, b)

# Copied from https://github.com/CliMA/ClimaLSM.jl/blob/f41c497a12f91725ff23a9cd7ba8d563285f3bd8/examples/richards_implicit.jl#L152
# Checked against soiltest https://github.com/CliMA/ClimaLSM.jl/blob/e7eaf2e6fffaf64b2d824e9c5755d2f60fa17a69/test/Soil/soiltest.jl#L459
function dψdθ(θ, ν, θ_r, vg_α, vg_n, vg_m, S_s)
    S = (θ - θ_r) / (ν - θ_r) # TODO use effective_saturation
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
Used for the implicit tendency when running a mixed implicit/explicit solver.
Adapted from https://github.com/CliMA/ClimaLSM.jl/blob/f41c497a12f91725ff23a9cd7ba8d563285f3bd8/examples/richards_implicit.jl#L173
"""
function make_implicit_tendency(model::Soil.RichardsModel)
    update_aux! = ClimaLSM.make_update_aux(model)
    function implicit_tendency!(dY, Y, p, t)
        update_aux!(p, Y, t)

        z = ClimaCore.Fields.coordinate_field(model.domain.space).z
        Δz_top, Δz_bottom = get_Δz(z)

        top_flux_bc = ClimaLSM.boundary_flux(
            model.boundary_conditions.top.water,
            ClimaLSM.TopBoundary(),
            Δz_top,
            p,
            t,
            model.parameters,
        )
        bot_flux_bc = ClimaLSM.boundary_flux(
            model.boundary_conditions.bottom.water,
            ClimaLSM.BottomBoundary(),
            Δz_bottom,
            p,
            t,
            model.parameters,
        )

        divf2c_op = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector.(top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector.(bot_flux_bc)),
        )
        gradc2f_op = Operators.GradientC2F()
        interpc2f_op = Operators.InterpolateC2F()

        @. dY.soil.ϑ_l =
            -(divf2c_op(-interpc2f_op(p.soil.K) * gradc2f_op(p.soil.ψ + z)))
    end
    return implicit_tendency!
end

"""
    make_explicit_tendency(model::Soil.RichardsModel)

Construct the tendency function for the explicit terms of the RHS.
Used for the explicit tendency when running a mixed implicit/explicit solver.
Adapted from https://github.com/CliMA/ClimaLSM.jl/blob/f41c497a12f91725ff23a9cd7ba8d563285f3bd8/examples/richards_implicit.jl#L204
"""
function make_explicit_tendency(model::Soil.RichardsModel)
    update_aux! = ClimaLSM.make_update_aux(model)
    function explicit_tendency!(dY, Y, p, t)
        # set dY before updating it
        dY .= FT(0)

        update_aux!(p, Y, t)

        z = ClimaCore.Fields.coordinate_field(model.domain.space).z
        ClimaLSM.Soil.horizontal_components!(dY, model.domain, model, p, z)

        # Source terms
        for src in model.sources
            ClimaLSM.source!(dY, src, Y, p, model.parameters)
        end
    end
    return explicit_tendency!
end

end
