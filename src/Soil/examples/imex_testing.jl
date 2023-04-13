"""
Usage:
    Run this file using julia with the following args after the file name:
1. "imp" or "exp" - use the implicit or explicit solver
2. "true_picard" or "mod_picard" - specify the method for the solver to use
3. "iters_n" - use NewtonsMethod with `max_iters = n`
4. "dt_n" - use a timestep of `n` seconds

Ex: to run the implicit solver using modified Picard with 1 iteration at each
timestep, and a timestep of 1 second, run:
julia --project src/Soil/examples/imex_testing.jl imp mod_picard iters_1 dt_1

Note that when using the explicit solver, the second and third
arguments are irrelevant.
"""

using ClimaCore
using ClimaCore: Operators, Spaces, Fields, Geometry
using UnPack
using LinearAlgebra
using DocStringExtensions
using NVTX
using Colors
using DiffEqBase
using BenchmarkTools
using JLD2
using Plots

import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Soil
using ClimaLSM.Domains: HybridBox, SphericalShell, Column


FT = Float64

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

"""
    TridiagonalW{R, J, W, T}

Jacobian representation used for the modified Picard method.
"""
struct TridiagonalW{R, J, W, T}
    # reference to dtγ, which is specified by the ODE solver
    dtγ_ref::R

    # derivative of tendency - used to fill Jacobian
    ∂ϑₜ∂ϑ::J

    # array of tridiagonal matrices containing W for each column
    W_column_arrays::W

    # caches used to evaluate ldiv!
    temp1::T
    temp2::T

    # whether this struct is used to compute Wfact_t or Wfact
    transform::Bool
end

"""
    TridiagonalW(Y, transform::Bool)

Constructor for the TridiagonalW struct, used for modified Picard.
Uses the space information from Y to allocate space for the
necessary fields.
"""
function TridiagonalW(Y, transform::Bool)
    FT = eltype(Y.soil.ϑ_l)
    space = axes(Y.soil.ϑ_l)
    N = Spaces.nlevels(space)

    tridiag_type = Operators.StencilCoefs{-1, 1, NTuple{3, FT}}
    ∂ϑₜ∂ϑ = Fields.Field(tridiag_type, space)

    W_column_arrays = [
        LinearAlgebra.Tridiagonal(
            Array{FT}(undef, N - 1),
            Array{FT}(undef, N),
            Array{FT}(undef, N - 1),
        ) for _ in 1:Threads.nthreads()
    ]
    dtγ_ref = Ref(zero(FT))

    return TridiagonalW(
        dtγ_ref,
        ∂ϑₜ∂ϑ,
        W_column_arrays,
        similar(Y),
        similar(Y),
        transform,
    )
end

# We only use Wfact, but the implicit/IMEX solvers require us to pass
# jac_prototype, then call similar(jac_prototype) to obtain J and Wfact. Here
# is a temporary workaround to avoid unnecessary allocations.
Base.similar(w::IdentityW) = w
Base.similar(w::TridiagonalW) = w

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
    Wfact!(W::TridiagonalW, Y, p, dtγ, t)

Compute the entries of the Jacobian and overwrite W with them.
Used for the modified Picard method.
See overleaf for Jacobian derivation: https://www.overleaf.com/project/63be02f455f84a77642ef485
"""
function Wfact!(W::TridiagonalW, Y, p, dtγ, t)
    (; dtγ_ref, ∂ϑₜ∂ϑ) = W
    dtγ_ref[] = dtγ
    # TODO add parameters and BCs to TridiagonalW so we don't access global vars here
    (; ν, vg_α, vg_n, vg_m, S_s, θ_r) = soil.parameters

    if axes(Y.soil.ϑ_l) isa Spaces.CenterFiniteDifferenceSpace
        face_space = Spaces.FaceFiniteDifferenceSpace(axes(Y.soil.ϑ_l))
    elseif axes(Y.soil.ϑ_l) isa Spaces.CenterExtrudedFiniteDifferenceSpace
        face_space = Spaces.FaceExtrudedFiniteDifferenceSpace(axes(Y.soil.ϑ_l))
    else
        error("invalid model space")
    end

    z = ClimaCore.Fields.coordinate_field(axes(Y.soil.ϑ_l)).z
    Δz_top, Δz_bottom = get_Δz(z)

    top_flux_bc = ClimaLSM.boundary_flux(
        soil.boundary_conditions.top.water,
        ClimaLSM.TopBoundary(),
        Δz_top,
        p,
        t,
        soil.parameters,
    )
    bot_flux_bc = ClimaLSM.boundary_flux(
        soil.boundary_conditions.bottom.water,
        ClimaLSM.BottomBoundary(),
        Δz_bottom,
        p,
        t,
        soil.parameters,
    )

    divf2c_op = Operators.DivergenceF2C(
        top = Operators.SetValue(Geometry.WVector.(top_flux_bc)),
        bottom = Operators.SetValue(Geometry.WVector.(bot_flux_bc)),
    )
    divf2c_stencil = Operators.Operator2Stencil(divf2c_op)
    gradc2f_op = Operators.GradientC2F(
        top = Operators.SetGradient(Geometry.WVector.(FT(0))),
        bottom = Operators.SetGradient(Geometry.WVector.(FT(0))),
    )
    gradc2f_stencil = Operators.Operator2Stencil(gradc2f_op)
    interpc2f_op = Operators.InterpolateC2F()
    compose = Operators.ComposeStencils()

    # TODO create field of ones on faces once and store in W to reduce allocations
    ones_face_space = ones(face_space)
    @. ∂ϑₜ∂ϑ = compose(
        divf2c_stencil(Geometry.WVector(ones_face_space)),
        (
            interpc2f_op(p.soil.K) * to_scalar_coefs(
                gradc2f_stencil(
                    dψdθ(Y.soil.ϑ_l, ν, θ_r, vg_α, vg_n, vg_m, S_s),
                ),
            )
        ),
    )
end

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

# Keep linsolve! function for compatibility with OrdinaryDiffEq
linsolve!(::Type{Val{:init}}, f, u0; kwargs...) = _linsolve!
_linsolve!(x, A, b, update_matrix = false; kwargs...) =
    LinearAlgebra.ldiv!(x, A, b)

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

# Function required by Krylov.jl (x and b can be AbstractVectors)
# See https://github.com/JuliaSmoothOptimizers/Krylov.jl/issues/605 for a
# related issue that requires the same workaround.
function LinearAlgebra.ldiv!(x, A::TridiagonalW, b)
    A.temp1 .= b
    LinearAlgebra.ldiv!(A.temp2, A, A.temp1)
    x .= A.temp2
end

function LinearAlgebra.ldiv!(
    x::Fields.FieldVector,
    A::TridiagonalW,
    b::Fields.FieldVector,
)
    _ldiv!(x, A, b, axes(x.soil.ϑ_l))
end

# 2D space - only one column
function _ldiv!(
    x::Fields.FieldVector,
    A::TridiagonalW,
    b::Fields.FieldVector,
    space::ClimaCore.Spaces.FiniteDifferenceSpace,
)
    (; dtγ_ref, ∂ϑₜ∂ϑ, W_column_arrays, transform) = A
    dtγ = dtγ_ref[]

    _ldiv_serial!(
        x.soil.ϑ_l,
        b.soil.ϑ_l,
        dtγ,
        transform,
        ∂ϑₜ∂ϑ,
        W_column_arrays[Threads.threadid()], # can / should this be colidx?
    )
end

# 3D space - iterate over columns
function _ldiv!(
    x::Fields.FieldVector,
    A::TridiagonalW,
    b::Fields.FieldVector,
    space::ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace,
)
    (; dtγ_ref, ∂ϑₜ∂ϑ, W_column_arrays, transform) = A
    dtγ = dtγ_ref[]

    NVTX.@range "linsolve" color = Colors.colorant"lime" begin
        Fields.bycolumn(axes(x.soil.ϑ_l)) do colidx
            _ldiv_serial!(
                x.soil.ϑ_l[colidx],
                b.soil.ϑ_l[colidx],
                dtγ,
                transform,
                ∂ϑₜ∂ϑ[colidx],
                W_column_arrays[Threads.threadid()], # can / should this be colidx?
            )
        end
    end
end


function _ldiv_serial!(
    x_column,
    b_column,
    dtγ,
    transform,
    ∂ϑₜ∂ϑ_column,
    W_column_array,
)
    x_column .= b_column

    x_column_view = parent(x_column)

    @views W_column_array.dl .= dtγ .* parent(∂ϑₜ∂ϑ_column.coefs.:1)[2:end]
    W_column_array.d .= -1 .+ dtγ .* parent(∂ϑₜ∂ϑ_column.coefs.:2)
    @views W_column_array.du .=
        dtγ .* parent(∂ϑₜ∂ϑ_column.coefs.:3)[1:(end - 1)]

    thomas_algorithm!(W_column_array, x_column_view)

    # Apply transform (if needed)
    if transform
        x_column .*= dtγ
    end
    return nothing
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


is_imex_CTS_algo(::CTS.IMEXAlgorithm) = true
is_imex_CTS_algo(::DiffEqBase.AbstractODEAlgorithm) = false

is_implicit(::ODE.OrdinaryDiffEqImplicitAlgorithm) = true
is_implicit(::ODE.OrdinaryDiffEqAdaptiveImplicitAlgorithm) = true
is_implicit(ode_algo) = is_imex_CTS_algo(ode_algo)

is_rosenbrock(::ODE.Rosenbrock23) = true
is_rosenbrock(::ODE.Rosenbrock32) = true
is_rosenbrock(::DiffEqBase.AbstractODEAlgorithm) = false
use_transform(ode_algo) =
    !(is_imex_CTS_algo(ode_algo) || is_rosenbrock(ode_algo))


# Setup largely taken from ClimaLSM.jl/test/Soil/soiltest.jl
if isinteractive()
    is_true_picard = false
else
    is_true_picard = (ARGS[2] == "true_picard")
end

# parameters for clay from Bonan 2019 supplemental program 8.2
ν = FT(0.495)
K_sat = FT(0.0443 / 3600 / 100) # m/s
vg_n = FT(1.43)
vg_α = FT(0.026 * 100) # inverse meters
vg_m = FT(1) - FT(1) / vg_n
θ_r = FT(0.124)
S_s = FT(1e-3) #inverse meters

zmax = FT(0)
zmin = FT(-1.5)
nelems = 150

# soil_domain = HybridBox(;
#     zlim = (-10.0, 0.0),
#     xlim = (0.0, 100.0),
#     ylim = (0.0, 100.0),
#     nelements = (10, 10, 10),
#     npolynomial = 1,
#     periodic = (true, true),
# )
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems);

top_bc = Soil.MoistureStateBC((p, t) -> eltype(t)(ν - 1e-3))
bot_bc = Soil.FreeDrainage()
sources = ()
boundary_fluxes = (; top = (water = top_bc,), bottom = (water = bot_bc,))
params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)

soil = Soil.RichardsModel{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = sources,
)

Y, p, coords = initialize(soil)
Y.soil.ϑ_l = FT(0.24)

if isinteractive()
    dt = FT(64)
else
    dt = parse(FT, replace(ARGS[4][4:end], "p" => "."))
end

t_start = FT(0)
t_end = FT(1e6) #FT(256)

if isinteractive()
    max_iters = 1
else
    max_iters = parse(Int64, ARGS[3][7:end])
end

# if isinteractive()
#     convergence_cond = CTS.MaximumRelativeError(1e-10)
# else
#     convergence_cond = CTS.MaximumRelativeError(parse(FT, ARGS[5][13:end]))
# end

# TODO when conv_checker was used before, the solver thought it had converged
#  when it hadn't. For now, don't use a conv checker, but maybe try again later
# conv_checker = CTS.ConvergenceChecker(component_condition = convergence_cond)
conv_checker = nothing
ode_algo = CTS.IMEXAlgorithm(
    CTS.ARS343(),
    CTS.NewtonsMethod(
        max_iters = max_iters,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)

if isinteractive()
    is_imp = true
else
    is_imp = ARGS[1] == "imp"
end

if is_imp # mixed implicit/explicit case
    transform = use_transform(ode_algo)

    W =
        is_true_picard ? IdentityW(Ref(zero(FT)), transform) :
        TridiagonalW(Y, transform)

    implicit_tendency! = make_implicit_tendency(soil)
    explicit_tendency! = make_explicit_tendency(soil)

    jac_kwargs = if use_transform(ode_algo)
        (; jac_prototype = W, Wfact_t = Wfact!)
    else
        (; jac_prototype = W, Wfact = Wfact!)
    end

    problem = ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = explicit_tendency!,
            T_imp! = ODEFunction(implicit_tendency!; jac_kwargs...),
            dss! = ClimaLSM.Soil.dss!,
        ),
        Y,
        (t_start, t_end),
        p,
    )
else # explicit-only case
    ode_function! = make_ode_function(soil)

    problem = ODE.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = ode_function!,
            T_imp! = nothing,
            dss! = ClimaLSM.Soil.dss!,
        ),
        Y,
        (t_start, t_end),
        p,
    )
end

integrator = init(
    problem,
    ode_algo;
    dt = dt,
    adaptive = false,
    progress = true,
    saveat = t_start:1:t_end,
)

sol = ODE.solve!(integrator)

sol_vals = parent(sol.u[end].soil.ϑ_l)
zs = parent(Fields.coordinate_field(axes(sol.u[end].soil.ϑ_l)))

isdir("imex-res") ? nothing : mkpath("imex-res")
filename = "imex-res/rre_sol"
for a in ARGS
    a = string(a)
    replace(a, "." => "p")
    global filename *= "-" * string(a)
end

plot(
    parent(sol.u[1].soil.ϑ_l)[51:150],
    zs[51:150],
    xlims = (0.2, 0.5),
    xticks = 0.2:0.05:0.5,
    yticks = -1.0:0.2:0.0,
    title = "implicit sol'n with dt=$dt at t=0",
    xlabel = "vol. water content (m^3 m^-3)",
    ylabel = "depth (m)",
)
savefig(filename * "-start")

t_end = Int(t_end)
plot(
    parent(sol.u[end].soil.ϑ_l)[51:150],
    zs[51:150],
    xlims = (0.2, 0.5),
    xticks = 0.2:0.05:0.5,
    yticks = -1.0:0.2:0.0,
    title = "implicit sol'n with dt=$dt at t=$t_end s",
    xlabel = "vol. water content (m^3 m^-3)",
    ylabel = "depth (m)",
)
savefig(filename * "-$t_end")


b = @benchmark ODE.solve!(integrator)

filename *= ".jld2"

jldsave(filename; ϑ_l = sol_vals, z = zs, benchmark = b, dt = dt)
@show filename

@assert all(x -> x > θ_r, sol_vals)
@assert all(x -> x < ν * 1.03, sol_vals)