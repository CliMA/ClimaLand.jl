module TridiagonalJacobian
using ClimaCore: Operators, Spaces, Fields, Geometry
using LinearAlgebra
using NVTX
using Colors
using ClimaLSM
using ClimaLSM.Soil
using ClimaLSM.Domains: HybridBox, SphericalShell, Column

export TridiagonalW, make_Wfact, make_implicit_tendency, make_explicit_tendency
to_scalar_coefs(vector_coefs) = map(vector_coef -> vector_coef.u₃, vector_coefs)

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

function TridiagonalW(Y, transform::Bool)
    FT = eltype(Y.soil.ϑ_l)
    space = axes(Y.soil.ϑ_l)
    N = Spaces.nlevels(space)

    tridiag_type = Operators.StencilCoefs{-1, 1, NTuple{3, FT}}
    ∂ϑₜ∂ϑ = Fields.Field(tridiag_type, space)

    ∂ϑₜ∂ϑ.coefs[1] .= FT(NaN)
    ∂ϑₜ∂ϑ.coefs[2] .= FT(NaN)
    ∂ϑₜ∂ϑ.coefs[3] .= FT(NaN)

    arr1 = Array{FT}(undef, N - 1)

    W_column_arrays = [
        LinearAlgebra.Tridiagonal(
            zeros(FT, N - 1) .+ FT(NaN),
            zeros(FT, N) .+ FT(NaN),
            zeros(FT, N - 1) .+ FT(NaN),
        ) for _ in 1:Threads.nthreads()
    ]
    # dtγ_ref = Ref(zero(FT))
    dtγ_ref = Ref(FT(NaN))

    temp1 = similar(Y)
    @. temp1.soil.ϑ_l = FT(NaN)
    temp2 = similar(Y)
    @. temp2.soil.ϑ_l = FT(NaN)

    return TridiagonalW(
        dtγ_ref,
        ∂ϑₜ∂ϑ,
        W_column_arrays,
        temp1,
        temp2,
        transform,
    )
end

Base.similar(w::TridiagonalW) = w



function make_Wfact(model::RichardsModel)
    update_aux! = make_update_aux(model)
    function Wfact!(W::TridiagonalW, Y, p, dtγ, t)

        FT = eltype(Y.soil.ϑ_l)
        (; dtγ_ref, ∂ϑₜ∂ϑ) = W
        dtγ_ref[] = dtγ

        (; ν, vg_α, vg_n, vg_m, S_s, θ_r) = model.parameters
        update_aux!(p, Y, t)

        if axes(Y.soil.ϑ_l) isa Spaces.CenterFiniteDifferenceSpace
            face_space = Spaces.FaceFiniteDifferenceSpace(axes(Y.soil.ϑ_l))
        elseif axes(Y.soil.ϑ_l) isa Spaces.CenterExtrudedFiniteDifferenceSpace
            face_space =
                Spaces.FaceExtrudedFiniteDifferenceSpace(axes(Y.soil.ϑ_l))
        else
            error("invalid model space")
        end

        z = Fields.coordinate_field(axes(Y.soil.ϑ_l)).z
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

        ∂T_bc∂ϑN = ClimaLSM.Soil.∂tendency_bc_∂ϑN(model.boundary_conditions.top.water, ClimaLSM.TopBoundary(), Δz_top, Y, p, t, model.parameters)
            
        divf2c_op = Operators.DivergenceF2C(
        top = Operators.SetValue(Geometry.WVector.(FT(0))),
        bottom = Operators.SetValue(Geometry.WVector.(FT(0)))
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
            divf2c_stencil(Geometry.Covariant3Vector(ones_face_space)),
            (
                interpc2f_op(p.soil.K) * to_scalar_coefs(
                    gradc2f_stencil(
                        ClimaLSM.Soil.dψdθ(Y.soil.ϑ_l, ν, θ_r, vg_α, vg_n, vg_m, S_s),
                    ),
                )
            ),
        )
        #TODO: fix, right now it is hardcoded for a single column

        parent(∂ϑₜ∂ϑ.coefs.:2)[end] = parent(ClimaLSM.Domains.top_center_to_surface(∂ϑₜ∂ϑ.coefs.:2) .+ ∂T_bc∂ϑN)[1]

    end
    return Wfact!
end

# Function required by OrdinaryDiffEq.jl
linsolve!(::Type{Val{:init}}, f, u0; kwargs...) = _linsolve!
_linsolve!(x, A, b, update_matrix = false; kwargs...) =
    LinearAlgebra.ldiv!(x, A, b)

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
    x::Fields.FieldVector, # Δx
    A::TridiagonalW, # f'(x)
    b::Fields.FieldVector, # -f(x)
    space::Spaces.FiniteDifferenceSpace,
)
    (; dtγ_ref, ∂ϑₜ∂ϑ, W_column_arrays, transform) = A
    dtγ = dtγ_ref[]

    _ldiv_serial!(
        x.soil.ϑ_l,
        b.soil.ϑ_l,
        dtγ,
        transform,
        ∂ϑₜ∂ϑ,
        W_column_arrays[Threads.threadid()], # can / should this be colidx? TODO do we need this in single-column case
    )
end

# 3D space - iterate over columns
function _ldiv!(
    x::Fields.FieldVector,
    A::TridiagonalW,
    b::Fields.FieldVector,
    space::Spaces.ExtrudedFiniteDifferenceSpace,
)
    (; dtγ_ref, ∂ϑₜ∂ϑ, W_column_arrays, transform) = A
    dtγ = dtγ_ref[]

    NVTX.@range "linsolve" color = Colors.colorant"lime" begin
        Fields.bycolumn(space) do colidx
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


function make_implicit_tendency(model::Soil.RichardsModel)
    update_aux! = ClimaLSM.make_update_aux(model)
    function implicit_tendency!(dY, Y, p, t)
        update_aux!(p, Y, t)

        z = Fields.coordinate_field(model.domain.space).z
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

function explicit_tendency!(dY, Y, p, t) end

end
