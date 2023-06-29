### These functions are here for now, but linsolve!, ldiv!, and
### thomas_algorithm! will soon be moved to ClimaCore.

import LinearAlgebra
import ClimaCore: Spaces

export linsolve!, thomas_algorithm!

# Function required by OrdinaryDiffEq.jl
linsolve!(::Type{Val{:init}}, f, u0; kwargs...) = _linsolve!
_linsolve!(x, A, b, update_matrix = false; kwargs...) =
    LinearAlgebra.ldiv!(x, A, b)

# Function required by Krylov.jl (x and b can be AbstractVectors)
# See https://github.com/JuliaSmoothOptimizers/Krylov.jl/issues/605 for a
# related issue that requires the same workaround.
function LinearAlgebra.ldiv!(x, A::AbstractTridiagonalW, b)
    A.temp1 .= b
    LinearAlgebra.ldiv!(A.temp2, A, A.temp1)
    x .= A.temp2
end

function LinearAlgebra.ldiv!(
    x::Fields.FieldVector,
    A::AbstractTridiagonalW,
    b::Fields.FieldVector,
)
    _ldiv!(x, A, b, axes(x.soil.ϑ_l))
end

# 2D space - only one column
function _ldiv!(
    x::Fields.FieldVector, # Δx
    A::AbstractTridiagonalW, # f'(x)
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
    A::AbstractTridiagonalW,
    b::Fields.FieldVector,
    space::Spaces.ExtrudedFiniteDifferenceSpace,
)
    (; dtγ_ref, ∂ϑₜ∂ϑ, W_column_arrays, transform) = A
    dtγ = dtγ_ref[]

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
