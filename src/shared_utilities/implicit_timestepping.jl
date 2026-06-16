using ClimaCore.MatrixFields
import ClimaCore.MatrixFields: @name, ⋅
using ClimaCore: Spaces
import LinearAlgebra
import LinearAlgebra: I
using NVTX

export make_compute_jacobian, set_dfluxBCdY!, initialize_jacobian

"""
    make_compute_jacobian(model::AbstractModel)

Creates and returns a function which computes the entries
of the Jacobian matrix `W` in place.

If the implicit tendency function is given by
`T!(dY, Y, p, t) = make_implicit_tendency(model)`, the Jacobian
should be given by `W_{i,j}! = ∂T!_i/∂Y_j`, where `Y_j` is the
`j-th` state variable
and `T!_i` is the implicit tendency of the `i-th` state variable.

The default is that no updates are required, but this function
must be extended for models that use implicit timestepping.
"""
function make_compute_jacobian(model::AbstractModel)
    function compute_jacobian!(W, Y, p, dtγ, t) end
    return compute_jacobian!
end

"""
    set_dfluxBCdY!(::AbstractModel,
                  ::AbstractBC,
                  ::AbstractBoundary,
                  _...)::Union{ClimaCore.Fields.FieldVector, Nothing}

A function stub which returns the derivative of the implicit tendency
term of the `model` arising from the boundary condition,
with respect to the state Y.
"""
function set_dfluxBCdY!(
    ::AbstractModel,
    ::AbstractBC,
    ::AbstractBoundary,
    _...,
)::Union{ClimaCore.Fields.FieldVector, Nothing} end

"""
    initialize_jacobian(Y::ClimaCore.Fields.FieldVector)

Constructs a `MatrixFields.FieldMatrixWithSolver` residual Jacobian matrix
populated with ClimaLand-specific values based on the state vector `Y`;
`Δt ∂Ẏ_i/∂Y_j - δ_{i,j}.

For variables `Y_i` that have zero values of
`∂Ẏ_i∂Y_i`, the corresponding block in the residual matrix
is a negative identity. We refer to these variables ``identity block"
variables. Note that their contribution to the residual matrix is
*not* zero.

Variables `Y_i` that have a nonzero `∂Ẏ_i∂Y_i` will have a term plus
the negative of identity matrix. We refer to these variables as
``non identity block" variables.

Offdiagonal blocks are assumed to be zero by ClimaTimeSteppers, and are only added
as needed.
"""
function initialize_jacobian(Y::ClimaCore.Fields.FieldVector)
    FT = eltype(Y)
    # Only add jacobian blocks for fields that are in Y for this model
    is_in_Y(var) = MatrixFields.has_field(Y, var)

    nonid_block_vars = (
        @name(soil.ϑ_l),
        @name(soil.ρe_int),
        @name(canopy.energy.T),
        @name(soilco2.CO2),
        @name(soilco2.O2),
    )
    id_block_vars = (# can replace scalar_field_names(Y), setdiff(names, non_identity_blocks_in_field_name_set)
        @name(∫F_vol_e_dt),
        @name(∫F_vol_liq_water_dt),
        @name(soilco2.SOC),
        @name(soil.θ_i),
        @name(canopy.hydraulics.ϑ_l),
        @name(snow.S),
        @name(snow.S_l),
        @name(snow.U),
        @name(snow.Z),
        @name(snow.P_avg),
        @name(snow.T_avg),
        @name(snow.R_avg),
        @name(snow.Qrel_avg),
        @name(snow.u_avg),
        @name(snow.A),
        @name(lake.U),
    )

    # Filter out the variables that are not in this model's state, `Y`
    available_nonid_vars =
        MatrixFields.unrolled_filter(is_in_Y, nonid_block_vars)
    available_id_vars = MatrixFields.unrolled_filter(is_in_Y, id_block_vars)

    # Helper functions for getting the non-identity blocks
    # intialized correctly
    get_jac_type(
        space::Union{
            Spaces.FiniteDifferenceSpace,
            Spaces.ExtrudedFiniteDifferenceSpace,
        },
        FT,
    ) = MatrixFields.TridiagonalMatrixRow{FT}
    get_jac_type(
        space::Union{Spaces.PointSpace, Spaces.SpectralElementSpace2D},
        FT,
    ) = MatrixFields.DiagonalMatrixRow{FT}

    get_j_field(space, FT) = zeros(get_jac_type(space, FT), space)

    # Diagonal blocks
    # Non identity
    nonid_blocks = MatrixFields.unrolled_map(
        var ->
            (var, var) =>
                get_j_field(axes(MatrixFields.get_field(Y, var)), FT),
        available_nonid_vars,
    )

    # Identity
    # Note: We have to use FT(-1) * I instead of -I because inv(-1) == -1.0,
    # which means that multiplying inv(-1) by a Float32 will yield a Float64.
    id_blocks = MatrixFields.unrolled_map(
        var -> (var, var) => FT(-1) * I,
        available_id_vars,
    )

    # Offdiagonal blocks
    # Here, we follow the convention that each pair has order `Ẏ_i, Y_j` to produce `∂Ẏ_i∂Y_j`
    off_diagonal_pairs = (
        (@name(soil.ρe_int), @name(soil.ϑ_l)),
        (@name(∫F_vol_e_dt), @name(canopy.energy.T)),
    )
    available_off_diagonal_pairs = MatrixFields.unrolled_filter(
        pair -> all(is_in_Y.(pair)),
        off_diagonal_pairs,
    )
    offdiagonal_blocks = MatrixFields.unrolled_map(
        pair ->
            (pair[1], pair[2]) =>
                get_j_field(axes(MatrixFields.get_field(Y, pair[1])), FT),
        available_off_diagonal_pairs,
    )
    matrix = MatrixFields.FieldMatrix(
        nonid_blocks...,
        offdiagonal_blocks...,
        id_blocks...,
    )

    # Choose algorithm based on whether off-diagonal blocks are present
    if is_in_Y(@name(soil.ρe_int)) && is_in_Y(@name(soil.ϑ_l))
        if is_in_Y(@name(canopy.energy.T)) && is_in_Y(@name(∫F_vol_e_dt))
            # Set up lower triangular solver for block Jacobian with off-diagonal blocks
            # Specify which variable to compute ∂T_x/∂x for first
            alg₁ = MatrixFields.BlockLowerTriangularSolve(@name(soil.ϑ_l))
            alg = MatrixFields.BlockLowerTriangularSolve(
                (
                    @name(soil.ϑ_l),
                    @name(soil.ρe_int),
                    @name(canopy.energy.T)
                )...;
                alg₁,
            )
        else
            alg = MatrixFields.BlockLowerTriangularSolve(@name(soil.ϑ_l))
        end
    else
        # Set up block diagonal solver for block Jacobian with no off-diagonal blocks
        alg = MatrixFields.BlockDiagonalSolve()
    end
    solver = MatrixFields.FieldMatrixSolver(alg, matrix, Y)

    return MatrixFields.FieldMatrixWithSolver(matrix, solver)
end
