export make_update_jacobian, make_jacobian, AbstractTridiagonalW, ∂tendencyBC∂Y

"""
    make_update_jacobian(model::AbstractModel)

Creates and returns a function which updates the entries
of the Jacobian matrix `W` in place.

If the implicit tendency function is given by
`T!(dY, Y, p, t) = make_implicit_tendency(model)`, the Jacobian
should be given by `W_{i,j}! = ∂T!_i/∂Y_j`, where `Y_j` is the
`j-th` state variable
and `T!_i` is the implicit tendency of the `i-th` state variable.

This is a stub function to be extended for concrete instances of
AbstractImExModels.
"""
function make_update_jacobian(model::AbstractModel)
    function update_jacobian!(W, Y, p, dtγ, t) end
    return update_jacobian!
end

"""
    ∂tendencyBC∂Y(::AbstractModel,
                  ::AbstractBC,
                  ::AbstractBoundary,
                  _...)::Union{ClimaCore.Fields.FieldVector, Nothing}

A function stub which returns the derivative of the implicit tendency
term of the `model` arising from the boundary condition,
with respect to the state Y.
"""
function ∂tendencyBC∂Y(
    ::AbstractModel,
    ::AbstractBC,
    ::AbstractBoundary,
    _...,
)::Union{ClimaCore.Fields.FieldVector, Nothing} end

"""
    AbstractTridiagonalW

An abstract type for tridiagonal Jacobian matrices.
"""
abstract type AbstractTridiagonalW end

Base.similar(w::AbstractTridiagonalW) = w

"""
    make_jacobian(model::AbstractModel, Y::ClimaCore.Fields.FieldVector;
                  transform::Bool = false)::Union{Nothing,AbstractTridiagonalW}

Creates and returns a struct with allocated memory
for storing the Jacobian entries.

If the implicit tendency function is given by
`T!(dY, Y, p, t) = make_implicit_tendency(model)`, the Jacobian
should be given by `W_{i,j}! = ∂T!_i/∂Y_j`, where `Y_j` is the
`j-th` state variable
and `T!_i` is the implicit tendency of the `i-th` state variable.

This is a stub function to be extended for concrete instances of
AbstractImExModels.
"""
make_jacobian(
    model::AbstractModel,
    Y::ClimaCore.Fields.FieldVector;
    transform::Bool = false,
)::Union{Nothing, AbstractTridiagonalW} = nothing
