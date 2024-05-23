export make_tendency_jacobian, make_update_jacobian, set_dfluxBCdY!


"""
   make_tendency_jacobian(model::AbstractModel)

Creates and returns a function which updates the auxiliary
variables `p` in place and then updates the entries of the
Jacobian matrix `W` for the `model` in place.

The default is that no updates are required, no implicit tendency is
present, and hence the timestepping is entirely explicit.

Note that the returned function `tendency_jacobian!` should be
used as `Wfact!` in `ClimaTimeSteppers.jl` and `SciMLBase.jl`.
"""
function make_tendency_jacobian(model::AbstractModel)
    update_aux! = make_update_aux(model)
    update_boundary_fluxes! = make_update_boundary_fluxes(model)
    update_jacobian! = make_update_jacobian(model)
    function tendency_jacobian!(W, Y, p, dtγ, t)
        update_aux!(p, Y, t)
        update_boundary_fluxes!(p, Y, t)
        update_jacobian!(W, Y, p, dtγ, t)
    end
    return tendency_jacobian!
end

"""
    make_update_jacobian(model::AbstractModel)

Creates and returns a function which updates the entries
of the Jacobian matrix `W` in place.

If the implicit tendency function is given by
`T!(dY, Y, p, t) = make_implicit_tendency(model)`, the Jacobian
should be given by `W_{i,j}! = ∂T!_i/∂Y_j`, where `Y_j` is the
`j-th` state variable
and `T!_i` is the implicit tendency of the `i-th` state variable.

The default is that no updates are required, no implicit tendency is
present, and hence the timestepping is entirely explicit.
"""
function make_update_jacobian(model::AbstractModel)
    function update_jacobian!(W, Y, p, dtγ, t) end
    return update_jacobian
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
