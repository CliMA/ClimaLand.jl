export AbstractModel,
    AbstractImExModel,
    AbstractExpModel,
    make_imp_tendency,
    make_exp_tendency,
    make_compute_imp_tendency,
    make_compute_exp_tendency,
    make_update_aux,
    initialize_prognostic,
    initialize_auxiliary,
    make_update_boundary_fluxes,
    initialize_vars,
    initialize,
    prognostic_vars,
    auxiliary_vars,
    prognostic_types,
    auxiliary_types,
    prognostic_domain_names,
    auxiliary_domain_names,
    make_set_initial_cache,
    make_update_cache,
    add_drivers_to_cache,
    get_drivers,
    name,
    total_liq_water_vol_per_area!,
    total_energy_per_area!

import ClimaComms
import .Domains: coordinates
## Default methods for all models - to be in a seperate module at some point.
"""
    abstract type AbstractModel{FT <: AbstractFloat}

An abstract type for all models.
"""
abstract type AbstractModel{FT <: AbstractFloat} end

"""
    AbstractImExModel{FT} <: AbstractModel{FT}

An abstract type for models which must be treated implicitly (and which may
also have tendency terms that can be treated explicitly).
This inherits all the default function definitions from AbstractModel, as well
as `make_imp_tendency` and `make_compute_imp_tendency` defaults.
"""
abstract type AbstractImExModel{FT} <: AbstractModel{FT} end

"""
    AbstractExpModel{FT} <: AbstractModel{FT}

An abstract type for models which must be treated explicitly.
This inherits all the default function definitions from AbstractModel, as well
as a `make_imp_tendency` default.
"""
abstract type AbstractExpModel{FT} <: AbstractModel{FT} end

"""
    name(model::AbstractModel)

Returns a symbol of the model component name, e.g. :soil or :vegetation.
"""
name(model::AbstractModel) =
    error("`name` not implemented for $(Base.typename(typeof(model)).wrapper)")

"""
   prognostic_vars(m::AbstractModel)

Returns the prognostic variable symbols for the model in the form of a tuple.

Note that this default suggests that a model has no prognostic variables,
which is an invalid model setup. This function is meant to be extended for
all models.
"""
prognostic_vars(m::AbstractModel) = ()

"""
   prognostic_domain_names(m::AbstractModel)

Returns the domain names for the prognostic variables in the form of a tuple.

Examples: (:surface, :surface, :subsurface).

Note that this default suggests that a model has no prognostic variables,
which is an invalid model setup. This function is meant to be extended for
all models.
"""
prognostic_domain_names(m::AbstractModel) = ()

"""
   prognostic_types(m::AbstractModel{FT}) where {FT}

Returns the prognostic variable types for the model in the form of a tuple.

Types provided must have `ClimaCore.RecursiveApply.rzero(T::DataType)`
 defined. Common examples
 include
- Float64, Float32 for scalar variables (a scalar value at each
coordinate point)
- SVector{k,Float64} for a mutable but statically sized array of
 length `k` at each coordinate point.

Here, the coordinate points are those returned by coordinates(model).

Note that this default suggests that a model has no prognostic variables,
which is an invalid model setup. This function is meant to be extended for
all models.
"""
prognostic_types(m::AbstractModel) = ()

"""
   auxiliary_vars(m::AbstractModel)

Returns the auxiliary variable symbols for the model in the form of a tuple.
"""
auxiliary_vars(m::AbstractModel) = ()

"""
   auxiliary_domain_names(m::AbstractModel)

Returns the domain names for the auxiliary variables in the form of a tuple.

Examples: (:surface, :surface, :subsurface).
"""
auxiliary_domain_names(m::AbstractModel) = ()

"""
   auxiliary_types(m::AbstractModel{FT}) where {FT}

Returns the auxiliary variable types for the model in the form of a tuple.

Types provided must have `ClimaCore.RecursiveApply.rzero(T::DataType)`
defined. Common examples
 include
- Float64, Float32 for scalar variables (a scalar value at each
coordinate point)
- SVector{k,Float64} for a mutable but statically sized array of
 length `k` at each coordinate point.
- Note that Arrays, MVectors are not isbits and cannot be used.

Here, the coordinate points are those returned by coordinates(model).
"""
auxiliary_types(m::AbstractModel) = ()

"""
    make_update_aux(model::AbstractModel)

Return an `update_aux!` function that updates auxiliary parameters `p`.
"""
function make_update_aux(model::AbstractModel)
    function update_aux!(p, Y, t) end
    return update_aux!
end

"""
    make_update_boundary_fluxes(model::AbstractModel)

Return an `update_boundary_fluxes!` function that updates the auxiliary parameters in `p`
corresponding to boundary fluxes or interactions between componets..
"""
function make_update_boundary_fluxes(model::AbstractModel)
    function update_boundary_fluxes!(p, Y, t) end
    return update_boundary_fluxes!
end


"""
     make_update_cache(model::AbstractModel)

A helper function which updates all cache variables of a model;
currently only used in `set_initial_cache` since not all
cache variables are updated at the same time.
"""
function make_update_cache(model::AbstractModel)
    # if not forced using atmospheric/radiatiave drivers
    # this return (nothing, nothing)
    drivers = get_drivers(model)
    update_drivers! = make_update_drivers(drivers)
    update_aux! = make_update_aux(model)
    update_boundary_fluxes! = make_update_boundary_fluxes(model)
    function update_cache!(p, Y, t)
        update_drivers!(p, t)
        update_aux!(p, Y, t)
        update_boundary_fluxes!(p, Y, t)
    end
    return update_cache!
end

"""
    add_drivers_to_cache(p, model::AbstractModel)

Adds driver variables to the cache; the default is not
add anything, consistent with the default of no additional
driver variables in the cache.
"""
add_drivers_to_cache(p, model::AbstractModel) = p

"""
    make_imp_tendency(model::AbstractImExModel)

Returns an `imp_tendency` that updates auxiliary variables and
updates the prognostic state of variables that are stepped implicitly.

`compute_imp_tendency!` should be compatible with SciMLBase.jl solvers.
"""
function make_imp_tendency(model::AbstractImExModel)
    compute_imp_tendency! = make_compute_imp_tendency(model)
    update_aux! = make_update_aux(model)
    update_boundary_fluxes! = make_update_boundary_fluxes(model)
    function imp_tendency!(dY, Y, p, t)
        update_aux!(p, Y, t)
        update_boundary_fluxes!(p, Y, t)
        compute_imp_tendency!(dY, Y, p, t)
    end
    return imp_tendency!
end

"""
    make_imp_tendency(model::AbstractModel)

Returns an `imp_tendency` that does nothing. This model type is not
stepped explicity.
"""
function make_imp_tendency(model::AbstractModel)
    compute_imp_tendency! = make_compute_imp_tendency(model)
    function imp_tendency!(dY, Y, p, t)
        compute_imp_tendency!(dY, Y, p, t)
    end
end

"""
    make_exp_tendency(model::AbstractModel)

Returns an `exp_tendency` that updates auxiliary variables and
updates the prognostic state of variables that are stepped explicitly.

`compute_exp_tendency!` should be compatible with SciMLBase.jl solvers.
"""
function make_exp_tendency(model::AbstractModel)
    compute_exp_tendency! = make_compute_exp_tendency(model)
    update_aux! = make_update_aux(model)
    update_boundary_fluxes! = make_update_boundary_fluxes(model)
    function exp_tendency!(dY, Y, p, t)
        update_aux!(p, Y, t)
        update_boundary_fluxes!(p, Y, t)
        compute_exp_tendency!(dY, Y, p, t)
    end
    return exp_tendency!
end

"""
    make_compute_imp_tendency(model::AbstractModel)

Return a `compute_imp_tendency!` function that updates state variables
that we will be stepped implicitly.
This fallback sets all tendencies of this model to zero, which is appropriate
for models that do not have any implicit tendencies to update.
Note that we cannot set `dY .= 0` here because this would overwrite the
tendencies of all models in the case of an integrated LSM.

`compute_imp_tendency!` should be compatible with SciMLBase.jl solvers.
"""
function make_compute_imp_tendency(model::AbstractModel)
    function compute_imp_tendency!(dY, Y, p, t)
        getproperty(dY, name(model)) .= 0
    end
    return compute_imp_tendency!
end

"""
    make_compute_exp_tendency(model::AbstractModel)

Return a `compute_exp_tendency!` function that updates state variables
that we will be stepped explicitly.
This fallback sets all tendencies of this model to zero, which is appropriate
for models that do not have any explicit tendencies to update.
Note that we cannot set `dY .= 0` here because this would overwrite the
tendencies of all models in the case of an integrated LSM.

`compute_exp_tendency!` should be compatible with SciMLBase.jl solvers.
"""
function make_compute_exp_tendency(model::AbstractModel)
    function compute_exp_tendency!(dY, Y, p, t)
        getproperty(dY, name(model)) .= 0
    end
    return compute_exp_tendency!
end

"""
    make_set_initial_cache(model::AbstractModel)

Returns the set_initial_cache! function, which updates the auxiliary
state `p` in place with the initial values corresponding to Y(t=t0) = Y0.

In principle, this function is not needed, because in the very first evaluation
of either `explicit_tendency` or `implicit_tendency`, at t=t0, the auxiliary
state is updated using the initial conditions for Y=Y0. However,
without setting the initial `p` state prior to running the simulation,
the value of `p` in the saved output at t=t0 will be unset.

Furthermore, specific methods of this function may be useful for models
which store time indepedent spatially varying parameter fields in
the auxiliary state. In this case, `update_aux!` does not need to do
anything, but they do need to be set with the initial (constant) values
before the simulation can be carried out.
"""
function make_set_initial_cache(model::AbstractModel)
    update_cache! = make_update_cache(model)
    function set_initial_cache!(p, Y0, t0)
        update_cache!(p, Y0, t0)
    end
    return set_initial_cache!
end

"""
    initialize_prognostic(model::AbstractModel, state::NamedTuple)

Returns a FieldVector of prognostic variables for `model` with the required
structure, with values equal to `similar(state)`. This assumes that all
prognostic variables are defined over the entire domain,
and that all prognostic variables have the same dimension and type.

If a model has no prognostic variables,
the returned FieldVector contains only an empty array.

The input `state` is an array-like object, usually a
ClimaCore Field or a Vector{FT}.

Adjustments to this - for example because different prognostic variables
have different dimensions - require defining a new method.
"""
function initialize_prognostic(model::AbstractModel{FT}, state) where {FT}
    if length(prognostic_vars(model)) == 0
        throw(
            AssertionError("Model must have at least one prognostic variable."),
        )
    end
    state_nt = initialize_vars(
        prognostic_vars(model),
        prognostic_types(model),
        prognostic_domain_names(model),
        state,
        name(model),
    )
    return ClimaCore.Fields.FieldVector(; state_nt...)
end

"""
    initialize_auxiliary(model::AbstractModel, state::NamedTuple)

Returns a NamedTuple of auxiliary variables for `model` with the required
structure, with values equal to `similar(state)`. This assumes that all
 auxiliary variables are defined over the entire domain,
and that all auxiliary variables have the same dimension and type. The auxiliary
variables NamedTuple can also hold preallocated objects which are not Fields.

If a model has no auxiliary variables,
the returned NamedTuple contains only an empty array.

The input `state` is an array-like object, usually a
ClimaCore Field or a Vector{FT}.

Adjustments to this - for example because different auxiliary variables
have different dimensions - require defining a new method.
"""
function initialize_auxiliary(model::AbstractModel{FT}, state) where {FT}
    p = initialize_vars(
        auxiliary_vars(model),
        auxiliary_types(model),
        auxiliary_domain_names(model),
        state,
        name(model),
    )
    if :domain ∈ propertynames(model)
        p = add_dss_buffer_to_aux(p, model.domain)
    else
        error(
            "Your model does not contain a domain. If this is intended, you will need a new method of initialize_auxiliary.",
        )
    end
    return p
end

function initialize_vars(keys, types, domain_names, state, model_name)
    FT = eltype(state)
    if length(keys) == 0
        return (; model_name => nothing)
    else
        zero_states = map(zip(types, domain_names)) do (T, D)
            zero_instance = ClimaCore.RecursiveApply.rzero(T)
            map(_ -> zero_instance, getproperty(state, D))
        end
        return (; model_name => (; zip(keys, zero_states)...))
    end
end

Domains.coordinates(model::AbstractModel) = Domains.coordinates(model.domain)

"""
    add_drivers_to_cache(p::NamedTuple, model::AbstractModel, coords)

Creates the driver variable NamedTuple (atmospheric and radiative forcing, etc),
and merges it into `p` under the key `drivers`. If no driver variables
are required, `p` is returned unchanged.
"""
function add_drivers_to_cache(p::NamedTuple, model::AbstractModel, coords)
    drivers = get_drivers(model)
    if hasproperty(model, :parameters) &&
       hasproperty(model.parameters, :earth_param_set) &&
       hasproperty(drivers, :atmos) &&
       drivers.atmos isa ClimaLand.PrescribedAtmosphere
        if LP.thermodynamic_parameters(model.parameters.earth_param_set) !=
           drivers.atmos.thermo_params
            error(
                "earth_param_set is inconsistent between the model and the atmosphere",
            )
        end
    end
    driver_nt = initialize_drivers(drivers, coords)
    if driver_nt == (;)
        return p
    else
        return merge(p, (; drivers = driver_nt))
    end
end

"""
    get_drivers(model::AbstractModel)

Returns the `driver` objects for the model - atmospheric and radiative forcing, etc - as a tuple (atmos, radiation, ...). If no drivers are needed
by a model, an empty tuple should be returned
"""
function get_drivers(model::AbstractModel)
    return ()
end

"""
    initialize(model::AbstractModel)

Creates the prognostic and auxiliary states structures, but with unset
values; constructs and returns the coordinates for the `model` domain.
We may need to consider this default more as we add diverse components and
`Simulations`.
"""
function initialize(model::AbstractModel{FT}) where {FT}
    coords = Domains.coordinates(model)
    Y = initialize_prognostic(model, coords)
    p = initialize_auxiliary(model, coords)
    p = add_drivers_to_cache(p, model, coords)
    return Y, p, coords
end

function ClimaComms.context(model::AbstractModel)
    if :domain ∈ propertynames(model)
        return ClimaComms.context(model.domain)
    else
        error(
            "Your model does not contain a domain. If this is intended, you will need a new method of ClimaComms.context.",
        )
    end
end

function ClimaComms.device(model::AbstractModel)
    if :domain ∈ propertynames(model)
        return ClimaComms.device(model.domain)
    else
        error(
            "Your model does not contain a domain. If this is intended, you will need a new method of ClimaComms.device.",
        )
    end
end

"""
    total_liq_water_vol_per_area!(cache, model::AbstractModel, Y, p, t)

A function which updates `cache` in place with the total liquid water volume
per unit ground area for the `model`, computed from `Y`, `p`, and `t`.

While water mass is the fundamentally conserved quantity, soil modelling represents
water by an equivalent water volume using the density of water and ice at standard
temperature and pressure. Because of that, we report here the total volume of water
present in a model (per unit area) that would arise if all the water was in liquid phase.
This can be converted to a mass using the density of liquid water.

This includes the water in multiple phases. For example, if ice is present, the water
volume is computed using ratio of the density of ice to the density of liquid water.
"""
function total_liq_water_vol_per_area!(cache, model::AbstractModel, Y, p, t) end

"""
    total_energy_per_area!(cache, model::AbstractModel, Y, p, t)

A function which updates `cache` in place with the total energy
per unit ground area for the `model`, computed from `Y`, `p`, and `t`.
"""
function total_energy_per_area!(cache, model::AbstractModel, Y, p, t) end
