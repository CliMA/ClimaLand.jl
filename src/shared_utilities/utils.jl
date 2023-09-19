import SciMLBase

"""
     heaviside(x::FT)::FT where {FT}

Computes the heaviside function.
"""
function heaviside(x::FT)::FT where {FT}
    if x > eps(FT)
        return FT(1.0)
    else
        return FT(0.0)
    end
end

"""
    to_scalar_coefs(vector_coefs)

Helper function used in computing tendencies of vertical diffusion terms.
"""
to_scalar_coefs(vector_coefs) = map(vector_coef -> vector_coef.u₃, vector_coefs)


"""
    add_dss_buffer_to_aux(p::NamedTuple, domain::Domains.AbstractDomain)

Fallback method for `add_dss_buffer_to_aux` which does not add a dss buffer.
"""
add_dss_buffer_to_aux(p::NamedTuple, domain::Domains.AbstractDomain) = p

"""
    add_dss_buffer_to_aux(
        p::NamedTuple,
        domain::Union{Domains.Plane, Domains.SphericalSurface},
    )

Adds a dss buffer corresponding to `domain.space` to `p` with the name `dss_buffer_2d`,
appropriate for a 2D domain.

This buffer is added so that we preallocate memory for the dss step and do not allocate it
at every timestep. We use a name which specifically denotes that
the buffer is on a 2d space. This is because some models
require both a buffer on the 3d space as well as on the surface
2d space, e.g. in the case when they have prognostic variables that are only
defined on the surface space.
"""
function add_dss_buffer_to_aux(
    p::NamedTuple,
    domain::Union{Domains.Plane, Domains.SphericalSurface},
)
    buffer =
        ClimaCore.Spaces.create_dss_buffer(ClimaCore.Fields.zeros(domain.space))
    return merge(p, (; dss_buffer_2d = buffer))
end

"""
    add_dss_buffer_to_aux(
        p::NamedTuple,
        domain::Union{Domains.HybridBox, Domains.SphericalShell},
    )

Adds a 3d dss buffer corresponding to `domain.space` to `p` with the name `dss_buffer_3d`,
appropriate for a 3D domain.

This buffer is added so that we preallocate memory for the dss step and do not allocate it
at every timestep. We use a name which specifically denotes that
the buffer is on a 3d space. This is because some models
require both a buffer on the 3d space as well as on the surface
2d space, e.g. in the case when they have prognostic variables that are only
defined on the surface space.
"""
function add_dss_buffer_to_aux(
    p::NamedTuple,
    domain::Union{Domains.HybridBox, Domains.SphericalShell},
)
    buffer =
        ClimaCore.Spaces.create_dss_buffer(ClimaCore.Fields.zeros(domain.space))
    return merge(p, (; dss_buffer_3d = buffer))
end


"""
      dss!(Y::ClimaCore.Fields.FieldVector, p::NamedTuple, t::FT)

Computes the weighted direct stiffness summation and updates `Y` in place.
In the case of a column domain, no dss operations are performed.
"""
function dss!(
    Y::ClimaCore.Fields.FieldVector{FT},
    p::NamedTuple,
    t::FT,
) where {FT}
    for key in propertynames(Y)
        property = getproperty(Y, key)
        dss_helper!(property, axes(property), p)
    end
end

"""
    dss_helper!(field_vec::ClimaCore.Fields.FieldVector, space, p::NamedTuple)

Method of `dss_helper!` which unpacks properties of Y when on a
domain that is 2-dimensional in the horizontal.

The assumption is that Y contains FieldVectors which themselves contain either
FieldVectors or Fields, and that the final unpacked variable is a Field.
This method is invoked when the current property itself contains additional
property(ies).
"""
function dss_helper!(
    field_vec::ClimaCore.Fields.FieldVector,
    space,
    p::NamedTuple,
)
    for key in propertynames(field_vec)
        property = getproperty(field_vec, key)
        dss_helper!(property, axes(property), p)
    end
end

"""
    dss_helper!(
        field::ClimaCore.Fields.Field,
        space::ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace,
        p::NamedTuple)

Method of `dss_helper!` which performs dss on a Field which is defined
on a 3-dimensional domain.

The assumption is that Y contains FieldVectors which themselves contain either
FieldVectors or Fields, and that the final unpacked variable is a Field.
This method is invoked when the element cannot be unpacked further.
We further assume that all fields in `Y` are defined on cell centers.
"""
function dss_helper!(
    field::ClimaCore.Fields.Field,
    space::ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace,
    p::NamedTuple,
)
    buffer = p.dss_buffer_3d
    ClimaCore.Spaces.weighted_dss!(field, buffer)
end

"""
    dss_helper!(
        field::ClimaCore.Fields.Field,
        space::ClimaCore.Spaces.AbstractSpectralElementSpace,
        p::NamedTuple)

Method of `dss_helper!` which performs dss on a Field which is defined
on a 2-dimensional domain.

The assumption is that Y contains FieldVectors which themselves contain either
FieldVectors or Fields, and that the final unpacked variable is a Field.
This method is invoked when the element cannot be unpacked further.
We further assume that all fields in `Y` are defined on cell centers.
"""
function dss_helper!(
    field::ClimaCore.Fields.Field,
    space::ClimaCore.Spaces.AbstractSpectralElementSpace,
    p::NamedTuple,
)
    buffer = p.dss_buffer_2d
    ClimaCore.Spaces.weighted_dss!(field, buffer)
end

"""
    dss_helper!(
        field::Union{ClimaCore.Fields.Field, Vector},
        space::Union{
            ClimaCore.Spaces.FiniteDifferenceSpace,
            ClimaCore.Spaces.PointSpace,
            Tuple,
        }, _)

Method of `dss_helper!` which does not perform dss.

This is intended for spaces that don't use spectral
elements (FiniteDifferenceSpace, PointSpace, etc).
Model components with no prognostic variables appear in Y as empty
Vectors, and also do not need dss.
"""
function dss_helper!(
    field::Union{ClimaCore.Fields.Field, Vector},
    space::Union{
        ClimaCore.Spaces.FiniteDifferenceSpace,
        ClimaCore.Spaces.PointSpace,
        Tuple,
    },
    _,
) end

"""
    SavingAffect{saveatType}

This struct is used by `NonInterpSavingCallback` to fill `saved_values` with
values of `p` at various timesteps. The `saveiter` field allows us to
allocate `saved_values` before the simulation and fill it during the run,
rather than pushing to an initially empty structure.
"""
mutable struct SavingAffect{saveatType}
    saved_values::NamedTuple
    saveat::saveatType
    saveiter::Int
end

"""
    condition(saveat)

This function returns a function with the type signature expected by
`SciMLBase.DiscreteCallback`, and determines whether `affect!` gets
called in the callback. This implementation simply checks if the current time
is contained in the list of save times used for the callback.
"""
condition(saveat) = (_, t, _) -> t in saveat

"""
    (affect!::SavingAffect)(integrator)

This function is used by `NonInterpSavingCallback` to perform the saving.
"""
function (affect!::SavingAffect)(integrator)
    while !isempty(affect!.saveat) && first(affect!.saveat) <= integrator.t
        affect!.saveiter += 1
        curr_t = popfirst!(affect!.saveat)

        # @assert curr_t == integrator.t
        if curr_t == integrator.t
            affect!.saved_values.t[affect!.saveiter] = curr_t
            affect!.saved_values.saveval[affect!.saveiter] =
                deepcopy(integrator.p)
        end
    end
end

"""
    saving_initialize(cb, u, t, integrator)

This function saves t and p at the start of the simulation, as long as the
initial time is in `saveat`. To run the simulation without saving these
initial values, don't pass the `initialize` argument to the `DiscreteCallback`
constructor.
"""
function saving_initialize(cb, u, t, integrator)
    (t in cb.affect!.saveat) && cb.affect!(integrator)
end

"""
    NonInterpSavingCallback(saved_values, saveat::Vector{FT})

Constructs a DiscreteCallback which saves the time and cache `p` at each time
specified by `saveat`. The first argument must be a named
tuple containing `t` and `saveval`, each having the same length as `saveat`.

Important: The times in `saveat` must be times the simulation is
evaluated at for this function to work.

Note that unlike SciMLBase's SavingCallback, this version does not
interpolate if a time in saveat is not a multiple of our timestep. This
function also doesn't work with adaptive timestepping.
"""
function NonInterpSavingCallback(saved_values, saveat::Vector{FT}) where {FT}
    # This assumes that saveat contains multiples of the timestep
    cond = condition(saveat)
    saveiter = 0
    affect! = SavingAffect(saved_values, saveat, saveiter)

    SciMLBase.DiscreteCallback(
        cond,
        affect!;
        initialize = saving_initialize,
        save_positions = (false, false),
    )
end
