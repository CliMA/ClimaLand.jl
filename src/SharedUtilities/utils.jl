using DiffEqCallbacks

"""
     heaviside(x::FT)::FT where {FT}

Computes the heaviside function.
"""
function heaviside(x::FT)::FT where {FT}
    if x >= FT(0.0)
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
   dss!(Y::ClimaCore.Fields.FieldVector, p::ClimaCore.Fields.FieldVector, T::FT)

Computes the weighted direct stiffness summation and updates `Y` in place.
In the case of a column domain, no dss operations are performed.
"""
function dss!(Y::ClimaCore.Fields.FieldVector{FT}, _, _) where {FT}
    for key in propertynames(Y)
        property = getproperty(Y, key)
        dss_helper!(property, axes(property))
    end
end

"""
    dss_helper!(field_vec::ClimaCore.Fields.FieldVector, _)

Method of `dss_helper!` which unpacks properties of Y when on a
domain that is 2-dimensional in the horizontal.

The assumption is that Y contains FieldVectors which themselves contain either
FieldVectors or Fields, and that the final unpacked variable is a Field.
This method is invoked when the current property itself contains additional
property(ies).
"""
function dss_helper!(field_vec::ClimaCore.Fields.FieldVector, _)
    for key in propertynames(field_vec)
        property = getproperty(field_vec, key)
        dss_helper!(property, axes(property))
    end
end

"""
    dss_helper!(
        field::ClimaCore.Fields.Field,
        domain::Union{
            ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace,
            ClimaCore.Spaces.AbstractSpectralElementSpace,
        })

Method of `dss_helper!` which performs dss on fields of Y when on a
domain that is 2-dimensional in the horizontal.

The assumption is that Y contains FieldVectors which themselves contain either
FieldVectors or Fields, and that the final unpacked variable is a Field.
This method is invoked when the element cannot be unpacked further.
"""
function dss_helper!(
    field::ClimaCore.Fields.Field,
    domain::Union{
        ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace,
        ClimaCore.Spaces.AbstractSpectralElementSpace,
    },
)
    Spaces.weighted_dss!(field)
end

"""
    dss_helper!(
        field::Union{ClimaCore.Fields.Field, Vector},
        domain::Union{
            ClimaCore.Spaces.FiniteDifferenceSpace,
            ClimaCore.Spaces.PointSpace,
            Tuple,
        })

Computes the appropriate weighted direct stiffness summation based on
the domain type, updates `Y` in place.

For spaces that don't use spectral elements (FiniteDifferenceSpace, PointSpace,
etc), no dss is needed.
Model components with no prognostic variables appear in Y as empty Vectors, and
also do not need dss.
"""
function dss_helper!(
    field::Union{ClimaCore.Fields.Field, Vector},
    domain::Union{
        ClimaCore.Spaces.FiniteDifferenceSpace,
        ClimaCore.Spaces.PointSpace,
        Tuple,
    },
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
`DiffEqCallbacks.DiscreteCallback`, and determines whether `affect!` gets
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
            affect!.saved_values.saveval[affect!.saveiter] = copy(integrator.p)
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

Note that unlike DiffEqCallbacks's SavingCallback, this version does not
interpolate if a time in saveat is not a multiple of our timestep. This
function also doesn't work with adaptive timestepping.
"""
function NonInterpSavingCallback(saved_values, saveat::Vector{FT}) where {FT}
    # This assumes that saveat contains multiples of the timestep
    cond = condition(saveat)
    saveiter = 0
    affect! = SavingAffect(saved_values, saveat, saveiter)

    DiffEqCallbacks.DiscreteCallback(
        cond,
        affect!;
        initialize = saving_initialize,
        save_positions = (false, false),
    )
end
