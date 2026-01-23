import ClimaCore
import SciMLBase
import ClimaDiagnostics.Schedules:
    EveryCalendarDtSchedule, EveryDtSchedule, DivisorSchedule
import Dates

using ClimaComms
using ClimaUtilities.ClimaArtifacts
import Interpolations
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.OnlineLogging: WallTimeInfo, report_walltime
import ClimaUtilities.TimeManager: ITime, date

export FTfromY

"""
    partial_q_sat_partial(q_sat::FT, T::FT, phase, earth_param_set) where {FT}

Computes the quantity ∂q_sat∂T at temperature T assuming q_sat has already been 
computed from T and ρ over the correct phase `phase`.. 
"""
function partial_q_sat_partial_T(
    q_sat::FT,
    T::FT,
    phase,
    earth_param_set,
) where {FT}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    L = ifelse(
        phase isa Thermodynamics.Liquid,
        Thermodynamics.latent_heat_vapor(thermo_params, T),
        Thermodynamics.latent_heat_sublim(thermo_params, T),
    )
    return Thermodynamics.∂q_vap_sat_∂T(thermo_params, nothing, T, q_sat, L)
end

"""
    check_land_equality(field1, field2)

Check that two fields are equal for masked domains;
this assumes that the underlying space of field1 is the same
as for field2, and this specific test is only required if
the two fields have masked spaces, where the masked values
are uninitialized. In this case, the masked values may be
anything and a standard equality field1 == field2 will fail.

This allocates two fields equal in size to field1/field2.
"""
function check_land_equality(
    field1::ClimaCore.Fields.Field,
    field2::ClimaCore.Fields.Field,
)
    tmp1 = ClimaCore.Fields.zeros(axes(field1)) # inherits the mask
    tmp2 = ClimaCore.Fields.zeros(axes(field2)) # inherits the mask
    FT = eltype(tmp1)
    fill!(ClimaCore.Fields.field_values(tmp1), FT(0)) # fills all values to 0
    fill!(ClimaCore.Fields.field_values(tmp2), FT(0)) # fills all values to 0
    tmp1 .= field1 # only sets unmasked values/mask aware operation
    tmp2 .= field2 # only sets unmasked values/mask aware operation
    @assert tmp1 == tmp2 # checks all values
end

"""
    check_land_equality(f1::Union{FT, Vector{FT}}, f2::Union{FT, Vector{FT}})
        where {FT <: AbstractFloat}

Check that two floats or vectors of floats are equal.
"""
function check_land_equality(
    f1::Union{FT, Vector{FT}},
    f2::Union{FT, Vector{FT}},
) where {FT <: AbstractFloat}
    @assert f1 == f2 # checks all values
end


"""
     heaviside(x::FT)::FT where {FT}

Computes the heaviside function H(x) = 1 (x ≥0), 0 (x < 0)
"""
function heaviside(x::FT)::FT where {FT}
    heaviside(x, FT(0))
end

"""
     heaviside(x::FT, a::FT)::FT where {FT}

Computes the heaviside function H(y) = 1 (y ≥0), 0 (y < 0),
where y = x-a.
"""
function heaviside(x::FT, a::FT)::FT where {FT}
    if x - a > eps(FT)
        return FT(1.0)
    else
        return FT(0.0)
    end
end

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
    # With npolynomial = 0, we don't need DSS (and DSS will fail with MPI)
    if domain.npolynomial == 0
        return p
    else
        buffer = ClimaCore.Spaces.create_dss_buffer(
            ClimaCore.Fields.zeros(domain.space.surface),
        )
        return merge(p, (; dss_buffer_2d = buffer))
    end
end

"""
    add_dss_buffer_to_aux(
        p::NamedTuple,
        domain::Union{Domains.HybridBox, Domains.SphericalShell},
    )

Adds a 2d and 3d dss buffer corresponding to `domain.space` to `p` with the names
`dss_buffer_3d`, and `dss_buffer_2d`.

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
    # With npolynomial = 0, we don't need DSS (and DSS will fail with MPI)
    if domain.npolynomial == 0
        return p
    else
        buffer_2d = ClimaCore.Spaces.create_dss_buffer(
            ClimaCore.Fields.zeros(domain.space.surface),
        )
        buffer_3d = ClimaCore.Spaces.create_dss_buffer(
            ClimaCore.Fields.zeros(domain.space.subsurface),
        )
        return merge(
            p,
            (; dss_buffer_3d = buffer_3d, dss_buffer_2d = buffer_2d),
        )
    end
end


"""
      dss!(Y::ClimaCore.Fields.FieldVector, p::NamedTuple, t)

Computes the weighted direct stiffness summation and updates `Y` in place.
In the case of a column domain, no dss operations are performed.
"""
function dss!(Y::ClimaCore.Fields.FieldVector, p::NamedTuple, t)
    #for key in propertynames(Y)
    #    property = getproperty(Y, key)
    #    dss_helper!(property, axes(property), p)
    #end
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
    # The buffer is set up expecting a single scalar at each point in
    # the space. if the field element type is a scalar, we are ready to
    # compute the dss. However, if the field is Tuple-valued,
    # we need to carry out the dss for each element of the tuple (each
    # of which is equivalent to a scalar valued field)
    field_element_type = eltype(field)
    if field_element_type <: AbstractFloat
        ClimaCore.Spaces.weighted_dss!(field, buffer)
    elseif field_element_type <: Tuple
        n = length(field_element_type.types)
        for i in 1:n
            ClimaCore.Spaces.weighted_dss!(field.:($i), buffer)
        end
    else
        @error("No DSS method defined for your field type.")
    end

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
    DriverUpdateCallback(
        updatefunc,
        update_period,
        t0;
        dt = nothing,
    )

Constructs a DiscreteCallback which updates the cache `p.drivers` at each time
specified by `updateat`, using the function `updatefunc` which takes as arguments (p,t).
If `update_period` is `nothing`, then `updatefunc` will not be called during the simulation
or during intialization. This can be used to disable the driver update for a simulation.
"""
function DriverUpdateCallback(updatefunc, update_period, t0; dt = nothing)
    affect! = (integrator) -> updatefunc(integrator.p, integrator.t)
    # when update_period is nothing, disable driver updates
    disable_update = isnothing(update_period)
    initialize =
        disable_update ? (_, _, _, _) -> nothing :
        (cb, u, t, integrator) -> affect!(integrator)
    update_period = disable_update ? zero(t0) : update_period


    return IntervalBasedCallback(update_period, t0, dt, affect!; initialize)
end

"""
    CheckpointCallback(
        checkpoint_period::Union{Number, Dates.Period, ITime},
        output_dir,
        t0;
        model,
        dt = nothing,
    )

Constructs a DiscreteCallback which saves the state to disk with the
`save_checkpoint` function.

# Arguments
- `checkpoint_period`: The interval between times where checkpoints are saved. Can be
  specified as a float (in seconds) `Dates.Period`, or `ITime`.
- `output_dir`: The directory where the checkpoint files will be saved.
- `t0`: The start of the simulation.
- `model`: The ClimaLand model object.
- `dt`: The timestep of the model (optional), used to check for consistency.

The callback uses `ClimaDiagnostics.EveryCalendarDtSchedule` to determine when
to save checkpoints based on the `checkpoint_period`. The schedule is
initialized with the `t0` to ensure that the first
checkpoint is saved at the correct time.

The `save_checkpoint` function is called with the current state vector `u`, the
current time `t`, and the `output_dir` to save the checkpoint to disk.
"""
function CheckpointCallback(
    checkpoint_period::Union{Number, Dates.Period, ITime},
    output_dir,
    t0;
    model,
    dt = nothing,
)
    affect! = let output_dir = output_dir, model = model
        (integrator) ->
            save_checkpoint(integrator.u, integrator.t, output_dir; model)
    end
    return IntervalBasedCallback(checkpoint_period, t0, dt, affect!)
end



"""
    SavingAffect{NT}

This struct is used by `NonInterpSavingCallback` to fill `saved_values` with
values of `p` at various timesteps. The `saveiter` field allows us to
allocate `saved_values` before the simulation and fill it during the run,
rather than pushing to an initially empty structure.
"""
mutable struct SavingAffect{NT <: NamedTuple}
    saved_values::NT
    saveiter::Int
end

"""
    (affect!::SavingAffect)(integrator)

This function is used by `NonInterpSavingCallback` to perform the saving.
"""
function (affect!::SavingAffect)(integrator)
    T_saved_t = Base.typesplit(eltype(affect!.saved_values.t), Nothing)
    affect!.saveiter += 1
    if integrator.t isa T_saved_t
        # needs a special case because ITime(::ITime) is not defined
        affect!.saved_values.t[affect!.saveiter] = integrator.t
    else
        affect!.saved_values.t[affect!.saveiter] = T_saved_t(integrator.t)
    end
    affect!.saved_values.saveval[affect!.saveiter] = deepcopy(integrator.p)
end


"""
    NonInterpSavingCallback(
        start_date::Dates.DateTime,
        stop_date::Dates.DateTime,
        callback_period::Dates.Period;
        first_save_date = start_date,
    )

Constructs a `DiscreteCallback` which saves the time and cache `p` at `first_save_date` and
every `callback_period` after. Times are saved as DateTimes, and the saved times and caches
are each saved in a vector. They are stored in a `NamedTuple`, which is the `saved_values`
property of the `affect!` property of the `DiscreteCallback`. For example:

`cb = NonInterpSavingCallback(start_date, stop_date, callback_period)`
`saved_values = cb.affect!.saved_values`
`saved_times = saved_values.t` - Vector of DateTimes
`saved_cache = saved_values.saveval` - A vector of caches

Note: This method constructs a callback that uses `ITime` for scheduling, and should only
be used with simulations that are run with `ITime`.
"""
function NonInterpSavingCallback(
    start_date::Dates.DateTime,
    stop_date::Dates.DateTime,
    callback_period::Dates.Period;
    first_save_date = start_date,
)
    first_save_during_init = first_save_date == start_date
    callback_start_date =
        first_save_during_init ? first_save_date :
        first_save_date - callback_period
    n_saves = length(first_save_date:callback_period:stop_date)
    sv = (;
        t = Vector{Union{Dates.DateTime, Nothing}}(nothing, n_saves),
        saveval = Vector{Any}(undef, n_saves),
    )
    affect! = SavingAffect(sv, 0)
    initialize =
        first_save_during_init ? (_, _, _, x) -> affect!(x) :
        (_, _, _, _) -> nothing
    t0 = ITime(0, epoch = start_date)
    callback_start = ITime(
        Dates.value(convert(Dates.Second, callback_start_date - start_date)),
        epoch = start_date,
    )
    return IntervalBasedCallback(
        callback_period,
        t0,
        nothing,
        affect!;
        initialize,
        callback_start,
    )
end

"""
    NonInterpSavingCallback(t0, tf, callback_period; t_first_save = t0)

Constructs a `DiscreteCallback` which saves the time and cache `p` at `t_first_save`
and every `callback_period` after. Times are saved as the same type as `t0`,
and the saved times and caches are each saved in a vector. They are stored in a `NamedTuple`,
which is the `saved_values` property of the `affect!` property of the `DiscreteCallback`.
For example:

`cb = NonInterpSavingCallback(t0, tf, callback_period)`
`saved_values = cb.affect!.saved_values`
`saved_times = saved_values.t` - Vector of typeof(t0)
`saved_cache = saved_values.saveval` - A vector of caches

Note: This method constructs a callback that uses `ITime` for scheduling, and should only
be used with simulations that are run with `ITime`.
"""
function NonInterpSavingCallback(t0, tf, callback_period; t_first_save = t0)
    t0_itime = ITime(t0)
    tf_itime = ITime(tf)
    callback_period_itime = ITime(callback_period)
    sv_itime = ITime(t_first_save)
    t0_itime, tf_itime, callback_period_itime, sv_itime =
        promote(t0_itime, tf_itime, callback_period_itime, sv_itime)

    first_save_during_init = t_first_save == t0
    callback_start =
        first_save_during_init ? sv_itime : sv_itime - callback_period_itime
    n_saves = length(sv_itime:callback_period_itime:tf_itime)
    sv = (;
        t = Vector{Union{typeof(t0), Nothing}}(nothing, n_saves),
        saveval = Vector{Any}(undef, n_saves),
    )
    affect! = SavingAffect(sv, 0)
    initialize =
        first_save_during_init ? (_, _, _, x) -> affect!(x) :
        (_, _, _, _) -> nothing
    return IntervalBasedCallback(
        callback_period_itime,
        t0_itime,
        nothing,
        affect!;
        initialize,
        callback_start,
    )
end


function FTfromY(Y::ClimaCore.Fields.FieldVector)
    return eltype(Y)
end


"""
    isdivisible(dt_large::Dates.Period, dt_small::Dates.Period)

Check if two periods are evenly divisible, i.e., if the larger period can be
expressed as an integer multiple of the smaller period.

In this, take into account the case when periods do not have fixed size, e.g.,
one month is a variable number of days.

# Examples
```
julia> isdivisible(Dates.Year(1), Dates.Month(1))
true

julia> isdivisible(Dates.Month(1), Dates.Day(1))
true

julia> isdivisible(Dates.Month(1), Dates.Week(1))
false
```

## Notes

Not all the combinations are fully implemented. If something is missing, please
consider adding it.
"""
function isdivisible(dt_large::Dates.Period, dt_small::Dates.Period)
    @warn "The combination $(typeof(dt_large)) and $(dt_small) was not covered. Please add a method to handle this case."
    return false
end

# For FixedPeriod and OtherPeriod, it is easy, we can directly divide the two
# (as long as they are both the same)
function isdivisible(dt_large::Dates.FixedPeriod, dt_small::Dates.FixedPeriod)
    return isinteger(dt_large / dt_small)
end

function isdivisible(dt_large::Dates.OtherPeriod, dt_small::Dates.OtherPeriod)
    return isinteger(dt_large / dt_small)
end

function isdivisible(dt_large::ITime, dt_small::ITime)
    return isinteger(dt_large / dt_small)
end

function isdivisible(dt_large::Number, dt_small::Number)
    return isinteger(dt_large / dt_small)
end

"""
    call_count_nans_state(state, mask::ClimaCore.Fields.Field = nothing, verbose = false)

This calls the function which counts the number of NaNs in the state, e.g. the FieldVector given by
`sol.u[end]` after calling `solve`. This function is useful for
debugging simulations to determine quantitatively if a simulation is stable; this is only
for use in non-column simulations.

If this function is called on a FieldVector, it will recursively call itself
on each Field in the FieldVector. If it is called on a Field, it will count
the number of NaNs in the Field and produce a warning if any are found; this behavior
is different for 2d and 3d fields, which is why we dispatch on the `space`.

If a ClimaCore Field is provided as `mask`, the function will only count NaNs
in the state variables where the mask is 1. This is intended to be used with
the land/sea mask, to avoid counting NaNs over the ocean. Note this assumes
the mask is 1 over land and 0 over ocean, and that the mask is a 2D field.

The `verbose` argument toggles whether the function produces output when no
NaNs are found.
"""
function call_count_nans_state(
    state::ClimaCore.Fields.FieldVector;
    mask = nothing,
    verbose = false,
)
    for pn in propertynames(state)
        state_new = getproperty(state, pn)
        @info "Checking NaNs in $pn"
        call_count_nans_state(state_new; mask, verbose)
    end
    return nothing
end

function call_count_nans_state(
    state::ClimaCore.Fields.Field;
    mask = nothing,
    verbose = false,
)
    call_count_nans_state(state, axes(state); mask, verbose)
end

function call_count_nans_state(
    state::ClimaCore.Fields.Field,
    space::Union{
        ClimaCore.Spaces.AbstractSpectralElementSpace,
        ClimaCore.Spaces.AbstractPointSpace,
    };
    mask = nothing,
    verbose = false,
)
    return count_nans_state(state; mask = mask, verbose = verbose)
end

function call_count_nans_state(
    state::ClimaCore.Fields.Field,
    space::Union{
        ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace,
        ClimaCore.Spaces.FiniteDifferenceSpace,
    };
    mask = nothing,
    verbose = false,
)
    return count_nans_state(
        ClimaLand.Domains.top_center_to_surface(state);
        mask = mask,
        verbose = verbose,
    )
end

"""
    count_nans_state(
        state::ClimaCore.Fields.Field,
        mask = nothing,
        verbose = false,
    )

Counts the NaNs in the field `state`.
"""
function count_nans_state(
    state::ClimaCore.Fields.Field;
    mask = nothing,
    verbose = false,
)
    if isnothing(mask)
        num_nans = count(isnan, parent(state))
    else
        num_nans = mapreduce(
            (s, m) -> m != 0 && isnan(s),
            Base.add_sum,
            parent(state),
            parent(mask),
        )
    end
    if isapprox(num_nans, 0)
        verbose && @info "No NaNs found"
    else
        @warn "$num_nans NaNs found"
    end
    return nothing
end


"""
    NaNCheckCallback(
        nancheck_period::Union{Number, Dates.Period, ITime},
        t0;
        dt = nothing,
        mask = nothing,
    )

Constructs a DiscreteCallback which counts the number of NaNs in the state
and produces a warning if any are found.

# Arguments
- `nancheck_period`: The interrval between times when the state is checked for NaNs.
  Can be specified as an `ITime`, a `Number`, or a `Dates.Period`.
- `t0`: The start of the simulation.
- `dt`: The timestep of the model (optional), used to check for consistency.
- `mask`: NaNs will not be counted in areas where `mask` is zero

The callback uses `ClimaDiagnostics` schedules to determine when
to check for NaNs based on the `nancheck_period`.
"""
function NaNCheckCallback(
    nancheck_period::Union{Number, Dates.Period, ITime},
    t0;
    dt = nothing,
    mask = nothing,
)

    affect! = (integrator) -> call_count_nans_state(integrator.u; mask)

    return IntervalBasedCallback(nancheck_period, t0, dt, affect!)
end

"""
    IntervalBasedCallback(
        period::Dates.Period,
        t0::ITime,
        dt::Union{ITime, Nothing},
        affect!;
        initialize = (_, _, _, _) -> nothing,
        callback_start::ITime = t0,
    )

Returns a SciML DiscreteCallback that calls `affect!` every `period`. This method is used
when the type of `period` is a `Dates.Period`

The returned callback has a condition function that is built on a ClimaDiagnostics
`EveryCalendarDtSchedule`. The `initialize` argument is passed
to the DiscreteCallback as keyword arguments. Excluding initialization, the first time
the callback will be called is when the simulation is at or past `period` + `callback_start`.
"""
function IntervalBasedCallback(
    period::Dates.Period,
    t0::ITime,
    dt::Union{ITime, Nothing},
    affect!;
    initialize = (_, _, _, _) -> nothing,
    callback_start::ITime = t0,
)
    schedule = EveryCalendarDtSchedule(
        period;
        start_date = t0,
        date_last = callback_start,
    )

    if !isnothing(dt)
        dt_period = Dates.Millisecond(1000float(dt))
        if period isa Dates.FixedPeriod && !isdivisible(period, dt_period)
            @warn "Callback period ($period) is not an integer multiple of dt $dt_period"
        end
    end

    cond = let schedule = schedule
        (_, _, integrator) -> schedule(integrator)
    end
    return SciMLBase.DiscreteCallback(cond, affect!; initialize)
end

"""
    IntervalBasedCallback(
        period::Union{Number, ITime},
        t0::Union{Number, ITime},
        dt::Union{Number, ITime, Nothing},
        affect!;
        initialize = (_, _, _, _) -> nothing,
        callback_start::Union{Number, ITime} = t0,
    )
Returns a SciML DiscreteCallback that calls `affect!` every `period`. This method is
used when `period` is not a `Dates.Period`.

The returned callback has a condition function that is built on a ClimaDiagnostics
`EveryDtSchedule`. The `initialize` argument is passed
to the DiscreteCallback as keyword arguments. Excluding initialization, the first time
the callback will be called is when the simulation is at or past `period` + `callback_start`.
"""
function IntervalBasedCallback(
    period::Union{Number, ITime},
    t0::Union{Number, ITime},
    dt::Union{Number, ITime, Nothing},
    affect!;
    initialize = (_, _, _, _) -> nothing,
    callback_start::Union{Number, ITime} = t0,
)
    if Any <: typejoin(typeof(period), typeof(callback_start)) # checks if they can be promoted to same type
        error("Unsupported combination of `period` and `callback_start` types")
    end
    period, callback_start = promote(period, callback_start)
    schedule = EveryDtSchedule(period; t_last = callback_start)

    if !isnothing(dt)
        if !isdivisible(period, dt)
            @warn "Callback period ($(period)) is not an integer multiple of dt $(dt)"
        end
    end

    cond = let schedule = schedule
        (u, t, integrator) -> schedule(integrator)
    end
    return SciMLBase.DiscreteCallback(cond, affect!; initialize)
end

"""
    ReportCallback(period, t0; dt = nothing)

Return a callback that prints performance and progress summaries every `period`
"""
function ReportCallback(period, t0; dt = nothing)
    walltime_info = WallTimeInfo()
    report = let wt = walltime_info
        (integrator) -> report_walltime(wt, integrator)
    end
    report_cb = IntervalBasedCallback(period, t0, dt, report)
    return report_cb
end
