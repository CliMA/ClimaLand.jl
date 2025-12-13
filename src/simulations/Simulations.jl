module Simulations
using ClimaTimeSteppers
using ClimaComms
import ClimaComms: context, device
using SciMLBase
using Dates
import ClimaUtilities.TimeManager: ITime, date
import ClimaDiagnostics
using ClimaLand
import ..Parameters as LP
export step!, solve!, LandSimulation
using NCDatasets
include("initial_conditions.jl")

import ..Diagnostics: close_output_writers

"""
    LandSimulation{
        M <: ClimaLand.AbstractModel,
        T <: ClimaTimeSteppers.DistributedODEAlgorithm,
        UC,
        DI,
        RC,
        CA <: SciMLBase.CallbackSet,
        I <: SciMLBase.DEIntegrator,
    }

the ClimaLand LandSimulation struct, which specifies
- the discrete set of equations to solve (defined by the `model`);
- the timestepping algorithm;
- user callbacks (passed as a tuple) to be executed at specific times in the simulations;
- the diagnostics to output (optional).

User callbacks are optional: examples currently include callbacks that estimate the time
to solution and SYPD of the simulation as it runs, checkpoint the state, or check the solution
for NaNs. Others can be added here.

Diagnostics are implemented as callbacks, and are also optional.
However, a default is provided. `diagnostics` is expected to be a
list of `ClimaDiagnostics.ScheduledDiagnostics`.

Finally, the private field _required_callbacks consists of callbacks that are required for the
simulation to run correctly. Currently, this includes the callbacks which update the atmospheric
forcing and update the LAI using prescribed data.
"""
struct LandSimulation{
    M <: ClimaLand.AbstractModel,
    T <: ClimaTimeSteppers.DistributedODEAlgorithm,
    UC,
    DI,
    RC,
    CA <: SciMLBase.CallbackSet,
    I <: SciMLBase.DEIntegrator,
}
    model::M
    timestepper::T
    start_date::Union{DateTime, Nothing}
    user_callbacks::UC
    diagnostics::DI
    required_callbacks::RC
    callbacks::CA
    _integrator::I
end

"""
    LandSimulation(
        t0::ITime,
        tf::ITime,
        Δt::ITime,
        model;
        outdir = ".",
        set_ic! = make_set_initial_state_from_file(
            ClimaLand.Artifacts.soil_ic_2008_50m_path(;
                context = ClimaComms.context(model),
            ),
            model,
        ),
        timestepper = ClimaTimeSteppers.IMEXAlgorithm(
            ClimaTimeSteppers.ARS111(),
            ClimaTimeSteppers.NewtonsMethod(
                max_iters = 3,
                update_j = ClimaTimeSteppers.UpdateEvery(
                    ClimaTimeSteppers.NewNewtonIteration,
                ),
            ),
        ),
        user_callbacks = (
            ClimaLand.NaNCheckCallback(
                isnothing(t0.epoch) ? div((tf - t0), 10) : Dates.Month(1),
                t0;
                dt = Δt,
                mask = ClimaLand.Domains.landsea_mask(ClimaLand.get_domain(model)),
            ),
            ClimaLand.ReportCallback(div((tf - t0), 10), t0),
        ),
        diagnostics = ClimaLand.default_diagnostics(model, t0, outdir),
        updateat = ITime(3600 * 3),
        solver_kwargs = (;),
    )

Creates a `LandSimulation` object for `model` with the  `_integrator` field initialized
with a problem from `t0` to `tf`, a time step of `Δt`, and with `timestepper` as the
time-stepping algorithm. During the construction process, the initial conditions
are set using `set_ic!`, the initial cache is set using `make_set_initial_cache(model)`,
and driver and diagnostic callbacks are created. If `updateat` is a vector, the drivers
will be updated at those times; if it is a single time, the drivers will be updated
at that interval. The order the callbacks are called while stepping the simulation is:
1. user_callbacks
2. driver update callback
3. diagnostics callback
"""
function LandSimulation(
    t0::ITime,
    tf::ITime,
    Δt::ITime,
    model;
    outdir = ".",
    set_ic! = make_set_initial_state_from_file(
        ClimaLand.Artifacts.soil_ic_2008_50m_path(;
            context = ClimaComms.context(model),
        ),
        model,
    ),
    timestepper = ClimaTimeSteppers.IMEXAlgorithm(
        ClimaTimeSteppers.ARS111(),
        ClimaTimeSteppers.NewtonsMethod(
            max_iters = 3,
            update_j = ClimaTimeSteppers.UpdateEvery(
                ClimaTimeSteppers.NewNewtonIteration,
            ),
        ),
    ),
    user_callbacks = (
        ClimaLand.NaNCheckCallback(
            isnothing(t0.epoch) ? div((tf - t0), 10) : Dates.Month(1),
            t0;
            dt = Δt,
            mask = ClimaLand.Domains.landsea_mask(ClimaLand.get_domain(model)),
        ),
        ClimaLand.ReportCallback(div((tf - t0), 10), t0),
    ),
    diagnostics = ClimaLand.default_diagnostics(model, t0, outdir),
    updateat = ITime(Δt),
    solver_kwargs = (;),
)
    # Enforce `h_atmos >= h_canopy` for models with canopy
    if hasproperty(model, :canopy)
        h_atmos_min = first(extrema(model.canopy.boundary_conditions.atmos.h))
        h_canopy = model.canopy.hydraulics.compartment_surfaces[end]
        @assert h_atmos_min >= h_canopy "Atmospheric height must be greater than or equal to canopy height. Got min atmos height $h_atmos_min and canopy height $h_canopy"
    elseif model isa ClimaLand.Canopy.CanopyModel
        h_atmos_min = first(extrema(model.boundary_conditions.atmos.h))
        h_canopy = model.hydraulics.compartment_surfaces[end]
        @assert h_atmos_min >= h_canopy "Atmospheric height must be greater than or equal to canopy height. Got min atmos height $h_atmos_min and canopy height $h_canopy"
    end


    if !isnothing(diagnostics) &&
       !isempty(diagnostics) &&
       !(
           first(diagnostics).output_writer isa
           ClimaDiagnostics.Writers.DictWriter
       ) &&
       first(diagnostics).output_writer.output_dir != outdir
        @warn "Note that the kwarg outdir and outdir used in diagnostics are inconsistent; using $(first(diagnostics).output_writer.output_dir)"
    end

    # set initial conditions
    Y, p, cds = initialize(model)
    set_ic!(Y, p, t0, model)
    set_initial_cache! = make_set_initial_cache(model)
    set_initial_cache!(p, Y, t0)

    # Create tendencies and jacobian update function
    exp_tendency! = make_exp_tendency(model)
    if model isa ClimaLand.AbstractExpModel
        T_imp! = nothing
    else
        imp_tendency! = ClimaLand.make_imp_tendency(model)
        jacobian! = ClimaLand.make_jacobian(model)
        jac_kwargs = (;
            jac_prototype = ClimaLand.FieldMatrixWithSolver(Y),
            Wfact = jacobian!,
        )
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...)
    end

    # Create SciML ODE Problem
    problem = SciMLBase.ODEProblem(
        ClimaTimeSteppers.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = T_imp!,
            dss! = ClimaLand.dss!,
        ),
        Y,
        (t0, tf),
        p,
    )

    # Required callbacks
    drivers = ClimaLand.get_drivers(model)
    updatefunc = ClimaLand.make_update_drivers(drivers)
    driver_cb =
        ClimaLand.DriverUpdateCallback(updatefunc, updateat, t0; dt = Δt)
    model_callbacks = ClimaLand.get_model_callbacks(model; t0, Δt)# everything else you need should be in the model!

    required_callbacks = (driver_cb, model_callbacks...) # TBD: can we update each step?

    diagnostics = isnothing(diagnostics) ? () : diagnostics
    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(diagnostics, Y, p, t0; dt = Δt)
    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)


    # Collect all callbacks #TODO: ordering can be confusing as the state can be saved
    # in both user_cbs and diag_cbs, and the driver update happens between them
    callbacks =
        SciMLBase.CallbackSet(user_callbacks..., required_callbacks..., diag_cb)
    if haskey(solver_kwargs, :saveat) && !(solver_kwargs[:saveat] isa Array)
        solver_kwargs = merge(
            (; solver_kwargs...),
            (; :saveat => collect(t0:solver_kwargs[:saveat]:tf)),
        )
    end
    _integrator = SciMLBase.init(
        problem,
        timestepper;
        dt = Δt,
        callback = callbacks,
        merge((; adaptive = false), solver_kwargs)...,
    )
    return LandSimulation(
        model,
        timestepper,
        t0.epoch,
        user_callbacks,
        diagnostics,
        required_callbacks,
        callbacks,
        _integrator,
    )
end

"""
    LandSimulation(
        t0::AbstractFloat,
        tf::AbstractFloat,
        Δt::AbstractFloat,
        args...;
        kwargs...,
    )

A convenience constructor for `LandSimulation` that converts `t0`, `tf`, `Δt` into `ITime`(s), setting the epoch of the ITime to nothing.
If the `kwargs` contain `updateat`, it will
convert the update frequency to an `ITime` as well. The same applies to values in `saveat` in `solver_kwargs`.

If your simulation has a notion of a calendar date, please use one of the
other constructors.
"""
function LandSimulation(
    t0::AbstractFloat,
    tf::AbstractFloat,
    Δt::AbstractFloat,
    args...;
    kwargs...,
)
    t0_itime, tf_itime, Δt_itime = promote(ITime.((t0, tf, Δt))...)
    kwargs = convert_kwarg_updates(t0_itime, kwargs)
    LandSimulation(t0_itime, tf_itime, Δt_itime, args...; kwargs...)
end

"""
    LandSimulation(
        start_date::Dates.DateTime,
        stop_date::Dates.DateTime,
        Δt::Union{AbstractFloat, Dates.Second},
        args...;
        kwargs...,
    )

A convenience constructor for `LandSimulation` that converts `start_date`, `stop_date`, and `Δt`
into `ITime`(s), using `start_date` as the epoch. If the `kwargs` contain `updateat`, it will
convert the update times to `ITime`(s) as well. The same applies to `saveat` in `solver_kwargs`.
"""
function LandSimulation(
    start_date::Dates.DateTime,
    stop_date::Dates.DateTime,
    Δt::Union{AbstractFloat, Dates.Second},
    args...;
    kwargs...,
)
    t0 = ITime(0, Dates.Second(1), start_date)
    tf = ITime(
        Dates.value(convert(Dates.Second, stop_date - start_date)),
        epoch = start_date,
    )
    Δt = ITime(Δt isa Dates.Second ? Δt.value : Δt, epoch = start_date)
    t0_itime, tf_itime, Δt_itime = promote(t0, tf, Δt)
    kwargs = convert_kwarg_updates(t0_itime, kwargs)
    LandSimulation(t0_itime, tf_itime, Δt_itime, args...; kwargs...)
end

"""
    step!(landsim::LandSimulation)

Advances the land simulation `landsim` forward in time by one step,
updating `landsim` in place.
"""
function step!(landsim::LandSimulation)
    SciMLBase.step!(landsim._integrator)
end

"""
    solve!(landsim::LandSimulation)

Advances the land simulation `landsim` forward from the initial to final time,
updating `landsim` in place.
"""
function solve!(landsim::LandSimulation)
    try
        SciMLBase.solve!(landsim._integrator)
    catch ret_code
        @error "ClimaLand simulation crashed. Stacktrace for failed simulation" exception =
            (ret_code, catch_backtrace())
    finally
        close_output_writers(landsim.diagnostics)
    end
end

"""
    ClimaComms.context(landsim::LandSimulation)

Returns the context of the land simulation (distributed or
on a single node).
"""
function ClimaComms.context(landsim::LandSimulation)
    return ClimaComms.context(landsim.model)
end

"""
    ClimaComms.device(landsim::LandSimulation)

Returns the device of the land simulation (indicating
if the simulation is on CPU or GPU).
"""
function ClimaComms.device(landsim::LandSimulation)
    return ClimaComms.device(landsim.model)
end

"""
    Base.show(io::IO, landsim::LandSimulation)

An extension of Base.show which prints the type of device,
model, and current date of the simulation `landsim`.
"""
function Base.show(io::IO, landsim::LandSimulation)
    device_type = nameof(typeof(ClimaComms.device(landsim)))
    model_type = nameof(typeof(landsim.model))
    t = landsim._integrator.t
    final_line = if isnothing(t.epoch)
        "└── Current simulation time: $(t)\n"
    else
        "└── Current date: $(date(t))\n"
    end
    return print(
        io,
        "$(model_type) Simulation\n",
        "├── Running on: $(device_type)\n",
        final_line,
    )
end

"""
    convert_kwarg_updates(t0::ITime, kwargs)

Converts the below kwargs passed to the LandSimulation constructor, which are
update vectors/intervals, into the same type as `t0`.
- `updateat`: a vector of update times, or an update interval
- `solver_kwargs[:updateat]`: a vector of update times, or an update interval
"""
function convert_kwarg_updates(t0::ITime, kwargs)
    if haskey(kwargs, :solver_kwargs)
        solver_kwargs =
            convert_kwarg_updates(t0, kwargs[:solver_kwargs], :saveat)
        kwargs = merge((; kwargs...), (; solver_kwargs,))
    end
    kwargs = convert_kwarg_updates(t0, kwargs, :updateat)
    return kwargs
end

function convert_kwarg_updates(t0::ITime, kwargs, key::Symbol)
    haskey(kwargs, key) || return kwargs # if the key is not present, return the original kwargs
    updates = kwargs[key]
    converted_updates = convert_updates.(t0, updates)
    if converted_updates isa Tuple
        converted_updates = collect(converted_updates)
    end
    return merge((; kwargs...), (; key => converted_updates))
end

convert_updates(t0::ITime, update_time::ITime) = promote(t0, update_time)[2]
convert_updates(t0::ITime, update_time::Dates.DateTime) = promote(
    t0,
    ITime(
        Dates.value(convert(Dates.Second, update_time - t0.epoch));
        epoch = t0.epoch,
    ),
)[2]
convert_updates(t0::ITime, update_time::Dates.ConvertiblePeriod) =
    promote(t0, ITime(Dates.value(convert(Dates.Second, update_time));))[2]
convert_updates(t0::ITime, update_time::AbstractFloat) =
    promote(t0, ITime(update_time, epoch = t0.epoch))[2]
convert_updates(t0, update_time) = t0 # fallback (used for non convertible Dates.Period(s))
end#module
