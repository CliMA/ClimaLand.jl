module Simulations
using ClimaTimeSteppers
using ClimaComms
import ClimaComms: context, device
using SciMLBase
using Dates
import ClimaUtilities.TimeManager: ITime, date
import ClimaDiagnostics
using ClimaLand
import ClimaLand.Domains: SphericalShell, HybridBox, SphericalSurface
export step!, solve!, LandSimulation
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
    user_callbacks::UC
    diagnostics::DI
    required_callbacks::RC
    callbacks::CA
    _integrator::I
end
function LandSimulation(
    t0::AbstractFloat,
    tf::AbstractFloat,
    Δt::AbstractFloat,
    args...;
    start_date = nothing,
    kwargs...,
)
    t0_itime, tf_itime, Δt_itime =
        promote(ITime.((t0, tf, Δt); epoch = start_date)...)
    if haskey(kwargs, :solver_kwargs) && haskey(kwargs[:solver_kwargs], :saveat)
        _, saveat_itime... =
            promote(t0_itime, ITime.(kwargs[:solver_kwargs][:saveat])...)
        kwargs = merge(
            (; kwargs...),
            (;
                solver_kwargs = merge(
                    kwargs[:solver_kwargs],
                    (; saveat = [saveat_itime...]),
                ),
            ),
        )
    end
    if haskey(kwargs, :updateat)
        _, updateat_itime... = promote(t0_itime, ITime.(kwargs[:updateat])...)
        kwargs = merge((; kwargs...), (; updateat = [updateat_itime...]))
    end
    LandSimulation(t0_itime, tf_itime, Δt_itime, args...; kwargs...)
end

function LandSimulation(
    start_date::Dates.DateTime,
    stop_date::Dates.DateTime,
    Δt::AbstractFloat,
    args...;
    kwargs...,
)
    t0 = ITime(0, Dates.Second(1), start_date)
    tf = ITime(
        Dates.value(convert(Dates.Second, stop_date - start_date)),
        epoch = start_date,
    )
    Δt = ITime(Δt, epoch = start_date)
    t0_itime, tf_itime, Δt_itime = promote(t0, tf, Δt)
    if haskey(kwargs, :solver_kwargs) && haskey(kwargs[:solver_kwargs], :saveat)
        saveat_converted = map(
            d -> ITime(
                Dates.value(convert(Dates.Second, d - start_date));
                epoch = start_date,
            ),
            kwargs[:solver_kwargs][:saveat],
        )
        _, saveat_itime... = promote(t0_itime, saveat_converted...)
        kwargs = merge(
            (; kwargs...),
            (;
                solver_kwargs = merge(
                    kwargs[:solver_kwargs],
                    (; saveat = [saveat_itime...]),
                ),
            ),
        )
    end
    if haskey(kwargs, :updateat)
        updateat_converted = map(
            d -> ITime(
                Dates.value(convert(Dates.Second, d - start_date));
                epoch = start_date,
            ),
            kwargs[:updateat],
        )
        _, updateat_itime... = promote(t0_itime, updateat_converted...)
        kwargs = merge((; kwargs...), (; updateat = [updateat_itime...]))
    end
    LandSimulation(t0_itime, tf_itime, Δt_itime, args...; kwargs...)
end
# TODO: Add doc string
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
            isnothing(t0.epoch) ? t0 * 10000 : Dates.Month(1),
            t0,
            Δt;
            mask = ClimaLand.Domains.landsea_mask(ClimaLand.get_domain(model)),
        ),
        ClimaLand.ReportCallback(1000),
    ),
    diagnostics = ClimaLand.default_diagnostics(model, t0, outdir),
    updateat = [promote(t0:(ITime(3600 * 3)):tf...)...],
    solver_kwargs = (;),
)
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
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    required_callbacks = (driver_cb,) # TBD: can we update each step?

    diagnostics = isnothing(diagnostics) ? () : diagnostics
    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(diagnostics, Y, p, t0; dt = Δt)
    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)


    # Collect all callbacks #TODO: ordering can be confusing as the state can be saved
    # in both user_cbs and diag_cbs, and the driver update happens between them
    callbacks =
        SciMLBase.CallbackSet(user_callbacks..., required_callbacks..., diag_cb)
    _integrator = SciMLBase.init(
        problem,
        timestepper;
        dt = Δt,
        callback = callbacks,
        adaptive = false,
        solver_kwargs...,
    )
    return LandSimulation(
        model,
        timestepper,
        user_callbacks,
        diagnostics,
        required_callbacks,
        callbacks,
        _integrator,
    )
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
end#module
