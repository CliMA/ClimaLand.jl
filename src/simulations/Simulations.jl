module Simulations
using ClimaTimeSteppers
using ClimaComms
import ClimaComms: context, device
using Dates
import ClimaUtilities.TimeManager: ITime, date
import ClimaDiagnostics
using ClimaLand
export step!, solve!, LandSimulation
using ClimaLand.ModelSetup
include("initial_conditions.jl")

"""
    LandSimulation

The ClimaLand LandSimulation struct, which specifies
- the discrete set of equations to solve (defined by the `model`);
- the timestepping algorithm (in `timestepper.algo`);
- user callbacks (passed as a tuple) to be executed at specific times in the simulations;
- the diagnostics to output (optional).

User callbacks are optional: examples currently include callbacks that estimate the time
to solution and SYPD of the simulation as it runs, checkpoint the state, or check the solution
for NaNs. Others can be added here.

`diagnostics` are provided as a list of `ClimaDiagnostics.ScheduledDiagnostics`,
with default specified by `default_diagnostics`.

The private field `_required_callbacks` consists of callbacks that are required
for the simulation to run correctly. Currently, this includes the callbacks
which update the atmospheric forcing and update the LAI using prescribed data.

Quick tips
==========

1. Accessing the state
```julia
sim.u
```
2. Accessing the current time
```julia
sim.t
```
3. Accessing the current date
```julia
import ClimaUtilities.TimeManager: date
date(sim.t)
```
4. Providing a specific new monthly diagnostics
```julia
diagnostic = ClimaLand.Diagnostics.monthly_average.(["lhf", "shf"])
```
"""
mutable struct LandSimulation{
    STATE,
    CACHE,
    TIME,
    TUP_TIME,
    M <: ClimaLand.AbstractModel,
    T,
    DI,
    RC,
    UC,
}
    u::STATE
    p::CACHE
    t::TIME
    dt::TIME
    tspan::TUP_TIME
    model::M
    timestepper::T
    diagnostic_handler::DI
    required_callbacks::RC
    user_callbacks::UC
end

function LandSimulation(
    FT,
    start_date::Dates.DateTime,
    stop_date::Dates.DateTime,
    Δt::AbstractFloat,
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
            Dates.Month(1),
            start_date,
            ITime(Δt, epoch = start_date),
            mask = ClimaLand.landsea_mask(ClimaLand.get_domain(model)),
        ),
        ClimaLand.ReportCallback(1000),
    ),
    diagnostics = ClimaLand.default_diagnostics(
        model,
        start_date;
        output_writer = ClimaDiagnostics.Writers.NetCDFWriter(
            ClimaLand.get_domain(model).space.subsurface,
            outdir;
            start_date,
        ),
    ),
)
    if !isnothing(diagnostics) &&
       !isempty(diagnostics) &&
       first(diagnostics).output_writer.output_dir != outdir
        @warn "Note that the kwarg outdir and outdir used in diagnostics are inconsistent; using $(first(diagnostics).output_writer.output_dir)"
    end

    domain = ClimaLand.get_domain(model)
    t0 = ITime(0, Dates.Second(1), start_date)
    tf = ITime(Dates.seconds(stop_date - start_date), epoch = start_date)
    Δt = ITime(Δt, epoch = start_date)
    t0, tf, Δt = promote(t0, tf, Δt)

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
        T_imp! = (; f = imp_tendency!, jac_kwargs...)
    end

    # Create SciML ODE Problem
    func = ClimaTimeSteppers.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = T_imp!,
        dss! = ClimaLand.dss!,
    )

    # Required callbacks
    updateat = [promote(t0:(ITime(3600 * 3)):tf...)...]
    drivers = ClimaLand.get_drivers(model)
    updatefunc = ClimaLand.make_update_drivers(drivers)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    _required_callbacks = (driver_cb,) # TBD: can we update each step?

    diagnostics = isnothing(diagnostics) ? () : diagnostics
    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(diagnostics, Y, p, t0; dt = Δt)

    algo = timestepper
    # u0 is used as prototype
    prob = (; u0 = Y, f = func)
    timestepper_cache = ClimaTimeSteppers.init_cache(prob, algo)
    isnothing(func.cache!) || func.cache!(Y, p, t0)
    return LandSimulation(
        Y,
        p,
        t0,
        Δt,
        (t0, tf),
        model,
        (; algo, func, cache = timestepper_cache),
        diagnostic_handler,
        _required_callbacks,
        user_callbacks,
    )
end

"""
   step!(landsim::LandSimulation)

Advances the land simulation `landsim` forward in time by one step,
updating `landsim` in place.
"""
function step!(landsim::LandSimulation)
    landsim.t += landsim.dt

    ClimaTimeSteppers.step_u!(landsim, landsim.timestepper.cache)
    for callback in landsim._required_callbacks
        callback.condition(landsim.u, landsim.t, landsim) &&
            callback.affect!(landsim)
    end
    for callback in landsim.user_callbacks
        callback.condition(landsim.u, landsim.t, landsim) &&
            callback.affect!(landsim)
    end
    ClimaDiagnostics.orchestrate_diagnostics(
        landsim,
        landsim.diagnostic_handler,
    )
end

# Compatibility with SciML
function Base.getproperty(landsim::LandSimulation, symbol::Symbol)
    if symbol === :alg
        return landsim.timestepper.algo
    elseif symbol === :step
        return landsim.t / landsim.dt
    elseif symbol === :sol
        return (;
            prob = (; f = landsim.timestepper.func, tspan = landsim.tspan)
        )
    else
        return Base.getfield(landsim, symbol)
    end
end

"""
   solve!(landsim::LandSimulation)

Advances the land simulation `landsim` forward from the initial to final time,
updating `landsim` in place.
"""
function solve!(landsim::LandSimulation)
    while landsim.t < last(landsim.tspan)
        step!(landsim)
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
    return print(
        io,
        "$(model_type) Simulation\n",
        "├── Running on: $(device_type)\n",
        "└── Current date: $(date(landsim.t))\n",
    )
end
end#module
