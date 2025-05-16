module Simulations
using ClimaTimeSteppers
using ClimaComms
import ClimaComms: context, device
using SciMLBase
using Dates
import ClimaUtilities.TimeManager: ITime, date
import ClimaDiagnostics
using ClimaLand
export step!, solve!, LandSimulation
using ClimaLand.ModelSetup
include("initial_conditions.jl")

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
    updateat = [promote(t0:(ITime(3600 * 3)):tf...)...]
    drivers = ClimaLand.get_drivers(model)
    updatefunc = ClimaLand.make_update_drivers(drivers)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    required_callbacks = (driver_cb,) # TBD: can we update each step?

    diagnostics = isnothing(diagnostics) ? () : diagnostics
    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(diagnostics, Y, p, t0; dt = Δt)
    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)


    # Collect all callbacks
    callbacks =
        SciMLBase.CallbackSet(user_callbacks..., required_callbacks..., diag_cb)

    _integrator = SciMLBase.init(
        problem,
        timestepper;
        dt = Δt,
        callback = callbacks,
        adaptive = false,
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
    SciMLBase.solve!(landsim._integrator)
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
        "└── Current date: $(date(landsim._integrator.t))\n",
    )
end
end#module
