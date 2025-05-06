module Simulations
using ClimaTimeSteppers
using SciMLBase
using Dates
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaUtilities.TimeManager: ITime, date
import ClimaDiagnostics
import ClimaParams as CP
using ClimaLand
import ClimaLand.Parameters as LP

include("domains.jl")
include("initial_conditions.jl")
include("spatial_parameters.jl")
include("model_setup.jl")

struct LandSimulation{
    C,
    P,
    M <: ClimaLand.AbstractLandModel,
    D <: ClimaLand.Domains.AbstractDomain,
    T <: ClimaTimeSteppers.IMEXAlgorithm,
    UC,
    DI,
    RC,
    CA <: SciMLBase.CallbackSet,
    P <: SciMLBase.ODEProblem,
}
    context::C
    params::P
    model::M
    domain::D
    timestepper::T
    user_callbacks::UC
    diagnostics::DI
    required_callbacks::RC
    callbacks::CA
    problem::P
end

function GlobalLandSimulation(
    FT,
    context,
    start_date,
    t0,
    tf,
    Δt;
    params = LP.LandParameters(FT),
    domain = ClimaLand.global_domain(FT; comms_ctx = context),
    model = LandModel(domain, start_date, params, FT), # add assert
    set_ic! = ClimaLand.set_ic_from_file(
        ClimaLand.Artifacts.soil_ic_2008_50m_path(; context = context),
    ), # could also support functions
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
            t0,# needs tobe ITime? promote
            Δt; # needs to be ITime? promote
            mask = ClimaLand.landsea_mask(domain),
        ),
        ClimaLand.ReportCallback(1000),
    ),
    diagnostics = (;
        output_vars = :short,
        average_period = :monthly,
        outdir = "",
    ), # need to generalize
)
    # convert times to Itime
    t0 = ITime(t0, epoch = start_date)
    tf = ITime(tf, epoch = start_date)
    Δt = ITime(Δt, epoch = start_date)
    t0, tf, Δt = promote(t0, tf, Δt)

    # set initial conditions
    Y, p, cds = initialize(model)
    set_ic!(Y, p, model, t0)
    set_initial_cache! = make_set_initial_cache(model)
    set_initial_cache!(p, Y, t0)

    # Create tendencies and jacobian update function
    exp_tendency! = make_exp_tendency(model)
    imp_tendency! = ClimaLand.make_imp_tendency(model)
    jacobian! = ClimaLand.make_jacobian(model)
    jac_kwargs = (;
        jac_prototype = ClimaLand.FieldMatrixWithSolver(Y),
        Wfact = jacobian!,
    )

    # Create SciML ODE Problem
    problem = SciMLBase.ODEProblem(
        ClimaTimeSteppers.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
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
    required_callbacks = (driver_cb,) # can we update each step?

    # Diagnostics callbacks - can be generalized in the future
    if !(diagnostics isa Nothing)
        nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
            domain.space.subsurface,
            diagnostics.outdir;
            start_date,
        )

        diags = ClimaLand.default_diagnostics(
            model,
            start_date;
            output_writer = nc_writer,
            output_vars = diagnostics.output_vars,
            average_period = diagnostics.average_period,
        )

        diagnostic_handler =
            ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)
        diag_cb = (ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler),)
    else
        diag_cb = ()
    end


    # Collect all callbacks
    callbacks = SciMLBase.CallbackSet(
        user_callbacks...,
        diag_cb...,
        required_callbacks...,
    )

    #_integrator = SciMLBase.init(problem, callbacks)
    
    return LandSimulation(
        context,
        params,
        model,
        domain,
        timestepper,
        user_callbacks,
        diagnostics,
        _required_callbacks,
        callbacks,
        problem,
        t
        Y,
        __p
        #_integrator
    )
end

function step!(landsimulation)
    integrator = SciMLBase.integrator(Y, p, t, callbacks)
    SciMLBase.step!(integrator)
    simulation.t[] = integrator.t
end

function solve!(simulation)
    SciMLBase.solve!(problem, ode_algo, callbacks)
end

    function plot(landsimulation)

        end
    
end#module
