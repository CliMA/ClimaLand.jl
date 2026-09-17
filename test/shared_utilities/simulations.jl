using Test
using Dates
using ClimaLand
using ClimaLand.Simulations: LandSimulation, step!, solve!
import ClimaLand.Parameters as LP
import ClimaComms
ClimaComms.@import_required_backends
import ClimaDiagnostics
import ClimaTimeSteppers
import ClimaUtilities.TimeManager: ITime, date

FT = Float32
toml_dict = LP.create_toml_dict(FT)
bucket_domain = ClimaLand.SphericalShell(;
    radius = FT(100),
    depth = FT(3.5),
    nelements = (1, 10),
)
bucket_atmos, bucket_rad = ClimaLand.prescribed_analytic_forcing(FT; toml_dict)
albedo = ClimaLand.Bucket.PrescribedBaregroundAlbedo{FT}(
    FT(0.8),
    (coordinate_point) -> 0.2,
    bucket_domain.space.surface,
)
parameters = ClimaLand.Bucket.BucketModelParameters(
    toml_dict;
    albedo,
    z_0m = FT(1e-2),
    z_0b = FT(1e-3),
    τc = FT(1),
)
model = ClimaLand.Bucket.BucketModel(;
    parameters,
    domain = bucket_domain,
    atmosphere = bucket_atmos,
    radiation = bucket_rad,
)
function set_ic!(Y, p, t0, model)
    Y.bucket.T .= FT(280)
    Y.bucket.W .= FT(0.5)
    Y.bucket.Ws .= FT(0.5)
    Y.bucket.σS .= FT(0)
end
Δt = 3600.0
tf = 4Δt
timestepper = ClimaTimeSteppers.ExplicitAlgorithm(ClimaTimeSteppers.RK4())

@testset "LandSimulation without a calendar, FT = $FT" begin
    simulation = LandSimulation(
        0.0,
        tf,
        Δt,
        model;
        set_ic!,
        timestepper,
        diagnostics = nothing,
        user_callbacks = (),
        updateat = Δt,
        solver_kwargs = (; saveat = 2Δt),
    )
    @test ClimaComms.context(simulation) == ClimaComms.context(model)
    @test ClimaComms.device(simulation) == ClimaComms.device(model)
    # A scalar saveat is expanded to the save times
    @test simulation._integrator._saveat ==
          collect(ITime(0.0):ITime(2Δt):ITime(tf))

    out = sprint(show, simulation)
    @test occursin("BucketModel Simulation", out)
    @test occursin("Current simulation time", out)

    step!(simulation)
    @test float(simulation._integrator.t) == Δt
    solve!(simulation)
    @test float(simulation._integrator.t) == tf
end

@testset "LandSimulation with a calendar, FT = $FT" begin
    start_date = DateTime(2005)
    stop_date = start_date + Second(round(Int, tf))
    # Update and save times may be given as periods, dates, floats or ITimes
    for (updateat, saveat) in (
        (Hour(1), [start_date, stop_date]),
        (1800.0, 2Δt),
        (ITime(1800), [ITime(0), ITime(7200)]),
    )
        simulation = LandSimulation(
            start_date,
            stop_date,
            Second(round(Int, Δt)),
            model;
            set_ic!,
            timestepper,
            diagnostics = nothing,
            user_callbacks = (),
            updateat,
            solver_kwargs = (; saveat),
        )
        @test simulation.start_date == start_date
        out = sprint(show, simulation)
        @test occursin("Current date: $(start_date)", out)
        solve!(simulation)
        @test date(simulation._integrator.t) == stop_date
    end
end

@testset "LandSimulation output directory warning, FT = $FT" begin
    start_date = DateTime(2005)
    stop_date = start_date + Second(round(Int, tf))
    ClimaLand.Diagnostics.define_diagnostics!(model, ["rn"])
    tmpdir = mktempdir(".")
    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
        bucket_domain.space.surface,
        tmpdir,
    )
    diagnostics = ClimaLand.Diagnostics.common_diagnostics(
        Val(:hourly),
        Val(:instantaneous),
        nc_writer,
        start_date,
        "rn",
    )
    @test_logs (:warn, r"inconsistent") match_mode = :any LandSimulation(
        start_date,
        stop_date,
        Δt,
        model;
        set_ic!,
        timestepper,
        diagnostics,
        user_callbacks = (),
        outdir = ".",
    )
end

@testset "LandSimulation crash reporting, FT = $FT" begin
    t0 = ITime(0.0)
    dt = ITime(Δt)
    crash_cb = ClimaLand.IntervalBasedCallback(
        dt,
        t0,
        dt,
        (integrator) -> error("Simulated crash"),
    )
    simulation = LandSimulation(
        0.0,
        tf,
        Δt,
        model;
        set_ic!,
        timestepper,
        diagnostics = nothing,
        user_callbacks = (crash_cb,),
    )
    @test_logs (:error, r"crashed") match_mode = :any solve!(simulation)
end
