using Test
using ClimaLand
using ClimaLand.Diagnostics: @with_error
import ClimaLand.Parameters as LP
import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams
import SciMLBase
import ClimaTimeSteppers
import ClimaDiagnostics
using Dates
using Statistics

@test isdefined(ClimaLand.Diagnostics, :compute_sw_albedo!)

# Define some diagnostics for a DummyModel

@test ClimaLand.Diagnostics.ALL_DIAGNOSTICS isa Dict
@test length(ClimaLand.Diagnostics.ALL_DIAGNOSTICS) == 0
struct DummyModel end
struct DummyModel2 end
ClimaLand.Diagnostics.@diagnostic_compute "sw_albedo" Union{
    DummyModel,
    DummyModel2,
} p.foo.bar

ClimaLand.Diagnostics.add_diagnostic_variable!(
    short_name = "alpha",
    long_name = "Albedo",
    standard_name = "albedo",
    units = "",
    compute! = (out, Y, p, t) -> compute_sw_albedo!(out, Y, p, t, land_model),
)

@test length(ClimaLand.Diagnostics.ALL_DIAGNOSTICS) == 1

# First, run a simulation for 1 hour
FT = Float32
seconds = 1.0;
minutes = 60seconds;
hours = 60minutes;
t0 = 0.0
Δt = 1minutes
tf = 2hours
bucket_domain = ClimaLand.SphericalShell(;
    radius = FT(100),
    depth = FT(3.5),
    nelements = (1, 10),
)

bucket_atmos, bucket_rad = ClimaLand.prescribed_analytic_forcing(FT)
τc = FT(1.0)
α_bareground_func = (coordinate_point) -> 0.2
α_snow = FT(0.8)
z_0m = FT(1e-2)
z_0b = FT(1e-3)
albedo = ClimaLand.Bucket.PrescribedBaregroundAlbedo{FT}(
    α_snow,
    α_bareground_func,
    bucket_domain.space.surface,
)
bucket_parameters =
    ClimaLand.Bucket.BucketModelParameters(FT; albedo, z_0m, z_0b, τc)

model = ClimaLand.Bucket.BucketModel(
    parameters = bucket_parameters,
    domain = bucket_domain,
    atmosphere = bucket_atmos,
    radiation = bucket_rad,
)

Y, p, coords = ClimaLand.initialize(model)
Y.bucket.T .= 280.0
Y.bucket.W .= 0.5 # no moisture
Y.bucket.Ws .= 0.5 # no runoff
Y.bucket.σS .= 0.0

exp_tendency! = ClimaLand.make_exp_tendency(model)
set_initial_cache! = ClimaLand.make_set_initial_cache(model)
set_initial_cache!(p, Y, t0)

prob = SciMLBase.ODEProblem(
    ClimaTimeSteppers.ClimaODEFunction((T_exp!) = exp_tendency!),
    Y,
    (t0, tf),
    p,
)

# ClimaDiagnostics

ClimaLand.Diagnostics.define_diagnostics!(model)
diags = ["rn", "lhf"]

tmpdir = mktempdir(".")
nc_writer =
    ClimaDiagnostics.Writers.NetCDFWriter(bucket_domain.space.surface, tmpdir)

out = ClimaLand.Diagnostics.hourly_averages(
    FT,
    diags...;
    output_writer = nc_writer,
    start_date = DateTime(2005),
)

diagnostic_handler = ClimaDiagnostics.DiagnosticsHandler(out, Y, p, t0; dt = Δt)

diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

updateat = collect(t0:Δt:tf);
drivers = ClimaLand.get_drivers(model)
updatefunc = ClimaLand.make_update_drivers(drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, diag_cb)
timestepper = ClimaTimeSteppers.RK4()
ode_algo = ClimaTimeSteppers.ExplicitAlgorithm(timestepper)
SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb)

using ClimaAnalysis
simdir = ClimaAnalysis.SimDir(tmpdir)
rn = get(simdir; short_name = "rn")
@test readdir(tmpdir) == ["lhf_1h_average.nc", "rn_1h_average.nc"]
@test length(rn.dims) == 3
@test mean(rn.data) != 0.0


@testset "runoff compute methods" begin
    using ClimaLand
    import ClimaParams
    import ClimaLand.Parameters as LP
    using ClimaDiagnostics
    using Dates

    FT = Float32
    domain = ClimaLand.Domains.global_domain(FT)
    atmos, radiation = ClimaLand.prescribed_analytic_forcing(FT)
    forcing = (; atmos, radiation)
    earth_param_set = LP.LandParameters(FT)
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

    # Soil model
    soil = ClimaLand.Soil.EnergyHydrology{FT}(
        domain,
        forcing,
        earth_param_set;
        prognostic_land_components,
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
        runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(), #NoRunoff(),
    )

    start_date = DateTime(2008, 1, 1)
    stop_date = DateTime(2008, 3, 1)
    dt = 450.0

    modis_lai_ncdata_path = ClimaLand.Artifacts.modis_lai_multiyear_paths(;
        start_date,
        end_date = stop_date,
    )
    LAI = ClimaLand.prescribed_lai_modis(
        modis_lai_ncdata_path,
        domain.space.surface,
        start_date,
    )

    model =
        ClimaLand.LandModel{FT}(forcing, LAI, earth_param_set, domain, dt; soil)

    # Test that diagnostics don't error when trying to produce nonexistent variables
    output_vars = ["sr", "ssr"]
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    average_period = :hourly

    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        model,
        start_date;
        output_writer,
        output_vars,
        average_period,
    )

    simulation = ClimaLand.Simulations.LandSimulation(
        start_date,
        stop_date,
        dt,
        model;
        diagnostics,
    )
end
