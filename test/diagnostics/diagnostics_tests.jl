using Test
using ClimaLand
using ClimaLand: Domains, Soil
using ClimaLand.Simulations: LandSimulation, step!
import ClimaLand.Parameters as LP
import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams
import SciMLBase
import ClimaTimeSteppers
import ClimaDiagnostics
using Dates
using Statistics

# Define some variables to reuse for each model
FT = Float32
default_params_filepath =
    joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
toml_dict = LP.create_toml_dict(FT, default_params_filepath);
earth_param_set = LP.LandParameters(toml_dict);

zmax = FT(0)
zmin = FT(-1.0)
longlat = FT.((-118.1, 34.1))
domain = Domains.Column(; zlim = (zmin, zmax), nelements = 10, longlat);
surface_space = domain.space.surface;

start_date = DateTime(2008);
stop_date = start_date + Second(60 * 60 * 72);
dt = 1000.0;

era5_ncdata_path =
    ClimaLand.Artifacts.era5_land_forcing_data2008_path(; lowres = true);
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    earth_param_set,
    FT,
);

@testset "EnergyHydrology diagnostics" begin
    model = Soil.EnergyHydrology{FT}(domain, (; atmos, radiation), toml_dict)

    function set_ic!(Y, p, t0, model)
        Y.soil.ϑ_l .= FT(0.24)
        Y.soil.θ_i .= FT(0.0)
        T = FT(290.15)
        ρc_s =
            Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Y.soil.θ_i,
                model.parameters.ρc_ds,
                model.parameters.earth_param_set,
            )
        Y.soil.ρe_int .=
            Soil.volumetric_internal_energy.(
                Y.soil.θ_i,
                ρc_s,
                T,
                model.parameters.earth_param_set,
            )
    end

    output_writer = ClimaDiagnostics.Writers.DictWriter()
    output_vars = ["swc", "sie", "swp"]
    # TODO compute and output at every dt instead of half hourly (this is the ClimaDiagnostics default)
    average_period = :halfhourly

    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        model,
        start_date;
        output_writer,
        output_vars,
        average_period,
    )

    simulation = LandSimulation(
        start_date,
        stop_date,
        dt,
        model;
        set_ic!,
        diagnostics,
        user_callbacks = (),
    )

    (; p, t) = simulation._integrator
    Y = simulation._integrator.u

    # Construct a DiagnosticsHandler object so we can manually compute diagnostics
    diagnostics_handler_init =
        ClimaDiagnostics.DiagnosticsHandler(simulation.diagnostics, Y, p, t; dt)
    ClimaDiagnostics.orchestrate_diagnostics(
        simulation._integrator,
        diagnostics_handler,
    )

    step!(simulation)
end

@testset "Diagnostics writing and reading" begin
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
        compute! = (out, Y, p, t) ->
            compute_sw_albedo!(out, Y, p, t, land_model),
    )

    @test length(ClimaLand.Diagnostics.ALL_DIAGNOSTICS) == 1

    # First, run a simulation for 1 hour
    seconds = 1.0
    minutes = 60seconds
    hours = 60minutes
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
    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
        bucket_domain.space.surface,
        tmpdir,
    )

    out = ClimaLand.Diagnostics.hourly_averages(
        FT,
        diags...;
        output_writer = nc_writer,
        start_date = DateTime(2005),
    )

    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(out, Y, p, t0; dt = Δt)

    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

    updateat = collect(t0:Δt:tf)
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
end
