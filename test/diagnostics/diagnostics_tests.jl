using Test
using ClimaLand
using ClimaLand: Domains, Soil, Canopy
using ClimaLand.Simulations: LandSimulation, step!
using ClimaLand.Diagnostics: @with_error
import ClimaLand.Parameters as LP
import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams
import ClimaLand.Parameters as LP
import SciMLBase
import ClimaTimeSteppers
import ClimaDiagnostics
import ClimaCore
using Dates
using Statistics

# This test must go first because it checks the global `ALL_DIAGNOSTICS`
@testset "Diagnostics writing and reading" begin
    FT = Float32
    @test isdefined(ClimaLand.Diagnostics, :compute_sw_albedo!)

    # Define some diagnostics for a DummyModel

    @test ClimaLand.Diagnostics.ALL_DIAGNOSTICS isa Dict
    # @test length(ClimaLand.Diagnostics.ALL_DIAGNOSTICS) == 0
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

    # @test length(ClimaLand.Diagnostics.ALL_DIAGNOSTICS) == 1

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

    toml_dict = LP.create_toml_dict(FT)
    bucket_atmos, bucket_rad =
        ClimaLand.prescribed_analytic_forcing(FT; toml_dict)
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
    bucket_parameters = ClimaLand.Bucket.BucketModelParameters(
        toml_dict;
        albedo,
        z_0m,
        z_0b,
        τc,
    )
    toml_dict = LP.create_toml_dict(FT)
    bucket_parameters = ClimaLand.Bucket.BucketModelParameters(
        toml_dict;
        albedo,
        z_0m,
        z_0b,
        τc,
    )

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

    start_date = DateTime(2005)
    out = ClimaLand.Diagnostics.common_diagnostics(
        Val(:hourly),
        Val(:average),
        nc_writer,
        start_date,
        diags...;
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
    reduction_period = :every_dt
    reduction_type = :instantaneous

    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        model,
        start_date;
        output_writer,
        output_vars,
        reduction_period,
        reduction_type,
        dt,
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

    # Test that the diagnostics were correctly created and computed once at initialization
    @test keys(simulation.diagnostics[1].output_writer.dict) ==
          Set(["swp_1000s_inst", "swc_1000s_inst", "sie_1000s_inst"])
    @test length(
        simulation.diagnostics[1].output_writer.dict["swp_1000s_inst"].keys,
    ) == 1 # number of diagnostic computations so far
    @test all(
        ClimaCore.Fields.field2array(
            simulation.diagnostics[1].output_writer.dict["swc_1000s_inst"].vals[1],
        ) .== FT(0.24),
    )

    # Step the simulation once to compute diagnostics again
    step!(simulation)

    @test length(
        simulation.diagnostics[1].output_writer.dict["swp_1000s_inst"].keys,
    ) == 2 # number of diagnostic computations so far
    # Check that the SWC values have changed after the second computation
    @test simulation.diagnostics[1].output_writer.dict["swc_1000s_inst"].vals[1] !=
          simulation.diagnostics[1].output_writer.dict["swc_1000s_inst"].vals[2]
end

@testset "LandModel diagnostics" begin
    LAI = TimeVaryingInput((t) -> FT(1.0))
    model = Soil.LandModel{FT}((; atmos, radiation), LAI, toml_dict, domain, dt)

    function set_ic!(Y, p, t0, model)
        Y.soil.ϑ_l .= FT(0.24)
        Y.soil.θ_i .= FT(0.0)
        T = FT(290.15)
        ρc_s =
            Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Y.soil.θ_i,
                model.soil.parameters.ρc_ds,
                model.soil.parameters.earth_param_set,
            )
        Y.soil.ρe_int .=
            Soil.volumetric_internal_energy.(
                Y.soil.θ_i,
                ρc_s,
                T,
                model.soil.parameters.earth_param_set,
            )

        Y.soilco2.C = FT(0.000412) # set to atmospheric co2, mol co2 per mol air

        Y.canopy.hydraulics.ϑ_l.:1 .= model.canopy.hydraulics.parameters.ν
        Y.canopy.energy.T = FT(297.5)
    end

    output_writer = ClimaDiagnostics.Writers.DictWriter()
    output_vars = ["swc", "ct", "sco2"]
    reduction_period = :every_dt
    reduction_type = :instantaneous

    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        model,
        start_date;
        output_writer,
        output_vars,
        reduction_period,
        reduction_type,
        dt,
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

    # Test that the diagnostics were correctly created and computed once at initialization
    @test keys(simulation.diagnostics[1].output_writer.dict) ==
          Set(["swc_1000s_inst", "ct_1000s_inst", "sco2_1000s_inst"])
    @test all(
        ClimaCore.Fields.field2array(
            simulation.diagnostics[1].output_writer.dict["swc_1000s_inst"].vals[1],
        ) .== FT(0.24),
    )
    @test all(
        ClimaCore.Fields.field2array(
            simulation.diagnostics[1].output_writer.dict["ct_1000s_inst"].vals[1],
        ) .== FT(297.5),
    )
    @test all(
        ClimaCore.Fields.field2array(
            simulation.diagnostics[1].output_writer.dict["sco2_1000s_inst"].vals[1],
        ) .== FT(0.000412),
    )
end

@testset "CanopyModel invalid windspeed diagnostic (coupled)" begin
    # surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface
    atmos_h = ClimaCore.Fields.ones(surface_space) .* FT(50)
    atmos = CoupledAtmosphere{FT}(surface_space, atmos_h)
    radiation = CoupledRadiativeFluxes{FT}()
    ground = ClimaLand.PrognosticGroundConditions{FT}()
    LAI = TimeVaryingInput((t) -> FT(1.0))
    model = ClimaLand.SoilCanopyModel{FT}(
        (; atmos, radiation, ground),
        LAI,
        toml_dict,
        domain,
    )

    output_writer = ClimaDiagnostics.Writers.DictWriter()
    output_vars = ["ws"]
    reduction_period = :every_dt
    reduction_type = :instantaneous

    # This will fail because wind speed is not an available diagnostic in this setup
    @test_throws AssertionError diagnostics =
        ClimaLand.Diagnostics.default_diagnostics(
            model,
            start_date;
            output_writer,
            output_vars,
            reduction_period,
            reduction_type,
            dt,
        )
end

@testset "Invalid diagnostic variable" begin
    ground = ClimaLand.PrognosticGroundConditions{FT}()
    LAI = TimeVaryingInput((t) -> FT(1.0))
    model = ClimaLand.SoilCanopyModel{FT}(
        (; atmos, radiation, ground),
        LAI,
        toml_dict,
        domain,
    )

    output_writer = ClimaDiagnostics.Writers.DictWriter()
    output_vars = ["swc", "ct", "invalid_diagnostic"]
    reduction_period = :every_dt
    reduction_type = :instantaneous

    # This will fail because "invalid_diagnostic" is not an available diagnostic for this model
    @test_throws AssertionError ClimaLand.Diagnostics.default_diagnostics(
        model,
        start_date;
        output_writer,
        output_vars,
        reduction_period,
        reduction_type,
        dt,
    )
end

@testset "Test runoff diagnostics" begin
    # TOPMODELRunoff by default
    model_topmodelrunoff =
        Soil.EnergyHydrology{FT}(domain, (; atmos, radiation), toml_dict)

    # Set up diagnostics with variables "sr" and "ssr"
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    output_vars = ["sr", "ssr"]
    reduction_period = :every_dt
    reduction_type = :instantaneous
    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        model_topmodelrunoff,
        start_date;
        output_writer,
        output_vars,
        reduction_period,
        reduction_type,
        dt,
    )

    simulation = LandSimulation(
        start_date,
        stop_date,
        dt,
        model_topmodelrunoff;
        diagnostics,
        user_callbacks = (),
    )

    # Test that the diagnostics were correctly created
    @test keys(simulation.diagnostics[1].output_writer.dict) ==
          Set(["sr_1000s_inst", "ssr_1000s_inst"])

    # SurfaceRunoff
    model_surfacerunoff = Soil.EnergyHydrology{FT}(
        domain,
        (; atmos, radiation),
        toml_dict;
        runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(),
    )

    # Set up diagnostics with variables "sr" and "ssr"
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    output_vars = ["sr"]
    reduction_period = :every_dt
    reduction_type = :instantaneous
    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        model_surfacerunoff,
        start_date;
        output_writer,
        output_vars,
        reduction_period,
        reduction_type,
        dt,
    )

    simulation = LandSimulation(
        start_date,
        stop_date,
        dt,
        model_surfacerunoff;
        diagnostics,
        user_callbacks = (),
    )

    # Test that the diagnostics were correctly created
    @test keys(simulation.diagnostics[1].output_writer.dict) ==
          Set(["sr_1000s_inst"])

    # NoRunoff
    model_norunoff = Soil.EnergyHydrology{FT}(
        domain,
        (; atmos, radiation),
        toml_dict;
        runoff = ClimaLand.Soil.Runoff.NoRunoff(),
    )

    # Set up diagnostics with variables "sr" and "ssr"
    output_writer = ClimaDiagnostics.Writers.DictWriter()
    output_vars = ["sr"]
    reduction_period = :every_dt
    reduction_type = :instantaneous

    # This will fail because "invalid_diagnostic" is not an available diagnostic for this model
    @test_throws AssertionError ClimaLand.Diagnostics.default_diagnostics(
        model_norunoff,
        start_date;
        output_writer,
        output_vars,
        reduction_period,
        reduction_type,
        dt,
    )
end
