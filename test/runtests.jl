using ClimaLand
using ParallelTestRunner

# The environment variable `CLIMALAND_TEST_GROUP` determines what set of tests
# to run. The available options are "all", "core", and "simulations".
const TEST_GROUP = get(ENV, "CLIMALAND_TEST_GROUP", "all")
TEST_GROUP in ("all", "core", "simulations") || error(
    "CLIMALAND_TEST_GROUP must be \"core\" or \"simulations\", got \"$TEST_GROUP\"",
)

# Simplify adding tests to the test suite; paths are absolute because workers
# start in the caller's working directory
macro include(files...)
    return Expr(:quote, :(foreach(include, $(joinpath.(@__DIR__, files)))))
end

#! format: off
core_tests = Dict(
    # Performance and code quality tests

    "Aqua tests" => @include("aqua.jl"),

    # Shared ClimaLand utilities tests
    "Richards model implicit timestepping tests" => @include("shared_utilities/implicit_timestepping/richards_model.jl"),
    "Full soil model implicit timestepping tests" => @include("shared_utilities/implicit_timestepping/energy_hydrology_model.jl"),
    "Domains module tests" => @include("shared_utilities/domains.jl"),
    "General utilities tests" => @include("shared_utilities/utilities.jl"),
    "Show methods tests" => @include("shared_utilities/show_methods.jl"),
    "Land simulation interface tests" => @include("shared_utilities/simulations.jl"),
    "Variable types tests" => @include("shared_utilities/variable_types.jl"),
    "Time integrated variables tests" => @include("shared_utilities/time_integrated_variables.jl"),

    "Driver tests" => @include("shared_utilities/drivers.jl", "shared_utilities/coupled_fluxes.jl"),
    "Stability cap tests" => @include("shared_utilities/stability_cap.jl"),

    # Standalone Bucket model tests
    "Bucket albedo types tests" => @include("standalone/Bucket/albedo_types.jl"),

    "Bucket snow tests" => @include("standalone/Bucket/snow_bucket_tests.jl"),
    "Bucket soil tests" => @include("standalone/Bucket/soil_bucket_tests.jl"),
    "Bucket restart tests" => @include("standalone/Bucket/restart.jl"),

    "Inland water" => @include("standalone/InlandWater/unit_tests.jl"),

    # Standalone Snow model tests
    "Snow parameterization tests" => @include("standalone/Snow/parameterizations.jl", "standalone/Snow/snow.jl"),
    "Neural Snow model tools tests" => @include("standalone/Snow/tool_tests.jl"),

    "Snow integrated water and energy content" => @include("standalone/Snow/conservation.jl"),

    "Snow parameter constructors" => @include("standalone/Snow/parameters.jl"),


    # Standalone Soil model tests
    "Soil Biogeochemistry module tests" => @include("standalone/Soil/Biogeochemistry/biogeochemistry_module.jl"),
    "Soil CO2 parameterization tests" => @include("standalone/Soil/Biogeochemistry/co2_parameterizations.jl"),

    "Soil CO2 saturation stability tests" => @include("standalone/Soil/Biogeochemistry/saturation_stability_test.jl"),

    "Soil climate drivers tests" => @include("standalone/Soil/climate_drivers.jl"),
    "Soil runoff tests" => @include("standalone/Soil/runoff.jl"),
    "Soil boundary condition tests" => @include("standalone/Soil/soil_bc.jl"),
    "Soil parameterization tests" => @include("standalone/Soil/soil_parameterizations.jl"),
    "Soil 3D domain tests" => @include("standalone/Soil/soil_test_3d.jl", "standalone/Soil/mask_test.jl"),
    "Soil integration tests" => @include("standalone/Soil/soiltest.jl"),

    "Soil spatial parameters and parameter constructors" => @include("standalone/Soil/parameters.jl"),

    # Standalone Surface Water model tests
    "Pond module tests" => @include("standalone/SurfaceWater/pond_test.jl"),

    # Standalone Vegetation model tests
    "Canopy module tests" => @include("standalone/Vegetation/canopy_model.jl"),

    "Canopy soil moisture stress tests" => @include("standalone/Vegetation/test_soil_moisture_stress.jl"),
    "Canopy PlantHydraulics tests" => @include("standalone/Vegetation/plant_hydraulics_test.jl"),
    "Big Leaf model tests" => @include("standalone/Vegetation/test_bigleaf_parameterizations.jl"),

    "Canopy helper function tests" => @include("standalone/Vegetation/test_canopy_helpers.jl"),
    "Two Stream model tests" => @include("standalone/Vegetation/test_two_stream.jl"),
    "Canopy integrated water and energy content" => @include("standalone/Vegetation/conservation.jl"),
    "Canopy prognostic carbon pools" => @include("standalone/Vegetation/test_prognostic_carbon.jl"),
    "Canopy spatial parameters" => @include("standalone/Vegetation/spatial_parameters.jl"),
    "Canopy spatially varying height" => @include("standalone/Vegetation/test_spatially_varying_canopy_height.jl"),
    "P model tests" => @include("standalone/Vegetation/test_pmodel.jl"),
    "Optimal LAI tests" => @include("standalone/Vegetation/test_optimal_lai.jl"),
    "PFT tests" => @include("standalone/Vegetation/test_pfts.jl"),

    # Integrated LSM tests
    "Integrated LSM unit tests" => @include("integrated/lsms.jl"),

    "Integrated soil/canopy unit tests" => @include("integrated/soil_canopy_lsm.jl"),
    "Integrated pond/soil LSM tests" => @include("integrated/pond_soil_lsm.jl"),
    "Integrated soil energy/hydrology/biogeochem LSM tests" => @include("integrated/soil_energy_hydrology_biogeochemistry.jl"),
    "Integrated soil and snow" => @include("integrated/soil_snow.jl"),

    # Test compatibility with automatic differentiation
    "Autodiff compatibility" => @include("shared_utilities/ad_compatibility.jl"),
)

simulation_tests = Dict(
    "Full land" => @include("integrated/full_land.jl", "integrated/conservation.jl"),
    "Full land restart tests" => @include("integrated/restart.jl"),
    "Soil carbon spin-up" => @include("integrated/soil_carbon_spinup.jl"),

    # FluxnetSimulations extension tests
    "Fluxnet simulation setup unit tests" => @include("integrated/fluxnet_sim.jl"),

    "Perturbed driver tests" => @include("shared_utilities/perturbed_drivers.jl"),

    # Diagnostics
    "Diagnostics" => @include("diagnostics/diagnostics_tests.jl"),
    "Diagnostics compute methods" => @include("diagnostics/compute_methods.jl"),

    "ILAMB setup" => @include("../experiments/ilamb/tests/test_ilamb_setup.jl"),
)
#! format: on

testsuite = Dict{String, Expr}()
TEST_GROUP in ("all", "core") && merge!(testsuite, core_tests)
TEST_GROUP in ("all", "simulations") && merge!(testsuite, simulation_tests)
# Increase max_worker_rss so that a worker doesn't restart as often
gpu_options =
    get(ENV, "CLIMACOMMS_DEVICE", "CPU") == "CUDA" ?
    (; max_worker_rss = 16 * 2^30, max_default_jobs = 3) : (;)
# Workers skip startup files, where some clusters add the ClimaArtifacts depots
env = ["JULIA_DEPOT_PATH" => join(DEPOT_PATH, Sys.iswindows() ? ';' : ':')]
# Run only as many workers as fit in memory at about 4 GiB each
runtests(
    ClimaLand,
    ARGS;
    testsuite,
    env,
    memory_per_worker = 4 * 2^30,
    gpu_options...,
)
