# # Running the tests
#
# Tests run in parallel across worker processes via ParallelTestRunner.jl. Each
# test file runs in a fresh sandbox module on a worker, so files must be
# self-contained.
#
# Run the whole suite with
#
#     julia --project=. -e 'using Pkg; Pkg.test("ClimaLand")'
#
# or, from the test environment,
#
#     julia --project=test test/runtests.jl
#
# Positional arguments select the tests whose name (the file path under `test/`
# without the `.jl` extension) starts with the given prefix:
#
#     Pkg.test("ClimaLand"; test_args = ["standalone/Soil"])
#     Pkg.test("ClimaLand"; test_args = ["integrated/restart", "aqua"])
#     julia --project=test test/runtests.jl shared_utilities
#
# Useful flags (pass them via `test_args` as well):
#
#   --list        print the test names and exit
#   --jobs=N      number of workers (default: from CPU count and free memory)
#   --quickfail   stop the whole run at the first failure
#   --verbose     print per-test start times
#
# `--jobs=1` is handy when a failure looks like a parallelism artifact, or when
# the interleaved output of several workers is hard to read.

using ClimaLand
using ParallelTestRunner

# Only the files listed here run: `Artifacts.jl` is a helper and
# `standalone/Snow/tool_tests` is currently disabled.
const TESTS = [
    # Performance and code quality
    "aqua",

    # Shared ClimaLand utilities
    "shared_utilities/implicit_timestepping/richards_model",
    "shared_utilities/implicit_timestepping/energy_hydrology_model",
    "shared_utilities/domains",
    "shared_utilities/utilities",
    "shared_utilities/variable_types",
    "shared_utilities/time_integrated_variables",
    "shared_utilities/drivers",
    "shared_utilities/coupled_fluxes",

    # Standalone Bucket model
    "standalone/Bucket/albedo_types",
    "standalone/Bucket/snow_bucket_tests",
    "standalone/Bucket/soil_bucket_tests",
    "standalone/Bucket/restart",

    # Standalone InlandWater model
    "standalone/InlandWater/unit_tests",

    # Standalone Snow model
    "standalone/Snow/parameterizations",
    "standalone/Snow/snow",
    "standalone/Snow/conservation",
    "standalone/Snow/parameters",

    # Standalone Soil model
    "standalone/Soil/Biogeochemistry/biogeochemistry_module",
    "standalone/Soil/Biogeochemistry/co2_parameterizations",
    "standalone/Soil/Biogeochemistry/saturation_stability_test",
    "standalone/Soil/climate_drivers",
    "standalone/Soil/runoff",
    "standalone/Soil/soil_bc",
    "standalone/Soil/soil_parameterizations",
    "standalone/Soil/soil_test_3d",
    "standalone/Soil/mask_test",
    "standalone/Soil/soiltest",
    "standalone/Soil/conservation",
    "standalone/Soil/parameters",

    # Standalone SurfaceWater model
    "standalone/SurfaceWater/pond_test",

    # Standalone Vegetation model
    "standalone/Vegetation/canopy_model",
    "standalone/Vegetation/test_soil_moisture_stress",
    "standalone/Vegetation/plant_hydraulics_test",
    "standalone/Vegetation/test_bigleaf_parameterizations",
    "standalone/Vegetation/test_two_stream",
    "standalone/Vegetation/conservation",
    "standalone/Vegetation/spatial_parameters",
    "standalone/Vegetation/test_spatially_varying_canopy_height",
    "standalone/Vegetation/test_pmodel",
    "standalone/Vegetation/test_optimal_lai",
    "standalone/Vegetation/test_pfts",

    # Integrated LSMs
    "integrated/lsms",
    "integrated/soil_canopy_lsm",
    "integrated/pond_soil_lsm",
    "integrated/soil_energy_hydrology_biogeochemistry",
    "integrated/soil_snow",
    "integrated/full_land",
    "integrated/restart",

    # FluxnetSimulations extension
    "integrated/fluxnet_sim",

    # Diagnostics
    "diagnostics/diagnostics_tests",
]

# Tests that build global domains at the default resolution use enough memory
# that running them alongside other workers risks out-of-memory failures.
const SERIAL = String[]

testsuite = find_tests(@__DIR__)
stale = setdiff([TESTS; SERIAL], keys(testsuite))
isempty(stale) || error("No such test file(s): $(join(sort!(stale), ", "))")
filter!(((name, _),) -> name in TESTS, testsuite)

# The ILAMB setup test lives next to the ILAMB experiment instead of in `test/`.
ilamb_test = abspath(
    joinpath(
        @__DIR__,
        "..",
        "experiments",
        "ilamb",
        "tests",
        "test_ilamb_setup.jl",
    ),
)
testsuite["ilamb_setup"] = :(include($ilamb_test))

# ClimaCore objects are heavily parametrized, so non-abbreviated stacktraces are
# hard to read. Julia only abbreviates them when running interactively, so force
# it on the workers as well.
init_worker_code = quote
    redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))
end

runtests(ClimaLand, ARGS; testsuite, init_worker_code, serial = SERIAL)
