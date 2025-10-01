# Calibration End-to-End Test Suite

# This file contains tests designed to be run after a test calibration and may
# not necessarily work with a different calibration configuration.

using Test
import ClimaLand
import Statistics: std

# Needed to access CalibrateConfig
include(
    joinpath(pkgdir(ClimaLand), "experiments/calibration/run_calibration.jl"),
)

@testset "Observations" begin
    (; obs_vec_filepath, sample_date_ranges, short_names) = CALIBRATE_CONFIG
    isfile(obs_vec_filepath) || error("Cannot find the observations generated")
    obs_vec = JLD2.load_object(obs_vec_filepath)

    @test length(sample_date_ranges) == length(obs_vec)

    metadata = first(EKP.get_metadata(first(obs_vec)))
    @test metadata.attributes["short_name"] == first(short_names)
    @test Dates.Second(first(metadata.dims["time"])) +
          Dates.DateTime(metadata.attributes["start_date"]) ==
          first(first(sample_date_ranges))
    @test Dates.Second(last(metadata.dims["time"])) +
          Dates.DateTime(metadata.attributes["start_date"]) ==
          first(last(sample_date_ranges))
end

@testset "Calibration results" begin
    isdir(CALIBRATE_CONFIG.output_dir) || error(
        "Cannot find output directory of calibration. These tests will not work without first running a test calibration",
    )

    (; output_dir, n_iterations) = CALIBRATE_CONFIG
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, n_iterations))

    @test EKP.get_N_iterations(ekp) == n_iterations

    # Spread should decrease
    ekp_u = EKP.get_u(ekp)
    spreads = map(std, ekp_u)
    @test 0.1 * first(spreads) > last(spreads)
end

nothing
