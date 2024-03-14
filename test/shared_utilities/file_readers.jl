#=
    Unit tests for ClimaLand FileReaders module
=#

using Dates
using Test

import ClimaLand
import ClimaLand: FileReaders
using NCDatasets

# Include testing artifacts from Bucket
albedo_temporal_path = ClimaLand.Bucket.cesm2_albedo_dataset_path()
albedo_bareground_path = ClimaLand.Bucket.bareground_albedo_dataset_path()
include(joinpath(pkgdir(ClimaLand), "test", "artifacts", "artifacts.jl"))
multiple_variables_path = era5_t2m_sp_u10n_dataset_path()

@testset "NCFileReader with time" begin
    PATH = albedo_temporal_path
    NCDataset(PATH) do nc
        available_dates = FileReaders.read_available_dates(nc)
        @test length(available_dates) == 1980
        @test available_dates[2] == DateTime(1850, 02, 14)

        ncreader = FileReaders.NCFileReader(PATH, "sw_alb")

        @test ncreader.dimensions[1] == nc["lon"][:]
        @test ncreader.dimensions[2] == nc["lat"][:]

        @test FileReaders.read(ncreader, DateTime(1850, 02, 14)) ==
              nc["sw_alb"][:, :, 2]

        # Read it a second time to check that the cache works
        @test FileReaders.read(ncreader, DateTime(1850, 02, 14)) ==
              nc["sw_alb"][:, :, 2]

        close(ncreader)

        @test isempty(FileReaders.OPEN_NCFILES)
    end

    PATH = multiple_variables_path
    NCDataset(PATH) do nc
        ncreader_sp = FileReaders.NCFileReader(PATH, "sp")
        ncreader_u = FileReaders.NCFileReader(PATH, "u10n")

        # Test that the underlying dataset is the same
        @test ncreader_u.dataset === ncreader_sp.dataset

        # Test that we need to close all the variables to close the file
        close(ncreader_sp)
        @test !isempty(FileReaders.OPEN_NCFILES)
        close(ncreader_u)
        @test isempty(FileReaders.OPEN_NCFILES)
    end
end

@testset "NCFileReader without time" begin
    PATH = albedo_bareground_path
    NCDataset(albedo_bareground_path) do nc
        available_dates = FileReaders.read_available_dates(nc)
        @test isempty(available_dates)

        ncreader = FileReaders.NCFileReader(PATH, "sw_alb")

        @test ncreader.dimensions[1] == nc["lon"][:]
        @test ncreader.dimensions[2] == nc["lat"][:]

        @test FileReaders.read(ncreader) == nc["sw_alb"][:, :]

        close(ncreader)
    end
end
