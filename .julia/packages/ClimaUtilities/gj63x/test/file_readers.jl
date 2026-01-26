using Artifacts
using Dates
using Test

import ClimaUtilities
import ClimaUtilities.FileReaders
using NCDatasets

@testset "NCFileReader with time" begin
    PATH = joinpath(artifact"era5_example", "era5_t2m_sp_u10n_20210101.nc")
    NCDataset(PATH) do nc
        ncreader_sp = FileReaders.NCFileReader(PATH, "sp")
        ncreader_u = FileReaders.NCFileReader(PATH, "u10n")

        # Test that the underlying dataset is the same
        @test ncreader_u.dataset === ncreader_sp.dataset

        @test length(ncreader_u.available_dates) == 24
        @test length(ncreader_sp.available_dates) == 24

        @test FileReaders.available_dates(ncreader_u) ==
              ncreader_u.available_dates

        available_dates = ncreader_sp.available_dates
        @test available_dates[2] == DateTime(2021, 01, 01, 01)

        @test ncreader_sp.dimensions[1] == nc["lon"][:]
        @test ncreader_sp.dimensions[2] == nc["lat"][:]

        @test FileReaders.read(ncreader_u, DateTime(2021, 01, 01, 01)) ==
              nc["u10n"][:, :, 2]

        @test FileReaders.read(ncreader_sp, DateTime(2021, 01, 01, 01)) ==
              nc["sp"][:, :, 2]

        # Read it a second time to check that the cache works
        @test FileReaders.read(ncreader_u, DateTime(2021, 01, 01, 01)) ==
              nc["u10n"][:, :, 2]

        # Test read!
        dest = copy(nc["u10n"][:, :, 2])
        fill!(dest, 0)
        FileReaders.read!(dest, ncreader_u, DateTime(2021, 01, 01, 01))
        @test dest == nc["u10n"][:, :, 2]

        # Test that we need to close all the variables to close the file
        open_ncfiles =
            Base.get_extension(ClimaUtilities, :ClimaUtilitiesNCDatasetsExt).NCFileReaderExt.OPEN_NCFILES

        close(ncreader_sp)
        @test !isempty(open_ncfiles)
        close(ncreader_u)
        @test isempty(open_ncfiles)
    end

    # Test times split across multiple files
    PATHS = [
        joinpath(@__DIR__, "test_data", "era5_1979_1.0x1.0_lai.nc"),
        joinpath(@__DIR__, "test_data", "era5_1980_1.0x1.0_lai.nc"),
    ]
    NCDataset(PATHS, aggdim = "time") do nc
        ncreader_agg = FileReaders.NCFileReader(PATHS, "lai_lv")
        FileReaders.available_dates(ncreader_agg) == nc["time"][:]
        length(FileReaders.available_dates(ncreader_agg)) == 104
    end
end

@testset "NCFileReader without time" begin
    PATH = joinpath(
        artifact"era5_static_example",
        "era5_t2m_sp_u10n_20210101_static.nc",
    )
    NCDataset(PATH) do nc
        read_dates_func =
            Base.get_extension(ClimaUtilities, :ClimaUtilitiesNCDatasetsExt).NCFileReaderExt.read_available_dates

        available_dates = read_dates_func(nc)
        @test isempty(available_dates)

        ncreader = FileReaders.NCFileReader(PATH, "u10n")

        @test ncreader.dimensions[1] == nc["lon"][:]
        @test ncreader.dimensions[2] == nc["lat"][:]

        @test FileReaders.read(ncreader) == nc["u10n"][:, :]

        # Test read!
        dest = copy(nc["u10n"][:, :])
        fill!(dest, 0)
        FileReaders.read!(dest, ncreader)
        @test dest == nc["u10n"][:, :]

        @test isempty(FileReaders.available_dates(ncreader))

        FileReaders.close_all_ncfiles()
        open_ncfiles =
            Base.get_extension(ClimaUtilities, :ClimaUtilitiesNCDatasetsExt).NCFileReaderExt.OPEN_NCFILES
        @test isempty(open_ncfiles)
    end
end

@testset "read_available_dates" begin
    read_dates_func =
        Base.get_extension(ClimaUtilities, :ClimaUtilitiesNCDatasetsExt).NCFileReaderExt.read_available_dates

    data_dir = mktempdir()
    NCDataset(joinpath(data_dir, "test_time_1.nc"), "c") do nc
        defDim(nc, "time", 2)
        times = [DateTime(2022), DateTime(2023)]
        defVar(nc, "time", times, ("time",))
        @test read_dates_func(nc) == times
    end
    NCDataset(joinpath(data_dir, "test_date_1.nc"), "c") do nc
        defDim(nc, "date", 2)
        times = [20220101, 20230101]
        defVar(nc, "date", times, ("date",))
        @test read_dates_func(nc) == DateTime.(string.(times), "yyyymmdd")
    end
end

@testset "read_missing_dims" begin
    @test_throws contains(
        "missing_dim.nc\"] does not contain information about dimensions (\"missing_dim\",)",
    ) FileReaders.NCFileReader(
        joinpath(@__DIR__, "test_data", "missing_dim.nc"),
        "test_var",
    )
end
