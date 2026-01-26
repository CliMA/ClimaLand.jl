using Dates
using Artifacts
using Test

import ClimaUtilities
import ClimaUtilities.DataHandling

import ClimaCore
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import Interpolations as Intp
import ClimaCoreTempestRemap
using NCDatasets
include("TestTools.jl")
const context = ClimaComms.context()
ClimaComms.init(context)

@testset "DataHandler, static" begin
    PATH = joinpath(
        artifact"era5_static_example",
        "era5_t2m_sp_u10n_20210101_static.nc",
    )
    for regridder_type in (:InterpolationsRegridder, :TempestRegridder)
        if regridder_type == :TempestRegridder && Sys.iswindows()
            continue
        end
        for FT in (Float32, Float64),
            use_spacefillingcurve in (true, false),
            varnames in (["sp"], ["sp", "t2m"])

            # TempestRegridder does not support multiple input variables
            if regridder_type == :TempestRegridder && length(varnames) > 1
                continue
            end

            radius = FT(6731e3)
            helem = 40
            Nq = 4

            horzdomain = ClimaCore.Domains.SphereDomain(radius)
            horzmesh =
                ClimaCore.Meshes.EquiangularCubedSphere(horzdomain, helem)
            if use_spacefillingcurve
                horztopology = ClimaCore.Topologies.Topology2D(
                    context,
                    horzmesh,
                    ClimaCore.Topologies.spacefillingcurve(horzmesh),
                )
            else
                horztopology =
                    ClimaCore.Topologies.Topology2D(context, horzmesh)
            end
            quad = ClimaCore.Spaces.Quadratures.GLL{Nq}()
            target_space =
                ClimaCore.Spaces.SpectralElementSpace2D(horztopology, quad)

            compose_function =
                length(varnames) == 1 ? identity : (x, y) -> x + y
            data_handler = DataHandling.DataHandler(
                PATH,
                varnames,
                target_space;
                regridder_type,
                compose_function,
            )

            @test DataHandling.available_dates(data_handler) == DateTime[]

            field = DataHandling.regridded_snapshot(data_handler)
            # Most basic testing, to get a sense that we are doing things right
            @test field isa ClimaCore.Fields.Field
            @test axes(field) == target_space
            @test eltype(field) == FT

            # For one variable, check that the regridded data is the same as the original data
            if length(varnames) == 1
                @test minimum(field) >= minimum(
                    data_handler.file_readers[first(varnames)].dataset[first(
                        varnames,
                    )][
                        :,
                        :,
                    ],
                )
                @test maximum(field) <= maximum(
                    data_handler.file_readers[first(varnames)].dataset[first(
                        varnames,
                    )][
                        :,
                        :,
                    ],
                )
            else
                # For more than one variable, check that both point to the same file path
                @test data_handler.file_readers[varnames[1]].file_paths ==
                      data_handler.file_readers[varnames[2]].file_paths
                # For more than one variable, check that the compose function was applied correctly
                @test data_handler.compose_function == compose_function

                @test minimum(field) >= minimum(
                    compose_function.(
                        data_handler.file_readers[varnames[1]].dataset[varnames[1]][
                            :,
                            :,
                        ],
                        data_handler.file_readers[varnames[2]].dataset[varnames[2]][
                            :,
                            :,
                        ],
                    ),
                )
                @test maximum(field) <= maximum(
                    compose_function.(
                        data_handler.file_readers[varnames[1]].dataset[varnames[1]][
                            :,
                            :,
                        ],
                        data_handler.file_readers[varnames[2]].dataset[varnames[2]][
                            :,
                            :,
                        ],
                    ),
                )
            end
            close(data_handler)
        end
    end

    # Test passing arguments down
    varnames = ["sp"]
    radius = 6731e3
    helem = 40
    Nq = 4

    horzdomain = ClimaCore.Domains.SphereDomain(radius)
    horzmesh = ClimaCore.Meshes.EquiangularCubedSphere(horzdomain, helem)
    horztopology = ClimaCore.Topologies.Topology2D(context, horzmesh)
    quad = ClimaCore.Spaces.Quadratures.GLL{Nq}()
    target_space = ClimaCore.Spaces.SpectralElementSpace2D(horztopology, quad)

    data_handler = DataHandling.DataHandler(
        PATH,
        varnames,
        target_space;
        regridder_type = :InterpolationsRegridder,
        file_reader_kwargs = (; preprocess_func = (data) -> 0.0 * data),
        regridder_kwargs = (;
            extrapolation_bc = (Intp.Flat(), Intp.Flat(), Intp.Flat())
        ),
    )

    @test data_handler.regridder.extrapolation_bc ==
          (Intp.Flat(), Intp.Flat(), Intp.Flat())
    field = DataHandling.regridded_snapshot(data_handler)
    @test extrema(field) == (0.0, 0.0)
end

@testset "Datahandler incorrect dims" begin
    target_space = make_spherical_space(Float64; context)
    # netcdf file with dimensions as expected
    sorted_dims_path = joinpath(@__DIR__, "test_data", "sorted_dims.nc")
    # netcdf file with descending z and reverse lat
    reversed_dims_path = joinpath(@__DIR__, "test_data", "reversed_dims.nc")
    data_handler_sorted = DataHandling.DataHandler(
        sorted_dims_path,
        ["test_var"],
        target_space.hybrid;
        regridder_type = :InterpolationsRegridder,
        regridder_kwargs = (;
            extrapolation_bc = (Intp.Periodic(), Intp.Flat(), Intp.Flat())
        ),
    )
    data_handler_reversed = DataHandling.DataHandler(
        reversed_dims_path,
        ["test_var"],
        target_space.hybrid;
        regridder_type = :InterpolationsRegridder,
        regridder_kwargs = (;
            extrapolation_bc = (Intp.Periodic(), Intp.Flat(), Intp.Flat())
        ),
    )
    # check that the dimension transforms are correct
    @test data_handler_reversed.regridder.dim_increasing == (true, false, false)
    @test data_handler_sorted.regridder.dim_increasing == (true, true, true)
    regridded_reversed = DataHandling.regridded_snapshot(data_handler_reversed)
    regridded_sorted = DataHandling.regridded_snapshot(data_handler_sorted)
    @test regridded_reversed == regridded_sorted
    close(data_handler_reversed)
    close(data_handler_sorted)
    # netcdf file not easily correctable lon dimension (swapped two points in lon)
    incorrect_lon_path = joinpath(@__DIR__, "test_data", "incorrect_dims.nc")
    # test that informative error is thrown when lon is not monotonically increasing or decreasing
    @test_throws AssertionError DataHandling.DataHandler(
        incorrect_lon_path,
        ["test_var"],
        target_space.hybrid;
        regridder_type = :InterpolationsRegridder,
        regridder_kwargs = (;
            extrapolation_bc = (Intp.Periodic(), Intp.Flat(), Intp.Flat())
        ),
    )

end

@testset "DataHandler errors" begin
    # Create dummy file paths and variable names
    path1 = "abc"
    path2 = "def"

    var1 = "var1"
    var2 = "var2"
    var3 = "var3"

    # Construct space
    FT = Float32
    radius = FT(6731e3)
    helem = 40
    Nq = 4

    horzdomain = ClimaCore.Domains.SphereDomain(radius)
    horzmesh = ClimaCore.Meshes.EquiangularCubedSphere(horzdomain, helem)
    horztopology = ClimaCore.Topologies.Topology2D(
        context,
        horzmesh,
        ClimaCore.Topologies.spacefillingcurve(horzmesh),
    )
    quad = ClimaCore.Spaces.Quadratures.GLL{Nq}()
    target_space = ClimaCore.Spaces.SpectralElementSpace2D(horztopology, quad)

    regridder_type = :InterpolationsRegridder
    compose_function = (x, y) -> x + y

    # Providing multiple variables without a compose function is not allowed
    @test_throws ErrorException DataHandling.DataHandler(
        path1,
        [var1, var2],
        target_space;
        regridder_type,
    )

    # TempestRegridder does not support multiple input variables
    regridder_type_tr = :TempestRegridder
    @test_throws ErrorException DataHandling.DataHandler(
        path1,
        [var1, var2],
        target_space;
        regridder_type = regridder_type_tr,
        compose_function,
    )

    # TempestRegridder does not support multiple input files
    regridder_type_tr = :TempestRegridder
    @test_throws ErrorException DataHandling.DataHandler(
        [path1, path2],
        var1,
        target_space;
        regridder_type = regridder_type_tr,
        compose_function,
    )
end

@testset "DataHandler, time data" begin
    if !Sys.iswindows()
        PATH = joinpath(artifact"era5_example", "era5_t2m_sp_u10n_20210101.nc")
        varname = "sp"
        for regridder_type in (:InterpolationsRegridder, :TempestRegridder)
            for FT in (Float32, Float64)
                radius = FT(6731e3)
                helem = 40
                Nq = 4

                horzdomain = ClimaCore.Domains.SphereDomain(radius)
                horzmesh =
                    ClimaCore.Meshes.EquiangularCubedSphere(horzdomain, helem)
                horztopology =
                    ClimaCore.Topologies.Topology2D(context, horzmesh)
                quad = ClimaCore.Spaces.Quadratures.GLL{Nq}()
                target_space =
                    ClimaCore.Spaces.SpectralElementSpace2D(horztopology, quad)

                data_handler = DataHandling.DataHandler(
                    PATH,
                    varname,
                    target_space;
                    start_date = Dates.DateTime(2000, 1, 1),
                    regridder_type,
                )

                @test DataHandling.available_dates(data_handler) ==
                      data_handler.file_readers[varname].available_dates
                @test data_handler.start_date .+
                      Second.(DataHandling.available_times(data_handler)) ==
                      DataHandling.available_dates(data_handler)

                available_times = DataHandling.available_times(data_handler)
                available_dates = DataHandling.available_dates(data_handler)

                @test DataHandling.dt(data_handler) ==
                      available_times[2] - available_times[1]

                @test DataHandling.time_to_date(data_handler, 0.0) ==
                      data_handler.start_date
                @test DataHandling.time_to_date(data_handler, 1.0) ==
                      data_handler.start_date + Second(1)

                @test DataHandling.date_to_time(
                    data_handler,
                    data_handler.start_date,
                ) == 0.0
                @test DataHandling.date_to_time(
                    data_handler,
                    data_handler.start_date + Second(1),
                ) == 1.0

                # Previous time with time
                @test_throws ErrorException DataHandling.previous_time(
                    data_handler,
                    available_times[1] - 1,
                )

                @test DataHandling.previous_time(
                    data_handler,
                    available_times[10] + 1,
                ) == available_times[10]
                @test DataHandling.previous_time(
                    data_handler,
                    available_times[1] + 1,
                ) == available_times[1]

                # Previous time with time, boundaries (return the node)
                @test DataHandling.previous_time(
                    data_handler,
                    available_times[1],
                ) == available_times[1]
                @test DataHandling.previous_time(
                    data_handler,
                    available_times[end],
                ) == available_times[end]

                # Previous time with date
                @test DataHandling.previous_time(
                    data_handler,
                    available_dates[10] + Second(1),
                ) == available_times[10]
                @test DataHandling.previous_time(
                    data_handler,
                    available_dates[10],
                ) == available_times[10]

                # On node
                @test DataHandling.previous_time(
                    data_handler,
                    available_times[10],
                ) == available_times[10]
                @test DataHandling.previous_time(
                    data_handler,
                    available_times[1],
                ) == available_times[1]
                @test DataHandling.previous_time(
                    data_handler,
                    available_dates[10],
                ) == available_times[10]

                # Previous date
                @test_throws ErrorException DataHandling.previous_date(
                    data_handler,
                    available_dates[1] - Second(1),
                )

                @test DataHandling.previous_date(
                    data_handler,
                    available_dates[10] + Second(1),
                ) == available_dates[10]
                @test DataHandling.previous_date(
                    data_handler,
                    available_dates[1] + Second(1),
                ) == available_dates[1]

                # Previous time with time, boundaries (return the node)
                @test DataHandling.previous_date(
                    data_handler,
                    available_dates[1],
                ) == available_dates[1]
                @test DataHandling.previous_date(
                    data_handler,
                    available_dates[end],
                ) == available_dates[end]

                # Next time with time
                @test_throws ErrorException DataHandling.next_time(
                    data_handler,
                    available_times[end] + 1,
                )

                @test DataHandling.next_time(
                    data_handler,
                    available_times[10] + 1,
                ) == available_times[11]
                @test DataHandling.next_time(
                    data_handler,
                    available_times[1] + 1,
                ) == available_times[2]
                # Next time with date
                @test DataHandling.next_time(
                    data_handler,
                    available_dates[10] + Second(1),
                ) == available_times[11]

                # On node
                @test DataHandling.next_time(
                    data_handler,
                    available_times[10],
                ) == available_times[11]
                @test DataHandling.next_time(
                    data_handler,
                    available_dates[10],
                ) == available_times[11]

                # Next date
                @test_throws ErrorException DataHandling.next_date(
                    data_handler,
                    available_dates[end] + Second(1),
                )

                @test DataHandling.next_date(
                    data_handler,
                    available_dates[10] + Second(1),
                ) == available_dates[11]
                @test DataHandling.next_date(
                    data_handler,
                    available_dates[1] + Second(1),
                ) == available_dates[2]

                # On node
                @test DataHandling.next_date(
                    data_handler,
                    available_dates[10],
                ) == available_dates[11]
                @test DataHandling.next_date(
                    data_handler,
                    available_dates[10],
                ) == available_dates[11]

                # Asking for a regridded_snapshot without specifying the time
                @test_throws ErrorException DataHandling.regridded_snapshot(
                    data_handler,
                )

                # Asking for a regridded_snapshot with time that does not exist
                @test_throws ErrorException DataHandling.regridded_snapshot(
                    data_handler,
                    -1234.0,
                )

                field = DataHandling.regridded_snapshot(
                    data_handler,
                    available_times[10],
                )

                # Most basic testing, to get a sense that we are doing things right
                @test field isa ClimaCore.Fields.Field
                @test axes(field) == target_space
                @test eltype(field) == FT
                @test minimum(field) >= minimum(
                    data_handler.file_readers[varname].dataset[varname][
                        :,
                        :,
                        10,
                    ],
                )
                @test maximum(field) <= maximum(
                    data_handler.file_readers[varname].dataset[varname][
                        :,
                        :,
                        10,
                    ],
                )

                close(data_handler)
            end
        end



        # Testing the logic on handling multiple files/variables. We don't test
        # the numeric values here because we assume that they are correctly
        # tested above.

        radius = 6731e3
        helem = 40
        Nq = 4

        horzdomain = ClimaCore.Domains.SphereDomain(radius)
        horzmesh = ClimaCore.Meshes.EquiangularCubedSphere(horzdomain, helem)
        horztopology = ClimaCore.Topologies.Topology2D(context, horzmesh)
        quad = ClimaCore.Spaces.Quadratures.GLL{Nq}()
        target_space =
            ClimaCore.Spaces.SpectralElementSpace2D(horztopology, quad)

        # Test reading one variable with time split across multiple files

        # TODO: Move this to test/Artifacts.toml
        file_paths = [
            joinpath(@__DIR__, "test_data", "era5_1979_1.0x1.0_lai.nc"),
            joinpath(@__DIR__, "test_data", "era5_1980_1.0x1.0_lai.nc"),
        ]

        data_handler_time = DataHandling.DataHandler(
            file_paths,
            "lai_hv",
            target_space;
            regridder_type = :InterpolationsRegridder,
        )
        @test length(data_handler_time.available_dates) == 104

        # Test reading three variables with time split across multiple files
        data_handler_three = DataHandling.DataHandler(
            file_paths,
            ["lai_hv", "lai_lv", "lai_hv"],
            target_space;
            regridder_type = :InterpolationsRegridder,
            compose_function = (x, y, z) -> x + y + z,
        )

        # Test the list-of-list constructor
        data_handler_lol = DataHandling.DataHandler(
            [file_paths, file_paths, file_paths],
            ["lai_hv", "lai_lv", "lai_hv"],
            target_space;
            regridder_type = :InterpolationsRegridder,
            compose_function = (x, y, z) -> x + y + z,
        )

        @test DataHandling.regridded_snapshot(
            data_handler_lol,
            first(DataHandling.available_times(data_handler_lol)),
        ) == DataHandling.regridded_snapshot(
            data_handler_three,
            first(DataHandling.available_times(data_handler_three)),
        )

    end
end
