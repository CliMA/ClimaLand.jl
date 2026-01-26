using Test
import Dates: DateTime, Day, Date, Year

import ClimaUtilities
import ClimaUtilities.Utils: period_to_seconds_float, linear_interpolation
import ClimaUtilities: DataHandling
import ClimaUtilities: TimeVaryingInputs
import ClimaUtilities.TimeVaryingInputs: LinearPeriodFillingInterpolation

import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

import ClimaCore
import Interpolations
import NCDatasets

const context = ClimaComms.context()
ClimaComms.init(context)

include("TestTools.jl")

@testset "InterpolatingTimeVaryingInput23D with LinearPeriodFillingInterpolation" begin

    # First, test the helper functions
    _interpolable_range = Base.get_extension(
        ClimaUtilities,
        :ClimaUtilitiesClimaCoreNCDatasetsExt,
    ).TimeVaryingInputsExt._interpolable_range
    _move_date_to_period = Base.get_extension(
        ClimaUtilities,
        :ClimaUtilitiesClimaCoreNCDatasetsExt,
    ).TimeVaryingInputsExt._move_date_to_period

    dates_period_left =
        [DateTime(1985, 1, 15), DateTime(1985, 2, 17), DateTime(1985, 12, 11)]
    dates_period_right =
        [DateTime(1995, 1, 8), DateTime(1995, 2, 18), DateTime(1995, 12, 18)]
    period = Year(1)

    @test _interpolable_range(dates_period_left, dates_period_right, period) ==
          (DateTime(1985, 1, 15), DateTime(1985, 12, 11))
    @test _interpolable_range(
        [DateTime(1985, 1, 1), DateTime(1985, 12, 11)],
        dates_period_right,
        period,
    ) == (DateTime(1985, 1, 8), DateTime(1985, 12, 11))
    @test _interpolable_range(
        dates_period_left,
        [DateTime(1995, 1, 3), DateTime(1995, 12, 21)],
        period,
    ) == (DateTime(1985, 1, 15), DateTime(1985, 12, 11))

    @test_throws ErrorException _interpolable_range(
        [DateTime(1985, 1, 3), DateTime(1995, 12, 21)],
        dates_period_right,
        period,
    )
    @test_throws ErrorException _interpolable_range(
        dates_period_left,
        [DateTime(1985, 1, 3), DateTime(1995, 12, 21)],
        period,
    )

    @test _move_date_to_period(
        DateTime(1987, 3, 7),
        DateTime(1985, 1, 1),
        Year(1),
    ) == DateTime(1985, 3, 7)
    @test _move_date_to_period(
        DateTime(1987, 3, 7),
        DateTime(1995, 1, 1),
        Year(1),
    ) == DateTime(1995, 3, 7)

    # First, we create a test NetCDF files
    mktempdir() do tmpdir
        data_path = joinpath(tmpdir, "somedata.nc")

        lon_array = [-120.0, 0.0, 120]
        lat_array = [-60.0, 0.0, 60]

        # The dates are chosen so that we cover an non trivial overlap between the periods
        # of the years 1985 and 1995.
        dates = [
            DateTime(1985, 1, 10),
            DateTime(1985, 12, 13),
            DateTime(1995, 1, 6),
            DateTime(1995, 12, 11),
        ]
        start_date = DateTime(1985, 1, 1)

        times = map(d -> period_to_seconds_float(d - start_date), dates)

        ones_array = ones(Float64, (3, 3))

        # data is a Matrix (3, 3, length(times)), each time slice has values of times
        #
        # The input data is linear so linear interpolation will be perfect
        data = stack(map(t -> ones_array * t, times))

        NCDatasets.NCDataset(data_path, "c") do nc
            NCDatasets.defDim(nc, "lon", length(lon_array))
            NCDatasets.defDim(nc, "lat", length(lat_array))
            NCDatasets.defDim(nc, "time", length(times))
            NCDatasets.defVar(nc, "lon", lon_array, ("lon",))
            NCDatasets.defVar(nc, "lat", lon_array, ("lat",))
            NCDatasets.defVar(nc, "time", dates, ("time",))
            NCDatasets.defVar(nc, "data", data, ("lon", "lat", "time"))
        end

        regridder_type = :InterpolationsRegridder
        FT = Float64
        target_space = make_spherical_space(FT; context).horizontal
        dest = ClimaCore.Fields.zeros(target_space)
        one_field = ClimaCore.Fields.ones(target_space)

        data_handler = DataHandling.DataHandler(
            data_path,
            "data",
            target_space;
            start_date,
            regridder_type,
        )

        function totime(date)
            return DataHandling.date_to_time(data_handler, date)
        end

        method = LinearPeriodFillingInterpolation()
        input =
            TimeVaryingInputs.TimeVaryingInput(data_handler; method, context)

        function test_date(target_date)
            target_time = totime(target_date)
            expected = target_time .* one_field
            TimeVaryingInputs.evaluate!(dest, input, target_time)
            @test parent(dest) ≈ parent(expected)
        end

        # Outside of the domain of definition of the input data
        @test_throws ErrorException TimeVaryingInputs.evaluate!(
            dest,
            input,
            times[begin] - 1,
        )

        @test_throws ErrorException TimeVaryingInputs.evaluate!(
            dest,
            input,
            times[end] + 1,
        )

        # Target time within the leftmost period
        test_date(dates[begin] + Day(1))

        # Target time within the leftmost period, but after the last date
        test_date(DateTime(1985, 12, 14))

        # Target time within the rightmost period
        test_date(dates[end] - Day(1))

        # Target time within the rightmost period, but before the first date
        test_date(DateTime(1995, 1, 5))

        # Target time falling into the interpolable region between two periods
        test_date(DateTime(1987, 6, 18))

        # Target time falling at the left of the interpolable region for both
        test_date(DateTime(1987, 1, 1))

        # Target time falling at the right of the interpolable region for both
        test_date(DateTime(1987, 12, 31))

        # Target time falling at the left of the interpolable region for one
        test_date(DateTime(1987, 1, 8))

        # Target time falling at the right of the interpolable region for one
        test_date(DateTime(1987, 12, 12))
    end
end
