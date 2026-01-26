using Test

using Dates

import ClimaUtilities.Utils:
    searchsortednearest,
    linear_interpolation,
    isequispaced,
    wrap_time,
    beginningofperiod,
    endofperiod,
    bounding_dates,
    period_to_seconds_float,
    unique_periods,
    sort_by_creation_time

@testset "searchsortednearest" begin
    A = 10 * collect(range(1, 10))

    @test searchsortednearest(A, 0) == 1
    @test searchsortednearest(A, 1000) == 10
    @test searchsortednearest(A, 20) == 2
    @test searchsortednearest(A, 21) == 2
    @test searchsortednearest(A, 29) == 3
end

@testset "linearinterpolation" begin
    @test linear_interpolation([1.0, 2.0, 3.0], [2.0, 4.0, 6.0], 1.5) ≈ 3.0
    @test linear_interpolation([1.0, 2.0, 3.0], [2.0, 4.0, 6.0], 2.5) ≈ 5.0

    # Test first element
    @test linear_interpolation([1.0, 2.0, 3.0], [2.0, 4.0, 6.0], 1.0) ≈ 2.0
    # Test last element
    @test linear_interpolation([1.0, 2.0, 3.0], [2.0, 4.0, 6.0], 3.0) ≈ 6.0
    # Test with different starting value
    @test linear_interpolation([0.0, 1.0, 2.0], [1.0, 4.0, 9.0], 0.5) ≈ 2.5
end

@testset "isequispaced" begin
    @test isequispaced([1, 2, 3, 4, 5])
    @test !isequispaced([1, 2, 4, 8, 16])
    @test isequispaced([1.0, 2.0, 3.0, 4.0, 5.0])
    @test !isequispaced([1.0, 2.0, 3.1, 4.0, 5.0])
end

@testset "wrap_time" begin
    t_init = 10
    t_end = 20
    dt = 1

    @test wrap_time(15, t_init, t_end) == 15
    @test wrap_time(25, t_init, t_end) == 15
    @test wrap_time(5, t_init, t_end) == 15

    # Wrapping when time equals t_init or t_end
    @test wrap_time(10, t_init, t_end) == 10
    @test wrap_time(20, t_init, t_end) == 10

    date_init = Dates.DateTime(2010) + Dates.Second(10)
    date_end = Dates.DateTime(2010) + Dates.Second(20)
    dt = Dates.Second(1)

    @test wrap_time(date_init + 5dt, date_init, date_end) == date_init + 5dt
    @test wrap_time(date_end + 5dt, date_init, date_end) == date_init + 5dt
    @test wrap_time(date_init - 5dt, date_init, date_end) == date_init + 5dt

    # Wrapping when time equals date_init or date_end
    @test wrap_time(date_init, date_init, date_end) == date_init
    @test wrap_time(date_end, date_init, date_end) == date_init
end

@testset "bounding_dates" begin
    @test beginningofperiod(Date(1993, 11, 19), Year(1)) == DateTime(1993, 1, 1)
    @test beginningofperiod(Date(1993, 11, 19), Month(1)) ==
          DateTime(1993, 11, 1)
    @test beginningofperiod(Date(1993, 11, 19), Week(1)) ==
          DateTime(1993, 11, 15)
    @test beginningofperiod(Date(1993, 11, 19), Day(1)) ==
          DateTime(1993, 11, 19)
    @test beginningofperiod(DateTime(1993, 11, 19, 1), Day(1)) ==
          DateTime(1993, 11, 19)

    @test endofperiod(Date(1993, 11, 19), Year(1)) ==
          DateTime(1993, 12, 31, 23, 59, 59)
    @test endofperiod(Date(1993, 11, 19), Month(1)) ==
          DateTime(1993, 11, 30, 23, 59, 59)
    @test endofperiod(Date(1993, 11, 19), Week(1)) ==
          DateTime(1993, 11, 21, 23, 59, 59)
    @test endofperiod(Date(1993, 11, 19), Day(1)) ==
          DateTime(1993, 11, 19, 23, 59, 59)

    dates = [
        Date(1993, 8, 13),
        Date(1993, 8, 18),
        Date(1993, 11, 19),
        Date(1994, 1, 1),
        Date(1998, 1, 17),
    ]

    @test_throws ErrorException bounding_dates(dates, Date(2000, 1, 1), Year(1))
    @test bounding_dates(dates, Date(1993, 8, 15), Year(1)) ==
          (Date(1993, 08, 13), Date(1993, 11, 19))
    @test bounding_dates(dates, Date(1993, 8, 15), Month(1)) ==
          (Date(1993, 08, 13), Date(1993, 08, 18))
end

@testset "period_to_seconds_float" begin
    @test period_to_seconds_float(Millisecond(1)) == 0.001
    @test period_to_seconds_float(Second(1)) == 1.0
    @test period_to_seconds_float(Minute(1)) == 60.0
    @test period_to_seconds_float(Hour(1)) == 3600.0
    @test period_to_seconds_float(Day(1)) == 86400.0
    @test period_to_seconds_float(Week(1)) == 604800.0
    @test period_to_seconds_float(Month(1)) == 2.629746e6
    @test period_to_seconds_float(Year(1)) == 3.1556952e7
end

@testset "unique_periods" begin
    dates = [
        DateTime(2022, 1, 1),
        DateTime(2022, 5, 15),
        DateTime(2023, 1, 1),
        DateTime(2023, 5, 5),
        DateTime(2024, 10, 10),
    ]

    # Test with non simple period
    @test_throws ErrorException unique_periods(dates, Year(2))

    @test unique_periods(dates, Year(1)) ==
          [DateTime(2022, 1, 1), DateTime(2023, 1, 1), DateTime(2024, 1, 1)]

    @test unique_periods(dates, Month(1)) == [
        DateTime(2022, 1, 1),
        DateTime(2022, 5, 1),
        DateTime(2023, 1, 1),
        DateTime(2023, 5, 1),
        DateTime(2024, 10, 1),
    ]

    @test unique_periods(dates, Week(1)) == [
        DateTime(2021, 12, 27),
        DateTime(2022, 5, 9),
        DateTime(2022, 12, 26),
        DateTime(2023, 5, 1),
        DateTime(2024, 10, 7),
    ]

    @test unique_periods(dates, Day(1)) == [
        DateTime(2022, 1, 1),
        DateTime(2022, 5, 15),
        DateTime(2023, 1, 1),
        DateTime(2023, 5, 5),
        DateTime(2024, 10, 10),
    ]
end

@testset "sort_by_creation_time" begin
    mktempdir() do basedir
        files = map(
            f -> joinpath(basedir, f),
            ["file3.txt", "file1.txt", "file2.txt"];
        )
        touch(files[2])
        sleep(0.1)
        touch(files[3])
        sleep(0.1)
        touch(files[1])

        sorted_files = sort_by_creation_time(files)
        @test sorted_files == [files[2], files[3], files[1]]
    end
end
