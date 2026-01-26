using CFTime


# time units
@testset "CF time units" begin
    t0, plength = CFTime.timeunits("days since 1950-01-02T03:04:05Z")
    @test t0 == DateTimeStandard(1950, 1, 2, 3, 4, 5)
    @test plength == 86400000


    t0, plength = CFTime.timeunits("days since -4713-01-01T00:00:00Z")
    @test t0 == DateTimeStandard(-4713, 1, 1)
    @test plength == 86400000


    t0, plength = CFTime.timeunits("days since -4713-01-01")
    @test t0 == DateTimeStandard(-4713, 1, 1)
    @test plength == 86400000


    t0, plength = CFTime.timeunits("days since 2000-01-01 0:0:0")
    @test t0 == DateTimeStandard(2000, 1, 1)
    @test plength == 86400000

    t0, plength = CFTime.timeunits("days since 2000-01-01 00:00")
    @test t0 == DateTimeStandard(2000, 1, 1)
    @test plength == 86400000

    # issue 24
    t0, plength = CFTime.timeunits("hours since 1900-01-01 00:00:00.0")
    @test t0 == DateTimeStandard(1900, 1, 1)
    @test plength == 86400000 ÷ 24


    t0, plength = CFTime.timeunits("seconds since 1992-10-8 15:15:42.5")
    @test t0 == DateTimeStandard(1992, 10, 8, 15, 15, 42, 500)
    @test plength == 1000

    units = "microseconds since 2000-01-01T23:59:59.12345678"
    origintuple, ratio = timeunits(Tuple, units)
    @test origintuple == (2000, 1, 1, 23, 59, 59, 123, 456, 780)


    for (calendar, DT) in [
            ("standard", DateTimeStandard),
            ("gregorian", DateTimeStandard),
            ("proleptic_gregorian", DateTimeProlepticGregorian),
            ("julian", DateTimeJulian),
            ("noleap", DateTimeNoLeap),
            ("365_day", DateTimeNoLeap),
            ("all_leap", DateTimeAllLeap),
            ("366_day", DateTimeAllLeap),
            ("360_day", DateTime360Day),
        ]

        calendart0, calendarplength = CFTime.timeunits("days since 2000-1-1 0:0:0", calendar)
        @test calendart0 == DT(2000, 1, 1)
        @test calendarplength == 86400000
    end

    @test_throws ErrorException CFTime.timeunits("fortnights since 2000-01-01")
    @test_throws ErrorException CFTime.timeunits("days since 2000-1-1 0:0:0", "foo")
end

@testset "Decoding/decoding CF time data" begin

    # value from python's cftime
    # print(cftime.DatetimeJulian(-4713,1,1) + datetime.timedelta(2455512,.375 * 24*60*60))
    # 2010-10-29 09:00:00

    #@test timedecode([2455512.375],"days since -4713-01-01T00:00:00","julian")
    #   == DateTimeJulian(2010,10,29,09)

    # values from
    # https://web.archive.org/web/20180212214229/https://en.wikipedia.org/wiki/Julian_day

    # Modified JD
    @test CFTime.timedecode([58160.6875], "days since 1858-11-17", "standard") ==
        [DateTime(2018, 2, 11, 16, 30, 0)]

    # CNES JD
    @test CFTime.timedecode([24878.6875], "days since 1950-01-01", "standard") ==
        [DateTime(2018, 2, 11, 16, 30, 0)]

    # Unix time
    # wikipedia pages reports 1518366603 but it should be 1518366600
    @test CFTime.timedecode([1518366600], "seconds since 1970-01-01", "standard") ==
        [DateTime(2018, 2, 11, 16, 30, 0)]

    # The Julian Day Number (JDN) is the integer assigned to a whole solar day in
    # the Julian day count starting from noon Universal time, with Julian day
    # number 0 assigned to the day starting at noon on Monday, January 1, 4713 BC,
    # proleptic Julian calendar (November 24, 4714 BC, in the proleptic Gregorian
    # calendar),

    # Julian Day Number of 12:00 UT on January 1, 2000, is 2 451 545
    # https://web.archive.org/web/20180613200023/https://en.wikipedia.org/wiki/Julian_day


    @test CFTime.timedecode(DateTimeStandard, 2_451_545, "days since -4713-01-01T12:00:00") ==
        DateTimeStandard(2000, 01, 01, 12, 00, 00)

    # Note for DateTime, 1 BC is the year 0!
    # DateTime(1,1,1)-Dates.Day(1)
    # 0000-12-31T00:00:00

    @test CFTime.timedecode(DateTime, 2_451_545, "days since -4713-11-24T12:00:00") ==
        DateTime(2000, 01, 01, 12, 00, 00)

    if CFTime._hasyear0(CFTime.DateTimeProlepticGregorian)
        units = "days since -4713-11-24T12:00:00"
    else
        units = "days since -4714-11-24T12:00:00"
    end
    @test CFTime.timedecode(DateTimeProlepticGregorian, 2_451_545, units) ==
        DateTimeProlepticGregorian(2000, 01, 01, 12, 00, 00)


    @test CFTime.timedecode([2455512.375], "days since -4713-01-01T00:00:00", "julian", prefer_datetime = false) ==
        [DateTimeJulian(2010, 10, 29, 9, 0, 0)]

    @test CFTime.timeencode([DateTimeJulian(2010, 10, 29, 9, 0, 0)], "days since -4713-01-01T00:00:00", "julian") ==
        [2455512.375]


    @test CFTime.timedecode(DateTime, [22280.0f0], "days since 1950-01-01 00:00:00") == [DateTime(2011, 1, 1)]

    @test_throws ErrorException CFTime.timeencode(
        [DateTimeJulian(2010, 10, 29, 9, 0, 0)],
        "days since -4713-01-01T00:00:00", "360_day"
    )

    # Transition between Julian and Gregorian Calendar

    #=
    In [11]: cftime.DatetimeGregorian(1582,10,4) + datetime.timedelta(1)
    Out[11]: cftime.DatetimeGregorian(1582, 10, 15, 0, 0, 0, 0, -1, 1)

    In [12]: cftime.DatetimeProlepticGregorian(1582,10,4) + datetime.timedelta(1)
    Out[12]: cftime.DatetimeProlepticGregorian(1582, 10, 5, 0, 0, 0, 0, -1, 1)

    In [13]: cftime.DatetimeJulian(1582,10,4) + datetime.timedelta(1)
    Out[13]: cftime.DatetimeJulian(1582, 10, 5, 0, 0, 0, 0, -1, 1)
    =#

    @test DateTimeStandard(1582, 10, 4) + Dates.Day(1) == DateTimeStandard(1582, 10, 15)
    @test DateTimeProlepticGregorian(1582, 10, 4) + Dates.Day(1) == DateTimeProlepticGregorian(1582, 10, 5)
    @test DateTimeJulian(1582, 10, 4) + Dates.Day(1) == DateTimeJulian(1582, 10, 5)


    @test CFTime.datetuple(CFTime.timedecode(0, "days since -4713-01-01T12:00:00", "julian", prefer_datetime = false)) ==
        (-4713, 1, 1, 12, 0, 0, 0)


    # test sub-milliseconds

    t = CFTime.timedecode(0, "seconds since 2000-01-01", "proleptic_gregorian", prefer_datetime = true)
    @test typeof(t) <: DateTime
    @test CFTime.datetuple(t)[1:3] == (2000, 1, 1)

    t = CFTime.timedecode(1, "microseconds since 2000-01-01 00:00:00.000001", "proleptic_gregorian", prefer_datetime = false)
    @test CFTime.datetuple(t)[8] == 2

    t = CFTime.timedecode(1, "nanoseconds since 2000-01-01", "proleptic_gregorian", prefer_datetime = false)
    @test typeof(t) <: DateTimeProlepticGregorian
    @test CFTime.datetuple(t)[9] == 1

    t = CFTime.timedecode(1.0e9, "nanoseconds since 2000-01-01 00:00:00.001", "proleptic_gregorian", prefer_datetime = false)
    @test CFTime.second(t) == 1

    t = CFTime.timedecode(1.0e9, "nanoseconds since 2000-01-01 00:00:00.001", "proleptic_gregorian")
    @test CFTime.second(t) == 1
end
