# x is date relative to Jan 1, with 1.00 representing Jan 1 00:00
function xtojandate(x, year)
    basedate = Dates.DateTime(year, 1, 1)
    deltat = Dates.Second(round((x - 1) * IP.day(param_set)))
    date = basedate + deltat
    return date
end

# Earth-Sun distance
function edist(x, year, od)
    date = xtojandate(x, year)
    _, dist = daily_zenith_angle(date, od, FT(0), param_set)
    return dist / IP.orbit_semimaj(param_set)
end

years = 1900:2100
days = zeros(length(years))
od = Insolation.OrbitalData()
for (i, year) in enumerate(years)
    f = (x -> edist(x, year, od))
    res = optimize(f, 1.0, 30)
    days[i] = Optim.minimizer(res)[1]
end

# test mean is about Jan 3.5
@test mean(days) ≈ 3.5 atol = 1

# test increasing
@test mean(days[:100]) < mean(days[100:end])
