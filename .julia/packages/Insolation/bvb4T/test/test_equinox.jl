# Difference in NH and SH zenith angles at time x in given year
function zdiff(x, year, od)
    date = xtomarchdate(x, year)
    theta_s, dist = daily_zenith_angle(date, od, FT(-45), param_set)
    theta_n, dist = daily_zenith_angle(date, od, FT(45), param_set)
    return theta_n - theta_s
end

# x is date relative to March 1, with 1.00 representing March 1 00:00
function xtomarchdate(x, year)
    basedate = Dates.DateTime(year, 3, 1)
    deltat = Dates.Second(round((x - 1) * IP.day(param_set)))
    return basedate + deltat
end

od = Insolation.OrbitalData()
days = zeros(length(1900:2100))
for (i, year) in enumerate(1900:2100)
    f = (x -> zdiff(x, year, od))
    days[i] = find_zeros(f, 1.0, 30)[1]
end

# test mean is about March 21
@test mean(days) ≈ 21 atol = 1

# test decreasing
@test mean(days[:100]) > mean(days[100:end])
