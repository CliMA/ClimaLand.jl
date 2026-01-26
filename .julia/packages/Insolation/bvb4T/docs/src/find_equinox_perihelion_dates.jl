using Plots
using Dates

using Insolation
import Insolation.Parameters as IP
import ClimaParams as CP

FT = Float64
param_set = IP.InsolationParameters(FT)

# Difference in NH and SH zenith angles at time x in given year
function zdiff(x, year, od)
    date = xtomarchdate(x, year)
    theta_s, dist =
        daily_zenith_angle(date, od, -45.0, param_set, milankovitch = true)
    theta_n, dist =
        daily_zenith_angle(date, od, 45.0, param_set, milankovitch = true)
    return theta_n - theta_s
end

# x is date relative to March 1, with 1.00 representing March 1 00:00
function xtomarchdate(x, year)
    basedate = Dates.DateTime(year, 3, 1)
    deltat = Dates.Second(round((x - 1) * IP.day(param_set)))
    return basedate + deltat
end

# Earth-Sun distance
function edist(x, year, od)
    date = xtojandate(x, year)
    _, dist = daily_zenith_angle(date, od, 0.0, param_set, milankovitch = true)
    return dist / IP.orbit_semimaj(param_set)
end

# x is date relative to Jan 1, with 1.00 representing Jan 1 00:00
function xtojandate(x, year)
    basedate = Dates.DateTime(year, 1, 1)
    deltat = Dates.Second(round((x - 1) * IP.day(param_set)))
    date = basedate + deltat
    return date
end
