date = Dates.DateTime(2020, 2, 20, 11, 11, 0)
lon, lat = [FT(80.0), FT(20.0)]

od = Insolation.OrbitalData()

args = (
    Insolation.helper_instantaneous_zenith_angle(
        date,
        param_set,
        eot_correction = false,
    )...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)

F = insolation(sza, d, param_set)
@test typeof(sza) == FT
@test typeof(azi) == FT
@test typeof(d) == FT
@test typeof(F) == FT

args = (
    Insolation.helper_instantaneous_zenith_angle(
        date,
        od,
        param_set,
        eot_correction = false,
    )...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)

F = insolation(sza, d, param_set)
@test typeof(sza) == FT
@test typeof(azi) == FT
@test typeof(d) == FT
@test typeof(F) == FT

args = (
    Insolation.helper_instantaneous_zenith_angle(
        date,
        param_set,
        eot_correction = true,
    )...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)

F = insolation(sza, d, param_set)
@test typeof(sza) == FT
@test typeof(azi) == FT
@test typeof(d) == FT
@test typeof(F) == FT

args = (
    Insolation.helper_instantaneous_zenith_angle(
        date,
        od,
        param_set,
        eot_correction = true,
    )...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)

F = insolation(sza, d, param_set)
@test typeof(sza) == FT
@test typeof(azi) == FT
@test typeof(d) == FT
@test typeof(F) == FT
