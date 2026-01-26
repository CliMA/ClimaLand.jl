rtol = 1e-2
od = Insolation.OrbitalData()

# sunrise at equator
date = Dates.DateTime(2020, 2, 20, 6, 11, 0)
lon, lat = [FT(0.0), FT(0.0)]
args = (
    Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)
@test sza ≈ π / 2 rtol = rtol

# solar noon at equator
date = Dates.DateTime(2020, 2, 20, 12, 14, 0)
args = (
    Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)
@test azi ≈ 3π / 2 rtol = rtol

# sunset at equator
date = Dates.DateTime(2020, 2, 20, 18, 17, 0)
args = (
    Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)
@test sza ≈ π / 2 rtol = rtol

# solar noon at equator, eot correction = false
date = Dates.DateTime(2020, 2, 20, 12, 0, 0)
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
@test azi ≈ 3π / 2 rtol = rtol

# sunset at equator, eot correction = false
date = Dates.DateTime(2020, 2, 20, 18, 0, 0)
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
@test sza ≈ π / 2 rtol = rtol

## Test Polar Night
# polar night NH 1
date = Dates.DateTime(2020, 12, 20, 11, 0, 0)
lon, lat = [FT(0.0), FT(80.0)]
args = (
    Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)
@test sza > π / 2

# polar night NH 2
date = Dates.DateTime(2020, 12, 20, 23, 0, 0)
args = (
    Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)
@test sza > π / 2

# polar night SH 1
date = Dates.DateTime(2020, 6, 20, 11, 0, 0)
lon, lat = [FT(0.0), FT(-80.0)]
args = (
    Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)
@test sza > π / 2

# polar night SH 2
date = Dates.DateTime(2020, 6, 20, 23, 0, 0)
args = (
    Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)
@test sza > π / 2

## Test Distance
date = Dates.DateTime(2000, 3, 22, 0, 0, 0)
args = (
    Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)
@test d ≈ IP.orbit_semimaj(param_set) rtol = rtol
