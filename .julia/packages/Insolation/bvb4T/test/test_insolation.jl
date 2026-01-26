atol = 1e-4
rtol = 1e-2
rtol_insol = 0.1

od = Insolation.OrbitalData()
## Test zero insolation at night
# sunrise at equator
date = Dates.DateTime(2020, 1, 1, 6, 0, 0)
lon, lat = [FT(0.0), FT(0.0)]
args = (
    Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)
F = insolation(sza, d, param_set)
@test F ≈ 0.0 atol = atol

S, mu = solar_flux_and_cos_sza(date, od, lon, lat, param_set)
@test S ≈ IP.tot_solar_irrad(param_set) rtol = rtol_insol
@test mu ≈ 0.0 atol = atol

# polar night NH 1
date = Dates.DateTime(2020, 12, 20, 11, 0, 0)
lon, lat = [FT(0.0), FT(80.0)]
args = (
    Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)
F = insolation(sza, d, param_set)
@test F ≈ 0.0 atol = atol

S, mu = solar_flux_and_cos_sza(date, od, lon, lat, param_set)
@test S ≈ IP.tot_solar_irrad(param_set) rtol = rtol_insol
@test mu ≈ 0.0 atol = atol

# polar night NH 2
date = Dates.DateTime(2020, 12, 20, 23, 0, 0)
args = (
    Insolation.helper_instantaneous_zenith_angle(date, od, param_set)...,
    lon,
    lat,
)
sza, azi, d = instantaneous_zenith_angle(args...)
F = insolation(sza, d, param_set)
@test F ≈ 0.0 atol = atol

S, mu = solar_flux_and_cos_sza(date, od, lon, lat, param_set)
@test S ≈ IP.tot_solar_irrad(param_set) rtol = rtol_insol
@test mu ≈ 0.0 atol = atol

## Test symmetry of insolation at equinox
nlats = 181
date = Dates.DateTime(2021, 3, 20, 9, 37, 0) # vernal equinox 2021
l_arr = FT.(collect(range(-90, stop = 90, length = nlats)))
F_arr = zeros(nlats)
for (i, lat) in enumerate(l_arr)
    θ, dist = daily_zenith_angle(date, od, lat, param_set)
    F_arr[i] = insolation(θ, dist, param_set)
end
F_NH = sort(F_arr[l_arr .>= 0])
F_SH = sort(F_arr[l_arr .<= 0])
@test F_NH ≈ F_SH rtol = rtol

## Test globally averaged insolation ≈ TSI
ndays, nlats = [365, 361]
d_arr = Array{Int}(round.(collect(range(0, stop = 365, length = ndays))))
l_arr = FT.(collect(range(-90, stop = 90, length = nlats)))
F_arr = zeros(ndays, nlats)

for (i, d) in enumerate(d_arr)
    for (j, lat) in enumerate(l_arr)
        datei = Dates.DateTime(2000, 1, 1) + Dates.Day(d)
        θ, dist = daily_zenith_angle(datei, od, lat, param_set)
        F_arr[i, j] = insolation(θ, dist, param_set)
    end
end

zonal_mean_insol = mean(F_arr, dims = 1)
area_fac = abs.(cosd.(l_arr))
global_mean_insol = sum(zonal_mean_insol * area_fac) / sum(area_fac)
@test global_mean_insol ≈ IP.tot_solar_irrad(param_set) / 4 rtol = rtol

## Test invariance of zonal-mean insolation under rotation of ϖ
ϖ0 = IP.lon_perihelion_epoch(param_set)
param_set = IP.InsolationParameters(FT, (; lon_perihelion_epoch = ϖ0 + π))

for (i, d) in enumerate(d_arr)
    for (j, lat) in enumerate(l_arr)
        datei = Dates.DateTime(2000, 1, 1) + Dates.Day(d)
        θ, dist = daily_zenith_angle(datei, od, lat, param_set)
        F_arr[i, j] = insolation(θ, dist, param_set)
    end
end

zonal_mean_insol_rotate = mean(F_arr, dims = 1)
@test zonal_mean_insol_rotate ≈ zonal_mean_insol rtol = rtol

param_set = IP.InsolationParameters(FT, (; lon_perihelion_epoch = ϖ0))
