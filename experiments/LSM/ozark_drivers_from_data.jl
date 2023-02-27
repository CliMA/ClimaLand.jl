using ArtifactWrappers
using DelimitedFiles
using Dierckx
using Thermodynamics

# Fluxnet driver data (precipitation, ET) and test data for comparison (soil moisture, leaf water potential)
af = ArtifactFile(
    url = "https://caltech.box.com/shared/static/3d03opwwczvxia0fks38uxr62qne57f0.csv",
    filename = "Ozark_test_param_set.csv",
)
dataset = ArtifactWrapper(@__DIR__, "driver", ArtifactFile[af]);
dataset_path = get_data_folder(dataset);
data = joinpath(dataset_path, "Ozark_test_param_set.csv")
driver_data = readdlm(data, ',', skipstart = 1)

# Natan's model output data for comparison
af2 = ArtifactFile(
    url = "https://caltech.box.com/shared/static/qkm7lrqaphlehe6p5hvccldhie6gxo0v.csv",
    filename = "holtzman_clima_output_april1.csv",
)
dataset2 = ArtifactWrapper(@__DIR__, "Natan", ArtifactFile[af2]);
dataset_path2 = get_data_folder(dataset2);
data2 = joinpath(dataset_path2, "holtzman_clima_output_april1.csv")
Natan_data = readdlm(data2, ',', skipstart = 1)

# 2005 data
precip_flux = driver_data[17569:2:35088, 5] # m/s
et_flux = driver_data[17569:2:35088, 7] # m/s
tree_species = driver_data[125:440, 15] # names of species
observed_predawn_lwp_2005 = driver_data[125:440, 16] # MPa
white_oak = observed_predawn_lwp_2005[tree_species .== "white oak"] # MPa
black_oak = observed_predawn_lwp_2005[tree_species .== "black oak"] # MPa
eastern_redcedar =
    observed_predawn_lwp_2005[tree_species .== "eastern redcedar"] # MPa
shagbark_hickory =
    observed_predawn_lwp_2005[tree_species .== "shagbark hickory"] # MPa
sugar_maple = observed_predawn_lwp_2005[tree_species .== "sugar maple"] # MPa
white_ash = observed_predawn_lwp_2005[tree_species .== "white ash"] # MPa
day_of_year_of_observed_predawn_lwp_2005 = driver_data[125:440, 14] # unitless, number of day in year
tree_observation_date_white_oak =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "white oak"] # MPa
tree_observation_date_black_oak =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "black oak"] # MPa
tree_observation_date_eastern_redcedar =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "eastern redcedar"] # MPa
tree_observation_date_shagbark_hickory =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "shagbark hickory"] # MPa
tree_observation_date_sugar_maple =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "sugar maple"] # MPa
tree_observation_date_white_ash =
    day_of_year_of_observed_predawn_lwp_2005[tree_species .== "white ash"] # MPa
unique_tree_observation_date_white_oak = unique(tree_observation_date_white_oak) #  unitless, number of day in year
mean_white_oak = [
    mean(
        white_oak[tree_observation_date_white_oak .== unique_tree_observation_date_white_oak[k]],
    ) for k in 1:1:length(unique_tree_observation_date_white_oak)
]

t = Array(0:(60 * 60):((length(precip_flux) - 1) * 60 * 60)) # s
p_spline = Spline1D(t, -precip_flux) # m/s
et_spline = Spline1D(t, et_flux) # m/s

# Natan's data
natan_et = Natan_data[2:end, 15] .* 18 / 1e3 / 1e3 # m^3/m^2/s
swc_column = Natan_data[2:end, 20] # m^3/m^3
swc_surface = Natan_data[2:end, 19] # m^3/m^3
swc_obs_surface = Natan_data[2:end, 12] # m^3/m^3
lwp = Natan_data[2:end, 18] # MPa
dates = Natan_data[2:end, 1]
RH = Natan_data[2:end, 6]
T_air = Natan_data[2:end, 4] .+ 273.15
dates_julia = tryparse.(DateTime, dates)
our_year = dates_julia[Dates.year.(dates_julia) .== 2005]
seconds = Dates.value.(our_year .- our_year[1]) ./ 1000
our_year_swc_column = FT.(swc_column[Dates.year.(dates_julia) .== 2005])
our_year_swc_surface = FT.(swc_surface[Dates.year.(dates_julia) .== 2005])
our_year_swc_obs_surface =
    FT.(swc_obs_surface[Dates.year.(dates_julia) .== 2005])
our_year_lwp = FT.(lwp[Dates.year.(dates_julia) .== 2005])
our_year_et = FT.(natan_et[Dates.year.(dates_julia) .== 2005])
our_year_T = FT.(T_air[Dates.year.(dates_julia) .== 2005])
thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
esat =
    Thermodynamics.saturation_vapor_pressure.(
        Ref(thermo_params),
        our_year_T,
        Ref(Thermodynamics.Liquid()),
    )
e = @. RH[Dates.year.(dates_julia) .== 2005] * esat * 100
q = @. 0.622 * 2 ./ (101325 - 0.378 * e)
q_spline = Spline1D(seconds, q)
natan_et_spline = Spline1D(seconds, our_year_et)
T_air_spline = Spline1D(seconds, our_year_T)

precipitation_function(t::FT) where {FT} = p_spline(t) < 0.0 ? p_spline(t) : 0.0 # m/s
snow_precip(t) = eltype(t)(0)
atmos_T(t) = T_air_spline(t)
atmos_q(t) = q_spline(t)
transpiration_function(t) = et_spline(t)
# These are guesses
atmos_u(t) = eltype(t)(3)
atmos_p(t) = eltype(t)(101325)
atmos_h = FT(40)
atmos_co2 = (t) -> eltype(t)(4.11e-4) # mol/mol



lat = FT(38.7441) # degree
long = FT(-92.2000) # degree

function zenith_angle(
    t::FT;
    latitude = lat,
    longitude = long,
    insol_params = earth_param_set.insol_params,
) where {FT}
    dt = DateTime("2005-01-01", "yyyy-mm-dd") + Dates.Second(t)
    return FT(
        instantaneous_zenith_angle(dt, longitude, latitude, insol_params)[1],
    )
end

# We need to get this from the data!
function shortwave_radiation(
    t::FT;
    latitude = lat,
    longitude = long,
    insol_params = earth_param_set.insol_params,
) where {FT}
    dt = DateTime("2005-01-01", "yyyy-mm-dd") + Dates.Second(t)
    S, μ =
        Insolation.solar_flux_and_cos_sza(dt, longitude, latitude, insol_params)
    return S * μ #W/m^2
end

function longwave_radiation(
    t::FT;
    latitude = lat,
    longitude = long,
    insol_params = earth_param_set.insol_params,
) where {FT}
    dt = DateTime("2005-01-01", "yyyy-mm-dd") + Dates.Second(t)
    S, μ =
        Insolation.solar_flux_and_cos_sza(dt, longitude, latitude, insol_params)
    return S * μ * 1 / FT(4) # W/m^2
end


atmos = ClimaLSM.PrescribedAtmosphere(
    precipitation_function,
    snow_precip,
    atmos_T,
    atmos_u,
    atmos_q,
    atmos_p,
    atmos_h;
    c_co2 = atmos_co2,
)
radiation = ClimaLSM.PrescribedRadiativeFluxes(
    FT,
    shortwave_radiation,
    longwave_radiation;
    θs = zenith_angle,
)

transpiration = PrescribedTranspiration{FT}(transpiration_function)
