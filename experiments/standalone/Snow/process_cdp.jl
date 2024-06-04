using NCDatasets
using Dates
using ClimaLand: PrescribedAtmosphere, PrescribedRadiativeFluxes
import ClimaLand.Parameters as LP
using Thermodynamics
using Statistics
using Insolation

# Data
met_data_path, snow_data_path = ClimaLand.Artifacts.esm_snowmip_data_path()
met_data = NCDataset(met_data_path)
snow_data = NCDataset(snow_data_path)

# Process Met Data
timestamp = met_data["time"][:]
year = Dates.year.(timestamp)
# convert from milliseconds to seconds
mask = (year .== 2009) .| (year .== 2010) .| (year .== 2011)
seconds = Dates.value.(timestamp[mask] .- timestamp[mask][1]) / 1000


LWdown = met_data["LWdown"][:][mask]
SWdown = met_data["SWdown"][:][mask]
Psurf = met_data["Psurf"][:][mask]
q_air = met_data["Qair"][:][mask]
Tair = met_data["Tair"][:][mask]
u = met_data["Wind"][:][mask]
# convert mass flux into volume flux of liquid water
# additionally, our convention is that a negative flux is downwards
rainf = -met_data["Rainf"][:][mask] ./ 1000.0 # convert to dS/dt
snowf = -met_data["Snowf"][:][mask] ./ 1000.0 # convert to dS/dt

# "Radiation"
SW_d = TimeVaryingInput(seconds, SWdown; context)
LW_d = TimeVaryingInput(seconds, LWdown; context)
# double check!
lat = FT(45.28)
long = FT(5.77)
ref_time = DateTime("2009-01-01-00", "yyyy-mm-dd-HH")
function zenith_angle(
    t,
    ref_time;
    latitude = lat,
    longitude = long,
    insol_params::Insolation.Parameters.InsolationParameters{FT} = param_set.insol_params,
) where {FT}
    # This should be time in UTC
    current_datetime = ref_time + Dates.Second(round(t))

    # Orbital Data uses Float64, so we need to convert to our sim FT
    d, δ, η_UTC =
        FT.(
            Insolation.helper_instantaneous_zenith_angle(
                current_datetime,
                ref_time,
                insol_params,
            )
        )

    FT(
        Insolation.instantaneous_zenith_angle(
            d,
            δ,
            η_UTC,
            longitude,
            latitude,
        )[1],
    )
end

radiation = ClimaLand.PrescribedRadiativeFluxes(
    FT,
    SW_d,
    LW_d,
    ref_time;
    θs = zenith_angle,
)


"Atmos"
earth_param_set = LP.LandParameters(FT)
liquid_precip = TimeVaryingInput(seconds, rainf; context)
snow_precip = TimeVaryingInput(seconds, snowf; context)
T_atmos = TimeVaryingInput(seconds, Tair; context)
u_atmos = TimeVaryingInput(seconds, u; context)
q_atmos = TimeVaryingInput(seconds, q_air; context)
# What is this value?
h_atmos = FT(3)
P_atmos = TimeVaryingInput(seconds, Psurf; context)
atmos = PrescribedAtmosphere(
    liquid_precip,
    snow_precip,
    T_atmos,
    u_atmos,
    q_atmos,
    P_atmos,
    ref_time,
    h_atmos,
    earth_param_set,
)

# Process Snow Data
albedo = snow_data["albs"][:][mask]
z = snow_data["snd_man"][:][mask]
mass = snow_data["snw_man"][:][mask]

snow_data_avail = .!(typeof.(mass) .<: Missing)
T_snow = snow_data["ts"][:][mask][snow_data_avail]
ρ_snow = mass[snow_data_avail] ./ z[snow_data_avail]
SWE = z[snow_data_avail] .* ρ_snow ./ 1000.0
α = median(
    albedo[snow_data_avail][.!(typeof.(albedo[snow_data_avail]) .<: Missing)],
)
ρ = median(ρ_snow[.~isnan.(ρ_snow)])
