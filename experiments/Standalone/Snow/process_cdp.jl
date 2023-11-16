using NCDatasets
using ArtifactWrappers
using Dierckx
using Dates
using ClimaLSM: PrescribedAtmosphere, PrescribedRadiativeFluxes
using Thermodynamics
using Statistics
using Insolation

# Meteorological data
af = ArtifactFile(
    url = "https://caltech.box.com/shared/static/q8ju0yy782j3f3jjrqowakpm5s3qf40v.nc",
    filename = "met_insitu_cdp_1994_2014.nc",
)
dataset = ArtifactWrapper(@__DIR__, "driver", ArtifactFile[af]);
dataset_path = get_data_folder(dataset);
met_data = NCDataset(joinpath(dataset_path, af.filename))

timestamp = met_data["time"][:]
year = Dates.year.(timestamp)
# convert from milliseconds to seconds
mask = (year .== 2009) .| (year .== 2010) .| (year .== 2011)
seconds = Dates.value.(timestamp[mask] .- timestamp[mask][1]) / 1000


LWdown = met_data["LWdown"][mask]
SWdown = met_data["SWdown"][mask]
Psurf = met_data["Psurf"][mask]
q_air = met_data["Qair"][mask]
Tair = met_data["Tair"][mask]
u = met_data["Wind"][mask]
# convert mass flux into volume flux of liquid water
# additionally, our convention is that a negative flux is downwards
rainf = -met_data["Rainf"][mask] ./ 1000.0 # convert to dS/dt
snowf = -met_data["Snowf"][mask] ./ 1000.0 # convert to dS/dt

# convert to spline functions
"Radiation"
SW_d = (t) -> eltype(t)(Spline1D(seconds, SWdown)(t))
LW_d = (t) -> eltype(t)(Spline1D(seconds, LWdown)(t))
# double check!
lat = FT(45.28)
long = FT(5.77)
ref_time = DateTime("2009-01-01-00", "yyyy-mm-dd-HH")
function zenith_angle(
    t::FT,
    orbital_data,
    ref_time;
    latitude = lat,
    longitude = long,
    insol_params = earth_param_set.insol_params,
) where {FT}
    # This should be time in UTC - double check!
    dt = ref_time + Dates.Second(round(t))
    FT(
        instantaneous_zenith_angle(
            dt,
            orbital_data,
            longitude,
            latitude,
            insol_params,
        )[1],
    )
end

rad = PrescribedRadiativeFluxes(
    FT,
    SW_d,
    LW_d,
    ref_time;
    θs = zenith_angle,
    orbital_data = Insolation.OrbitalData(),
)

"Atmos"
liquid_precip = (t) -> eltype(t)(Spline1D(seconds, rainf)(t))
snow_precip = (t) -> eltype(t)(Spline1D(seconds, snowf)(t))
T_atmos = (t) -> eltype(t)(Spline1D(seconds, Tair)(t))
u_atmos = (t) -> eltype(t)(Spline1D(seconds, u)(t))
q_atmos = (t) -> eltype(t)(Spline1D(seconds, q_air)(t))
# What is this value?
h_atmos = FT(3)
P_atmos = (t) -> eltype(t)(Spline1D(seconds, Psurf)(t))
atmos = PrescribedAtmosphere(
    liquid_precip,
    snow_precip,
    T_atmos,
    u_atmos,
    q_atmos,
    P_atmos,
    ref_time,
    h_atmos,
)
# Meteorological data
snow_af = ArtifactFile(
    url = "https://caltech.box.com/shared/static/umwr6wqjhsz3zq6acht68m4d28a64ueb.nc",
    filename = "obs_insitu_cdp_1994_2014.nc",
)
snow_dataset = ArtifactWrapper(@__DIR__, "snow_driver", ArtifactFile[snow_af]);
snow_dataset_path = get_data_folder(snow_dataset);
snow_data = NCDataset(joinpath(snow_dataset_path, snow_af.filename))

albedo = snow_data["albs"][mask]
z = snow_data["snd_man"][mask]
mass = snow_data["snw_man"][mask]

snow_data_avail = .!(typeof.(mass) .<: Missing)
T_snow = snow_data["ts"][mask][snow_data_avail]
ρ_snow = mass[snow_data_avail] ./ z[snow_data_avail]
SWE = z[snow_data_avail] .* ρ_snow ./ 1000.0
α = median(
    albedo[snow_data_avail][.!(typeof.(albedo[snow_data_avail]) .<: Missing)],
)
ρ = median(ρ_snow[.~isnan.(ρ_snow)])
Δt = FT(60 * 60 * 6)
parameters =
    SnowParameters{FT}(Δt; α_snow = α, ρ_snow = ρ, earth_param_set = param_set)
S_0 = FT(SWE[1])
q_l_0 = FT(0)
U_0 = ClimaLSM.Snow.energy_from_q_l_and_swe(S_0, q_l_0, parameters)
