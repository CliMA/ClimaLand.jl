# Fetch site data
path = ClimaLand.Artifacts.esm_snowmip_data_path(; context = context)
metadata_path = joinpath(path, "site_metadata.txt")
site_metadata = readdlm(metadata_path, '\t', skipstart = 1)
# Get the latitude, longitude, start year, and end year for a given site from the metadata file
function get_loc_yr(metadata::Array{Any}, site::String)
    try
        i_site = findfirst(isequal(site), metadata[:, 1])

        if i_site === nothing
            throw(ErrorException("No metadata found for $site in $metadata_fp"))
        end

        lat, lon, start_yr, end_yr = metadata[i_site, 2],
        metadata[i_site, 3],
        metadata[i_site, 6],
        metadata[i_site, 7]

        return lat, lon, start_yr, end_yr
    catch e
        println(e)
        return nothing
    end
end
lat, long, start_year, end_year = get_loc_yr(site_metadata, SITE_NAME)
site_data_stem = "insitu_$(SITE_NAME)_$(start_year)_$(end_year).nc"
met_filename = join(["met_", site_data_stem])
obs_filename = join(["obs_", site_data_stem])


# Meteorological data
timestamp, mask, seconds, LWdown, SWdown, Psurf, q_air, Tair, u, rainf, snowf =
    NCDataset(joinpath(path, met_filename)) do met_data
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
        return timestamp,
        mask,
        seconds,
        LWdown,
        SWdown,
        Psurf,
        q_air,
        Tair,
        u,
        rainf,
        snowf
    end
# Required parameters for forcing
default_params_filepath =
    joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
toml_dict = LP.create_toml_dict(FT, default_params_filepath)
earth_param_set = LP.LandParameters(toml_dict)
thermo_params = LP.thermodynamic_parameters(earth_param_set)
# "Radiation"
SW_d = TimeVaryingInput(seconds, SWdown; context)
LW_d = TimeVaryingInput(seconds, LWdown; context)

start_date = timestamp[mask][1]

zenith_angle =
    (t, s) -> ClimaLand.default_zenith_angle(
        t,
        s;
        insol_params = earth_param_set.insol_params,
        latitude = FT(lat),
        longitude = FT(long),
    )

radiation = ClimaLand.PrescribedRadiativeFluxes(
    FT,
    SW_d,
    LW_d,
    start_date;
    θs = zenith_angle,
    earth_param_set = earth_param_set,
)


"Atmos"
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
    start_date,
    h_atmos,
    earth_param_set,
)
forcing = (; atmos, radiation)
# Snow data
albedo, z, mass, ts = NCDataset(joinpath(path, obs_filename)) do snow_data
    albedo = snow_data["albs"][:][mask]
    z = snow_data["snd_man"][:][mask]
    mass = snow_data["snw_man"][:][mask]
    ts = snow_data["ts"][:][mask]
    return albedo, z, mass, ts
end



mass_data_avail = .!(typeof.(mass) .<: Missing)
# Although snow mass data is present, other data may be missing
# We replace Missing with NaN, here, so that we can plot the data
# with gaps
fill_missing(x, FT) = typeof(x) <: Missing ? FT(NaN) : x

T_snow = fill_missing.(ts[mass_data_avail], FT)
ρ_snow = fill_missing.(mass[mass_data_avail] ./ z[mass_data_avail], FT)
depths = fill_missing.(z[mass_data_avail], FT)
SWE = depths .* ρ_snow ./ 1000.0
α = median(
    albedo[mass_data_avail][.!(typeof.(albedo[mass_data_avail]) .<: Missing)],
)
ρ = median(ρ_snow[.~isnan.(ρ_snow)])
