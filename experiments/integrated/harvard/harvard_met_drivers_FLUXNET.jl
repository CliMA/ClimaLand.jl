using ArtifactWrappers
using DelimitedFiles
using Dierckx
using Thermodynamics
using Dates

function replace_missing_with_mean!(field, flag)
    good_indices = (flag .== 0) .|| (flag .== 1)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end

# Fluxnet Harvard Forest US-Ha1 (CO2 and H2O fluxes and met drivers)

# Save the data file as an artifact
af = ArtifactFile(
    url = "https://caltech.box.com/shared/static/xixaod6511cutz51ag81k1mtvy05hbol.csv",
    filename = "combined_US-Ha1-2010.csv",
)
dataset = ArtifactWrapper(@__DIR__, "ameriflux_data_Ha1", ArtifactFile[af]);
dataset_path = get_data_folder(dataset);
data = joinpath(dataset_path, "combined_US-Ha1-2010.csv");
driver_data = readdlm(data, ',')

# Read the data in each column. For the Harvard site, we have hourly data
column_names = driver_data[1, :]
TA = driver_data[2:end, column_names .== "TA_F"] .+ 273.15; # convert C to K
VPD = driver_data[2:end, column_names .== "VPD_F"] .* 100; # convert hPa to Pa
PA = driver_data[2:end, column_names .== "PA_F"] .* 1000; # convert kPa to Pa
P = driver_data[2:end, column_names .== "P_F"] ./ (1000 * 3600); # convert mm/HH to m/s
WS = driver_data[2:end, column_names .== "WS_F"]; # already m/s
LW_IN = driver_data[2:end, column_names .== "LW_IN_F"]
SW_IN = driver_data[2:end, column_names .== "SW_IN_F"]
CO2_F = driver_data[2:end, column_names .== "CO2_F_MDS_QC"]
CO2 = driver_data[2:end, column_names .== "CO2_F_MDS"] .* 1e-6; # convert \mumol to mol
replace_missing_with_mean!(CO2, CO2_F)
# TODO: Find SWC data somewhere
TS_F = driver_data[2:end, column_names .== "TS_F_MDS_2_QC"] # Most likely 5cm depth TODO: (switched to 2, check for dependencies)
TS = driver_data[2:end, column_names .== "TS_F_MDS_2"] .+ 273.15;# convert C to K
replace_missing_with_mean!(TS, TS_F)
LE = driver_data[2:end, column_names .== "LE_F_MDS"]
GPP = driver_data[2:end, column_names .== "GPP_DT_VUT_REF"] .* 1e-6 # to convert from micromol to mol.
H = driver_data[2:end, column_names .== "H_F_MDS"]
LAI_data = driver_data[2:end, column_names .== "value_mean"]

# Get the data timestamps both in local and UTC times
LOCAL_DATETIME = DateTime.(string.(driver_data[2:end, 1]), "yyyymmddHHMM")
UTC_DATETIME = LOCAL_DATETIME .+ Dates.Hour(5)
thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
esat =
    Thermodynamics.saturation_vapor_pressure.(
        Ref(thermo_params),
        TA,
        Ref(Thermodynamics.Liquid()),
    )

# Compute atmospheric pressures
e = @. esat - VPD
q = @. 0.622 * e ./ (PA - 0.378 * e)

#Make a bunch of splines
seconds = FT.(0:3600:((length(UTC_DATETIME) - 1) * 3600));
p_spline = Spline1D(seconds, -P[:]) # m/s
atmos_q = Spline1D(seconds, q[:])
atmos_T = Spline1D(seconds, TA[:])
atmos_p = Spline1D(seconds, PA[:])
atmos_co2 = Spline1D(seconds, CO2[:])
atmos_u = Spline1D(seconds, WS[:])
LW_IN_spline = Spline1D(seconds, LW_IN[:])
SW_IN_spline = Spline1D(seconds, SW_IN[:])
atmos_h = FT(32)
precipitation_function(t::FT) where {FT} = p_spline(t) < 0.0 ? p_spline(t) : 0.0 # m/s
snow_precip(t) = eltype(t)(0) # this is likely not correct
LAIspline = Spline1D(seconds, LAI_data[:])
LAIfunction = (t) -> eltype(t)(LAIspline(t))


# Construct the drivers
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

# Site latitude and longitude data
lat = FT(42.5378) # degree
long = FT(-72.1715) # degree

# Find the function for solar zenith angle vs. time for this site
function zenith_angle(
    t::FT,
    orbital_data;
    latitude = lat,
    longitude = long,
    insol_params = earth_param_set.insol_params,
) where {FT}
    # This should be time in UTC
    dt = DateTime("2005-01-01-05", "yyyy-mm-dd-HH") + Dates.Second(round(t))
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

# Prescribed radiation function for this site
radiation = ClimaLSM.PrescribedRadiativeFluxes(
    FT,
    SW_IN_spline,
    LW_IN_spline;
    θs = zenith_angle,
    orbital_data = Insolation.OrbitalData(),
)
