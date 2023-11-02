## Dataset obtained from Ameriflux. Funding for the AmeriFlux data portal was provided by the U.S. Department of Energy Office of Science.
## Missouri Ozark site
## Citation: Jeffrey Wood, Lianhong Gu (2021),
## AmeriFlux FLUXNET-1F US-MOz Missouri Ozark Site, Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1854370


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

function replace_missing_with_mean_by_value!(field)
    good_indices = .~(field .== -9999)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end

# Fluxnet Ozark (CO2 and H2O fluxes and met drivers)
# 2005 data extracted from the full Fluxnet data for the site as follows:
#af = ArtifactFile(
#    url = "https://caltech.box.com/shared/static/cy4jlul43kx72r2pqthvm4isjwkatgxy.csv",
#    filename = "AMF_US-MOz_FLUXNET_FULLSET_HH_2004-2019_3-5.csv",
#)
#dataset = ArtifactWrapper(@__DIR__, "ameriflux_data_all_years", ArtifactFile[af]);
#dataset_path = get_data_folder(dataset);
#data = joinpath(dataset_path, "AMF_US-MOz_FLUXNET_FULLSET_HH_2004-2019_3-5.csv")
#driver_data = readdlm(data, ',')
#mask = Dates.year.(LOCAL_DATETIME) .== 2005
#data_from_2005 = vcat(driver_data[[1],:], driver_data[2:end, :][mask, :])
#open("AMF_US-MOz_FLUXNET_FULLSET_HH_2005.csv"), "w") do io
#             writedlm(io, data_from_2005, ',')
#         end

af = ArtifactFile(
    url = "https://caltech.box.com/shared/static/1uwg8rjg2wx7y0vp8j9kv2d44y3fajyk.csv",
    filename = "AMF_US-MOz_FLUXNET_FULLSET_HH_2005.csv",
)
dataset = ArtifactWrapper(@__DIR__, "ameriflux_data-US-MOz", ArtifactFile[af]);
dataset_path = get_data_folder(dataset);
data = joinpath(dataset_path, af.filename)
driver_data = readdlm(data, ',')

column_names = driver_data[1, :]
TA = driver_data[2:end, column_names .== "TA_F"] .+ 273.15; # convert C to K
VPD = driver_data[2:end, column_names .== "VPD_F"] .* 100; # convert hPa to Pa
PA = driver_data[2:end, column_names .== "PA_F"] .* 1000; # convert kPa to Pa
P = driver_data[2:end, column_names .== "P_F"] ./ (1000 * 1800); # convert mm/HH to m/s
WS = driver_data[2:end, column_names .== "WS_F"]; # already m/s
LW_IN = driver_data[2:end, column_names .== "LW_IN_F"]
SW_IN = driver_data[2:end, column_names .== "SW_IN_F"]
CO2_F = driver_data[2:end, column_names .== "CO2_F_MDS_QC"]
CO2 = driver_data[2:end, column_names .== "CO2_F_MDS"] .* 1e-6; # convert \mumol to mol
replace_missing_with_mean!(CO2, CO2_F)
SWC_F = driver_data[2:end, column_names .== "SWC_F_MDS_1_QC"] # Most likely 5cm depth
SWC = driver_data[2:end, column_names .== "SWC_F_MDS_1"] ./ 100; # to convert from % to m^3/m^3
replace_missing_with_mean!(SWC, SWC_F)
TS_F = driver_data[2:end, column_names .== "TS_F_MDS_1_QC"] # Most likely 5cm depth
TS = driver_data[2:end, column_names .== "TS_F_MDS_1"] .+ 273.15;# convert C to K
replace_missing_with_mean!(TS, TS_F)
GPP = driver_data[2:end, column_names .== "GPP_DT_VUT_REF"] .* 1e-6 # to convert from micromol to mol.
LE = driver_data[2:end, column_names .== "LE_CORR"]
H_CORR = driver_data[2:end, column_names .== "H_CORR"]
H = driver_data[2:end, column_names .== "H_F_MDS"]
H_F = driver_data[2:end, column_names .== "H_F_MDS_QC"]
replace_missing_with_mean!(H, H_F)
G = driver_data[2:end, column_names .== "G_F_MDS"]
G_F = driver_data[2:end, column_names .== "G_F_MDS_QC"]
replace_missing_with_mean!(G, G_F)

LW_OUT = driver_data[2:end, column_names .== "LW_OUT"]# This has missing data
SW_OUT = driver_data[2:end, column_names .== "SW_OUT"]# This has missing data
replace_missing_with_mean_by_value!(SW_OUT)
replace_missing_with_mean_by_value!(LW_OUT)

LOCAL_DATETIME = DateTime.(string.(driver_data[2:end, 1]), "yyyymmddHHMM")
UTC_DATETIME = LOCAL_DATETIME .+ Dates.Hour(6)
DATA_DT = Second(LOCAL_DATETIME[2] - LOCAL_DATETIME[1]).value # seconds

thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
esat =
    Thermodynamics.saturation_vapor_pressure.(
        Ref(thermo_params),
        TA,
        Ref(Thermodynamics.Liquid()),
    )
e = @. esat - VPD
q = @. 0.622 * e ./ (PA - 0.378 * e)

#Make a bunch of splines
seconds = FT.(0:DATA_DT:((length(UTC_DATETIME) - 1) * DATA_DT));
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


# Construct the drivers. For the reference time we will use the UTC time at the
# start of the simulation
atmos = ClimaLSM.PrescribedAtmosphere(
    precipitation_function,
    snow_precip,
    atmos_T,
    atmos_u,
    atmos_q,
    atmos_p,
    UTC_DATETIME[1],
    atmos_h;
    c_co2 = atmos_co2,
)

lat = FT(38.7441) # degree
long = FT(-92.2000) # degree

function zenith_angle(
    t,
    orbital_data,
    ref_time;
    latitude = lat,
    longitude = long,
    insol_params::Insolation.Parameters.InsolationParameters{FT} = earth_param_set.insol_params,
) where {FT}
    # This should be time in UTC
    dt = ref_time + Dates.Second(round(t))
    # Orbital Data uses Float64, so we need to convert to our sim
    # FT.
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


radiation = ClimaLSM.PrescribedRadiativeFluxes(
    FT,
    SW_IN_spline,
    LW_IN_spline,
    UTC_DATETIME[1];
    θs = zenith_angle,
    orbital_data = Insolation.OrbitalData(),
)


# LAI data
lai_af = ArtifactFile(
    url = "https://caltech.box.com/shared/static/yfqj0yqkx8yps7ydltsaixuy99083pon.csv",
    filename = "Ozark_MODIS_LAI_2005.csv",
)
lai_dataset = ArtifactWrapper(@__DIR__, "modis_data", ArtifactFile[lai_af]);
lai_dataset_path = get_data_folder(lai_dataset);
lai_raw_data = joinpath(lai_dataset_path, "Ozark_MODIS_LAI_2005.csv");
LAI_data = readdlm(lai_raw_data, ',') #m2.m-2
LAI_column_names = LAI_data[1, :]
LAI_data = LAI_data[2:end, LAI_column_names .== "lai_mean"] #m2.m-2
# This has the same timestamp as the driver data, so it's ok to use the time column from that file here
LAIspline = Spline1D(seconds, LAI_data[:])
LAIfunction = (t) -> eltype(t)(LAIspline(t))
