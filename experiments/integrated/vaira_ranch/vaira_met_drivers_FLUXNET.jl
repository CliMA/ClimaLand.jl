## Dataset obtained from Ameriflux. Funding for the AmeriFlux data portal was provided by the U.S. Department of Energy Office of Science.
## Vaira Ranch site
## Citation: Siyan Ma, Liukang Xu, Joseph Verfaillie, Dennis Baldocchi (2023),
## AmeriFlux FLUXNET-1F US-Var Vaira Ranch- Ione, Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1993904

## LAI data courtesy of Mitra Asadollahi

using DelimitedFiles
using Dierckx
using Thermodynamics
using Dates
using ArtifactWrappers
import Insolation

function replace_missing_with_mean!(field, flag)
    good_indices = (flag .== 0) .|| (flag .== 1)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end
1

function replace_missing_with_mean_by_value!(field)
    good_indices = .~(field .== -9999)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end

function replace_missing_with_zero_by_value!(field)
    good_indices = .~(field .== -9999)
    field[.~good_indices] .= 0.0
    return field
end

af = ArtifactFile(
    url = "https://caltech.box.com/shared/static/dx0p5idbsbrdebsda10t9pfv2lbdaz95.csv",
    filename = "AMF_US-Var_FLUXNET_FULLSET_HH_2003-2006.csv",
)
dataset = ArtifactWrapper(@__DIR__, "ameriflux_data-US-Var", ArtifactFile[af]);
dataset_path = get_data_folder(dataset);
data = joinpath(dataset_path, af.filename)
driver_data = readdlm(data, ',')
column_names = driver_data[1, :]

indices = 1:1:length(driver_data[2:end, 1])
LOCAL_DATETIME = DateTime.(string.(driver_data[2:end, 1]), "yyyymmddHHMM")
subset = indices[(LOCAL_DATETIME .>= Dates.DateTime(
    "200301010000",
    "yyyymmddHHMM",
)) .& (Dates.Year.(LOCAL_DATETIME) .<= Dates.Year(2003))]
LOCAL_DATETIME = LOCAL_DATETIME[subset]
UTC_DATETIME = LOCAL_DATETIME .+ Dates.Hour(8)
DATA_DT = Second(LOCAL_DATETIME[2] - LOCAL_DATETIME[1]).value # seconds

CO2_F = driver_data[2:end, column_names .== "CO2_F_MDS_QC"][subset]
CO2 = driver_data[2:end, column_names .== "CO2_F_MDS"][subset] .* 1e-6; # convert \mumol to mol
replace_missing_with_mean!(CO2, CO2_F)

RECO = driver_data[2:end, column_names .== "RECO_DT_VUT_REF"] .* 1e-6 # to convert from micromol to mol.
TA = driver_data[2:end, column_names .== "TA_F"][subset] .+ 273.15; # convert C to K
VPD = driver_data[2:end, column_names .== "VPD_F"][subset] .* 100; # convert hPa to Pa
PA = driver_data[2:end, column_names .== "PA_F"][subset] .* 1000; # convert kPa to Pa
P = driver_data[2:end, column_names .== "P_F"][subset] ./ (1000 * 1800); # convert mm/HH to m/s
WS = driver_data[2:end, column_names .== "WS_F"][subset]; # already m/s
LW_IN = driver_data[2:end, column_names .== "LW_IN_F"][subset]
SW_IN = driver_data[2:end, column_names .== "SW_IN_F"][subset]
thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
esat =
    Thermodynamics.saturation_vapor_pressure.(
        Ref(thermo_params),
        TA,
        Ref(Thermodynamics.Liquid()),
    )
e = @. esat - VPD
q = @. 0.622 * e ./ (PA - 0.378 * e)
SWC_F = driver_data[2:end, column_names .== "SWC_F_MDS_1_QC"][subset] # Most likely 5cm depth
SWC_1 = driver_data[2:end, column_names .== "SWC_F_MDS_1"][subset] ./ 100; # to convert from % to m^3/m^3
replace_missing_with_mean!(SWC_1, SWC_F)

SWC_F = driver_data[2:end, column_names .== "SWC_F_MDS_2_QC"][subset] # Most likely 20cm depth
SWC_2 = driver_data[2:end, column_names .== "SWC_F_MDS_2"][subset] ./ 100; # to convert from % to m^3/m^3
replace_missing_with_mean!(SWC_2, SWC_F)

SWC_F = driver_data[2:end, column_names .== "SWC_F_MDS_3_QC"][subset] # Most likely 50cm depth
SWC_3 = driver_data[2:end, column_names .== "SWC_F_MDS_3"][subset] ./ 100; # to convert from % to m^3/m^3
replace_missing_with_mean!(SWC_3, SWC_F)

TS_F = driver_data[2:end, column_names .== "TS_F_MDS_1_QC"][subset] # -20
TS_1 = driver_data[2:end, column_names .== "TS_F_MDS_1"][subset] .+ 273.15;# convert C to K
replace_missing_with_mean!(TS_1, TS_F)

TS_F = driver_data[2:end, column_names .== "TS_F_MDS_2_QC"][subset] # -40
TS_2 = driver_data[2:end, column_names .== "TS_F_MDS_2"][subset] .+ 273.15;# convert C to K
replace_missing_with_mean!(TS_2, TS_F)

TS_F = driver_data[2:end, column_names .== "TS_F_MDS_3_QC"][subset] # -80
TS_3 = driver_data[2:end, column_names .== "TS_F_MDS_3"][subset] .+ 273.15;# convert C to K
replace_missing_with_mean!(TS_3, TS_F)

TS_F = driver_data[2:end, column_names .== "TS_F_MDS_4_QC"][subset] # -160
TS_4 = driver_data[2:end, column_names .== "TS_F_MDS_4"][subset] .+ 273.15;# convert C to K
replace_missing_with_mean!(TS_4, TS_F)

TS_F = driver_data[2:end, column_names .== "TS_F_MDS_5_QC"][subset] # -320
TS_5 = driver_data[2:end, column_names .== "TS_F_MDS_5"][subset] .+ 273.15;# convert C to K
replace_missing_with_mean!(TS_5, TS_F)

GPP = driver_data[2:end, column_names .== "GPP_DT_VUT_REF"][subset] .* 1e-6 # to convert from micromol to mol.
LE = driver_data[2:end, column_names .== "LE_CORR"][subset]
H = driver_data[2:end, column_names .== "H_CORR"][subset]
G = driver_data[2:end, column_names .== "G_F_MDS"][subset]
G_F = driver_data[2:end, column_names .== "G_F_MDS_QC"][subset]
replace_missing_with_mean!(G, G_F)

LW_OUT = driver_data[2:end, column_names .== "LW_OUT"][subset]# This is all missing
SW_OUT = driver_data[2:end, column_names .== "SW_OUT"][subset]# This has missing data
replace_missing_with_zero_by_value!(SW_OUT)
replace_missing_with_mean_by_value!(LW_OUT)


#Make a bunch of splines
seconds = 0.0:DATA_DT:((length(UTC_DATETIME) - 1) * DATA_DT);
p_spline = Spline1D(seconds, -P[:]) # m/s
atmos_q = Spline1D(seconds, q[:])
atmos_T = Spline1D(seconds, TA[:])
atmos_p = Spline1D(seconds, PA[:])
atmos_co2 = Spline1D(seconds, CO2[:])
atmos_u = Spline1D(seconds, WS[:])
LW_IN_spline = Spline1D(seconds, LW_IN[:])
SW_IN_spline = Spline1D(seconds, SW_IN[:])
atmos_h = FT(2)
precipitation_function(t) = p_spline(t) < 0.0 ? p_spline(t) : 0.0 # m/s
snow_precip(t) = 0.0 # Not correct in winter



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

lat = FT(38.4133) # degree
long = FT(-120.9508) # degree

function zenith_angle(
    t,
    ref_time;
    latitude = lat,
    longitude = long,
    insol_params::Insolation.Parameters.InsolationParameters{FT} = earth_param_set.insol_params,
) where {FT}
    # This should be time in UTC
    current_datetime = ref_time + Dates.Second(round(t))

    # Orbital Data uses Float64, so we need to convert to our sim FT.
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


radiation = ClimaLSM.PrescribedRadiativeFluxes(
    FT,
    SW_IN_spline,
    LW_IN_spline,
    UTC_DATETIME[1];
    θs = zenith_angle,
)

# LAI data
LAI_af = ArtifactFile(
    url = "https://caltech.box.com/shared/static/y5vf8s9qkoogglc1bc2eyu1k95sbjsc3.csv",
    filename = "US-Var_fSun_Vcmax_DD.csv",
)
LAI_dataset = ArtifactWrapper(@__DIR__, "lai_data-US-Var", ArtifactFile[LAI_af]);
LAI_datapath = joinpath(get_data_folder(LAI_dataset), LAI_af.filename)
LAI_data = readdlm(LAI_datapath, ',') #m2.m-2
LAI_column_names = LAI_data[1, :]
indices = 1:1:length(LAI_data[2:end, 1])
LAI_DATETIME = DateTime.(string.(LAI_data[2:end, 1]), "yyyymmddHHMM")
subset = indices[(LAI_DATETIME .>= Dates.DateTime(
    "200301010000",
    "yyyymmddHHMM",
)) .& (Dates.Year.(LAI_DATETIME) .<= Dates.Year(2003))]
LAI_DATETIME = LAI_DATETIME[subset]
LAI_timeseries = LAI_data[2:end, LAI_column_names .== "LAI"][subset] #m2.m-2

LAI_DATA_DT = Second(LAI_DATETIME[2] - LAI_DATETIME[1]).value # seconds

LAI_seconds = 0.0:LAI_DATA_DT:((length(LAI_DATETIME) - 1) * LAI_DATA_DT);
LAIspline = Spline1D(LAI_seconds, LAI_timeseries)
LAIfunction = (t) -> LAIspline(t)
