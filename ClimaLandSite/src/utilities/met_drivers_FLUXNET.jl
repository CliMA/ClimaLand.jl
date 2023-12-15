using ArtifactWrappers
using DelimitedFiles
using Dierckx
using Thermodynamics
using Dates
using Formatting
using HTTP
using JSON
using Insolation

using ClimaLSM

# Methods for reading in the LAI data from MODIS data
include(
    joinpath(pkgdir(ClimaLSM), "experiments/integrated/fluxnet/pull_MODIS.jl"),
)

af = ArtifactFile(
    url = data_link,
    filename = "AMF_$(site_ID)_FLUXNET_FULLSET.csv",
)
dataset = ArtifactWrapper(
    "$(@__DIR__)" * "/$site_ID",
    "ameriflux_data",
    ArtifactFile[af],
);
dataset_path = get_data_folder(dataset);
data = joinpath(dataset_path, "AMF_$(site_ID)_FLUXNET_FULLSET.csv");
driver_data = readdlm(data, ',')

LOCAL_DATETIME = DateTime.(format.(driver_data[2:end, 1]), "yyyymmddHHMM")
UTC_DATETIME = LOCAL_DATETIME .+ Dates.Hour(time_offset)
DATA_DT = Second(LOCAL_DATETIME[2] - LOCAL_DATETIME[1]).value # seconds

# Label of every data column to be collected from the driver data file
labels = (
    :TA,
    :VPD,
    :PA,
    :P,
    :WS,
    :LW_IN,
    :SW_IN,
    :G,
    :GPP,
    :LE,
    :H,
    :LW_OUT,
    :SW_OUT,
    :SWC,
    :TS,
    :CO2,
    :RECO,
)

# For every data column to be collected, name of the column in the file, 
# transformation function to desired unit, and string of final unit
collect_args = [
    ("TA_F", (x) -> x .+ 273.15, "K")
    ("VPD_F", (x) -> x .* 100, "Pa")
    ("PA_F", (x) -> x .* 1000, "Pa")
    ("P_F", (x) -> x ./ (1000 * DATA_DT), "m/s")
    ("WS_F", (x) -> x, "m/s")
    ("LW_IN_F", (x) -> x, "W/m^2")
    ("SW_IN_F", (x) -> x, "W/m^2")
    ("G_F_MDS", (x) -> x, "W/m^2")
    ("GPP_DT_VUT_REF", (x) -> x .* 1e-6, "mol/m^2/s")
    ("LE_CORR", (x) -> x, "W/m^2")
    ("H_CORR", (x) -> x, "W/m^2")
    ("LW_OUT", (x) -> x, "W/m^2")
    ("SW_OUT", (x) -> x, "W/m^2")
    ("SWC_F_MDS_1", (x) -> x ./ 100, "m^3/m^3")
    ("TS_F_MDS_1", (x) -> x .+ 273.15, "K")
    ("CO2_F_MDS", (x) -> x .* 1e-6, "mol")
    ("RECO_DT_VUT_REF", (x) -> x .* 1e-6, "mol/m^-2/s")
]

# Named tuple mapping every label to a DataColumn with the correct transformed 
# unit and data status. Automatically fills gaps in the data 
drivers = (;
    zip(
        labels,
        [
            transform_column(
                VerifiedColumn(driver_data, name, ""),
                transform,
                unit,
            ) for (name, transform, unit) in collect_args
        ],
    )...
)
# Check that all required driver data is present, otherwise throw error
required = [
    "TA" => drivers.TA,
    "VPD" => drivers.VPD,
    "PA" => drivers.PA,
    "P" => drivers.P,
    "WS" => drivers.WS,
    "LW_IN" => drivers.LW_IN,
    "SW_IN" => drivers.SW_IN,
]
missing_drivers = []
for (name, column) in required
    if column.status == absent
        push!(missing_drivers, name)
    end
end
if length(missing_drivers) != 0
    error("Driver data missing for columns: $([missing_drivers[i] * " " for 
        i in 1:length(missing_drivers)]...)")
end

thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
esat =
    Thermodynamics.saturation_vapor_pressure.(
        Ref(thermo_params),
        drivers.TA.values,
        Ref(Thermodynamics.Liquid()),
    )
e = @. esat - drivers.VPD.values
q = @. 0.622 * e ./ (drivers.PA.values - 0.378 * e)

# Create splines for each atmospheric driver needed for PrescribedAtmosphere
# and for PrescribedRadiation
seconds = FT.(0:DATA_DT:((length(UTC_DATETIME) - 1) * DATA_DT));
p_spline = Spline1D(seconds, -drivers.P.values[:]) # m/s
atmos_q = Spline1D(seconds, q[:])
atmos_T = Spline1D(seconds, drivers.TA.values[:])
atmos_p = Spline1D(seconds, drivers.PA.values[:])
atmos_co2 = Spline1D(seconds, drivers.CO2.values[:])
atmos_u = Spline1D(seconds, drivers.WS.values[:])
LW_IN_spline = Spline1D(seconds, drivers.LW_IN.values[:])
SW_IN_spline = Spline1D(seconds, drivers.SW_IN.values[:])
precipitation_function(t::FT) where {FT} = p_spline(t) < 0.0 ? p_spline(t) : 0.0 # m/s
snow_precip(t) = FT(0) # this is likely not correct

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

function zenith_angle(
    t,
    ref_time;
    latitude = lat,
    longitude = long,
    insol_params::Insolation.Parameters.InsolationParameters{FT} = earth_param_set.insol_params,
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

radiation = ClimaLSM.PrescribedRadiativeFluxes(
    FT,
    SW_IN_spline,
    LW_IN_spline,
    UTC_DATETIME[1];
    θs = zenith_angle,
)

# Start and end dates of data in MODIS format
modis_start_date = "A$(Dates.year(UTC_DATETIME[1]))$(lpad(Dates.dayofyear(UTC_DATETIME[1]), 3, "0"))"
modis_end_date = "A$(Dates.year(UTC_DATETIME[end]))$(lpad(Dates.dayofyear(UTC_DATETIME[end]), 3, "0"))"

MODIS_LAI = single_col_data_matrix(
    parse_response(
        check_response(
            send_get_subset(
                "MCD15A2H",
                modis_start_date,
                modis_end_date,
                site_ID,
                band = "Lai_500m",
            ),
        ),
    ),
)

LAI_dt = Second(MODIS_LAI[2, 1] - MODIS_LAI[1, 1]).value
LAI_seconds = FT.(0:LAI_dt:((length(MODIS_LAI[:, 1]) - 1) * LAI_dt))

# LAI spline for radiative transfer
LAIspline = Spline1D(LAI_seconds, MODIS_LAI[:, 2])
LAIfunction = (t) -> eltype(t)(LAIspline(t))

# Necessary inputs from LAI Data
# Note that f_root_to_shoot, capacity, and h_leaf are site-specific parameters
# defined in the parameters file for each site
maxLAI = FT(maximum(MODIS_LAI[:, 2]))
RAI = maxLAI * f_root_to_shoot
plant_ν = capacity / (maxLAI * h_leaf) / FT(1000)
