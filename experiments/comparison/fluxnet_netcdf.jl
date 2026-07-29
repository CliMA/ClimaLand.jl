# Write a site's FLUXNET CSV forcing into a NetCDF laid out like the real site-level met
# files, so the CSV forcing path and the NetCDF forcing path can be driven from identical
# physical values.
#
# No site on this machine has both a CSV and a real NetCDF met file (the
# `callmip_phase1_forcing` artifact has no download URL and is not in the depot), so the
# NetCDF is generated here. Every conversion mirrors what `prescribed_forcing_fluxnet`
# applies to the CSV, which is what makes the two paths comparable: timestamps are LOCAL
# standard time at the CSV's START/END midpoint (`prescribed_forcing_netcdf` applies the
# local -> UTC shift itself), Tair in K, Psurf in Pa, Precip as a rate rather than an
# accumulation, and Qair precomputed with the same `compute_q`. What is left over is the
# loading machinery, which is the thing under comparison.

import ClimaLand
using Dates
using NCDatasets

# The CSV reader, column mapping, `compute_q`, and the fill-value predicate are internal to
# the extension, and the generator has to apply exactly what the CSV forcing path applies.
const FSExt = Base.get_extension(ClimaLand, :FluxnetSimulationsExt)

# CSV column backing each generated variable, used only to count fill values.
const _CSV_SOURCE = Dict(
    :Tair => "TA_F",
    :Wind => "WS_F",
    :Qair => "TA_F",
    :Psurf => "PA_F",
    :SWdown => "SW_IN_F",
    :LWdown => "LW_IN_F",
    :CO2air => "CO2_F_MDS",
    :Precip => "P_F",
    :VPD => "VPD_F",
)

const _TIMESTAMP_NAMES = ("TIMESTAMP_START", "TIMESTAMP_END")

# Local standard time at the midpoint of each averaging period, which is the timestamp the CSV
# forcing path uses and the one written to the NetCDF.
function fluxnet_local_dates(data, name_map)
    stamp(name) =
        DateTime.(string.(Int.(data[:, name_map[name]])), "yyyymmddHHMM")
    local_start = stamp("TIMESTAMP_START")
    local_end = stamp("TIMESTAMP_END")
    return local_start .+ (local_end .- local_start) ./ 2
end

"""
    fluxnet_data_dt(site_ID)

Return the spacing, in seconds, between `site_ID`'s forcing observations.

Sizing a warmup needs this: the NetCDF path only touches disk when the pair of dates
bracketing `t` changes, which happens once per `fluxnet_data_dt` of simulated time.
"""
function fluxnet_data_dt(site_ID)
    (data, columns) = FSExt.read_fluxnet_data(site_ID)
    name_map = FSExt.get_column_name_map(_TIMESTAMP_NAMES, columns)
    local_dates = fluxnet_local_dates(data, name_map)
    return Dates.value(Second(local_dates[2] - local_dates[1]))
end

"""
    write_fluxnet_netcdf(path, site_ID, lat, long, time_offset, thermo_params, FT)

Write `site_ID`'s FLUXNET CSV forcing to `path` as a FLUXNET-convention NetCDF and return
`(; path, first_precip_date)`.

`first_precip_date` is the UTC datetime of the first nonzero precipitation in the record, or
`nothing` if the site records none. A comparison needs it: until precipitation occurs, the
composed `Precip`/`VPD` forcing is identically zero in both paths and its agreement carries
no information.
"""
function write_fluxnet_netcdf(
    path,
    site_ID,
    lat,
    long,
    time_offset,
    thermo_params,
    ::Type{FT},
) where {FT}
    (data, columns) = FSExt.read_fluxnet_data(site_ID)
    varnames = (
        "TIMESTAMP_START",
        "TIMESTAMP_END",
        "TA_F",
        "VPD_F",
        "PA_F",
        "P_F",
        "WS_F",
        "LW_IN_F",
        "SW_IN_F",
        "CO2_F_MDS",
    )
    name_map = FSExt.get_column_name_map(varnames, columns)
    col(name) = FT.(data[:, name_map[name]])

    local_dates = fluxnet_local_dates(data, name_map)
    data_dt = Dates.value(Second(local_dates[2] - local_dates[1]))

    TA_F, VPD_F, PA_F, P_F = col("TA_F"), col("VPD_F"), col("PA_F"), col("P_F")
    # `Precip` is a rate in kg/m^2/s while the CSV holds an accumulation in mm over the
    # interval; 1 mm == 1 kg/m^2, so dividing by the interval converts exactly, and the two
    # paths' precipitation then agree without either sign or /1000 being duplicated.
    vars = (
        Tair = TA_F .+ 273.15,
        Wind = col("WS_F"),
        Qair = FSExt.compute_q.(PA_F, VPD_F, TA_F; thermo_params),
        Psurf = PA_F .* 1000,
        SWdown = col("SW_IN_F"),
        LWdown = col("LW_IN_F"),
        CO2air = col("CO2_F_MDS"),
        Precip = P_F ./ data_dt,
        VPD = VPD_F,
    )

    # The CSV path drops entries flagged -9999; the handler path interpolates between the
    # dates in the file and cannot drop them, so report what is present rather than hide it.
    for name in keys(vars)
        n_missing =
            count(FSExt.var_missing, data[:, name_map[_CSV_SOURCE[name]]])
        n_missing > 0 && @warn "Fill values in CSV source" name n_missing
    end

    # The real files' layout: (x = 1, y = 1, time) with scalar latitude and longitude, which
    # is what the handler matches this file to a column by.
    NCDataset(path, "c") do ds
        defDim(ds, "x", 1)
        defDim(ds, "y", 1)
        defDim(ds, "time", length(local_dates))
        defVar(ds, "x", [1.0], ("x",))
        defVar(ds, "y", [1.0], ("y",))
        defVar(
            ds,
            "time",
            local_dates,
            ("time",);
            attrib = ["units" => "seconds since 1970-01-01 00:00:00"],
        )
        for (name, values) in pairs(vars)
            defVar(
                ds,
                String(name),
                reshape(values, :, 1, 1),
                ("time", "y", "x"),
            )
        end
        defVar(ds, "latitude", fill(lat, 1, 1), ("y", "x"))
        defVar(ds, "longitude", fill(long, 1, 1), ("y", "x"))
    end

    first_precip = findfirst(>(0), P_F)
    first_precip_date =
        isnothing(first_precip) ? nothing :
        local_dates[first_precip] - FSExt.hour_offset_to_period(time_offset)
    return (; path, first_precip_date)
end
