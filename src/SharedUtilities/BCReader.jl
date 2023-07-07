"""
    BCReader

This module coordinates reading of boundary conditions from NetCDF files,
as well as regridding calls and temporal interpolations from
monthly to daily intervals.

This is based on ClimaCoupler.jl's BCReader and TimeManager modules.
"""
module BCReader

using ClimaComms
using ClimaCore: Fields
using Dates
using JLD2
import CFTime

using ClimaLSM.Regridder

export to_datetime,
    datetime_to_strdate,
    BCFileInfo,
    bcfile_info_init,
    float_type_bcf,
    update_midmonth_data!,
    next_date_in_file,
    interpolate_midmonth_data,
    no_scaling


"""
    to_datetime(date)

Convert a DateTime-like object (e.g. DateTimeNoLeap) to a DateTime,
using CFTime.jl. We need this since the CESM2 albedo file contains
DateTimeNoLeap objects for dates, which can't be used for math with DateTimes.

Note that this conversion may fail if the date to convert doesn't
exist in the DateTime calendar.

# Arguments
- `date`: object to be converted to DateTime
"""
to_datetime(date) = CFTime.reinterpret(DateTime, date)

"""
    datetime_to_strdate(datetime::DateTime)

Convert from Date to String ("YYYYMMDD") format.

# Arguments
- `datetime`: [DateTime] object to be converted to string
"""
datetime_to_strdate(datetime::DateTime) =
    string(lpad(Dates.year(datetime), 4, "0")) *
    string(string(lpad(Dates.month(datetime), 2, "0"))) *
    string(lpad(Dates.day(datetime), 2, "0"))

"""
    BCFileInfo

Stores information specific to each boundary condition from a file and each variable.

# Inputs:
- path::b                         # directory of the BC file
- comms::X                    # communication context used for MPI
- regrid_dirpath::S               # filename root for regridded data
- varname::V                      # name of the variable
- all_dates::D                    # vector of all dates contained in the original data file
- monthly_fields::C               # tuple of the two monthly fields, that will be used for the daily interpolation
- scaling_function::O             # function that scales, offsets or transforms the raw variable
- land_fraction::M                # fraction with 1 = 100% land, 0 = 100% ocean and/or sea-ice
- segment_idx::Vector{Int}        # index of the monthly data in the file
- segment_idx0::Vector{Int}       # `segment_idx` of the file data that is closest to date0
- segment_length::Vector{Int}     # length of each month segment (used in the daily interpolation)
- interpolate_daily::Bool         # switch to trigger daily interpolation
"""
struct BCFileInfo{FT <: Real, B, X, S, V, D, C, O, M, VI}
    path::B
    comms::X
    regrid_dirpath::S
    varname::V
    all_dates::D
    monthly_fields::C
    scaling_function::O
    land_fraction::M
    segment_idx::VI
    segment_idx0::VI
    segment_length::VI
    interpolate_daily::Bool
end

BCFileInfo{FT}(args...) where {FT} =
    BCFileInfo{FT, typeof.(args[1:9])...}(args...)
float_type_bcf(::BCFileInfo{FT}) where {FT} = FT


"""
    bcfile_info_init(
        FT,
        regrid_dirpath,
        path,
        varname,
        boundary_space,
        comms;
        interpolate_daily = false,
        segment_idx0 = nothing,
        scaling_function = no_scaling,
        land_fraction = nothing,
        date0 = nothing,
        mono = true,
    )

Regrids from lat-lon grid to cgll grid, saving the output in a new file,
and returns the info packaged in a single struct.

Note: For now we assume our simulations always start at date0.
This will need to be changed if we read in mulitple files over time that
start at different dates, or if we want to run our simulation from a
different date.
We could address this by adding a `date_start` field,
where the date at any time `t` in a simulation is given by `date_start + t`.
This is relevant since the land model keeps track of time
as the number of seconds elapsed since the starting time, rather
than as a standalone DateTime.

# Arguments
- `FT`: [DataType] Float type.
- `regrid_dirpath`: [String] directory the BC file is stored in.
- `path`: [String] file containing data to regrid.
- `varname`: [String] name of the variable to be regridded.
- `boundary_space`: [Spaces.AbstractSpace] the space to which we are mapping.
- `comms`: [ClimaComms.AbstractCommsContext] context used for this operation.
- `interpolate_daily`: [Bool] switch to trigger daily interpolation.
- `segment_idx0`: [Vector{Int}] index of the file data that is closest to date0.
- `scaling function`: [Function] scales, offsets or transforms `varname`.
- `land_fraction`: [Fields.field] fraction with 1 = land, 0 = ocean / sea-ice.
- `date0`: [Dates.DateTime] start date of the file data.
- `mono`: [Bool] flag for monotone remapping of `path`.

# Returns
- `BCFileInfo`
"""
function bcfile_info_init(
    FT,
    regrid_dirpath,
    path,
    varname,
    boundary_space,
    comms;
    interpolate_daily = false,
    segment_idx0 = nothing,
    scaling_function = no_scaling,
    land_fraction = nothing,
    date0 = nothing,
    mono = true,
)

    # regrid all times and save to hdf5 files
    regrid_dirpath = varname * "_cgll"
    if ClimaComms.iamroot(comms)
        Regridder.hdwrite_regridfile_rll_to_cgll(
            FT,
            regrid_dirpath,
            path,
            varname,
            boundary_space;
            hd_outfile_root = regrid_dirpath,
            mono = mono,
        )
    end
    ClimaComms.barrier(comms)
    data_dates = JLD2.load(
        joinpath(regrid_dirpath, regrid_dirpath * "_times.jld2"),
        "times",
    )

    # init time tracking info
    current_fields =
        Fields.zeros(FT, boundary_space), Fields.zeros(FT, boundary_space)
    segment_length = [Int(0)]

    # unless the start file date is specified, find the closest one to the start date
    date0 == nothing ? date0 = data_dates[1] : nothing
    segment_idx0 =
        segment_idx0 != nothing ? segment_idx0 :
        [
            argmin(
                abs.(
                    parse(FT, datetime_to_strdate(date0)) .-
                    parse.(FT, datetime_to_strdate.(data_dates[:]))
                ),
            ),
        ]

    return BCFileInfo{FT}(
        FT,
        comms,
        regrid_dirpath,
        varname,
        data_dates,
        current_fields,
        scaling_function,
        land_fraction,
        deepcopy(segment_idx0),
        segment_idx0,
        segment_length,
        interpolate_daily,
    )
end

"""
    update_midmonth_data!(date, bcf_info::BCFileInfo{FT}) where {FT}

Extracts boundary condition data from regridded (to model grid) NetCDF files.
The times for which data is extracted depends on the specifications in the
`bcf_info` struct).

# Arguments
- `date`: [Dates.DateTime] start date for data.
- `bcf_info`: [BCFileInfo] containing boundary condition data.
"""
function update_midmonth_data!(date, bcf_info::BCFileInfo{FT}) where {FT}
    # monthly count
    (;
        regrid_dirpath,
        comms,
        regrid_dirpath,
        varname,
        all_dates,
        scaling_function,
    ) = bcf_info
    midmonth_idx = bcf_info.segment_idx[1]
    midmonth_idx0 = bcf_info.segment_idx0[1]

    if (midmonth_idx == midmonth_idx0) &&
       (Dates.days(date - all_dates[midmonth_idx]) <= 0) # for init
        if date !== all_dates[midmonth_idx]
            @warn "this time period is before BC data - using file from $(all_dates[midmonth_idx0])"
        end
        midmonth_idx = bcf_info.segment_idx[1] -= Int(1)
        midmonth_idx =
            midmonth_idx < Int(1) ? midmonth_idx + Int(1) : midmonth_idx
        bcf_info.monthly_fields[1] .= scaling_function(
            Regridder.read_from_hdf5(
                regrid_dirpath,
                regrid_dirpath,
                all_dates[Int(midmonth_idx0)],
                varname,
                comms,
            ),
            bcf_info,
        )
        bcf_info.monthly_fields[2] .= deepcopy(bcf_info.monthly_fields[1])
        bcf_info.segment_length .= Int(0)

    elseif Dates.days(date - all_dates[end - 1]) > 0 # for fini
        @warn "this time period is after BC data - using file from $(all_dates[end - 1])"
        bcf_info.monthly_fields[1] .= scaling_function(
            Regridder.read_from_hdf5(
                regrid_dirpath,
                regrid_dirpath,
                all_dates[end],
                varname,
                comms,
            ),
            bcf_info,
        )
        bcf_info.monthly_fields[2] .= deepcopy(bcf_info.monthly_fields[1])
        bcf_info.segment_length .= Int(0)

        # throw error when there are closer initial indices for the bc file data that matches this date0
    elseif Dates.days(date - all_dates[Int(midmonth_idx + 1)]) > 2
        nearest_idx = argmin(
            abs.(
                parse(FT, datetime_to_strdate(date)) .-
                parse.(FT, datetime_to_strdate.(all_dates[:]))
            ),
        )
        bcf_info.segment_idx[1] = midmonth_idx = midmonth_idx0 = nearest_idx
        @warn "init data does not correspond to start date. Initializing with `SIC_info.segment_idx[1] = midmonth_idx = midmonth_idx0 = $nearest_idx` for this start date"

        # date crosses to the next month
    elseif Dates.days(date - all_dates[Int(midmonth_idx)]) > 0
        midmonth_idx = bcf_info.segment_idx[1] += Int(1)
        @warn "On $date updating monthly data files: mid-month dates = [ $(all_dates[Int(midmonth_idx)]) , $(all_dates[Int(midmonth_idx+1)]) ]"
        bcf_info.segment_length .=
            (
                all_dates[Int(midmonth_idx + 1)] - all_dates[Int(midmonth_idx)]
            ).value
        bcf_info.monthly_fields[1] .= scaling_function(
            Regridder.read_from_hdf5(
                regrid_dirpath,
                regrid_dirpath,
                all_dates[Int(midmonth_idx)],
                varname,
                comms,
            ),
            bcf_info,
        )
        bcf_info.monthly_fields[2] .= scaling_function(
            Regridder.read_from_hdf5(
                regrid_dirpath,
                regrid_dirpath,
                all_dates[Int(midmonth_idx + 1)],
                varname,
                comms,
            ),
            bcf_info,
        )

    else
        throw(ErrorException("Check boundary file specification"))
    end
end

"""
    next_date_in_file(bcf_info)

Returns the next date stored in the file `bcfile_info` struct after the
current date index given by `segment_idx`.
Note: this function does not update `segment_idx`, so repeated calls will
return the same value unless `segment_idx` is modified elsewhere in between.

# Arguments
- `bcf_info`: [BCFileInfo] containing the date information.

# Returns
- Dates.DateTime
"""
next_date_in_file(bcf_info::BCFileInfo{FT}) where {FT} =
    bcf_info.all_dates[bcf_info.segment_idx[1] + Int(1)]

# IO - daily
"""
    interpolate_midmonth_data(date, bcf_info::BCFileInfo{FT}) where {FT}

Interpolates linearly between two `Fields` in the `bcf_info` struct,
or returns the first Field if interpolation is switched off.

# Arguments
- `date`: [Dates.DateTime] start date for data.
- `bcf_info`: [BCFileInfo] contains fields to be interpolated.

# Returns
- Fields.field
"""
# TODO rename function
function interpolate_midmonth_data(date, bcf_info::BCFileInfo{FT}) where {FT}
    if bcf_info.interpolate_daily && bcf_info.segment_length[1] > FT(0)
        (; segment_length, segment_idx, all_dates, monthly_fields) = bcf_info

        return interpol.(
            monthly_fields[1],
            monthly_fields[2],
            FT((date - all_dates[Int(segment_idx[1])]).value),
            FT(segment_length[1]),
        )
    else
        return bcf_info.monthly_fields[1]
    end
end

"""
    interpol(f1::FT, f2::FT, Δt_tt1::FT, Δt_t2t1::FT)

Performs linear interpolation of `f` at time `t` within
a segment `Δt_t2t1 = (t2 - t1)`, of fields `f1` and `f2`, with `t2 > t1`.

# Arguments
- `f1`: [FT] first value to be interpolated (`f(t1) = f1`).
- `f2`: [FT] second value to be interpolated.
- `Δt_tt1`: [FT] time between `t1` and some `t` (`Δt_tt1 = (t - t1)`).
- `Δt_t2t1`: [FT] time between `t1` and `t2`.

# Returns
- FT
"""
function interpol(f1::FT, f2::FT, Δt_tt1::FT, Δt_t2t1::FT) where {FT}
    @assert Δt_t2t1 > FT(0) "t2 must be > t1, but `Δt_t2t1` = $Δt_t2t1"
    interp_fraction = Δt_tt1 / Δt_t2t1
    @assert abs(interp_fraction) <= FT(1) "time interpolation weights must be <= 1, but `interp_fraction` = $interp_fraction"
    return f1 * interp_fraction + f2 * (FT(1) - interp_fraction)
end

"""
    no_scaling(field, bcf_info)

Remap the values of a `field` onto the space of the
`bcf_info`'s land_fraction without scaling.

# Arguments
- `field`: [Fields.Field] contains the values to be remapped.
- `bcf_info`: [BCFileInfo] contains a land_fraction to remap onto the space of.
"""
no_scaling(field::Fields.Field, bcf_info::BCFileInfo{FT}) where {FT} =
    Regridder.swap_space!(field, axes(bcf_info.land_fraction))

end
