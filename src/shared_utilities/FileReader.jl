"""
    FileReader

This module coordinates reading, regridding, and interpolating
data from NetCDF files that is required for global land model
simulations, including globally varying parameters which may
or may not change in time. It also includes regridding and
temporal interpolations of this data.

This is based on ClimaCoupler.jl's BCReader and TimeManager modules.
"""
module FileReader

using ClimaComms
using ClimaCore: Fields, Spaces
using Dates
using JLD2
using CFTime

using ClimaLSM.Regridder

export AbstractPrescribedData,
    PrescribedDataTemporal,
    PrescribedDataStatic,
    FileInfo,
    FileState,
    SimInfo,
    read_data_fields!,
    next_date_in_file,
    interpolate_data,
    to_datetime


"""
    abstract type AbstractPrescribedData

An abstract type for storing prescribed data info. Subtypes
include temporally-varying prescribed data and static prescribed data.
"""
abstract type AbstractPrescribedData end

"""
    PrescribedDataStatic <: AbstractPrescribedData

Stores information to read in a prescribed variable from a file.
The data is read in once and stored without changing for the duration of a
simulation. This type is meant to be used with input data that does not have
a time dimension.
Each of the fields of this struct is itself a struct.

# Inputs:
- file_info::FI     # unchanging info about the input data file
"""
struct PrescribedDataStatic{FI} <: AbstractPrescribedData
    file_info::FI
end

"""
    PrescribedDataTemporal <: AbstractPrescribedData

Stores information to read in a prescribed variable from a file.
Contains sufficient information to read in the variable at various
timesteps, and to coordinate this reading between variables coming from
different files. This type is meant to be used with input data that has
a time dimension.
Each of the fields of this struct is itself a struct.

# Inputs:
- file_info::FI     # unchanging info about the input data file
- file_state::FS    # info about the dates currently being read from file
- sim_info::S       # unchanging info about the start date/time of the simulation
"""
struct PrescribedDataTemporal{FI, FS, S} <: AbstractPrescribedData
    file_info::FI
    file_state::FS
    sim_info::S
end

"""
    FileInfo

Stores information about the current data being read in from a file.

# Inputs:
- infile_path::String        # path to the input NetCDF data file
- regrid_dirpath::String     # directory for storing files used in regridding
- varname::String            # name of the variable we're reading from the input file, which we assume is a scalar
- outfile_root::String       # root for regridded data files generated when writing data at each time from input file
- all_dates::Vector          # vector containing all dates of the input file, which we assume are `DateTime`s or `DateTimeNoLeap`s
- date_idx0::Vector{Int}     # index of the first data in the file being used for this simulation
"""
struct FileInfo
    infile_path::String
    regrid_dirpath::String
    varname::String
    outfile_root::String
    all_dates::Vector
    date_idx0::Vector{Int}
end

"""
    FileState

Stores information about the current data being read in from a file.

# Inputs:
- data_fields::F                # tuple of two fields at consecutive dates, that will be used for interpolation
- date_idx::Vector{Int}         # index in the input file of the first data field currently being used
- segment_length::Vector{Int}   # length of the time interval between the two data field entries; used in temporal interpolation
"""
struct FileState{F}
    data_fields::F
    date_idx::Vector{Int}
    segment_length::Vector{Int}
end

"""
    SimInfo

Stores information about the simulation being run. We may want to store
multiple copies of an instance of this struct in multiple PrescribedDataTemporal
objects if we're reading in data over time for multiple variables.

# Inputs:
- date_ref::D    # a reference date before or at the start of the simulation
- t_start        # time in seconds since `date_ref`
"""
struct SimInfo{D}
    date_ref::D
    t_start::Any
end

"""
    PrescribedDataStatic(
        infile_path::String,
        regrid_dirpath::String,
        varname::String,
    )

Constructor for the `PrescribedDataStatic`` type.
Creates a `FileInfo` object containing all the information needed to read in
the data stored in the input file, which will later be regridded to our
simulation grid. Date-related args (last 3 to FileInfo) are unused for static
data maps.
"""
function PrescribedDataStatic(
    infile_path::String,
    regrid_dirpath::String,
    varname::String,
)
    file_info = FileInfo(infile_path, regrid_dirpath, varname, "", [], [])
    return PrescribedDataStatic(file_info)
end


"""
    PrescribedDataTemporal{FT}(
        regrid_dirpath,
        infile_path,
        varname,
        date_ref,
        t_start,
        surface_space;
        mono = true,
    ) where {FT <: AbstractFloat}

Constructor for the `PrescribedDataTemporal` type.
Regrids from the input lat-lon grid to the simulation cgll grid, saving
the regridded output in a new file found at `regrid_dirpath`, and
returns the info required to run the simulation using this prescribed
data packaged into a single `PrescribedDataTemporal` struct.

# Arguments
- `regrid_dirpath`   # directory the data file is stored in.
- `infile_path`      # NCDataset file containing data to regrid.
- `varname`          # name of the variable to be regridded.
- `date_ref`         # reference date to coordinate start of the simulation
- `t_start`          # start time of the simulation relative to `date_ref` (date_start = date_ref + t_start)
- `surface_space`    # the space to which we are mapping.
- `mono`             # flag for monotone remapping of `infile_path`.

# Returns
- `PrescribedDataTemporal` object
"""
function PrescribedDataTemporal{FT}(
    regrid_dirpath::String,
    infile_path::String,
    varname::String,
    date_ref::Union{DateTime, DateTimeNoLeap},
    t_start,
    surface_space::Spaces.AbstractSpace;
    mono::Bool = true,
) where {FT <: AbstractFloat}
    comms_ctx = surface_space.topology.context
    outfile_root = varname * "_cgll"

    # Regrid data at all times from lat/lon (RLL) to simulation grid (CGLL)
    if ClimaComms.iamroot(comms_ctx)
        Regridder.hdwrite_regridfile_rll_to_cgll(
            FT,
            regrid_dirpath,
            infile_path,
            varname,
            surface_space,
            outfile_root;
            mono = mono,
        )
    end
    ClimaComms.barrier(comms_ctx)
    all_dates = JLD2.load(
        joinpath(regrid_dirpath, outfile_root * "_times.jld2"),
        "times",
    )

    # Init time tracking info
    data_fields =
        Fields.zeros(FT, surface_space), Fields.zeros(FT, surface_space)
    # Store `segment_length` as an array so we can modify it as a field of a struct
    segment_length = Int[0]

    date_start = date_ref + Dates.Second(t_start)
    if date_start < all_dates[1]
        @warn "Simulation start date is before first file data"
    end

    # Find the index of the start file closest to this simulation's start date
    # Like `segment_length`, store in an array so we can modify in struct
    date_idx0 =
        [argmin(abs.(Dates.value(date_start) .- Dates.value.(all_dates[:])))]

    # Construct component structs of PrescribedDataTemporal object
    file_info = FileInfo(
        infile_path,
        regrid_dirpath,
        varname,
        outfile_root,
        all_dates,
        date_idx0,
    )
    file_state = FileState(data_fields, copy(date_idx0), segment_length)
    sim_info = SimInfo(date_ref, t_start)

    args = (file_info, file_state, sim_info)
    return PrescribedDataTemporal{typeof.(args)...}(args...)
end

"""
    read_data_fields!(prescribed_data::PrescribedDataTemporal, date::DateTime, space::Spaces.AbstractSpace)

Extracts data from regridded (to model grid) NetCDF files.
The times for which data is extracted depends on the specifications in the
`prescribed_data` struct).
Data at one point in time is stored in `prescribed_data.file_state.data_fields[1]`, and
data at the next time is stored in `prescribed_data.file_state.data_fields[2]`. With these
two data fields saved, we can interpolate between them for any dates
in this range of time.

# Arguments
- `prescribed_data`      # containing parameter value data.
- `date`                 # current date to read in data for.
- `space`                # space we're remapping the data onto.
"""
function read_data_fields!(
    prescribed_data::PrescribedDataTemporal,
    date::DateTime,
    space::Spaces.AbstractSpace,
)
    comms_ctx = space.topology.context
    pd_file_info = prescribed_data.file_info
    pd_file_state = prescribed_data.file_state

    (; regrid_dirpath, outfile_root, varname, all_dates) = pd_file_info

    date_idx0 = pd_file_info.date_idx0[1]
    date_idx = pd_file_state.date_idx[1]


    # Case 1: Current date is before or at first date in data file
    #  Load in data at first date for both `data_fields[1]` and `data_fields[2]`
    if (date_idx == date_idx0) && (Dates.days(date - all_dates[date_idx]) <= 0)
        if date !== all_dates[date_idx]
            @warn "this time period is before input data - using file from $(all_dates[date_idx0])"
        end

        pd_file_state.data_fields[1] .= Regridder.swap_space!(
            Regridder.read_from_hdf5(
                regrid_dirpath,
                outfile_root,
                all_dates[Int(date_idx0)],
                varname,
                comms_ctx,
            ),
            space,
        )
        pd_file_state.data_fields[2] .= pd_file_state.data_fields[1]
        pd_file_state.segment_length .= 0

        # Case 2: current date is at or after last date in input file
        #  Load in data at last date for both `data_fields[1]` and `data_fields[2]`
    elseif Dates.days(date - all_dates[end - 1]) >= 0
        @warn "this time period is after input data - using file from $(all_dates[end - 1])"
        pd_file_state.data_fields[1] .= Regridder.swap_space!(
            Regridder.read_from_hdf5(
                regrid_dirpath,
                outfile_root,
                all_dates[end],
                varname,
                comms_ctx,
            ),
            space,
        )
        pd_file_state.data_fields[2] .= pd_file_state.data_fields[1]
        pd_file_state.segment_length .= 0

        # Case 3: current date is later than date of data being read in
        #  Load in data at most recent past date in `data_fields[1]` and
        #  next date in `data_fields[2]`
    elseif Dates.days(date - all_dates[Int(date_idx)]) > 0
        # Increment `date_idx` to use next date
        date_idx = pd_file_state.date_idx[1] += Int(1)
        @warn "On $date updating data files: dates = [ $(all_dates[Int(date_idx)]) , $(all_dates[Int(date_idx+1)]) ]"
        # Time between consecutive dates being stored gives `segment_length`
        pd_file_state.segment_length .=
            (all_dates[Int(date_idx + 1)] - all_dates[Int(date_idx)]).value
        # Read in data fields at both dates
        pd_file_state.data_fields[1] .= Regridder.swap_space!(
            Regridder.read_from_hdf5(
                regrid_dirpath,
                outfile_root,
                all_dates[Int(date_idx)],
                varname,
                comms_ctx,
            ),
            space,
        )
        pd_file_state.data_fields[2] .= Regridder.swap_space!(
            Regridder.read_from_hdf5(
                regrid_dirpath,
                outfile_root,
                all_dates[Int(date_idx + 1)],
                varname,
                comms_ctx,
            ),
            space,
        )
        # Case 4: Nearest date index not being used for initial date field
        #  Throw warning and update date index
    elseif Dates.days(date - all_dates[Int(date_idx + 1)]) > 2
        nearest_idx =
            argmin(abs.(Dates.value(date) .- Dates.value.(all_dates[:])))
        pd_file_state.date_idx[1] = date_idx = date_idx0 = nearest_idx
        @warn "init data does not correspond to start date. Initializing with `file_state.date_idx[1] = date_idx = date_idx0 = $nearest_idx` for this start date"
    else
        throw(ErrorException("Check input file specification"))
    end
end

"""
    next_date_in_file(prescribed_data::PrescribedDataTemporal)

Returns the next date stored in the file `prescribed_data` struct after the
current date index given by `date_idx`.
Note: this function does not update `date_idx`, so repeated calls will
return the same value unless `date_idx` is modified elsewhere in between.

# Arguments
- `prescribed_data`   # contains all input file information needed for the simulation.

# Returns
- DateTime or DateTimeNoLeap
"""
next_date_in_file(prescribed_data::PrescribedDataTemporal) =
    prescribed_data.file_info.all_dates[prescribed_data.file_state.date_idx[1] + Int(
        1,
    )]

"""
    interpolate_data(prescribed_data::PrescribedDataTemporal, date::Union{DateTime, DateTimeNoLeap}, space::Spaces.AbstractSpace)

Interpolates linearly between two `Fields` in the `prescribed_data` struct.

# Arguments
- `prescribed_data`      # contains fields to be interpolated.
- `date`                 # start date for data.
- `space`                # the space of our simulation.

# Returns
- Fields.field
"""
function interpolate_data(
    prescribed_data::PrescribedDataTemporal,
    date::Union{DateTime, DateTimeNoLeap},
    space::Spaces.AbstractSpace,
)
    FT = Spaces.undertype(space)
    (; segment_length, date_idx, data_fields) = prescribed_data.file_state
    (; all_dates) = prescribed_data.file_info

    # Interpolate if the time period between dates is nonzero
    if segment_length[1] > FT(0) && date != all_dates[Int(date_idx[1])]
        Δt_tt1 = FT((date - all_dates[Int(date_idx[1])]).value)
        return interpol.(
            data_fields[1],
            data_fields[2],
            Δt_tt1,
            FT(segment_length[1]),
        )
        # Otherwise use the data at the first date
    else
        return data_fields[1]
    end
end

"""
    interpol(f1::FT, f2::FT, Δt_tt1::FT, Δt_t2t1::FT)

Performs linear interpolation of `f` at time `t` within
a segment `Δt_t2t1 = (t2 - t1)`, of fields `f1` and `f2`, with `t2 > t1`.

# Arguments
- `f1`::FT          # first value to be interpolated (`f(t1) = f1`).
- `f2`::FT          # second value to be interpolated.
- `Δt_tt1`::FT      # time between `t1` and some `t` (`Δt_tt1 = (t - t1)`).
- `Δt_t2t1`::FT     # time between `t1` and `t2`.

# Returns
- FT
"""
function interpol(f1::FT, f2::FT, Δt_tt1::FT, Δt_t2t1::FT) where {FT}
    interp_fraction = Δt_tt1 / Δt_t2t1
    @assert abs(interp_fraction) <= FT(1) "time interpolation weights must be <= 1, but `interp_fraction` = $interp_fraction"
    return f1 * (FT(1) - interp_fraction) + f2 * (interp_fraction)
end

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

end
