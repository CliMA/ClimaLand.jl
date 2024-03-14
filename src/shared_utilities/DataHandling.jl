"""
    DataHandling

The `DataHandling` module is responsible for read data from files and resample the data onto
the simulation grid.

This is no trivial task. Among the challenges:
- data can be large and cannot be read all in one go and/or held in memory
- regridding onto the simulation grid can be very expensive
- IO can be very expensive
- CPU/GPU communication can be a bottleneck

The `DataHandling` takes the divide and conquer approach: the various core tasks and
features and split into other independent modules (chiefly `FileReaders`, and `Regridders`).
Such modules can be developed, tested, and extended independently (as long as they maintain
a consistent interface). For instance, if need arises, the `DataHandler` can be used
(almost) directly to process files with a different format from NetCDF.

The key struct in `DataHandling` is the `DataHandler`. The `DataHandler` contains a
`FileReaders`, a `Regridders`, and other metadata necessary to perform its operations (e.g.,
target `ClimaCore.Space`). The `DataHandler` can be used for static or temporal data, and
exposes the following key functions:
- `regridded_snapshot(time)`: to obtain the regridded field at the given `time`. `time` has to be
                    available in the data.
- `available_times` (`available_dates`): to list all the `times` (`dates`) over which the
                    data is defined.
- `previous_time(time/date)` (`next_time(time/date)`): to obtain the time of the snapshot
                         before the given `time` or `date`. This can be used to compute the
                         interpolation weight for linear interpolation, or in combination
                         with `regridded_snapshot` to read a particular snapshot
Most `DataHandling` functions take either `time` or `date`, with the difference being that
`time` is intended as "simulation time" and is expected to be in seconds; `date` is a
calendar date (from `Dates.DateTime`). Conversion between time and date is performed using
the reference date and simulation starting time provided to the `DataHandler`.

The `DataHandler` has a caching mechanism in place: once a field is read and regridded, it
is stored in the local cache to be used again (if needed).

"""
module DataHandling

using Dates

import ClimaCore

import ..FileReaders: AbstractFileReader, NCFileReader, read
import ..Regridders
import ..Regridders:
    AbstractRegridder, TempestRegridder, regrid, AVAILABLE_REGRIDDERS

"""
    DataHandler{
         FR <: AbstractFileReader,
         REG <: AbstractRegridder,
         SPACE <: ClimaCore.Spaces.AbstractSpace,
         REF_DATE <: Dates.DateTime,
         TSTART <: AbstractFloat,
         DATES <: AbstractArray{<:Dates.DateTime},
         DIMS,
         TIMES <: AbstractArray{<:AbstractFloat},
         CACHE,
     }

Currently, the `DataHandler` works with one variable at the time. This might not be the most
efficiently way to tackle the problem: we should be able to reuse the interpolation weights if
multiple variables are defined on the same grid and have to be remapped to the same grid.

Assumptions:
- There is only one file with the entire time development of the given variable
- The file has well-defined physical dimensions (e.g., lat/lon)
- Currently, the time dimension has to be either "time" or "date", the spatial
  dimensions have to be lat and lon (restriction from TempestRegridder)

DataHandler is meant to live on the CPU, but the Fields can be on the GPU as well.
"""
struct DataHandler{
    FR <: AbstractFileReader,
    REG <: AbstractRegridder,
    SPACE <: ClimaCore.Spaces.AbstractSpace,
    REF_DATE <: Dates.DateTime,
    TSTART <: AbstractFloat,
    DATES <: AbstractArray{<:Dates.DateTime},
    DIMS,
    TIMES <: AbstractArray{<:AbstractFloat},
    CACHE,
}
    """Object responsible for getting the data from disk to memory"""
    file_reader::FR

    """Object responsible for resampling a rectangular grid to the simulation grid"""
    regridder::REG

    """ClimaCore Space over which the data has to be resampled"""
    target_space::SPACE

    """Tuple of linear arrays where the data is defined (typically long/lat)"""
    dimensions::DIMS

    """Calendar dates over which the data is defined"""
    available_dates::DATES

    """Simulation time at the beginning of the simulation in seconds (typically 0, but
    could be different, e.g., for restarted simulations)"""
    t_start::TSTART

    """Reference calendar date at the beginning of the simulation."""
    reference_date::REF_DATE

    """Timesteps over which the data is defined (in seconds)"""
    available_times::TIMES

    """Private field where cached data is stored"""
    _cached_regridded_fields::CACHE
end

"""
    DataHandler(file_path::AbstractString,
                varname::AbstractString,
                target_space::ClimaCore.Spaces.AbstractSpace;
                reference_date::Dates.DateTime = Dates.DateTime(1979, 1, 1),
                t_start::AbstractFloat = 0.0,
                regridder_type = :TempestRegridder)

Create a `DataHandler` to read `varname` from `file_path` and remap it to `target_space`.

The DataHandler maintains a cache of Fields that were previously computed.

TODO: Add function to clear cache, and/or CACHE_MAX_SIZE (this will probably require
developing a LRU cache scheme)

Positional arguments
=====================

- `file_path`: Path of the NetCDF file that contains the data.
- `varname`: Name of the dataset in the NetCDF that has to be read and processed.
- `target_space`: Space where the simulation is run, where the data has to be regridded to.

Keyword arguments
===================

Time/date information will be ignored for static input files. (They are still set to make
everything more type stable.)

- `reference_date`: Calendar date corresponding to the start of the simulation.
- `t_start`: Simulation time at the beginning of the simulation. Typically this is 0
             (seconds), but if might be different if the simulation was restarted.
- `regridder_type`: What type of regridding to perform. Currently, the only one implemented
                    is `:tempest` to use `TempestRemap`. `TempestRemap` regrids everything
                    ahead of time and saves the result to HDF5 files.
"""
function DataHandler(
    file_path::AbstractString,
    varname::AbstractString,
    target_space::ClimaCore.Spaces.AbstractSpace;
    reference_date::Dates.DateTime = Dates.DateTime(1979, 1, 1),
    t_start::AbstractFloat = 0.0,
    regridder_type = :TempestRegridder,
)

    # File reader, deals with ingesting data, possibly buffered/cached
    file_reader = NCFileReader(file_path, varname)

    # Regridder, deals with converting and interpolating the input data onto the simulation
    # grid
    regridder_type in AVAILABLE_REGRIDDERS ||
        error("Regridder $regridder_type not implemented")

    regridder_args = ()

    if regridder_type == :TempestRegridder
        regridder_args = (target_space, mktempdir(), varname, file_path)
    elseif regridder_type == :InterpolationsRegridder
        regridder_args = (target_space,)
    end

    RegridderConstructor = getfield(Regridders, regridder_type)
    regridder = RegridderConstructor(regridder_args...)

    # NOTE: this is not concretely typed
    _cached_regridded_fields = Dict{Dates.DateTime, ClimaCore.Fields.Field}()

    available_dates = file_reader.available_dates
    times_s = Second.(available_dates .- reference_date) ./ Second(1)
    available_times = times_s .- t_start

    return DataHandler(
        file_reader,
        regridder,
        target_space,
        file_reader.dimensions,
        available_dates,
        t_start,
        reference_date,
        available_times,
        _cached_regridded_fields,
    )
end

"""
    close(data_handler::DataHandler)

Close any file associated to the given `data_handler`.
"""
function Base.close(data_handler::DataHandler)
    close(data_handler.file_reader)
    return nothing
end

"""
    available_times(data_handler::DataHandler)

Return the time in seconds of the snapshots in the data, measured considering
the starting time of the simulation and the reference date
"""
function available_times(data_handler::DataHandler)
    return data_handler.available_times
end

"""
    available_dates(data_handler::DataHandler)

Return the dates of the snapshots in the data.
"""
function available_dates(data_handler::DataHandler)
    return data_handler.available_dates
end

"""
    previous_time(data_handler::DataHandler, time::AbstractFloat)
    previous_time(data_handler::DataHandler, date::Dates.DateTime)

Return the time in seconds of the snapshot before the given `time`.
If `time` is one of the snapshots, return itself.
"""
function previous_time(data_handler::DataHandler, time::AbstractFloat)
    time in data_handler.available_times && return time
    index = searchsortedfirst(data_handler.available_times, time) - 1
    return data_handler.available_times[index]
end

function previous_time(data_handler::DataHandler, date::Dates.DateTime)
    if date in data_handler.available_dates
        index = searchsortedfirst(data_handler.available_dates, date)
    else
        index = searchsortedfirst(data_handler.available_dates, date) - 1
    end
    return data_handler.available_times[index]
end

"""
    next_time(data_handler::DataHandler, time::AbstractFloat)
    next_time(data_handler::DataHandler, date::Dates.DateTime)

Return the time in seconds of the snapshot after the given `time`.
If `time` is one of the snapshots, return the next time.
"""
function next_time(data_handler::DataHandler, time::AbstractFloat)
    index = searchsortedfirst(data_handler.available_times, time)
    time in data_handler.available_times && (index += 1)
    return data_handler.available_times[index]
end

function next_time(data_handler::DataHandler, date::Dates.DateTime)
    if date in data_handler.available_dates
        index = searchsortedfirst(data_handler.available_dates, date) + 1
    else
        index = searchsortedfirst(data_handler.available_dates, date)
    end
    return data_handler.available_times[index]
end

"""
    regridded_snapshot(data_handler::DataHandler, date::Dates.DateTime)
    regridded_snapshot(data_handler::DataHandler, time::AbstractFloat)
    regridded_snapshot(data_handler::DataHandler)

Return the regridded snapshot from `data_handler` associated to the given `time` (if relevant).

The `time` has to be available in the `data_handler`.

`regridded_snapshot` potentially modifies the internal state of `data_handler` and it might be a very
expensive operation.

TODO: Add `regridded_snapshot!`
"""
function regridded_snapshot(data_handler::DataHandler, date::Dates.DateTime)

    # Dates.DateTime(0) is the cache key for static maps
    if date != Dates.DateTime(0)
        date in data_handler.available_dates || error("date not available")
    end

    regridder_type = nameof(typeof(data_handler.regridder))
    regrid_args = ()

    regridded_snapshot = get!(data_handler._cached_regridded_fields, date) do
        if regridder_type == :TempestRegridder
            regrid_args = (date,)
        elseif regridder_type == :InterpolationsRegridder
            regrid_args = (
                read(data_handler.file_reader, date),
                data_handler.dimensions,
            )
        else
            error("Uncaught case")
        end
        regrid(data_handler.regridder, regrid_args...)
    end

    return regridded_snapshot
end

function regridded_snapshot(data_handler::DataHandler, time::AbstractFloat)
    date =
        data_handler.reference_date +
        Second(round(data_handler.t_start)) +
        Second(round(time))
    return regridded_snapshot(data_handler, date)
end

function regridded_snapshot(data_handler::DataHandler)
    # This function can be called only when there are no dates (ie, the dataset is static)
    isempty(data_handler.available_dates) ||
        error("DataHandler is function of time")

    # In this case, we use as cache key Dates.DateTime(0)
    return regridded_snapshot(data_handler, Dates.DateTime(0))
end

end
