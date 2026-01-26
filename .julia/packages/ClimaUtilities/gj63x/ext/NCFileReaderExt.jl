module NCFileReaderExt

import ClimaUtilities.DataStructures
import ClimaUtilities.FileReaders

import Dates
import NCDatasets

include("nc_common.jl")

# We allow multiple NCFileReader to share the same underlying NCDataset. For this, we put
# all the NCDataset into a dictionary where we keep track of them. OPEN_NCFILES is
# dictionary that maps Vector of file paths (or a single string) to a Tuple with the first
# element being the NCDataset and the second element being a Set of Strings, the variables
# that are being read from that file. Every time a NCFileReader is created, this Set is
# modified by adding or removing a varname.
const OPEN_NCFILES =
    Dict{Union{String, Vector{String}}, Tuple{NetCDFDataset, Set{String}}}()

"""
    NCFileReader

A struct to read and process NetCDF files.

NCFileReader wants to be smart, e.g., caching reads, or spinning the I/O off to a different
thread (not implemented yet). Multiple NetCDF files can be read at the same time as long as
they can be aggregated along the time dimension.
"""
struct NCFileReader{
    VSTR <: Vector{STR} where {STR <: AbstractString},
    STR2 <: AbstractString,
    DIMS <: Tuple,
    NC <: NetCDFDataset,
    DATES <: AbstractArray{Dates.DateTime},
    PREP <: Function,
    CACHE <: DataStructures.LRUCache{Dates.DateTime, <:AbstractArray},
} <: FileReaders.AbstractFileReader
    """Path of the NetCDF file(s)"""
    file_paths::VSTR

    """Name of the dataset in the NetCDF files"""
    varname::STR2

    """A tuple of arrays with the various physical dimensions where the data is defined
    (e.g., lon/lat)"""
    dimensions::DIMS

    """A tuple with the names of the physcial dimensions, in the same order as `dimensions`"""
    dim_names::Tuple

    """A vector of DateTime collecting all the available dates in the files"""
    available_dates::DATES

    """NetCDF dataset opened by NCDataset. Don't forget to close the reader!"""
    dataset::NC

    """Optional function that is applied to the read dataset. Useful to do unit-conversion
    or remove NaNs. Not suitable for anything more complicated than that."""
    preprocess_func::PREP

    """A place where to store values that have been already read. Uses an LRU cache,
    which contains a dictionary mapping dates to arrays, and has a fixed maximum size.
    For static data sets, a sentinel data `DateTime(0)` is used as key."""
    _cached_reads::CACHE

    """Index of the time dimension in the array (typically first). -1 for static datasets"""
    time_index::Int
end

"""
    FileReaders.NCFileReader(
        file_paths,
        varname::AbstractString;
        preprocess_func = identity,
        cache_max_size:Int = 128,
    )

A struct to efficiently read and process NetCDF files.

When more than one file is passed, the files should contain the time development of one or
multiple variables. Files are joined along the time dimension.

## Argument

`file_paths` can be a string, or a collection of paths to files that contain the
same variables but at different times.

"""
function FileReaders.NCFileReader(
    file_paths,
    varname::AbstractString;
    preprocess_func = identity,
    cache_max_size::Int = 128,
)
    # file_paths could be a vector/tuple or a string. Let's start by standarizing to a
    # vector
    file_paths isa AbstractString && (file_paths = [file_paths])
    only_one_file = length(file_paths) == 1

    # If we have more than one file, we have to aggregate them
    aggtime_kwarg = ()
    if !only_one_file
        # Let's first try to identify the time dimension, if it exists. To do that, we open the
        # first dataset. We need this to aggregate multiple datasets, if data is split across
        # multiple files
        NCDatasets.NCDataset(first(file_paths)) do first_dataset
            is_time = x -> x == "time" || x == "date" || x == "t"
            time_dims = filter(is_time, NCDatasets.dimnames(first_dataset))
            if !isempty(time_dims)
                # When loading multifile dataset using NCDatasets.jl, the NetCDF files are
                # not kept open due to the common limitation of 1024 open files per user on
                # Linux. However, for our use case, we will not reach this limit. Hence, we
                # keep the files open with deferopen = false.
                # See: https://github.com/JuliaGeo/NCDatasets.jl/issues/277
                aggtime_kwarg =
                    (:aggdim => first(time_dims), :deferopen => false)
            else
                error(
                    "Multiple files given, but no temporal dimension found. Combining multiple files is only possible along the temporal dimension.",
                )
            end
        end
    end

    # When we have only no time data, we have to pass this as a string
    file_path_to_ncdataset =
        isempty(aggtime_kwarg) ? first(file_paths) : file_paths

    # Get dataset from global dictionary. If not available, open the new dataset and add
    # entry to global dictionary
    dataset, open_varnames = get!(
        OPEN_NCFILES,
        file_paths,  # We map the collection of files to the dataset
        (
            NCDatasets.NCDataset(file_path_to_ncdataset; aggtime_kwarg...),
            Set([varname]),
        ),
    )
    # push! will do nothing when file is opened for the first time
    push!(open_varnames, varname)

    available_dates = read_available_dates(dataset)

    time_index = -1

    dim_names = NCDatasets.dimnames(dataset[varname])

    if !isempty(available_dates)
        is_time = x -> x == "time" || x == "date" || x == "t"

        time_index_vec = findall(is_time, dim_names)
        length(time_index_vec) == 1 ||
            error("Could not find (unique) time dimension")
        time_index = time_index_vec[]

        issorted(available_dates) || error(
            "Cannot process files that are not sorted in time in ($file_paths)",
        )

        # Remove time from the dim names
        dim_names = filter(!is_time, dim_names)
    end

    if all(d in keys(dataset) for d in dim_names)
        dimensions =
            Tuple(NCDatasets.nomissing(Array(dataset[d])) for d in dim_names)
    else
        error(
            "$file_paths does not contain information about dimensions $(filter(!in(keys(dataset)), dim_names))",
        )
    end

    # Use an LRU cache to store regridded fields
    _cached_reads = DataStructures.LRUCache{Dates.DateTime, Array}(
        max_size = cache_max_size,
    )

    return NCFileReader(
        file_paths,
        varname,
        dimensions,
        dim_names,
        available_dates,
        dataset,
        preprocess_func,
        _cached_reads,
        time_index,
    )
end

"""
    close(file_reader::NCFileReader)

Close `NCFileReader`. If no other `NCFileReader` is using the same file, close the NetCDF file.
"""
function Base.close(file_reader::NCFileReader)
    # If we don't have the key, we don't have to do anything (we already closed
    # the file)
    files_are_not_open = !haskey(OPEN_NCFILES, file_reader.file_paths)
    files_are_not_open && return nothing

    open_variables = OPEN_NCFILES[file_reader.file_paths][end]
    pop!(open_variables, file_reader.varname)
    if isempty(open_variables)
        NCDatasets.close(file_reader.dataset)
        delete!(OPEN_NCFILES, file_reader.file_paths)
    end
    return nothing
end

"""
    close_all_ncfiles()

Close all the `NCFileReader` currently open.
"""
function FileReaders.close_all_ncfiles()
    foreach(OPEN_NCFILES) do (_, ds_vars)
        NCDatasets.close(ds_vars[begin])
    end
    empty!(OPEN_NCFILES)
    return nothing
end

"""
    read(file_reader::NCFileReader, date::Dates.DateTime)

Read and preprocess the data at the given `date`.
"""
function FileReaders.read(file_reader::NCFileReader, date::Dates.DateTime)
    # For cache hits, return a copy to give away ownership of the data (if we were to just
    # return _cached_reads[date], modifying the return value would modify the private state
    # of the file reader)

    if haskey(file_reader._cached_reads, date)
        return copy(file_reader._cached_reads[date])
    end

    # DateTime(0) is the sentinel value for static datasets
    if date == Dates.DateTime(0)
        return get!(file_reader._cached_reads, date) do
            file_reader.preprocess_func.(
                Array(file_reader.dataset[file_reader.varname]),
            )
        end
    end

    index = findall(d -> d == date, file_reader.available_dates)
    length(index) == 1 ||
        error("Problem with date $date in one of $(file_reader.file_paths)")
    index = index[]

    var = file_reader.dataset[file_reader.varname]
    slicer = [
        i == file_reader.time_index ? index : Colon() for
        i in 1:length(NCDatasets.dimnames(var))
    ]
    return file_reader.preprocess_func.(
        file_reader.dataset[file_reader.varname][slicer...],
    )
end

"""
    available_dates(file_reader::NCFileReader)

Returns the dates in the given file.
"""
function FileReaders.available_dates(file_reader::NCFileReader)
    return file_reader.available_dates
end

"""
    read(file_reader::NCFileReader)

Read and preprocess data (for static datasets).
"""
function FileReaders.read(file_reader::NCFileReader)
    isempty(file_reader.available_dates) ||
        error("File contains temporal data, date required")

    # When there's no dates, we use DateTime(0) as key
    return get!(file_reader._cached_reads, Dates.DateTime(0)) do
        return file_reader.preprocess_func.(
            Array(file_reader.dataset[file_reader.varname]),
        )
    end
end

"""
    read!(dest, file_reader::NCFileReader)

Read and preprocess data (for static datasets), saving the output to `dest`.
"""
function FileReaders.read!(dest, file_reader::NCFileReader)
    dest .= FileReaders.read(file_reader)
    return nothing
end

"""
    read!(dest, file_reader::NCFileReader, date::Dates.DateTime)

Read and preprocess the data at the given `date`, saving the output to `dest`.
"""
function FileReaders.read!(
    dest,
    file_reader::NCFileReader,
    date::Dates.DateTime,
)
    dest .= FileReaders.read(file_reader, date)
    return nothing
end

end
