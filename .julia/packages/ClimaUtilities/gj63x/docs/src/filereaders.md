# [`FileReaders`](@id file_reader_module)

Reading files is a common need for most scientific projects. This can come with
a series of problems that have to be solved, from performance (accessing can be
a very computationally expensive operation), to dealing with multiple files that
are logically connected. The `FileReaders` provides an abstraction layer to
decouple the scientific needs with the technical implementation so that file
processing can be optimized and extended independently of the rest of the model.

At this point, the implemented `FileReaders` are always linked to a specific
variable and they come with a caching system to avoid unnecessary reads.

Future extensions might include:
- doing chunked reads;
- async reads.

## [`NCFileReaders`](@id ncfilereaders)

> This extension is loaded when loading `NCDatasets`

The only file reader currently implemented is the `NCFileReader`, used to read
NetCDF files. Each `NCFileReader` is associated to a collection of files
(possibly just one) and one variable (but multiple `NCFileReader`s can share the
same file). When a `NCFileReader` is constructed with multiple files, the
various files should contain the time development of the given variable.

Once created, `NCFileReader` is accessed with the `read!(file_reader, date)`
function, which returns the `Array` associated to given `date` (if available).
The `date` can be omitted if the data is static. The data is stored in a
preallocated array so it can be accessed multiple times without reallocating.

`NCFileReader`s implement two additional features: (1) optional preprocessing,
and (2) cache reads. `NCFileReader`s can be created with a `preprocessing_func`
keyword argument, function is applied to the read datasets when `read`ing.
`preprocessing_func` should be a lightweight function, such as removing `NaN`s
or changing units. Every time `read(file_reader, date)` is called, the
`NCFileReader` checks if the `date` is currently stored in the cache. If yes, it
just returns the value (without accessing the disk). If not, it reads and
process the data and adds it to the cache. This uses a least-recently-used (LRU)
cache implemented in `DataStructures`, which removes the least-recently-used
data stored in the cache when its maximum size is reached (the default max size
is 128).

It is good practice to always close the `NCFileReader`s when they are no longer
needed. The function `close_all_ncfiles` closes all the ones that are currently
open.

!!! note

    Currently, the order does not matter when passing multiple files. However, it
    is good practice to pass them in order.

### Example

Assume you have a file `era5_2000.nc`, which contains two variables `u` and `v`,
defined for the year 2000.

```julia
import ClimaUtilities.FileReaders
import NCDatasets
# Loading NCDatasets automatically loads `NCFileReaders`

u_var = FileReaders.NCFileReader("era5_2000.nc", "u")
# Change units for v
v_var = FileReaders.NCFileReader("era5_2000.nc", "u", preprocess_func = x -> 1000x)

dates = FileReaders.available_dates(u_var)
# dates is a vector of Dates.DateTime

first_date = dates[begin]

# The first time we call read, the file is accessed and read
u_array = FileReaders.read(u_var, first_date)
# As the name suggests, u_array is an Array

# All the other times, we access the cache, so no IO operation is involved
u_array_again = FileReaders.read(u_var, first_date)

close(u_var)
close(v_var)
# Alternatively: FileReaders.close_all_ncfiles()
```

Suppose now that the data is split in multiple years, we can read them as with a
single `NCFileReader` simply by passing the list of files:
```julia
u_var = FileReaders.NCFileReader(["era5_2000.nc", "era5_2001.nc", "era5_2002.nc"], "u")
```
While the order is not strictly required, it is still good practice to pass the
files in the correct order.

## API

```@docs
ClimaUtilities.FileReaders.NCFileReader
ClimaUtilities.FileReaders.read
ClimaUtilities.FileReaders.read!
ClimaUtilities.FileReaders.available_dates
ClimaUtilities.FileReaders.close_all_ncfiles
Base.close
```
