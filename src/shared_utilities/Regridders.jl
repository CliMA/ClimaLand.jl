"""
    Regridders

The `Regridders` module implement structs and functions to remap datasets to simulation
grids.

Currently, the schemes implemented are `TempestRegridder`, which uses
`ClimaCoreTempestRemap`, and `InterpolationsRegridder`, which uses `Interpolations.jl`.

The key function exposed by `Regridder` is the `regrid` method.
"""
module Regridders

using Dates

import Adapt

import ClimaComms
import ClimaCore
import ClimaCore: Spaces, Fields, Geometry
import ClimaCoreTempestRemap
import Interpolations as Intp
import ..FileReaders: read_available_dates

using NCDatasets

# When adding a new regridder, you also have to change some functions in the DataHandler
# module. Find where :TempestRegridder is used.
abstract type AbstractRegridder end

const AVAILABLE_REGRIDDERS = [:TempestRegridder, :InterpolationsRegridder]

"""
    InterpolationsRegridder

An online regridder that uses Interpolations.jl

InterpolationsRegridder is only implemented for LatLong and LatLongZ spaces. It performs
linear interpolation along each of the directions (separately), while imposing periodic
boundary conditions for longitude, flat for latitude, and throwing errors when extrapolating
in z.

InterpolationsRegridder is GPU and MPI compatible in the simplest possible way. Each MPI
process has the entire input data and everything is copied to GPU.
"""
struct InterpolationsRegridder{
    SPACE <: ClimaCore.Spaces.AbstractSpace,
    FIELD <: ClimaCore.Fields.Field,
    BC,
} <: AbstractRegridder
    target_space::SPACE
    coordinates::FIELD
    extrapolation_bc::BC
end

# Note, we swap Lat and Long! This is because according to the CF conventions longitude
# should be first, so files will have longitude as first dimension.
totuple(pt::Geometry.LatLongZPoint) = pt.long, pt.lat, pt.z
totuple(pt::Geometry.LatLongPoint) = pt.long, pt.lat

function InterpolationsRegridder(target_space::ClimaCore.Spaces.AbstractSpace)
    coordinates = Fields.coordinate_field(target_space)

    extrapolation_bc = ()
    if eltype(coordinates) <: Geometry.LatLongPoint
        extrapolation_bc = (Intp.Periodic(), Intp.Flat())
    elseif eltype(coordinates) <: Geometry.LatLongZPoint
        extrapolation_bc = (Intp.Periodic(), Intp.Flat(), Intp.Throw())
    else
        error("Only lat-long, lat-long-z spaces are supported")
    end

    return InterpolationsRegridder(target_space, coordinates, extrapolation_bc)
end

"""
    regrid(regridder::InterpolationsRegridder, data, dimensions)::Field

Regrid the given data as defined on the given dimensions to the `target_space` in `regridder`.

This function is allocating.
"""
function regrid(regridder::InterpolationsRegridder, data, dimensions)
    FT = Spaces.undertype(regridder.target_space)
    dimensions_FT = map(d -> FT.(d), dimensions)

    # Make a linear spline
    itp = Intp.extrapolate(
        Intp.interpolate(dimensions_FT, FT.(data), Intp.Gridded(Intp.Linear())),
        regridder.extrapolation_bc,
    )

    # Move it to GPU (if needed)
    gpuitp = Adapt.adapt(ClimaComms.array_type(regridder.target_space), itp)

    return map(regridder.coordinates) do coord
        gpuitp(totuple(coord)...)
    end
end

"""
    TempestRegridder

A struct to resample a given NetCDF file to `the target_space`.

Currently, this work with ClimaCoreTempestRemap. `ClimaCoreTempestRemap` uses TempestRemap,
a command-line utility. Hence, This function has to process files, so, for this regridder,
we cannot really split file processing and remapping.

TempestRegridder only works on CPU and on a single process.

Implicit assumptions in the code:
- input_file contains "lat" and "lon"
- the time dimension, if available is called "time" or "date"
- the name of the variable in the NetCDF file is `varname`
- all the data is in one single file
"""
struct TempestRegridder{
    SPACE <: ClimaCore.Spaces.AbstractSpace,
    STR1 <: AbstractString,
    STR2 <: AbstractString,
    STR3 <: AbstractString,
} <: AbstractRegridder
    target_space::SPACE
    regrid_dir::STR1
    outfile_root::STR2
    varname::STR3
    mono::Bool
end

"""
    TempestRegridder(target_space::ClimaCore.Spaces.AbstractSpace,
                     input_file::AbstractString,
                     varname::AbstractString,
                     regrid_dir::AbstractString;
                     mono::Bool)

Set up a `TempestRegridder` object to regrid the variable in the given `input_file` to the
`target_space`. `TempestRegridder` works only on CPU and on a single process.

Positional arguments
=====================

- `target_space`: the ClimaCore Space where the simulation is being performed.
- `input_file`: the path of the NetCDF file that has to be read and processed.
- `varname`: the name of the variable in the NetCDF file.
- `regrid_dir`: a folder where regridded Fields are saved as HDF5 files.

Keyword arguments
==================

- `mono`: Whether remapping has to be monotonic or not.
"""
function TempestRegridder(
    target_space::ClimaCore.Spaces.AbstractSpace,
    regrid_dir::AbstractString,
    varname::AbstractString,
    input_file::AbstractString;
    mono::Bool = true,
)
    space = target_space
    comms_ctx = ClimaComms.context(space)

    if ClimaComms.iamroot(comms_ctx)
        @info "Saving TempestRegrid files to $regrid_dir"

        FT = Spaces.undertype(space)
        outfile_root = varname
        varnames = [varname]

        hdwrite_regridfile_rll_to_cgll(
            FT,
            regrid_dir,
            input_file,
            varnames,
            target_space,
            outfile_root;
            mono,
        )
    end
    # We have to make sure that all the processes wait on the (expensive) regridding,
    # otherwise they won't have access to the regridded fields
    ClimaComms.barrier(comms_ctx)

    return TempestRegridder{
        typeof(target_space),
        typeof(regrid_dir),
        typeof(outfile_root),
        typeof(varname),
    }(
        target_space,
        regrid_dir,
        outfile_root,
        varname,
        mono,
    )
end

"""
    regrid(regridder::TempestRegridder, date::Dates.DateTime)

Return the field associated to the `regridder` at the given `date`.
"""
function regrid(regridder::TempestRegridder, date::Dates.DateTime)
    return read_from_hdf5(
        regridder.regrid_dir,
        regridder.outfile_root,
        date,
        regridder.varname,
        regridder.target_space,
    )
end

export hdwrite_regridfile_rll_to_cgll

"""
    reshape_cgll_sparse_to_field!(field::Fields.Field, in_array::Array, R)

Reshapes a sparse vector array `in_array` (CGLL, raw output of the TempestRemap),
and uses its data to populate the input Field object `field`.
Redundant nodes are populated using `dss` operations.

Code taken from ClimaCoupler.Regridder.

# Arguments
- `field`: [Fields.Field] object populated with the input array.
- `in_array`: [Array] input used to fill `field`.
- `R`: [NamedTuple] containing `target_idxs` and `row_indices` used for indexing.
"""
function reshape_cgll_sparse_to_field!(
    field::ClimaCore.Fields.Field,
    in_array::Array,
    R,
)
    field_array = parent(field)

    fill!(field_array, zero(eltype(field_array)))
    Nf = size(field_array, 3)

    for (n, row) in enumerate(R.row_indices)
        it, jt, et = (
            view(R.target_idxs[1], n),
            view(R.target_idxs[2], n),
            view(R.target_idxs[3], n),
        )
        for f in 1:Nf
            field_array[it, jt, f, et] .= in_array[row]
        end
    end

    # broadcast to the redundant nodes using unweighted dss
    space = axes(field)
    topology = ClimaCore.Spaces.topology(space)
    hspace = ClimaCore.Spaces.horizontal_space(space)
    target = ClimaCore.Fields.field_values(field)

    ClimaCore.Topologies.dss!(target, topology)
end

"""
    swap_space(field, new_space)

Update the space of a ClimaCore.Fields.Field object. Returns a new Field
object with the same values as the original field, but on the new space.
This is needed to correctly read in regridded files that may be reused
between simulations.

# Arguments
- `field`: The ClimaCore.Fields.Field object to swap the space of.
- `new_space`: The new ClimaCore.Spaces.AbstractSpace to assign to the field.
"""
function swap_space(field, new_space)
    return ClimaCore.Fields.Field(
        ClimaCore.Fields.field_values(field),
        new_space,
    )
end

"""
    read_from_hdf5(REGIRD_DIR, hd_outfile_root, time, varname,
        space)

Read in a variable `varname` from an HDF5 file onto the provided space.
If a CommsContext other than SingletonCommsContext is used in the `space`,
the input HDF5 file must be readable by multiple MPI processes.

Code taken from ClimaCoupler.Regridder.

# Arguments
- `REGRID_DIR`: [String] directory to save output files in.
- `hd_outfile_root`: [String] root of the output file name.
- `time`: [Dates.DateTime] the timestamp of the data being written.
- `varname`: [String] variable name of data.
- `space`: [ClimaCore.Spaces.AbstractSpace] to read the HDF5 file onto.
# Returns
- Field or FieldVector
"""
function read_from_hdf5(REGRID_DIR, hd_outfile_root, time, varname, space)
    comms_ctx = ClimaComms.context(space)
    # Include time component in HDF5 reader name if it's a valid date
    if !(time == Dates.DateTime(0))
        hdfreader_path = joinpath(
            REGRID_DIR,
            hd_outfile_root * "_" * varname * "_" * string(time) * ".hdf5",
        )
    else
        hdfreader_path =
            joinpath(REGRID_DIR, hd_outfile_root * "_" * varname * ".hdf5")
    end
    hdfreader = ClimaCore.InputOutput.HDF5Reader(hdfreader_path, comms_ctx)

    field = ClimaCore.InputOutput.read_field(hdfreader, varname)
    Base.close(hdfreader)

    # Ensure the field is on the correct space when reusing regridded files
    #  between simulations
    return swap_space(field, space)
end


"""
    write_to_hdf5(REGRID_DIR, hd_outfile_root, time, field, varname,
        comms_ctx = ClimaComms.SingletonCommsContext())
Function to save individual HDF5 files after remapping.
If a CommsContext other than SingletonCommsContext is used for `comms_ctx`,
the HDF5 output is readable by multiple MPI processes.

Code taken from ClimaCoupler.Regridder.


# Arguments
- `REGRID_DIR`: [String] directory to save output files in.
- `hd_outfile_root`: [String] root of the output file name.
- `time`: [Dates.DateTime] the timestamp of the data being written.
- `field`: [Fields.Field] object to be written.
- `varname`: [String] variable name of data.
- `comms_ctx`: [ClimaComms.AbstractCommsContext] context used for this operation.
"""
function write_to_hdf5(
    REGRID_DIR,
    hd_outfile_root,
    time,
    field,
    varname,
    comms_ctx = ClimaComms.SingletonCommsContext(),
)
    # Include time component in HDF5 writer name, and write time to file if it's a valid date
    if !(time == Dates.DateTime(0))
        hdfwriter = ClimaCore.InputOutput.HDF5Writer(
            joinpath(
                REGRID_DIR,
                hd_outfile_root * "_" * varname * "_" * string(time) * ".hdf5",
            ),
            comms_ctx,
        )

        t = Dates.datetime2unix.(time)
        ClimaCore.InputOutput.HDF5.write_attribute(
            hdfwriter.file,
            "unix time",
            t,
        ) # TODO: a better way to write metadata, CMIP convention
    else
        hdfwriter = ClimaCore.InputOutput.HDF5Writer(
            joinpath(REGRID_DIR, hd_outfile_root * "_" * varname * ".hdf5"),
            comms_ctx,
        )
    end
    ClimaCore.InputOutput.write!(hdfwriter, field, string(varname))
    Base.close(hdfwriter)
end


"""
    function hdwrite_regridfile_rll_to_cgll(
        FT,
        REGRID_DIR,
        datafile_rll,
        varnames,
        space,
        outfile_root;
        mono = false,
    )
Reads and regrids data of all `varnames` variables from an input NetCDF file and
saves it as another NetCDF file using Tempest Remap.
The input NetCDF fileneeds to be `Exodus` formatted, and can contain
time-dependent data. The output NetCDF file is then read back, the output
arrays converted into Fields and saved as HDF5 files (one per time slice).
This function should be called by the root process.
The saved regridded HDF5 output is readable by multiple MPI processes.
Assumes that all variables specified by `varnames` have the same dates
and grid.

Code taken from ClimaCoupler.Regridder.


# Arguments
- `FT`: [DataType] Float type.
- `REGRID_DIR`: [String] directory to save output files in.
- `datafile_rll`: [String] filename of RLL dataset to be mapped to CGLL.
- `varnames`: [Vector{String}] the name of the variable to be remapped.
- `space`: [ClimaCore.Spaces.AbstractSpace] the space to which we are mapping.
- `outfile_root`: [String] root of the output file name.
- `mono`: [Bool] flag to specify monotone remapping.
"""
function hdwrite_regridfile_rll_to_cgll(
    FT,
    REGRID_DIR,
    datafile_rll,
    varnames::Vector{String},
    space,
    outfile_root;
    mono = false,
)
    out_type = "cgll"

    datafile_cgll = joinpath(REGRID_DIR, outfile_root * ".g")
    meshfile_rll = joinpath(REGRID_DIR, outfile_root * "_mesh_rll.g")
    meshfile_cgll = joinpath(REGRID_DIR, outfile_root * "_mesh_cgll.g")
    meshfile_overlap = joinpath(REGRID_DIR, outfile_root * "_mesh_overlap.g")
    weightfile = joinpath(REGRID_DIR, outfile_root * "_remap_weights.nc")

    # If doesn't make sense to regrid with GPUs/MPI processes
    cpu_context =
        ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())

    # Note: this topology gives us `space == space_undistributed` in the undistributed
    # case (as desired), which wouldn't hold if we used `spacefillingcurve` here.
    topology = ClimaCore.Topologies.Topology2D(
        cpu_context,
        ClimaCore.Spaces.topology(space).mesh,
    )
    Nq =
        ClimaCore.Spaces.Quadratures.polynomial_degree(
            ClimaCore.Spaces.quadrature_style(space),
        ) + 1
    space_undistributed = ClimaCore.Spaces.SpectralElementSpace2D(
        topology,
        ClimaCore.Spaces.Quadratures.GLL{Nq}(),
    )

    if isfile(datafile_cgll) == false
        isdir(REGRID_DIR) ? nothing : mkpath(REGRID_DIR)

        nlat, nlon = NCDatasets.NCDataset(datafile_rll) do ds
            (ds.dim["lat"], ds.dim["lon"])
        end
        # write lat-lon mesh
        ClimaCoreTempestRemap.rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

        # write cgll mesh, overlap mesh and weight file
        ClimaCoreTempestRemap.write_exodus(meshfile_cgll, topology)
        ClimaCoreTempestRemap.overlap_mesh(
            meshfile_overlap,
            meshfile_rll,
            meshfile_cgll,
        )

        # 'in_np = 1' and 'mono = true' arguments ensure mapping is conservative and monotone
        # Note: for a kwarg not followed by a value, set it to true here (i.e. pass 'mono = true' to produce '--mono')
        # Note: out_np = degrees of freedom = polynomial degree + 1
        kwargs = (; out_type = out_type, out_np = Nq)
        kwargs = mono ? (; (kwargs)..., in_np = 1, mono = mono) : kwargs
        ClimaCoreTempestRemap.remap_weights(
            weightfile,
            meshfile_rll,
            meshfile_cgll,
            meshfile_overlap;
            kwargs...,
        )

        ClimaCoreTempestRemap.apply_remap(
            datafile_cgll,
            datafile_rll,
            weightfile,
            varnames,
        )
    else
        @warn "Using the existing $datafile_cgll : check topology is consistent"
    end

    # weightfile info needed to populate all nodes and save into fields with
    #  sparse matrices
    _, _, row_indices = NCDataset(weightfile, "r") do ds_wt
        (Array(ds_wt["S"]), Array(ds_wt["col"]), Array(ds_wt["row"]))
    end

    target_unique_idxs =
        out_type == "cgll" ?
        collect(ClimaCore.Spaces.unique_nodes(space_undistributed)) :
        collect(ClimaCore.Spaces.all_nodes(space_undistributed))
    target_unique_idxs_i =
        map(row -> target_unique_idxs[row][1][1], row_indices)
    target_unique_idxs_j =
        map(row -> target_unique_idxs[row][1][2], row_indices)
    target_unique_idxs_e = map(row -> target_unique_idxs[row][2], row_indices)
    target_unique_idxs =
        (target_unique_idxs_i, target_unique_idxs_j, target_unique_idxs_e)

    R = (; target_idxs = target_unique_idxs, row_indices = row_indices)

    offline_field = ClimaCore.Fields.zeros(FT, space_undistributed)

    times = [DateTime(0)]
    # Save regridded HDF5 file for each variable in `varnames`
    for varname in varnames
        # read the remapped file with sparse matrices
        offline_outvector, times =
            NCDatasets.NCDataset(datafile_cgll, "r") do ds_wt
                (
                    # read the data in, and remove missing type (will error if missing data is present)
                    offline_outvector = nomissing(Array(ds_wt[varname])[:, :]), # ncol, times
                    times = read_available_dates(ds_wt),
                )
            end

        # Use dummy date when there are no real dates
        isempty(times) && (times = [DateTime(0)])

        # Convert input data float type if needed
        if eltype(offline_outvector) <: AbstractFloat &&
           eltype(offline_outvector) != FT
            @warn "Converting $varname data in $datafile_cgll from $(eltype(offline_outvector)) to $FT"
            offline_outvector = Array{FT}(offline_outvector)
        end
        length_times = length(times)
        size_outvector = size(offline_outvector, 2)
        @assert length_times == size_outvector "Inconsistent time dimension in $datafile_cgll for $varname ($length_times != $size_outvector)"

        offline_fields = ntuple(x -> similar(offline_field), length(times))
        ntuple(
            x -> reshape_cgll_sparse_to_field!(
                offline_fields[x],
                offline_outvector[:, x],
                R,
            ),
            length(times),
        )

        map(
            x -> write_to_hdf5(
                REGRID_DIR,
                outfile_root,
                times[x],
                offline_fields[x],
                varname,
                cpu_context,
            ),
            1:length(times),
        )
    end
    return times
end

end
