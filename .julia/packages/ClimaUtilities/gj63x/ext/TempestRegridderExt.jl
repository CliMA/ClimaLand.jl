module TempestRegridderExt

using Dates

import ClimaCoreTempestRemap
import ClimaCoreTempestRemap: ClimaCore
import ClimaCoreTempestRemap: ClimaComms
import ClimaCoreTempestRemap: NCDatasets

import ClimaUtilities.Regridders
import ClimaUtilities.FileReaders

include("nc_common.jl")

"""
    TempestRegridder

A struct to resample a given NetCDF file to `the target_space`.

Currently, this work with ClimaCoreTempestRemap. `ClimaCoreTempestRemap` uses TempestRemap,
a command-line utility. Hence, This function has to process files, so, for this regridder,
we cannot really split file processing and remapping.

TempestRegridder only works on CPU and on a single process, and for a non-3D target space.

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
} <: Regridders.AbstractRegridder
    target_space::SPACE
    regrid_dir::STR1
    varname::STR2
    mono::Bool
end

"""
    TempestRegridder(target_space::ClimaCore.Spaces.AbstractSpace,
                     varname::AbstractString,
                     input_file::AbstractString;
                     regrid_dir::AbstractString,
                     mono::Bool)

Set up a `TempestRegridder` object to regrid the variable in the given `input_file` to the
`target_space`. `TempestRegridder` works only on CPU and on a single process.

Positional arguments
=====================

- `target_space`: the ClimaCore Space where the simulation is being performed.
- `input_file`: the path of the NetCDF file that has to be read and processed.
- `regrid_dir`: the path where to save the regrid files created by TempestRemap.
- `varname`: the name of the variable in the NetCDF file.
- `regrid_dir`: a folder where regridded Fields are saved as HDF5 files.

Keyword arguments
==================

- `mono`: Whether remapping has to be monotonic or not.
"""
function Regridders.TempestRegridder(
    target_space::ClimaCore.Spaces.AbstractSpace,
    varname::AbstractString,
    input_file::AbstractString;
    regrid_dir::AbstractString,
    mono::Bool = true,
)
    space = target_space
    comms_ctx = ClimaComms.context(space)

    if ClimaComms.iamroot(comms_ctx)
        @info "Saving TempestRegrid files to $regrid_dir"

        FT = ClimaCore.Spaces.undertype(space)
        varnames = [varname]

        hdwrite_regridfile_rll_to_cgll(
            FT,
            regrid_dir,
            input_file,
            varnames,
            target_space;
            mono,
        )
    end
    # We have to make sure that all the processes wait on the (expensive) regridding,
    # otherwise they won't have access to the regridded fields
    ClimaComms.barrier(comms_ctx)

    return TempestRegridder{
        typeof(target_space),
        typeof(regrid_dir),
        typeof(varname),
    }(
        target_space,
        regrid_dir,
        varname,
        mono,
    )
end

"""
    regrid(regridder::TempestRegridder, date::Dates.DateTime)

Return the field associated to the `regridder` at the given `date`.
"""
function Regridders.regrid(regridder::TempestRegridder, date::Dates.DateTime)
    return read_from_hdf5(
        regridder.regrid_dir,
        date,
        regridder.varname,
        regridder.target_space,
    )
end

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
    read_from_hdf5(REGIRD_DIR, time, varname,
        space)

Read in a variable `varname` from an HDF5 file onto the provided space.
If a CommsContext other than SingletonCommsContext is used in the `space`,
the input HDF5 file must be readable by multiple MPI processes.

Code taken from ClimaCoupler.Regridder.

# Arguments
- `REGRID_DIR`: [String] directory to save output files in.
- `time`: [Dates.DateTime] the timestamp of the data being written.
- `varname`: [String] variable name of data.
- `space`: [ClimaCore.Spaces.AbstractSpace] to read the HDF5 file onto.
# Returns
- Field or FieldVector
"""
function read_from_hdf5(REGRID_DIR, time, varname, space)
    comms_ctx = ClimaComms.context(space)
    # Include time component in HDF5 reader name if it's a valid date
    if !(time == Dates.DateTime(0))
        hdfreader_path =
            joinpath(REGRID_DIR, varname * "_" * string(time) * ".hdf5")
    else
        hdfreader_path = joinpath(REGRID_DIR, varname * ".hdf5")
    end
    hdfreader = ClimaCore.InputOutput.HDF5Reader(hdfreader_path, comms_ctx)

    field = ClimaCore.InputOutput.read_field(hdfreader, varname)
    Base.close(hdfreader)

    # Ensure the field is on the correct space
    @assert ClimaCore.Spaces.topology(axes(field)) ==
            ClimaCore.Spaces.topology(space)
    @assert ClimaCore.Spaces.quadrature_style(axes(field)) ==
            ClimaCore.Spaces.quadrature_style(space)
    @assert ClimaCore.Spaces.global_geometry(axes(field)) ==
            ClimaCore.Spaces.global_geometry(space)

    # Since the spaces aren't the same, we need to copy the field values onto the space
    return swap_space(field, space)
end

"""
    write_to_hdf5(REGRID_DIR, time, field, varname,
        comms_ctx = ClimaComms.SingletonCommsContext())
Function to save individual HDF5 files after remapping.
If a CommsContext other than SingletonCommsContext is used for `comms_ctx`,
the HDF5 output is readable by multiple MPI processes.

Code taken from ClimaCoupler.Regridder.


# Arguments
- `REGRID_DIR`: [String] directory to save output files in.
- `time`: [Dates.DateTime] the timestamp of the data being written.
- `field`: [Fields.Field] object to be written.
- `varname`: [String] variable name of data.
- `comms_ctx`: [ClimaComms.AbstractCommsContext] context used for this operation.
"""
function write_to_hdf5(
    REGRID_DIR,
    time,
    field,
    varname,
    comms_ctx = ClimaComms.SingletonCommsContext(),
)
    # Include time component in HDF5 writer name, and write time to file if it's a valid date
    if !(time == Dates.DateTime(0))
        hdfwriter = ClimaCore.InputOutput.HDF5Writer(
            joinpath(REGRID_DIR, varname * "_" * string(time) * ".hdf5"),
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
            joinpath(REGRID_DIR, varname * ".hdf5"),
            comms_ctx,
        )
    end
    ClimaCore.InputOutput.write!(hdfwriter, field, string(varname))
    Base.close(hdfwriter)
end

"""
    construct_singleton_space(space)

Constructs a singleton space from a given space. This is not ideal, but
necessary to use TempestRemap regridding, which is not MPI or GPU compatible.
"""
function construct_singleton_space(space)
    # If doesn't make sense to regrid with GPUs/MPI processes
    cpu_context =
        ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())

    # Check if input space was constructed using `spacefillingcurve`
    use_spacefillingcurve =
        ClimaCore.Spaces.topology(space).elemorder isa CartesianIndices ?
        false : true

    mesh = ClimaCore.Spaces.topology(space).mesh
    topology = nothing
    if use_spacefillingcurve
        topology = ClimaCore.Topologies.Topology2D(
            cpu_context,
            mesh,
            ClimaCore.Topologies.spacefillingcurve(mesh),
        )
    else
        topology = ClimaCore.Topologies.Topology2D(cpu_context, mesh)
    end
    Nq =
        ClimaCore.Spaces.Quadratures.polynomial_degree(
            ClimaCore.Spaces.quadrature_style(space),
        ) + 1
    space_singleton = ClimaCore.Spaces.SpectralElementSpace2D(
        topology,
        ClimaCore.Spaces.Quadratures.GLL{Nq}(),
    )

    return space_singleton
end

"""
    function hdwrite_regridfile_rll_to_cgll(
        FT,
        REGRID_DIR,
        datafile_rll,
        varnames,
        space,
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
- `mono`: [Bool] flag to specify monotone remapping.
"""
function hdwrite_regridfile_rll_to_cgll(
    FT,
    REGRID_DIR,
    datafile_rll,
    varnames::Vector{String},
    space;
    mono = false,
)
    out_type = "cgll"

    datafile_cgll = joinpath(REGRID_DIR, "datafile.g")
    meshfile_rll = joinpath(REGRID_DIR, "mesh_rll.g")
    meshfile_cgll = joinpath(REGRID_DIR, "mesh_cgll.g")
    meshfile_overlap = joinpath(REGRID_DIR, "mesh_overlap.g")
    weightfile = joinpath(REGRID_DIR, "remap_weights.nc")

    # If doesn't make sense to regrid with GPUs/MPI processes
    space_singleton = construct_singleton_space(space)
    topology_singleton = ClimaCore.Spaces.topology(space_singleton)
    cpu_context = ClimaCore.Spaces.topology(space_singleton).context

    if isfile(datafile_cgll) == false
        isdir(REGRID_DIR) ? nothing : mkpath(REGRID_DIR)

        nlat, nlon = NCDatasets.NCDataset(datafile_rll) do ds
            (ds.dim["lat"], ds.dim["lon"])
        end
        # write lat-lon mesh
        ClimaCoreTempestRemap.rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

        # write cgll mesh, overlap mesh and weight file
        ClimaCoreTempestRemap.write_exodus(meshfile_cgll, topology_singleton)
        ClimaCoreTempestRemap.overlap_mesh(
            meshfile_overlap,
            meshfile_rll,
            meshfile_cgll,
        )

        # 'in_np = 1' and 'mono = true' arguments ensure mapping is conservative and monotone
        # Note: for a kwarg not followed by a value, set it to true here (i.e. pass 'mono = true' to produce '--mono')
        # Note: out_np = degrees of freedom = polynomial degree + 1
        Nq =
            ClimaCore.Spaces.Quadratures.polynomial_degree(
                ClimaCore.Spaces.quadrature_style(space),
            ) + 1
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
    _, _, row_indices = NCDatasets.NCDataset(weightfile, "r") do ds_wt
        (Array(ds_wt["S"]), Array(ds_wt["col"]), Array(ds_wt["row"]))
    end

    target_unique_idxs =
        out_type == "cgll" ?
        collect(ClimaCore.Spaces.unique_nodes(space_singleton)) :
        collect(ClimaCore.Spaces.all_nodes(space_singleton))
    target_unique_idxs_i =
        map(row -> target_unique_idxs[row][1][1], row_indices)
    target_unique_idxs_j =
        map(row -> target_unique_idxs[row][1][2], row_indices)
    target_unique_idxs_e = map(row -> target_unique_idxs[row][2], row_indices)
    target_unique_idxs =
        (target_unique_idxs_i, target_unique_idxs_j, target_unique_idxs_e)

    R = (; target_idxs = target_unique_idxs, row_indices = row_indices)

    offline_field = ClimaCore.Fields.zeros(FT, space_singleton)

    times = [DateTime(0)]
    # Save regridded HDF5 file for each variable in `varnames`
    for varname in varnames
        # read the remapped file with sparse matrices
        offline_outvector, times =
            NCDatasets.NCDataset(datafile_cgll, "r") do ds_wt
                (
                    # read the data in, and remove missing type (will error if missing data is present)
                    offline_outvector = NCDatasets.nomissing(
                        Array(ds_wt[varname])[:, :],
                    ), # ncol, times
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

        # Save regridded HDF5 file for each time slice
        map(
            x -> write_to_hdf5(
                REGRID_DIR,
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
