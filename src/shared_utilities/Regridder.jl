# Code from ClimaCoupler Regridder.jl
module Regridder
using ClimaCore
using ClimaComms
using NCDatasets
using ClimaCoreTempestRemap
using Dates
using JLD2
using DocStringExtensions

export MapInfo, regrid_netcdf_to_field

nans_to_zero(v::T) where {T} = isnan(v) ? T(0) : v


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

    ClimaCore.Spaces.dss2!(target, topology, hspace.quadrature_style)
end

"""
    read_from_hdf5(REGIRD_DIR, hd_outfile_root, time, varname,
        comms_ctx = ClimaComms.SingletonCommsContext())
Read in a variable `varname` from an HDF5 file.
If a CommsContext other than SingletonCommsContext is used for `comms_ctx`,
the input HDF5 file must be readable by multiple MPI processes.

Code taken from ClimaCoupler.Regridder.

# Arguments
- `REGRID_DIR`: [String] directory to save output files in.
- `hd_outfile_root`: [String] root of the output file name.
- `time`: [Dates.DateTime] the timestamp of the data being written.
- `varname`: [String] variable name of data.
- `comms_ctx`: [ClimaComms.AbstractCommsContext] context used for this operation.
# Returns
- Field or FieldVector
"""
function read_from_hdf5(
    REGRID_DIR,
    hd_outfile_root,
    time,
    varname,
    comms_ctx = ClimaComms.SingletonCommsContext(),
)
    hdfreader = ClimaCore.InputOutput.HDF5Reader(
        joinpath(REGRID_DIR, hd_outfile_root * "_" * string(time) * ".hdf5"),
        comms_ctx,
    )

    field = ClimaCore.InputOutput.read_field(hdfreader, varname)
    Base.close(hdfreader)
    return field
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
    t = Dates.datetime2unix.(time)
    hdfwriter = ClimaCore.InputOutput.HDF5Writer(
        joinpath(REGRID_DIR, hd_outfile_root * "_" * string(time) * ".hdf5"),
        comms_ctx,
    )

    ClimaCore.InputOutput.HDF5.write_attribute(hdfwriter.file, "unix time", t) # TODO: a better way to write metadata, CMIP convention
    ClimaCore.InputOutput.write!(hdfwriter, field, string(varname))
    Base.close(hdfwriter)
end


"""
    function hdwrite_regridfile_rll_to_cgll(
        FT,
        REGRID_DIR,
        datafile_rll,
        varname,
        space;
        hd_outfile_root = "data_cgll",
        mono = false,
    )
Reads and regrids data of the `varname` variable from an input NetCDF file and
saves it as another NetCDF file using Tempest Remap.
The input NetCDF fileneeds to be `Exodus` formatted, and can contain
time-dependent data. The output NetCDF file is then read back, the output
arrays converted into Fields and saved as HDF5 files (one per time slice).
This function should be called by the root process.
The saved regridded HDF5 output is readable by multiple MPI processes.

Code taken from ClimaCoupler.Regridder.


# Arguments
- `FT`: [DataType] Float type.
- `REGRID_DIR`: [String] directory to save output files in.
- `datafile_rll`: [String] filename of RLL dataset to be mapped to CGLL.
- `varname`: [String] the name of the variable to be remapped.
- `space`: [ClimaCore.Spaces.AbstractSpace] the space to which we are mapping.
- `hd_outfile_root`: [String] root of the output file name.
- `mono`: [Bool] flag to specify monotone remapping.
"""
function hdwrite_regridfile_rll_to_cgll(
    FT,
    REGRID_DIR,
    datafile_rll,
    varname,
    space,
    hd_outfile_root;
    mono = false,
)
    out_type = "cgll"

    outfile = hd_outfile_root * ".nc"
    outfile_root = mono ? outfile[1:(end - 3)] * "_mono" : outfile[1:(end - 3)]
    datafile_cgll = joinpath(REGRID_DIR, outfile_root * ".g")

    meshfile_rll = joinpath(REGRID_DIR, outfile_root * "_mesh_rll.g")
    meshfile_cgll = joinpath(REGRID_DIR, outfile_root * "_mesh_cgll.g")
    meshfile_overlap = joinpath(REGRID_DIR, outfile_root * "_mesh_overlap.g")
    weightfile = joinpath(REGRID_DIR, outfile_root * "_remap_weights.nc")

    topology = ClimaCore.Topologies.Topology2D(
        space.topology.mesh,
        ClimaCore.Topologies.spacefillingcurve(space.topology.mesh),
    )
    Nq =
        ClimaCore.Spaces.Quadratures.polynomial_degree(space.quadrature_style) +
        1
    space_undistributed = ClimaCore.Spaces.SpectralElementSpace2D(
        topology,
        ClimaCore.Spaces.Quadratures.GLL{Nq}(),
    )

    if isfile(datafile_cgll) == false
        isdir(REGRID_DIR) ? nothing : mkpath(REGRID_DIR)

        nlat, nlon = NCDataset(datafile_rll) do ds
            (ds.dim["lat"], ds.dim["lon"])
        end
        # write lat-lon mesh
        rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

        # write cgll mesh, overlap mesh and weight file
        write_exodus(meshfile_cgll, topology)
        overlap_mesh(meshfile_overlap, meshfile_rll, meshfile_cgll)

        # 'in_np = 1' and 'mono = true' arguments ensure mapping is conservative and monotone
        # Note: for a kwarg not followed by a value, set it to true here (i.e. pass 'mono = true' to produce '--mono')
        # Note: out_np = degrees of freedom = polynomial degree + 1
        kwargs = (; out_type = out_type, out_np = Nq)
        kwargs = mono ? (; (kwargs)..., in_np = 1, mono = mono) : kwargs
        remap_weights(
            weightfile,
            meshfile_rll,
            meshfile_cgll,
            meshfile_overlap;
            kwargs...,
        )
        apply_remap(datafile_cgll, datafile_rll, weightfile, [varname])
    else
        @warn "Using the existing $datafile_cgll : check topology is consistent"
    end

    function get_time(ds)
        if "time" in ds
            data_dates =
                Dates.DateTime.(
                    reinterpret.(
                        Ref(NCDatasets.DateTimeStandard),
                        ds["time"][:],
                    )
                )
        elseif "date" in ds
            data_dates = strdate_to_datetime.(string.(ds["date"][:]))
        else
            @warn "No dates available in file $datafile_rll"
            data_dates = [Dates.DateTime(0)]
        end
    end

    # read the remapped file with sparse matrices
    offline_outvector, times = NCDataset(datafile_cgll, "r") do ds_wt
        (
            offline_outvector = ds_wt[varname][:][:, :], # ncol, times
            times = get_time(ds_wt),
        )
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

    offline_fields = ntuple(x -> similar(offline_field), length(times))

    ntuple(
        x -> reshape_cgll_sparse_to_field!(
            offline_fields[x],
            offline_outvector[:, x],
            R,
        ),
        length(times),
    )

    # TODO: extend write! to handle time-dependent fields
    map(
        x -> write_to_hdf5(
            REGRID_DIR,
            hd_outfile_root,
            times[x],
            offline_fields[x],
            varname,
        ),
        1:length(times),
    )
    jldsave(
        joinpath(REGRID_DIR, hd_outfile_root * "_times.jld2");
        times = times,
    )
end


"""

This function is slightly modified from ClimaCoupler "land_sea_mask" function.
"""
function regrid_netcdf_to_field(
    FT,
    REGRID_DIR,
    comms_ctx::ClimaComms.AbstractCommsContext,
    infile,
    varname,
    boundary_space;
    outfile_root = string(varname, "_cgll"),
    mono = true,
)

    if ClimaComms.iamroot(comms_ctx)
        hdwrite_regridfile_rll_to_cgll(
            FT,
            REGRID_DIR,
            infile,
            varname,
            boundary_space,
            outfile_root;
            mono = mono,
        )
    end
    ClimaComms.barrier(comms_ctx)
    file_dates =
        load(joinpath(REGRID_DIR, outfile_root * "_times.jld2"), "times")
    field = read_from_hdf5(
        REGRID_DIR,
        outfile_root,
        file_dates[1],
        varname,
        comms_ctx,
    )
    field = swap_space!(field, boundary_space) # why do we need this?
    return nans_to_zero.(field)
end


function swap_space!(field, new_space)
    field_out = zeros(new_space)
    parent(field_out) .= parent(field)
    return field_out
end

"""
    MapInfo

A struct holding the information required to identify the
netcdf file where a global map of a particular parameter is stored,
and for carrying out regridding of that dataset to a ClimaCore.Domains.AbstractDomain.

$(DocStringExtensions.FIELDS)
"""
struct MapInfo
    "Local path to NetCDF file"
    path::String
    "Variable name of interest"
    varname::String
    "Path where temporary regrid files are stored."
    regrid_dirpath::String
    "Communication context. Not tested yet with MPI"
    comms::ClimaComms.AbstractCommsContext
end

end
