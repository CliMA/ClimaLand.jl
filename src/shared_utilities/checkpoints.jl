import ClimaCore: Fields, InputOutput
import ClimaUtilities

"""
    ClimaLand.find_restart(output_dir)

Find the most recent restart file in the specified output directory.

This function utilizes `ClimaUtilities.OutputPathGenerator.detect_restart_file`
to locate the latest restart file within the output directory structure,
assuming the `ActiveLinkStyle` is used for managing output folders.

# Arguments
- `output_dir`: The base output directory where the simulation results are stored.

# Returns
- The path to the most recent restart file found, or `nothing` if no restart
  file is found.
"""
function find_restart(output_dir)
    return ClimaUtilities.OutputPathGenerator.detect_restart_file(
        output_dir;
        style = ClimaUtilities.OutputPathGenerator.ActiveLinkStyle(),
    )
end

"""
    set_initial_conditions_from_checkpoint!(Y, restart_file; model::AbstractModel = nothing)

Read restart file in `restart_file` and write its content in `Y`.

The optional argument `model` is used to verify that the checkpoint contains the
same model (as verified by the hash).

See also [`ClimaLand.read_checkpoint`](@ref).
"""
function set_initial_conditions_from_checkpoint!(
    Y,
    restart_file;
    model = nothing,
)
    Y_restart, _ = read_checkpoint(restart_file; model)
    Y .= Y_restart
    return nothing
end

"""
    initial_time_from_checkpoint(restart_file; model::AbstractModel = nothing)

Read and return the restart time in `restart_file`.

The optional argument `model` is used to verify that the checkpoint contains the
same model (as verified by the hash).
"""
function initial_time_from_checkpoint(restart_file; model = nothing)
    _, t_restart = read_checkpoint(restart_file; model)
    return t_restart
end

"""
    _context_from_Y(Y)

Try extracting the context from the FieldVector Y.

Typically Y has a structure like:
```
Y
 .bucket
        .T
        .W
        .Ws
```

`_context_from_Y` tries to obtain the context from a Field in the hierarchy.
"""
function _context_from_Y(Y::Fields.FieldVector)
    for p in propertynames(Y)
        maybe_return = _context_from_Y(getproperty(Y, p))
        isnothing(maybe_return) || return maybe_return
    end
end

function _context_from_Y(Y::Fields.Field)
    return ClimaComms.context(Y)
end

function _context_from_Y(Y)
    return nothing
end

"""
    ClimaLand.save_checkpoint(Y, t, output_dir; model = nothing, context = ClimaComms.context(Y))

Save a simulation checkpoint to an HDF5 file.

This function saves the current state of the simulation, including the state
vector `Y` and the current simulation time `t`, to an HDF5 file within the
specified output directory.

# Arguments
- `Y`: The state of the simulation.
- `t`: The current simulation time.
- `output_dir`: The directory where the checkpoint file will be saved.
- `model` (Optional): The ClimaLand model object. If provided the hash of the model
  will be stored in the checkpoint file. Defaults to `nothing`. This is used
  to check for consistency.
- `context` (Optional): The ClimaComms context. This is used for distributed I/O
  operations. Defaults to the context extracted from the state vector `Y` or the `model`.
"""
function save_checkpoint(
    Y,
    t,
    output_dir;
    model = nothing,
    context = isnothing(model) ? _context_from_Y(Y) : ClimaComms.context(model),
)
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    output_file = joinpath(output_dir, "day$day.$sec.hdf5")
    hdfwriter = InputOutput.HDF5Writer(output_file, context)
    # If model was passed, add its hash, otherwise add nothing
    hash_model = isnothing(model) ? "nothing" : hash(model)
    InputOutput.write_attributes!(
        hdfwriter,
        "/",
        Dict("time" => t, "land_model_hash" => hash_model),
    )
    InputOutput.write!(hdfwriter, Y, "Y")
    Base.close(hdfwriter)
    return nothing
end

"""
    ClimaLand.read_checkpoint(file_path; model = nothing, context = ClimaComms.context())

Read a simulation checkpoint from an HDF5 file.

This function loads the simulation state from a previously saved checkpoint file.

# Arguments
- `file_path`: The path to the HDF5 checkpoint file.
- `model` (Optional): The ClimaLand model object. If provided the hash of the model
  stored in the checkpoint file will be compared with the hash of the provided
  model and a warning will be issued if they don't match. Defaults to `nothing`.
- `context` (Optional): The ClimaComms context. This is used for parallel I/O
  operations. Defaults to the default ClimaComms context.

# Returns
- `Y`: The state vector loaded from the checkpoint file.
- `t`: The simulation time loaded from the checkpoint file.
"""
function read_checkpoint(
    file_path;
    model = nothing,
    context = isnothing(model) ? ClimaComms.context() :
              ClimaComms.context(model),
)
    hdfreader = InputOutput.HDF5Reader(file_path, context)
    Y = InputOutput.read_field(hdfreader, "Y")
    attributes = InputOutput.read_attributes(hdfreader, "/")
    if !isnothing(model)
        if hash(model) != attributes["land_model_hash"]
            @warn "Restart file $(file_path) was constructed with a different land model"
        end
    end
    t = attributes["time"]
    Base.close(hdfreader)
    return Y, t
end
