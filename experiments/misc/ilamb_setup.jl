import Dates

"""
    setup_run_for_ilamb(
        ilamb_diagnostics_dir::String,
        ilamb_root::String,
        commit_id::String,
    )

Set up the ILAMB diagnostics produced from the run for ILAMB.

This involves making the `MODELS` directory in `ilamb_root` if it does not
exist, creating a model directory in the `MODELS` directory, and creating
symbolic links in the model directory.
"""
function setup_run_for_ilamb(
    ilamb_diagnostics_dir::String,
    ilamb_root::String,
    commit_id::String,
)
    # Verify directories exist
    isdir(ilamb_diagnostics_dir) ||
        error("$ilamb_diagnostics_dir is not a directory")
    isdir(ilamb_root) || error("$ilamb_root is not a directory")

    # The root ILAMB directory should contain MODELS (simulational data) and
    # DATA (observational data)
    # We can check if this is the correct directory by checking for the presence
    # of the MODELS directory. We expect the MODELS and DATA directories to be
    # at the same level
    "DATA" in readdir(ilamb_root) || error(
        "The DATA directory cannot be found in $ilamb_root. You cannot run ILAMB from this directory",
    )

    # Make the MODELS directory
    MODELS_dir = joinpath(ilamb_root, "MODELS")
    if !isdir(MODELS_dir)
        # This should not be called that often, since the MODELS directory should
        # persist from run to run
        mkdir(MODELS_dir)
    end

    # Make the directory with an appropriate model name
    model_name = create_model_name(commit_id)
    model_dir = joinpath(MODELS_dir, model_name)
    isdir(model_dir) || mkdir(model_dir)

    create_symlinks(ilamb_diagnostics_dir, model_dir)
    validate_symlinks(MODELS_dir)
end

"""
    create_model_name(commit_id)

Create a model name from the date and commit id.
"""
function create_model_name(commit_id::String)
    return "$(Dates.Date(now()))_$commit_id"
end

"""
    create_symlinks(ilamb_diagnostics_dir::String, model_dir::String)

Create symbolic links in `ilamb_data_dir` that point to the diagnostics produced
from the long runs in `ilamb_diagnostics_dir`.

Symbolic links are created rather than copied to avoid unnecessary data
duplication.
"""
function create_symlinks(ilamb_diagnostics_dir::String, model_dir::String)
    # Make the directory with the variable name if it doesn't exists
    # Then, check if a link exists or not
    # If it does, then skip it and return a warning
    # If it doesn't, make a new one (remember symlink(source, target))

    nc_filepaths = String[]
    # The directory should be flat, but we iterate recursively in case this
    # change for whatever reason
    for (root, _, files) in walkdir(ilamb_diagnostics_dir)
        for file in files
            if endswith(file, ".nc")
                push!(nc_filepaths, joinpath(root, file))
            end
        end
    end

    # The model_dir should contain symlinks to the NetCDF files in the
    # diagnostics directory for ILAMB
    for nc_filepath in nc_filepaths
        filename = basename(nc_filepath)
        target_path = joinpath(model_dir, filename)

        if islink(target_path)
            @warn "Symlink already exists at $target_path, removing and recreating"
            rm(target_path)
        elseif isfile(target_path)
            error("File (not symlink) already exists at $target_path")
        end

        # Create the symlink
        symlink(nc_filepath, target_path)
    end

    return nothing
end

"""
    validate_symlinks(model_dir)

Validate the symbolic links in `MODELS` directory where the ILAMB root is.

The function `create_symlinks` only create symbolic links for the current simulation. The
`MODELS` directory may contain the results of other longruns whose symbolic
links may or may not be valid. For example, if the diagnostics are deleted, then
the symbolic links will be invalid.

Note that this does not remove the directory if all the symlinks associated with
a particular run is deleted.
"""
function validate_symlinks(model_dir)
    model_dirs = String[]
    for maybe_dir in readdir(model_dir, join = true)
        isdir(maybe_dir) && push!(model_dirs, maybe_dir)
    end

    for single_model_dir in model_dirs
        for path in readdir(single_model_dir, join = true)
            if islink(path)
                # Check if the symlink target exists
                if !isfile(path) && !isdir(path)
                    @info "Invalid symlink found: $path -> $(readlink(path))"
                    rm(path)
                end
            end
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2
        error(
            "Usage: julia --project=.buildkite ilamb_create_symlinks.jl <ilamb_diagnostics_dir> <ilamb_root_dir> <commit_id>",
        )
    end
    ilamb_diagnostics_dir = ARGS[1]
    ilamb_root_dir = ARGS[2]
    # On Github, commit ids are shortened to 7 characters, so this script does
    # the same. Also, this should be gotten from the BUILDKITE_COMMIT
    # environment variable
    commit_id = first(ARGS[3], 7)
    setup_run_for_ilamb(ilamb_diagnostics_dir, ilamb_root_dir, commit_id)
end
