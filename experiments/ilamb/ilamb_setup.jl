import Dates

const NUM_MODELS_TO_KEEP = 3

"""
    setup_run_for_ilamb(
        ilamb_diagnostics_dir::String,
        ilamb_models_parent_dir::String,
        ilamb_build_dir::String,
        buildkite_id::String,
        commit_id::String,
        num_models_to_keep::Integer,
    )

Set up the ILAMB diagnostics produced from the run for ILAMB.

This involves making the `MODELS` directory in `ilamb_models_parent_dir` if it
does not exist, creating a model directory in the `MODELS` directory for this
run, and creating symbolic links in the model directory pointing to NetCDF files
in `ilamb_diagnostics_dir`.

If there are more models present than `num_models_to_keep`, models are
automatically removed starting with the oldest model.
"""
function setup_run_for_ilamb(
    ilamb_diagnostics_dir::String,
    ilamb_models_parent_dir::String,
    ilamb_build_dir::String,
    buildkite_id::String,
    commit_id::String,
    num_models_to_keep::Integer,
)
    # Verify directories exist
    isdir(ilamb_diagnostics_dir) ||
        error("$ilamb_diagnostics_dir is not a directory")
    isdir(ilamb_models_parent_dir) || error("$ilamb_root is not a directory")

    # Make the MODELS directory
    MODELS_dir = joinpath(ilamb_models_parent_dir, "MODELS")
    # This should not be called that often, because the MODELS directory
    # should persist from run to run
    isdir(MODELS_dir) || mkdir(MODELS_dir)

    # Make the directory with an appropriate model name
    model_name = create_model_name(buildkite_id, commit_id)
    model_dir = joinpath(MODELS_dir, model_name)
    isdir(model_dir) || mkdir(model_dir)

    create_symlinks(ilamb_diagnostics_dir, model_dir)
    validate_symlinks(MODELS_dir)

    delete_old_runs(MODELS_dir, num_models_to_keep)
    clean_build_dir(MODELS_dir, ilamb_build_dir)
end

"""
    create_model_name(buildkite_id, commit_id)

Create a model name from the current date, buildkite_id, and commit id.
"""
function create_model_name(buildkite_id::String, commit_id::String)
    # On GitHub, commit IDs are shortened to 7 characters, so this script does
    # the same. Also, `commit_id` should be passed by the BUILDKITE_COMMIT
    # environment variable. Similarly, `buildkite_id` should be passed by
    # the BUILDKITE_BUILD_NUMBER environment variable.
    short_commit_id = first(commit_id, 7)
    return "$(buildkite_id)_$(Dates.Date(Dates.now()))_$short_commit_id"
end

"""
    create_symlinks(ilamb_diagnostics_dir::String, model_dir::String)

Create symbolic links in `model_dir` that point to the diagnostics produced
from the long runs in `ilamb_diagnostics_dir`.

Creating symbolic links is preferred over copy files to avoid unnecessary data
duplication.
"""
function create_symlinks(ilamb_diagnostics_dir::String, model_dir::String)
    nc_filepaths = String[]
    # The directory should be flat, but we iterate recursively in case this
    # changes for whatever reason
    for (root, _, files) in walkdir(ilamb_diagnostics_dir)
        for file in files
            endswith(file, ".nc") && push!(nc_filepaths, joinpath(root, file))
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
    validate_symlinks(model_dir::String)

Validate the symbolic links in `MODELS` directory.

The function `create_symlinks` only creates symbolic links for the current
simulation. The `MODELS` directory may contain the results of other long runs
whose symbolic links may or may not be valid. For example, if the diagnostics
are deleted, then the symbolic links will be invalid.

Note that this does not remove the directory if all the symlinks associated with
a particular run is deleted.
"""
function validate_symlinks(model_dir::String)
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

"""
    delete_old_runs(MODELS_dir::String, N::Integer)

Delete model directories in `MODELS_dir` until there are only `N` model
directories remaining.

Only directories corresponding to ClimaLand simulation outputs (matching the
pattern created by `create_model_name`) can be deleted. Directories are sorted
by Buildkite ID, with lower IDs considered older and deleted first.
"""
function delete_old_runs(MODELS_dir::String, N::Integer)
    dirs = filter(isdir, readdir(MODELS_dir, join = true))
    length(dirs) <= N && return nothing
    num_to_delete = length(dirs) - N
    filter!(dir -> basename(dir) |> is_model_from_simulation, dirs)
    sort!(dirs, by = dir -> extract_buildkite_id(basename(dir)))
    dirs_to_remove = first(dirs, num_to_delete)
    @info "Removing directories: $dirs_to_remove"
    rm.(dirs_to_remove, recursive = true)
    return nothing
end

"""
    clean_build_dir(MODELS_dir::String, build_dir::String)

If the build directory exists, then clean the `_build` directory by deleting
all files that are not NetCDF files that correspond to models in `MODELS_dir`.

See this issue in the ILAMB repository for more information about the caching:
https://github.com/rubisco-sfa/ILAMB/issues/78. The ILAMB leaderboard is created
in two passes. In the first pass, NetCDF files are produced that contain the
analysis. The second pass produces the plots and figures that will eventually be
displayed on the websites. The second pass cannot be cached, because if a new
model is added, all plots need to be updated. However, the first pass can be
cached by not deleting the NetCDF files which this function accomplished.
"""
function clean_build_dir(MODELS_dir::String, build_dir::String)
    # If there isn't a build directory, then this is the first time ILAMB is
    # being run
    isdir(build_dir) || return nothing
    basename(build_dir) == "_build" ||
        error("This is not the build directory: $build_dir")
    dirs = filter(isdir, readdir(MODELS_dir, join = true))
    model_names = basename.(dirs)

    # To clean the build directory, all files that are not NetCDF files or
    # NetCDF files that do not correspond to any model directory will be removed
    for (root, _, files) in walkdir(build_dir)
        for file in files
            filepath = joinpath(root, file)
            if !endswith(filepath, ".nc") ||
               !any(model_name -> occursin(model_name, file), model_names)
                rm(filepath)
            end
        end
    end
    return nothing
end

"""
    is_model_from_simulation(model_name::String)

Return true if `model_name` matches the pattern produced by `create_model_name`.
"""
function is_model_from_simulation(model_name::String)
    # \d* - Match any number of digits (corresponds to buildkite ID)
    # \d{4}-\d{2}-\d{2} - Match the date format
    # [0-9a-f]{7} - Match the commit ID that is shortened to 7 characters
    return occursin(r"\d+_\d{4}-\d{2}-\d{2}_[0-9a-f]{7}", model_name)
end

"""
    extract_buildkite_id(model_name::String)

Extract the buildkite ID from `model_name` as an integer.
"""
function extract_buildkite_id(model_name::String)
    m = match(r"(\d+)_", model_name)
    return m === nothing ? nothing : parse(Int, m.captures[1])
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 4
        error(
            "Usage: julia --project=.buildkite ilamb_setup.jl <ilamb_diagnostics_dir> <ilamb_models_parent_dir> <buildkite_id> <commit_id>",
        )
    end
    ilamb_diagnostics_dir = ARGS[1]
    # This directory should store the MODELS directory
    ilamb_models_parent_dir = ARGS[2]
    ilamb_build_dir = joinpath(ilamb_models_parent_dir, "_build")
    # This can also be any unique identifier
    buildkite_id = ARGS[3]
    commit_id = ARGS[4]
    setup_run_for_ilamb(
        ilamb_diagnostics_dir,
        ilamb_models_parent_dir,
        ilamb_build_dir,
        buildkite_id,
        commit_id,
        NUM_MODELS_TO_KEEP,
    )
end
