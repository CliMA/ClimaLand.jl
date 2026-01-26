"""
The `OutputPathGenerator` module provides tools to prepare the directory structure to be
used as output for a simulation. This might entail creating folders, moving existing data,
et cetera.
"""
module OutputPathGenerator

import ClimaComms
import ClimaComms: barrier, iamroot, bcast
import ..Utils: sort_by_creation_time

import Base: rm

# Note, Styles have nothing to do with traits
abstract type OutputPathGeneratorStyle end

"""
    maybe_wait_filesystem(context,
                          path,
                          check_func = ispath,
                          sleep_time = 0.1,
                          max_attempts = 10)


Distributed filesystems might need some time to catch up a file/folder is created/removed.

This function watches the given `path` with `check_func` and returns when `check_func(path)`
returns true. This is done by trying up to `max_attempt` times and sleeping `sleep_time`
seconds in between. `sleep_time` is increased by 50 % after each attempt.

Example: when creating a file, we want to check that all the MPI processes see that new
file. In this case, `check_func` could be `ispath`. Another example is with removing files
and `check_func` would be `(f) -> !ispath(f)`.
"""
function maybe_wait_filesystem(
    context,
    path;
    check_func = ispath,
    sleep_time = 0.1,
    max_attempts = 10,
)
    barrier(context)
    attempt = 1
    while attempt < max_attempts
        check_func(path) && return nothing
        sleep(sleep_time)
        attempt = attempt + 1
        sleep_time = 1.5sleep_time
    end
    error(
        "Path $path not properly synced. On distributed systems, this is typically due to the slow response of the filesystem." *
        " We waited for a few seconds and tried multiple times, but different MPI processes still disagree on the state of the filesystem. Aborting.",
    )
    return nothing
end

"""
    RemovePreexistingStyle

With this option, the output directory is directly specified. If the directory already
exists, remove it. No confirmation is asked, so use at your own risk.
"""
struct RemovePreexistingStyle <: OutputPathGeneratorStyle end

"""
    ActiveLinkStyle

This style generates a unique output path within a base directory specified by
`output_path`. It ensures the base directory exists and creates it if necessary.
Additionally, it manages a sequence of subfolders and a symbolic link named "output_active"
for convenient access to the active output location.

This style is designed to:
- be non-destructive,
- provide a deterministic and fixed path for the latest available data,
- and have nearly zero runtime overhead.

`generate_output_path` returns path to the newly created folder with the next available
increment (of the form `output_1234`), and ensures that a valid `output_active` link points
to that folder.

# Examples:

Let us assume that `output_path = dormouse`.

- `dormouse` does not exist in the current working directory: `ActiveLinkStyle` will create
  it and return `dormouse/output_0000`. In the process, a symlink `dormouse/output_active`
  is also created. This symlink points to `dormouse/output_0000`.
- `dormouse` exists and contains a `output_active` link that points to
  `dormouse/output_0005`, `ActiveLinkStyle` will a create new directory
  `dormouse/output_0006`, return this path, and change the `output_active` to point to this
  directory.
- `dormouse` exists and does not contain a `output_active`, `ActiveLinkStyle` will check if
  any `dormouse/output_XXXX` exists. If not, it creates `dormouse/output_0000` and a link
  `dormouse/output_active` that points to this directory.

## A note for Windows users

Windows does not always allow the creation of symbolic links by unprivileged users. This
depends on the version of Windows, but also some of its settings. When the creation of
symbolic links is not possible, `OutputPathGenerator` will create NTFS junction points
instead. Junction points are similar to symbolic links, with the main difference that they
have to refer to directories and they have to be absolute paths. As a result, on systems
that do not allow unprivileged users to create symbolic links, moving the base output folder
results in breaking the `output_active` link.
"""
struct ActiveLinkStyle <: OutputPathGeneratorStyle end

"""
    generate_output_path(output_path,
                         context = nothing,
                         style::OutputPathGeneratorStyle = ActiveLinkStyle())

Process the `output_path` and return a string with the path where to write the output.

The `context` is a `ClimaComms` context and is required for MPI runs.

How the output should be structured (in terms of directory tree) is determined by the
`style`.

# Styles

- `RemovePreexistingStyle`: the `output_path` provided is the actual output path. If a directory
  already exists there, remove it without asking for confirmation.

- `ActiveLinkStyle`: the `output_path` returned is a new folder of the form
  `output_path/output_1234`, where the number is incremented every time this function is
  called. `ActiveLinkStyle` also creates a link `output_path/output_active` that ensures
  that the most recent output is always accessible at the `output_path/output_active` path.
  This is style is non-destructive.

(Note, "styles" have nothing to do with Julia traits.)
"""
function generate_output_path(
    output_path;
    context = ClimaComms.context(),
    style::OutputPathGeneratorStyle = ActiveLinkStyle(),
)
    output_path == "" && error("output_path cannot be empty")
    # Let's make sure we are synced before we do any filesystem operation
    ClimaComms.barrier(context)
    return generate_output_path(style, output_path; context)
end

"""
    generate_output_path(::RemovePreexistingStyle, output_path, context = nothing)

Documentation for this function is in the `RemovePreexistingStyle` struct.
"""
function generate_output_path(
    ::RemovePreexistingStyle,
    output_path;
    context = ClimaComms.context(),
)
    if iamroot(context)
        if isdir(output_path)
            @warn "Removing $output_path"
            rm(output_path, recursive = true)
        end
        mkpath(output_path)
    end
    maybe_wait_filesystem(context, output_path)
    return output_path
end


"""
    generate_output_path(::ActiveLinkStyle, output_path, context = nothing)

Documentation for this function is in the `ActiveLinkStyle` struct.
"""
function generate_output_path(::ActiveLinkStyle, output_path; context = nothing)
    # Only root needs to do file-system operations, we will broadcast the result
    # to the other MPI processes to the new_output_folder variable after
    # everything is done.


    # We have to be a little careful with errors: the block below is executed
    # only by root, but would like to throw errors. We cannot throw errors in
    # such a block because it is only executed by root, and the other MPI
    # processes would hang waiting for the barrier. So, we return an
    # error_message and broadcast it to the other processes

    # Filled below this block with a ClimaComms.bcast
    new_output_folder = ""
    error_message = ""

    if iamroot(context)
        # Ensure path ends with a trailing slash for consistency
        path_separator_str = Base.Filesystem.path_separator
        # We need the path_separator as a char to use it in rstrip
        path_separator_char = path_separator_str[1]
        output_path =
            rstrip(output_path, path_separator_char) * path_separator_char

        # Create folder if it does not exist
        isdir(output_path) || mkpath(output_path)

        name_rx = r"output_(\d\d\d\d)"

        # Look for a output_active link
        active_link = joinpath(output_path, "output_active")

        link_exists = islink(active_link)

        if link_exists
            target = readlink(active_link)
            counter_str = match(name_rx, target)
            if !isnothing(counter_str)
                # counter_str is the only capturing group
                counter = parse(Int, counter_str[1])
                next_counter = counter + 1

                # Remove old link
                rm(active_link)
            else
                error_message = "Link $target points to a folder with a name we do not handle"
            end
        else
            # The link does not exist, but maybe there are already output folders. We can try to
            # guess what was the last one by first filtering the folders that match the name,
            # and then taking the first one when sorted in reverse alphabetical order.
            existing_outputs =
                filter(x -> !isnothing(match(name_rx, x)), readdir(output_path))
            if length(existing_outputs) > 0
                @warn "$output_path already contains some output data, but no active link"
                latest_output = first(sort(existing_outputs, rev = true))
                counter_str = match(name_rx, latest_output)
                counter = parse(Int, counter_str[1])
                next_counter = counter + 1
                @warn "Restarting counter from $next_counter"
            else
                # This is our first output folder, initialize counter
                next_counter = 0
            end
        end

        if isempty(error_message)
            # Ensure that there are four digits
            next_counter_str = lpad(next_counter, 4, "0")

            # Create the next folder
            new_output_folder =
                joinpath(output_path, "output_$next_counter_str")

            # Make new folder
            mkpath(new_output_folder)

            # On Windows, creating symlinks might require admin privileges. This depends on the
            # version of Windows and some of its settings (e.g., if "Developer Mode" is enabled).
            # So, we first try creating a symlink. If this fails, we resort to creating a NTFS
            # junction. This is almost the same as a symlink, except it requires absolute paths.
            # In general, relative paths would be preferable because they make the output
            # completely relocatable (whereas absolute paths are not).
            try
                symlink(
                    "output_$next_counter_str",
                    active_link,
                    dir_target = true,
                )
            catch e
                if e isa Base.IOError && Base.uverrorname(e.code) == "EPERM"
                    active_link = abspath(active_link)
                    dest_active_link = abspath("output_$next_counter_str")
                    symlink(dest_active_link, active_link, dir_target = true)
                else
                    error_message = string(e)
                end
            end
        end
    end
    error_message = bcast(context, error_message)
    isempty(error_message) || error(error_message)
    new_output_folder = bcast(context, new_output_folder)
    maybe_wait_filesystem(context, new_output_folder)
    return new_output_folder
end


"""
    detect_restart_file(base_output_dir;
                        restart_file_rx = r"day\\d+\\.\\w+\\.hdf5",
                        sort_func = sort_by_creation_time,
                        style = ActiveLinkStyle()
                        )

Detects and returns the path to the most recent restart file within the directory structure
specified by `base_output_dir`.


Returns `nothing` if no suitable restart file is found.

This function searches for restart files within the directory structure organized according
to the provided `ActiveLinkStyle`. It identifies potential output directories based on the
style and then looks for files matching the `restart_file_rx` regular expression within
these directories.

By default, the function assumes restart files have names like "dayDDDD.SSSSS.hdf5", where
DDDD represents the day number and SSSSS represents the number of seconds.

If multiple restart files are found, the function uses the `sort_func` to determine the most
recent one. The default sorting function, `sort_by_creation_time`, sorts files based on
their creation timestamps, returning the file with the latest creation time. Users can
provide custom sorting functions to prioritize files based on other criteria, such as the
simulation time stored within the HDF5 file.

**Return Value:**

- If a suitable restart file is found, the function returns its full path as a string.
- If no output directory matching the `ActiveLinkStyle` or no restart file
  matching the `restart_file_rx` is found, the function returns `nothing`. This
  indicates that automatic restart is not possible.

**Arguments:**

- `output_dir_style`: An `ActiveLinkStyle` object defining the structure of output
  directories.
- `base_output_dir`: The base directory where the output directory structure is located.
- `restart_file_rx`: A regular expression used to identify restart files within output
  directories. Defaults to `r"day\\d+\\.\\w+\\.hdf5"`.
- `sort_func`: A function used to sort restart files and select the most recent one.
  Defaults to `sort_by_creation_time`.
"""
function detect_restart_file(
    base_output_dir;
    restart_file_rx = r"day\d+\.\w+\.hdf5",
    sort_func = sort_by_creation_time,
    style = ActiveLinkStyle(),
)
    style isa ActiveLinkStyle ||
        error("detect_restart_file works only with ActiveLinkStyle")

    # if base_output_dir does not exist, we return restart_file = nothing because there is
    # no restart file to be detected
    isdir(base_output_dir) || return nothing

    # output_dir will be something like ABC/DEF/output_1234
    name_rx = r"output_(\d\d\d\d)"
    restart_file = nothing

    existing_outputs =
        filter(x -> !isnothing(match(name_rx, x)), readdir(base_output_dir))

    isempty(existing_outputs) && return nothing

    # Walk directories backwards looking for restart files
    for output_folder in sort(existing_outputs, rev = true)
        previous_folder = joinpath(base_output_dir, output_folder)
        possible_restart_files =
            filter(f -> occursin(restart_file_rx, f), readdir(previous_folder))

        if !isempty(possible_restart_files)
            restart_file_name = last(sort_func(possible_restart_files))
            restart_file = joinpath(previous_folder, restart_file_name)
            return restart_file
        end
    end
    # Nothing was found
    return nothing
end

##############
# DEPRECATED #
##############
function detect_restart_file(
    style,
    base_output_dir;
    restart_file_rx = r"day\d+\.\w+\.hdf5",
    sort_func = sort_by_creation_time,
)
    Base.depwarn(
        "`style` as first argument was moved to a keyword argument",
        :detect_restart_file,
    )
    return detect_restart_file(
        base_output_dir;
        restart_file_rx,
        sort_func,
        style,
    )
end


end
