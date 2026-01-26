"""
    FileReaders

The `FileReaders` module implements backends to read and process input files.

Given that reading from disk can be an expensive operation, this module provides a pathway
to optimize the performance (if needed).

The FileReaders module contains a global cache of all the NCDatasets that are currently open.
This allows multiple NCFileReader to share the underlying file without overhead.
"""
module FileReaders

abstract type AbstractFileReader end

function NCFileReader end

function read end

function read! end

function available_dates end

function close_all_ncfiles end

extension_fns = [
    :NCDatasets => [
        :NCFileReader,
        :read,
        :read!,
        :available_dates,
        :close_all_ncfiles,
        :close,
    ],
]

"""
    is_pkg_loaded(pkg::Symbol)

Check if `pkg` is loaded or not.
"""
function is_pkg_loaded(pkg::Symbol)
    return any(k -> Symbol(k.name) == pkg, keys(Base.loaded_modules))
end

function __init__()
    # Register error hint if a package is not loaded
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(
            MethodError,
        ) do io, exc, _argtypes, _kwargs
            for (pkg, fns) in extension_fns
                if Symbol(exc.f) in fns && !is_pkg_loaded(pkg)
                    print(io, "\nImport $pkg to enable `$(exc.f)`.";)
                end
            end
        end
    end
end

end
