"""
    Regridders

The `Regridders` module implement structs and functions to remap datasets to simulation
grids.

Currently, the schemes implemented are `TempestRegridder`, which uses
`ClimaCoreTempestRemap`, and `InterpolationsRegridder`, which uses `Interpolations.jl`.

The key function exposed by `Regridders` is the `regrid` method.
"""
module Regridders

import ..ClimaUtilities

# When adding a new regridder, you also have to change some functions in the DataHandler
# module. Find where :TempestRegridder is used.
abstract type AbstractRegridder end

function TempestRegridder end

function InterpolationsRegridder end

function regrid end

function regrid! end

"""
    default_regridder_type()

Return the type of regridder to be used if the user doesn't specify one.
This function returns the first available regridder in the following order:
  - InterpolationsRegridder
  - TempestRegridder
based on which regridder(s) are currently loaded.
"""
function default_regridder_type()
    # Use InterpolationsRegridder if available
    if !isnothing(
        Base.get_extension(
            ClimaUtilities,
            :ClimaUtilitiesClimaCoreInterpolationsExt,
        ),
    )
        regridder_type = :InterpolationsRegridder
        # If InterpolationsRegridder isn't available, and TempestRegridder is, use TempestRegridder
    elseif !isnothing(
        Base.get_extension(
            ClimaUtilities,
            :ClimaUtilitiesClimaCoreTempestRemapExt,
        ),
    )
        regridder_type = :TempestRegridder
    else
        error("No regridder available")
    end
    return regridder_type
end

extension_fns = [
    :ClimaCoreTempestRemap => [:TempestRegridder, :regrid],
    :ClimaCore => [:InterpolationsRegridder, :regrid],
    :Interpolations => [:InterpolationsRegridder, :regrid],
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
