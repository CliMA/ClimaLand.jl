# SpaceVaryingInputs.jl
#
# This module contains methods to process external data, regrid it onto the
# model grid, and return the corresponding fields for use in the simulation.
# This module only concerns with external data which varies in space,
# and not time. For temporally varying input, we refer you to `TimeVaryingInputs.jl`.

# All spatially varying parameter fields are assumed to fit into memory,
# and on GPU runs, they have underlying CuArrays on the GPU.

# The planned parameter underlying arrays are:
# - one-dimensional (values prescribed as a function of depth at a site),
# - two-dimensional (values prescribed globally at each lat/lon),
# - three-dimensional (values prescribed as a function of depth globally)
# - analytic (functions of the coordinates of the space)

module SpaceVaryingInputs

function SpaceVaryingInput end

extension_fns =
    [:ClimaCore => [:SpaceVaryingInput], :NCDatasets => [:SpaceVaryingInput]]

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
                if Symbol(exc.f) == :SpaceVaryingInput
                    print(
                        io,
                        "\nYou might also need a regridder to use `$(exc.f)`.";
                    )
                end
            end
        end
    end
end


end
