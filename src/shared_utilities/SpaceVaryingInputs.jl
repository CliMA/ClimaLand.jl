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
using ClimaCore
using ClimaComms
using DocStringExtensions
import ..searchsortednearest
import ..linear_interpolation
using ClimaLand.FileReader
export SpaceVaryingInput

# Analytic case
"""
    SpaceVaryingInput(data_function::Function, space::ClimaCore.Spaces.AbstractSpace) 
  
Returns the parameter field to be used in the model; appropriate when
a parameter is defined using a function of the coordinates of the space.

Pass the ``data" as a function `data_function` which takes coordinates as arguments, 
and  the ClimaCore space of the model simulation.

This returns a scalar field.
Note that data_function is broadcasted over the coordinate field. Internally, inside
your function, this must be unpacked (coords.lat, coords.lon, e.g.) for 
use of the coordinate values directly.
"""
function SpaceVaryingInput(
    data_function::Function,
    space::ClimaCore.Spaces.AbstractSpace,
)
    model_value = ClimaCore.Fields.zeros(space)
    coords = ClimaCore.Fields.coordinate_field(space)
    return model_value .= data_function.(coords)
end

# 1-D Case
"""
    function SpaceVaryingInput(
        data_z::AbstractArray,
        data_values::AbstractArray,
        space::S,
    ) where {S <: ClimaCore.Spaces.CenterFiniteDifferenceSpace}

Given a set of depths `data_z` and the observed values `data_values` 
at those depths, create an interpolated field of values at each value
of z in the model grid - defined implicitly by `space`.

Returns a ClimaCore.Fields.Field of scalars.
"""
function SpaceVaryingInput(
    data_z::AbstractArray,
    data_values::AbstractArray,
    space::S,
) where {S <: ClimaCore.Spaces.CenterFiniteDifferenceSpace}
    model_value = ClimaCore.Fields.zeros(space)
    # convert the passed arrays to the appropriate type for the device
    device = ClimaComms.device(space)
    AT = ClimaComms.array_type(device)
    data_values = AT(data_values)
    data_z = AT(data_z)
    zvalues = ClimaCore.Fields.coordinate_field(space).z
    #now create the parameter field
    model_value .= map(zvalues) do z
        linear_interpolation(data_z, data_values, z)
    end
    return model_value
end


"""
    SpaceVaryingInputs.SpaceVaryingInput(
        data_z::AbstractArray,
        data_values::NamedTuple,
        space::S,
        dest_type::Type{DT},
    ) where {
        S <: ClimaCore.Spaces.CenterFiniteDifferenceSpace,
        DT,
    }

Returns a field of parameter structs to be used in the model; 
appropriate when the parameter struct values vary in depth; 
the `dest_type` argument is the struct type - we assumed that
your struct as a constructor which accepts the values of its arguments
by kwarg,
- `data_z` is where the measured values were obtained,
- `data_values` is a NamedTuple with keys equal to the argument names
of the struct, and with values equal to an array of measured values,
- `space` defines the model grid.

As an example, we can create a field of vanGenuchten structs as follows. This struct
requires two parameters, `α` and `n`. Let's assume that we have measurements of these
as a function of depth at the locations given by `data_z`, called `data_α` and `data_n`.
Then we can write
`vG_field = SpaceVaryingInput(data_z, (;α = data_α, n = data_n), space, vanGenuchten{Float32})`.
Under the hood, at each point in the model grid, we will create
`vanGenuchten{Float32}(;α = interp_α, n = interp_n)`, where `interp` indicates
the interpolated value at the model depth.

Returns a ClimaCore.Fields.Field of type DT.

"""
function SpaceVaryingInput(
    data_z::AbstractArray,
    data_values::NamedTuple,
    space::S,
    dest_type::Type{DT},
) where {S <: ClimaCore.Spaces.CenterFiniteDifferenceSpace, DT}
    zvalues = ClimaCore.Fields.coordinate_field(space).z
    # convert the passed arrays to the appropriate type for the device
    device = ClimaComms.device(space)
    AT = ClimaComms.array_type(device)
    data_z = AT(data_z)
    data_values_AT = map(AT, data_values)
    # now create the field of structs
    model_value = map(zvalues) do z
        args = map(data_values_AT) do value
            return linear_interpolation(data_z, value, z)
        end
        DT(; args...)
    end
    return model_value
end

"""
    SpaceVaryingInput(data::PDS, varname::N; space::S) 
    where {PDS <: FileReader.PrescribedDataStatic, N <: String, S <: ClimaCore.Spaces.SpectralElement2D}
  
Returns the parameter field to be used in the model; appropriate when
a parameter is defined on the surface of the Earth.

Pass the data as a `FileReader.PrescribedDataStatic` object
as well as the variable name of the variable in the data file, and the ClimaCore space
of the model simulation.

Returns a ClimaCore.Fields.Field of scalars; analogous to the 1D case which also
returns a ClimaCore.Fields.Field of scalars.
"""
function SpaceVaryingInput(
    data::PDS,
    varname::N,
    space::S,
) where {
    PDS <: FileReader.PrescribedDataStatic,
    N <: String,
    S <: ClimaCore.Spaces.SpectralElementSpace2D,
}
    return FileReader.get_data_at_date(data, space, varname)
end

end
