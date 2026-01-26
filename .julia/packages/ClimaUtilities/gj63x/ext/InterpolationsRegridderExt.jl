module InterpolationsRegridderExt

import Interpolations as Intp

import ClimaCore
import ClimaCore.Fields: Adapt
import ClimaCore.Fields: ClimaComms

import ClimaUtilities.Regridders

struct InterpolationsRegridder{
    SPACE <: ClimaCore.Spaces.AbstractSpace,
    FIELD <: ClimaCore.Fields.Field,
    IM,
    BC,
    DT <: Tuple,
} <: Regridders.AbstractRegridder

    """ClimaCore.Space where the output Field will be defined"""
    target_space::SPACE

    """ClimaCore.Field of physical coordinates over which the data will be interpolated"""
    coordinates::FIELD

    """Method of gridded interpolation as accepted by Interpolations.jl"""
    interpolation_method::IM

    """Tuple of extrapolation conditions as accepted by Interpolations.jl"""
    extrapolation_bc::BC

    """Tuple of booleans signifying if the dimension is monotonically increasing. True for
    dimensions that are monotonically increasing, false for dimensions that are monotonically decreasing."""
    dim_increasing::DT
end

# Note, we swap Lat and Long! This is because according to the CF conventions longitude
# should be first, so files will have longitude as first dimension.
totuple(pt::ClimaCore.Geometry.LatLongZPoint) = pt.long, pt.lat, pt.z
totuple(pt::ClimaCore.Geometry.LatLongPoint) = pt.long, pt.lat
totuple(pt::ClimaCore.Geometry.XYZPoint) = pt.x, pt.y, pt.z

"""
    InterpolationsRegridder(target_space::ClimaCore.Spaces.AbstractSpace
                            [; extrapolation_bc::Tuple,
                               dim_increasing::Union{Nothing, Tuple},
                               interpolation_method = Interpolations.Linear()])

An online regridder that uses Interpolations.jl

Currently, InterpolationsRegridder is only implemented for LatLong and LatLongZ spaces. It
performs linear interpolation along each of the directions (separately). By default, it
imposes periodic boundary conditions for longitude, flat for latitude, and throwing errors
when extrapolating in z. This can be customized by passing the `extrapolation_bc` keyword
argument.

InterpolationsRegridder is GPU and MPI compatible in the simplest possible way: each MPI
process has the entire input data and everything is copied to GPU.

Keyword arguments
=================

The optional keyword argument `extrapolation_bc` controls what should be done when the
interpolation point is not in the domain of definition. This has to be a tuple of N
elements, where N is the number of spatial dimensions. For 3D spaces, the default is
`(Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Throw())`.

The optional keyword argument `dim_increasing` controls which dimensions should
be reversed before performing interpolation. This must be a tuple of N booleans, where
N is the number of spatial dimensions. The default is the `true` for each
spatial dimension.
"""
function Regridders.InterpolationsRegridder(
    target_space::ClimaCore.Spaces.AbstractSpace;
    extrapolation_bc::Union{Nothing, Tuple} = nothing,
    dim_increasing::Union{Nothing, Tuple} = nothing,
    interpolation_method = Intp.Linear(),
)
    coordinates = ClimaCore.Fields.coordinate_field(target_space)
    # set default values for the extrapolation_bc and dim_increasing if they are not provided
    if eltype(coordinates) <: ClimaCore.Geometry.LatLongPoint
        isnothing(extrapolation_bc) &&
            (extrapolation_bc = (Intp.Periodic(), Intp.Flat()))
        isnothing(dim_increasing) && (dim_increasing = (true, true))
    elseif eltype(coordinates) <: ClimaCore.Geometry.LatLongZPoint
        isnothing(extrapolation_bc) &&
            (extrapolation_bc = (Intp.Periodic(), Intp.Flat(), Intp.Throw()))
        isnothing(dim_increasing) && (dim_increasing = (true, true, true))
    elseif eltype(coordinates) <: ClimaCore.Geometry.XYZPoint
        isnothing(extrapolation_bc) &&
            (extrapolation_bc = (Intp.Flat(), Intp.Flat(), Intp.Throw()))
        isnothing(dim_increasing) && (dim_increasing = (true, true, true))
    else
        error("Only lat-long, lat-long-z, and x-y-z spaces are supported")
    end

    return InterpolationsRegridder(
        target_space,
        coordinates,
        interpolation_method,
        extrapolation_bc,
        dim_increasing,
    )
end

"""
    regrid(regridder::InterpolationsRegridder, data, dimensions)::Field

Regrid the given data as defined on the given dimensions to the `target_space` in `regridder`.

This function is allocating.
"""
function Regridders.regrid(regridder::InterpolationsRegridder, data, dimensions)
    FT = ClimaCore.Spaces.undertype(regridder.target_space)

    # For a 2D space with LatLongZ coordinates, we need to drop the z dimension to regrid 2D data.
    if eltype(regridder.coordinates) <: ClimaCore.Geometry.LatLongZPoint &&
       !(
           regridder.target_space isa
           ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace
       ) &&
       length(dimensions) == 2

        coords = map(regridder.coordinates) do coord
            ClimaCore.Geometry.LatLongPoint(coord.lat, coord.long)
        end
    else
        coords = regridder.coordinates
    end

    dimensions_FT = map(dimensions, regridder.dim_increasing) do dim, increasing
        !increasing ? reverse(FT.(dim)) : FT.(dim)
    end

    data_transformed = data
    # Reverse the data if needed. This allocates, so ideally it should be done in preprocessing
    if !all(regridder.dim_increasing)
        decreasing_indices =
            Tuple([i for (i, d) in enumerate(regridder.dim_increasing) if !d])
        data_transformed = reverse(data, dims = decreasing_indices)
    end
    # Make a linear spline
    itp = Intp.extrapolate(
        Intp.interpolate(
            dimensions_FT,
            FT.(data_transformed),
            Intp.Gridded(regridder.interpolation_method),
        ),
        regridder.extrapolation_bc,
    )

    # Move it to GPU (if needed)
    gpuitp = Adapt.adapt(ClimaComms.array_type(regridder.target_space), itp)

    return map(coords) do coord
        gpuitp(totuple(coord)...)
    end
end

"""
     Regridders.regrid!(output, regridder::InterpolationsRegridder, data::AbstractArray{FT, N}, dimensions::NTuple{N, AbstractVector{FT}})

Regrid the given data as defined on the given dimensions to the `target_space` in `regridder`, and store the result in `output`.
This method does not automatically reverse dimensions, so the dimensions must be sorted before calling this method.
"""
function Regridders.regrid!(
    output,
    regridder::InterpolationsRegridder,
    data::AbstractArray{FT, N},
    dimensions::NTuple{N, AbstractVector{FT}},
) where {FT, N}
    all(regridder.dim_increasing) || error(
        "Dimensions must be monotonically increasing to use regrid!. Sort the dimensions first, or use regrid.",
    )
    itp = Intp.extrapolate(
        Intp.interpolate(
            dimensions,
            data,
            Intp.Gridded(regridder.interpolation_method),
        ),
        regridder.extrapolation_bc,
    )
    gpuitp = Adapt.adapt(ClimaComms.array_type(regridder.target_space), itp)
    output .= splat(gpuitp).(totuple.(regridder.coordinates))
    return nothing
end

end
