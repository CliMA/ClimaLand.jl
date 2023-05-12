export boundary_flux,
    AbstractBC,
    AbstractBoundary,
    TopBoundary,
    BottomBoundary,
    diffusive_flux,
    get_Δz

"""
    AbstractBC

An abstract type for types of boundary conditions, which will include
prescribed functions of space and time as Dirichlet conditions or
Neumann conditions, in addition to other 
convenient conditions.
"""
abstract type AbstractBC end

"""
    AbstractBoundary

An abstract type to indicate which boundary we are doing calculations for.
Currently, we support the top boundary (TopBoundary)
and bottom boundary (BottomBoundary).
"""
abstract type AbstractBoundary end


"""
    TopBoundary{} <: AbstractBoundary{}

A simple object which should be passed into a function to
indicate that we are considering the top boundary of the soil.
"""
struct TopBoundary <: AbstractBoundary end

"""
    BottomBoundary{} <: AbstractBoundary{}

A simple object which should be passed into a function to
indicate that we are considering the bottom boundary of the soil.
"""
struct BottomBoundary <: AbstractBoundary end

"""
    get_Δz(z::ClimaCore.Fields.Field)

A function to return a tuple containing the distance between the top boundary
and its closest center, and the bottom boundary and its closest center, 
both as Fields.
"""
function get_Δz(z::ClimaCore.Fields.Field)
    # Extract the differences between levels of the face space
    fs = ClimaLSM.Domains.obtain_face_space(axes(z))
    z_face = ClimaCore.Fields.coordinate_field(fs).z
    Δz = ClimaCore.Fields.Δz_field(z_face)

    Δz_top = Fields.level(
        Δz,
        ClimaCore.Utilities.PlusHalf(ClimaCore.Spaces.nlevels(fs) - 1),
    )
    Δz_bottom = Fields.level(Δz, ClimaCore.Utilities.PlusHalf(0))
    return Δz_top ./ 2, Δz_bottom ./ 2
end

"""
    diffusive_flux(K, x_2, x_1, Δz)

Calculates the diffusive flux of a quantity x (water content, temp, etc).
Here, x_2 = x(z + Δz) and x_1 = x(z), so x_2 is at a larger z by convention.
"""
function diffusive_flux(K, x_2, x_1, Δz)
    return @. -K * (x_2 - x_1) / Δz
end

"""
    boundary_flux(bc::AbstractBC, bound_type::AbstractBoundary, Δz, _...)::ClimaCore.Fields.Field
A function which returns the correct boundary flux  given
    any boundary condition (BC). 
"""
function boundary_flux(
    bc::AbstractBC,
    bound_type::AbstractBoundary,
    Δz,
    _...,
)::ClimaCore.Fields.Field end
