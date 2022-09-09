"""
    struct SphericalShell{FT} <: AbstractDomain{FT}
        radius::FT
        height::FT
        nelements::Tuple{Int, Int}
        npolynomial::Int
    end
A struct holding the necessary information to construct a domain, a mesh, 
a 2d spectral element space (non-radial directions) 
x a 1d finite difference space (radial direction),
 and the resulting coordinate field.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct SphericalShell{FT, S} <: AbstractDomain{FT}
    "The radius of the shell"
    radius::FT
    "The radial extent of the shell"
    height::FT
    "The number of elements to be used in the non-radial and radial directions"
    nelements::Tuple{Int, Int}
    "The polynomial order to be used in the non-radial directions"
    npolynomial::Int
    "The associated ClimaCore Space"
    space::S
end

"""
    SphericalShell(;
        radius::FT,
        height::FT,
        nelements::Tuple{Int, Int},
        npolynomial::Int,
    ) where {FT}
Outer constructor for the `SphericalShell` domain, using keyword arguments.
"""
function SphericalShell(;
    radius::FT,
    height::FT,
    nelements::Tuple{Int, Int},
    npolynomial::Int,
) where {FT}
    @assert 0 < radius
    @assert 0 < height
    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(FT(0)),
        ClimaCore.Geometry.ZPoint(FT(height));
        boundary_tags = (:bottom, :top),
    )

    vertmesh = ClimaCore.Meshes.IntervalMesh(
        vertdomain,
        ClimaCore.Meshes.Uniform(),
        nelems = nelements[2],
    )
    vert_center_space = ClimaCore.Spaces.CenterFiniteDifferenceSpace(vertmesh)

    horzdomain = ClimaCore.Domains.SphereDomain(radius)
    horzmesh = ClimaCore.Meshes.EquiangularCubedSphere(horzdomain, nelements[1])
    horztopology = ClimaCore.Topologies.Topology2D(horzmesh)
    quad = ClimaCore.Spaces.Quadratures.GLL{npolynomial + 1}()
    horzspace = ClimaCore.Spaces.SpectralElementSpace2D(horztopology, quad)

    hv_center_space = ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(
        horzspace,
        vert_center_space,
    )
    return SphericalShell{FT, typeof(hv_center_space)}(
        radius,
        height,
        nelements,
        npolynomial,
        hv_center_space,
    )
end


"""
    struct SphericalSurface{FT} <: AbstractDomain{FT}
        radius::FT
        nelements::Tuple{Int, Int}
        npolynomial::Int
    end
A struct holding the necessary information to construct a domain, a mesh, 
a 2d spectral element space (non-radial directions) and the resulting coordinate field.
# Fields
$(DocStringExtensions.FIELDS)
"""
struct SphericalSurface{FT, S} <: AbstractDomain{FT}
    "The radius of the surface"
    radius::FT
    "The number of elements to be used in the non-radial directions"
    nelements::Int
    "The polynomial order to be used in the non-radial directions"
    npolynomial::Int
    "The associated ClimaCore Space"
    space::S
end

"""
    SphericalSurface(;
        radius::FT,
        nelements::Int
        npolynomial::Int,
    ) where {FT}
Outer constructor for the `SphericalSurface` domain, using keyword arguments.
"""
function SphericalSurface(;
    radius::FT,
    nelements::Int,
    npolynomial::Int,
) where {FT}
    @assert 0 < radius
    horzdomain = ClimaCore.Domains.SphereDomain(radius)
    horzmesh = Meshes.EquiangularCubedSphere(horzdomain, nelements)
    horztopology = Topologies.Topology2D(horzmesh)
    quad = Spaces.Quadratures.GLL{npolynomial + 1}()
    horzspace = Spaces.SpectralElementSpace2D(horztopology, quad)
    return SphericalSurface{FT, typeof(horzspace)}(
        radius,
        nelements,
        npolynomial,
        horzspace,
    )
end



"""
    make_lsm_domain(atmos_boundary_space::ClimaCore.Spaces.SpectralElementSpace2D, zlim::Tuple{FT,FT}, nelements::Int) where {FT}

A helper function for coupled simulations where we want the LSM and the Atmosphere model to have the same resolution in the horizontal, assuming
that the atmosphere space is created first. 
"""
function make_lsm_domain(atmos_boundary_space::ClimaCore.Spaces.SpectralElementSpace2D, zlim::Tuple{FT,FT}, nelements::Int) where {FT}
    @assert zlim[1] < zlim[2]
    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(zlim[1]),
        ClimaCore.Geometry.ZPoint(zlim[2]);
        boundary_tags = (:bottom, :top),
    )

    vertmesh = ClimaCore.Meshes.IntervalMesh(
        vertdomain,
        ClimaCore.Meshes.Uniform(),
        nelems = nelements,
    )
    vert_center_space = ClimaCore.Spaces.CenterFiniteDifferenceSpace(vertmesh)
    subsurface_space = ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(
        atmos_boundary_space,
        vert_center_space,
    )
    horzmesh = compute mesh
    if typeof(horzmesh) <: Meshes.EquiangularCubedSphere
        radius = 0.0
        height = zlim[2]-zlim[1]
        nelements = (horz, nelements)
        npolynomial = 2
        surface_domain = SphericalSurface{FT, typeof(atmos_boundary_space)}(radius, nelements,npolynomial, atmos_boundary_space)
        subsurface_domain = SphericalShell{FT, typeof(subsurface_space)}(radius,height nelements,npolynomial, subsurface_space)
        return LSMSphericalShellDomain{FT, typeof(subsurface_domain), typeof(surface_domain)}(subsurface_domain, surface_domain)
    else
        xlim = [0,0]
        ylim = [0,0]
        nelements = (1,2,nelements)
        npolynomial = 2
        surface_domain = Plane{FT, typeof(atmos_boundary_space)}(xlim, ylim, nelements[1:2], (true,true), npolynomial,atmos_boundary_space)
        subsurface_domain = HybridBox{FT, typeof(subsurface_space)}(xlim, ylim, zlim, nelements, (true,true), npolynomial, subsurface_space)
        return LSMMultiColumnDomain{FT, typeof(subsurface_domain), typeof(surface_domain)}(subsurface_domain, surface_domain)
    end
end
