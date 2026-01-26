import Dates

import SciMLBase

import ClimaCore
import ClimaComms
@static if pkgversion(ClimaComms) >= v"0.6"
    ClimaComms.@import_required_backends
end
import ClimaTimeSteppers

function ColumnCenterFiniteDifferenceSpace(
    zelem = 10,
    context = ClimaComms.SingletonCommsContext();
    FT = Float64,
)
    return _column(
        zelem,
        ClimaCore.Spaces.CenterFiniteDifferenceSpace,
        context,
        FT,
    )
end

function ColumnFaceFiniteDifferenceSpace(
    zelem = 10,
    context = ClimaComms.SingletonCommsContext();
    FT = Float64,
)
    return _column(
        zelem,
        ClimaCore.Spaces.FaceFiniteDifferenceSpace,
        context,
        FT,
    )
end

function _column(zelem, constructor, context, FT)
    zlim = (FT(0.0), FT(1.0))
    domain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(zlim[1]),
        ClimaCore.Geometry.ZPoint(zlim[2]);
        boundary_names = (:bottom, :top),
    )
    mesh = ClimaCore.Meshes.IntervalMesh(domain, nelems = zelem)
    topology = ClimaCore.Topologies.IntervalTopology(context, mesh)
    return constructor(topology)
end


function SphericalShellSpace(;
    radius = 6371.0,
    height = 10.0,
    nelements = 10,
    zelem = 10,
    npolynomial = 4,
    context = ClimaComms.context(),
    FT = Float64,
)
    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(FT(0)),
        ClimaCore.Geometry.ZPoint(FT(height));
        boundary_names = (:bottom, :top),
    )
    vertmesh = ClimaCore.Meshes.IntervalMesh(vertdomain; nelems = zelem)
    if pkgversion(ClimaCore) >= v"0.14.10"
        vert_center_space = ClimaCore.Spaces.CenterFiniteDifferenceSpace(
            ClimaComms.device(context),
            vertmesh,
        )
    else
        vert_center_space =
            ClimaCore.Spaces.CenterFiniteDifferenceSpace(vertmesh)
    end

    horzdomain = ClimaCore.Domains.SphereDomain(FT(radius))
    horzmesh = ClimaCore.Meshes.EquiangularCubedSphere(horzdomain, nelements)
    horztopology = ClimaCore.Topologies.Topology2D(context, horzmesh)
    quad = ClimaCore.Spaces.Quadratures.GLL{npolynomial + 1}()
    horzspace = ClimaCore.Spaces.SpectralElementSpace2D(horztopology, quad)

    return ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(
        horzspace,
        vert_center_space,
    )
end

function BoxSpace(;
    FT = Float64,
    xlim = (-FT(1), FT(1)),
    ylim = (-FT(1), FT(1)),
    height = 10.0,
    lonlat = false,
    nelements = (10, 10),
    zelem = 10,
    npolynomial = 4,
    context = ClimaComms.context(),
)
    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(FT(0)),
        ClimaCore.Geometry.ZPoint(FT(height));
        boundary_names = (:bottom, :top),
    )
    vertmesh = ClimaCore.Meshes.IntervalMesh(vertdomain; nelems = zelem)
    if pkgversion(ClimaCore) >= v"0.14.10"
        vert_center_space = ClimaCore.Spaces.CenterFiniteDifferenceSpace(
            ClimaComms.device(context),
            vertmesh,
        )
    else
        vert_center_space =
            ClimaCore.Spaces.CenterFiniteDifferenceSpace(vertmesh)
    end

    # NOTE: here we assume LatPoint is first
    XPoint = lonlat ? ClimaCore.Geometry.LatPoint : ClimaCore.Geometry.XPoint
    YPoint = lonlat ? ClimaCore.Geometry.LongPoint : ClimaCore.Geometry.YPoint
    periodic = !lonlat
    boundary_names = (:first, :second)

    domain_x = ClimaCore.Domains.IntervalDomain(
        XPoint(xlim[1]),
        XPoint(xlim[2]);
        periodic,
        boundary_names,
    )

    domain_y = ClimaCore.Domains.IntervalDomain(
        YPoint(ylim[1]),
        YPoint(ylim[2]);
        periodic,
        boundary_names,
    )

    horzdomain = ClimaCore.Domains.RectangleDomain(domain_x, domain_y)
    horzmesh =
        ClimaCore.Meshes.RectilinearMesh(horzdomain, nelements[1], nelements[2])
    horztopology = ClimaCore.Topologies.Topology2D(context, horzmesh)
    quad = ClimaCore.Spaces.Quadratures.GLL{npolynomial + 1}()
    horzspace = ClimaCore.Spaces.SpectralElementSpace2D(horztopology, quad)

    return ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(
        horzspace,
        vert_center_space,
    )
end

"""
    create_problem()

An ODE problem for an exponential decay.
"""
function create_problem(space; t0 = 0.0, tf = 1.0, dt = 1e-3)
    # Let's solve an exponential decay

    Y = ClimaCore.Fields.FieldVector(; my_var = ones(space))
    p = (; tau = -0.1, start_date = Dates.DateTime(476, 9, 4))

    function exp_tendency!(dY, Y, p, t)
        @. dY.my_var = p.tau * Y.my_var
    end

    prob = SciMLBase.ODEProblem(
        ClimaTimeSteppers.ClimaODEFunction(T_exp! = exp_tendency!),
        Y,
        (t0, tf),
        p,
    )
    algo = ClimaTimeSteppers.ExplicitAlgorithm(ClimaTimeSteppers.RK4())

    args = prob, algo
    kwargs = Dict(:dt => dt)

    return args, kwargs
end
