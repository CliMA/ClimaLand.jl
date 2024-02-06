using Test
import ClimaLand
using ClimaLand.SpaceVaryingInputs: SpaceVaryingInput
using ClimaLand: Domains

using ClimaCore: Fields
using ClimaComms

device = ClimaComms.device()
context = ClimaComms.SingletonCommsContext(device)
AT = ClimaComms.array_type(device)
# This tests the analytic and 1d cases of the SpaceVaryingInput function
@testset "SpaceVaryingInput" begin
    FT = Float32
    zmin = FT(-1.0)
    zmax = FT(0.0)
    xlim = FT.((0.0, 10.0))
    ylim = FT.((0.0, 1.0))
    zlim = FT.((zmin, zmax))
    nelements = (1, 1, 10)
    radius = FT(100)
    depth = FT(30)
    n_elements_sphere = (6, 20)
    npoly_sphere = 3
    shell = Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = n_elements_sphere,
        npolynomial = npoly_sphere,
    )

    box = Domains.HybridBox(;
        xlim = xlim,
        ylim = ylim,
        zlim = zlim,
        nelements = nelements,
        npolynomial = 0,
    )

    column = Domains.Column(; zlim = zlim, nelements = nelements[3])

    domains = [shell, box, column]
    analytic_func = (coords) -> 2.0
    for domain in domains
        for space in domain.space
            coords = Fields.coordinate_field(space)
            @test SpaceVaryingInput(analytic_func, space) ==
                  FT.(analytic_func.(coords))
        end
    end

    # 1D cases
    data_z = collect(range(FT(-1.0), FT(0.0), 11))
    data_value = data_z .* 2
    space = column.space.subsurface
    field = SpaceVaryingInput(data_z, data_value, space)
    @test parent(field)[:] ≈ AT(collect(range(FT(-1.9), FT(-0.1), 10)))

    struct Tmp{FT}
        a::FT
        b::FT
        c::FT
        function Tmp{FT}(; a::FT, c::FT) where {FT}
            b = a * 2
            new{FT}(a, b, c)
        end
    end
    data_values = (; a = data_z .* 2, c = data_z .* 3)
    field_of_structs = SpaceVaryingInput(data_z, data_values, space, Tmp{FT})
    @test eltype(field_of_structs) == Tmp{FT}
    @test field_of_structs.a == field
    @test field_of_structs.b == 2 .* field_of_structs.a
    @test parent(field_of_structs.c)[:] ≈
          AT(collect(range(FT(-2.85), FT(-0.15), 10)))

    # 2D SpaceVaryingInput
    # the 2d from netcdf + regridding is tested now implicitly in the albedo_types test, since the BulkAlbedoStatic
    # uses SpaceVaryingInput now. 
end
