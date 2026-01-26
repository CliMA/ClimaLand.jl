using Test

import ClimaUtilities
using ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput

import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
using ClimaCore
using Interpolations

const context = ClimaComms.context()
ClimaComms.init(context)

include("TestTools.jl")

AT = ClimaComms.array_type(ClimaComms.device())

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

    spaces = make_spherical_space(FT; context)
    column = spaces.vertical

    analytic_func = (coords) -> 2.0
    for space in (spaces.horizontal, spaces.vertical)
        coords = ClimaCore.Fields.coordinate_field(space)
        @test SpaceVaryingInput(analytic_func, space) ==
              FT.(analytic_func.(coords))
    end

    # 1D cases
    data_z = collect(range(FT(0.0), FT(1.0), 11))
    data_value = data_z .* 2
    field = SpaceVaryingInput(data_z, data_value, column)
    @test parent(field)[:] ≈ AT(collect(range(FT(0.1), FT(1.9), 10)))

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
    field_of_structs = SpaceVaryingInput(data_z, data_values, column, Tmp{FT})
    @test eltype(field_of_structs) == Tmp{FT}
    @test field_of_structs.a == field
    @test field_of_structs.b == 2 .* field_of_structs.a
    @test parent(field_of_structs.c)[:] ≈
          AT(collect(range(FT(0.15), FT(2.85), 10)))

end
