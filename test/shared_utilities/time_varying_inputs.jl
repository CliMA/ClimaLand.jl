using Test
import ClimaLSM
import ClimaLSM: TimeVaryingInputs

import ClimaCore: Domains, Geometry, Fields, Meshes, Topologies, Spaces
import ClimaComms

device = ClimaComms.device()
context = ClimaComms.SingletonCommsContext(device)

@testset "Analytic TimeVaryingInput" begin
    fun = (x) -> 2x
    input = TimeVaryingInputs.TimeVaryingInput(fun)

    FT = Float32

    # Prepare a field
    domain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(0),
        Geometry.ZPoint{FT}(5),
        boundary_names = (:bottom, :top),
    )
    mesh = Meshes.IntervalMesh(domain; nelems = 10)
    topology = Topologies.IntervalTopology(context, mesh)

    column_space = Spaces.CenterFiniteDifferenceSpace(topology)
    column_field = Fields.zeros(column_space)

    TimeVaryingInputs.evaluate!(column_field, input, 10.0)
    @test Array(parent(column_field))[1] == fun(10.0)
end

@testset "Temporal TimeVaryingInput 1D" begin
    # Check times not sorted
    xs = [1.0, 0.0]
    ys = [1.0, 2.0]

    @test_throws ErrorException TimeVaryingInputs.TimeVaryingInput(xs, ys)

    for FT in (Float32, Float64)
        # Prepare spaces/fields

        domain = Domains.IntervalDomain(
            Geometry.ZPoint{FT}(0),
            Geometry.ZPoint{FT}(5),
            boundary_names = (:bottom, :top),
        )
        mesh = Meshes.IntervalMesh(domain; nelems = 10)
        topology = Topologies.IntervalTopology(context, mesh)

        column_space = Spaces.CenterFiniteDifferenceSpace(topology)
        point_space = Spaces.level(column_space, 1)
        column_field = Fields.zeros(column_space)
        point_field = Fields.zeros(point_space)

        times = collect(range(FT(0), FT(8π), 100))
        vals = sin.(times)

        # Nearest neighbor interpolation
        input = TimeVaryingInputs.TimeVaryingInput(
            times,
            vals;
            context,
            method = TimeVaryingInputs.NearestNeighbor(),
        )

        # Test in
        @test FT(3.0) in input
        @test !(FT(-3.0) in input)

        # Test with different types of spaces
        for dest in (point_field, column_field)
            # Time outside of range
            @test_throws ErrorException TimeVaryingInputs.evaluate!(
                dest,
                input,
                FT(-4),
            )

            TimeVaryingInputs.evaluate!(dest, input, times[10])

            @test Array(parent(dest))[1] == vals[10]

            # Linear interpolation
            input = TimeVaryingInputs.TimeVaryingInput(times, vals; context)

            TimeVaryingInputs.evaluate!(dest, input, 0.1)

            index = searchsortedfirst(times, 0.1)
            @test times[index - 1] <= 0.1 <= times[index]
            expected =
                vals[index - 1] +
                (vals[index] - vals[index - 1]) /
                (times[index] - times[index - 1]) * (0.1 - times[index - 1])

            @test Array(parent(dest))[1] ≈ expected

            # Check edge case
            TimeVaryingInputs.evaluate!(dest, input, 0.0)

            @test Array(parent(dest))[1] ≈ 0.0
        end
    end
end
