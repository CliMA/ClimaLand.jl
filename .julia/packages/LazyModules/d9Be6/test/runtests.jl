using LazyModules
using Test

module LazyOffsetArrays
    using LazyModules

    @lazy import OffsetArrays ="6fe1bfb0-de20-5000-8ca7-80f57d26f881"

    function zero_based(A)
        o = Base.invokelatest(OffsetArrays.Origin, 0)
        return Base.invokelatest(OffsetArrays.OffsetArray, A, o)
    end
    export zero_based
end

module LazyColors
    using LazyModules

    @lazy import Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    function fancy_color()
        LazyModules.require(Colors)
        return zero(Colors.RGB)
    end

    export fancy_color
end

@testset "LazyModules" begin
    using .LazyOffsetArrays

    A = rand(10, 10)
    AO = zero_based(A)
    @test axes(AO) == (0:9, 0:9)

    using OffsetArrays
    @test AO isa OffsetArray

    @static if VERSION >= v"1.6"
        @lazy import OffsetArrays as FOO = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
        AOO = FOO.OffsetArray(A, -1, -1)
        @test AOO == AO
    end

    using .LazyColors
    msg = "Colors is required to be loaded first, maybe `using Colors` or `import Colors` and try again."
    @test_throws ErrorException(msg) fancy_color()
    using Colors
    @test fancy_color() == zero(RGB)
end
