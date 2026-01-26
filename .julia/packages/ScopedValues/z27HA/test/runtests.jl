using Test
using ScopedValues

@testset "errors" begin
    @test ScopedValue{Float64}(1)[] == 1.0
    @test_throws InexactError ScopedValue{Int}(1.5)
    val = ScopedValue(1)
    @test_throws MethodError val[] = 2
    with() do
        @test_throws MethodError val[] = 2
    end
    val = ScopedValue{Int}()
    @test_throws KeyError val[]
    @test_throws MethodError ScopedValue()

    sval = ScopedValue(1)
    with(sval=>2.0) do
        @test sval[] == 2
    end
end

const sval = ScopedValue(1)
@testset "inheritance" begin
    @test sval[] == 1
    with() do
        @test sval[] == 1
        with() do
            @test sval[] == 1
        end
        with(sval => 2) do
            @test sval[] == 2
        end
        @test sval[] == 1
    end
    @test sval[] == 1
end

const sval_float = ScopedValue(1.0)

@testset "multiple scoped values" begin
    with(sval => 2, sval_float => 2.0) do
        @test sval[] == 2
        @test sval_float[] == 2.0
    end
    with(sval => 2, sval => 3) do
        @test sval[] == 3
    end
end

emptyf() = nothing

@testset "conversion" begin
    with(emptyf, sval_float=>2)
    @test_throws MethodError with(emptyf, sval_float=>"hello")
end

import Base.Threads: @spawn
@testset "tasks" begin
    @test fetch(@spawn begin
        sval[]
    end) == 1
    with(sval => 2) do
        @test fetch(@spawn begin
            sval[]
        end) == 2
    end
end

@testset "show" begin
    @test sprint(show, ScopedValue{Int}()) == "ScopedValue{$Int}(undefined)"
    @test sprint(show, sval) == "ScopedValue{$Int}(1)"
    @test sprint(show, ScopedValues.current_scope()) == "nothing"
    with(sval => 2.0) do
        @test sprint(show, sval) == "ScopedValue{$Int}(2)"
        objid = sprint(show, Base.objectid(sval))
        # Interpolate to handle `Base.ScopedValues.Scope` vs `ScopedValues.scope`
        @test sprint(show, ScopedValues.current_scope()) == "$(ScopedValues.Scope)(ScopedValue{$Int}@$objid => 2)"
    end
end

const depth = ScopedValue(0)
function nth_with(f, n)
    if n <= 0
        f()
    else
        with(depth => n) do
            nth_with(f, n-1)
        end
    end
end


@testset "nested with" begin
    @testset for depth in 1:16
        nth_with(depth) do
            @test sval_float[] == 1.0
        end
        with(sval_float=>2.0) do
            nth_with(depth) do
                @test sval_float[] == 2.0
            end
        end
        nth_with(depth) do
            with(sval_float=>2.0) do
                @test sval_float[] == 2.0
            end
        end
    end
    with(sval_float=>2.0) do
        nth_with(15) do
            @test sval_float[] == 2.0
            with(sval_float => 3.0) do
                @test sval_float[] == 3.0
            end
        end
    end
end

@testset "macro" begin
    @with sval=>2 sval_float=>2.0 begin
        @test sval[] == 2
        @test sval_float[] == 2.0
    end
    # Doesn't do much...
    ret = @with begin
        @test sval[] == 1
        @test sval_float[] == 1.0
        1.23 # return value
    end
    @test ret == 1.23
end

@testset "ScopedFunctor" begin
    function check_svals()
        @test sval[] == 8
        @test sval_float[] == 8.0
    end
    sf = nothing
    @with sval=>8 sval_float=>8.0 begin
        sf = ScopedFunctor(check_svals)
    end
    sf()
    @with sval=>8 sval_float=>8.0 begin
        sf2 = ScopedFunctor{Function}(check_svals)
    end
    sf2()
end
