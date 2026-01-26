using SymbolicIndexingInterface
using StaticArrays

sys = SymbolCache([:x, :y, :z], [:a, :b, :c], :t)

for (buf, newbuf, idxs, vals) in [
    # standard operation
    ([1.0, 2.0, 3.0], [2.0, 3.0, 4.0], [:x, :y, :z], [2.0, 3.0, 4.0]),
    # buffer type "demotion"
    ([1.0, 2.0, 3.0], [2, 2, 3], [:x], [2]),
    # buffer type promotion
    ([1, 2, 3], [2.0, 2.0, 3.0], [:x], [2.0]),
    # value type promotion
    ([1, 2, 3], [2.0, 3.0, 4.0], [:x, :y, :z], Real[2, 3.0, 4.0]),
    # standard operation
    ([1.0, 2.0, 3.0], [2.0, 3.0, 4.0], [:a, :b, :c], [2.0, 3.0, 4.0]),
    # buffer type "demotion"
    ([1.0, 2.0, 3.0], [2, 2, 3], [:a], [2]),
    # buffer type promotion
    ([1, 2, 3], [2.0, 2.0, 3.0], [:a], [2.0]),
    # value type promotion
    ([1, 2, 3], [2, 3.0, 4.0], [:a, :b, :c], Real[2, 3.0, 4.0]),
    # skip non-parameters
    ([1, 2, 3], [2.0, 3.0, 3.0], [:a, :b, :(a + b)], [2.0, 3.0, 5.0])
]
    @test_deprecated remake_buffer(sys, buf, Dict(idxs .=> vals))
    for varmap in [Dict(idxs .=> vals), Dict{Any, Any}(idxs .=> vals)]
        _newbuf = remake_buffer(sys, buf, keys(varmap), values(varmap))
        @test _newbuf != buf
        @test newbuf == _newbuf
        @test typeof(newbuf) == typeof(_newbuf)
    end
    for arrType in [Vector, SVector{3}, MVector{3}, SizedVector{3}]
        buf = arrType(buf)
        newbuf = arrType(newbuf)
        _newbuf = remake_buffer(sys, buf, idxs, vals)

        @test _newbuf != buf # should not alias
        @test newbuf == _newbuf # test values
        @test typeof(newbuf) == typeof(_newbuf) # ensure appropriate type
    end
end

for (buf, newbuf, idxs, vals) in [
    # standard operation
    ((1.0, 2.0, 3.0), (2.0, 3.0, 4.0), [:a, :b, :c], [2.0, 3.0, 4.0]),
    # buffer type "demotion"
    ((1.0, 2.0, 3.0), (2, 3, 4), [:a, :b, :c], [2, 3, 4]),
    # buffer type promotion
    ((1, 2, 3), (2.0, 3.0, 4.0), [:a, :b, :c], [2.0, 3.0, 4.0]),
    # value type promotion
    ((1, 2, 3), (2, 3.0, 4.0), [:a, :b, :c], Real[2, 3.0, 4.0]),
    # standard operation
    ((1.0, 2.0, 3.0), (2.0, 3.0, 4.0), [:x, :y, :z], [2.0, 3.0, 4.0]),
    # buffer type "demotion"
    ((1.0, 2.0, 3.0), (2, 3, 4), [:x, :y, :z], [2, 3, 4]),
    # buffer type promotion
    ((1, 2, 3), (2.0, 3.0, 4.0), [:x, :y, :z], [2.0, 3.0, 4.0]),
    # value type promotion
    ((1, 2, 3), (2, 3.0, 4.0), [:x, :y, :z], Real[2, 3.0, 4.0]),
    # skip non-variables
    ([1, 2, 3], [2.0, 3.0, 3.0], [:x, :y, :(x + y)], [2.0, 3.0, 5.0])
]
    @test_deprecated remake_buffer(sys, buf, Dict(idxs .=> vals))
    for varmap in [Dict(idxs .=> vals), Dict{Any, Any}(idxs .=> vals)]
        _newbuf = remake_buffer(sys, buf, keys(varmap), values(varmap))
        @test newbuf == _newbuf
        @test typeof(newbuf) == typeof(_newbuf)
    end
    _newbuf = remake_buffer(sys, buf, idxs, vals)
    @test newbuf == _newbuf # test values
    @test typeof(newbuf) == typeof(_newbuf) # ensure appropriate type
end

@test isnothing(remake_buffer(sys, nothing, [], []))

@testset "`remake_buffer` with `Dict`" begin
    sys = nothing
    buf = Dict("a" => 1, "b" => 2)
    buf2 = remake_buffer(sys, buf, collect(keys(buf)), collect(values(buf)))
    @test isequal(buf, buf2)
    buf2 = remake_buffer(sys, buf, buf)
    @test isequal(buf, buf2)
end
