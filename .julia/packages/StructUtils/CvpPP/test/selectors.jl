using Test, StructUtils
using StructUtils.Selectors: @selectors

# A simple wrapper struct for Vector{Pair{String, Any}}
struct StringPairs
    data::Vector{Pair{String, Any}}
end

# Apply the selectors macro to StringPairs
@selectors StringPairs

# Override applyeach to iterate over the wrapped Dict
function StructUtils.applyeach(st::StructUtils.StructStyle, f, x::StringPairs)
    for (k, v) in getfield(x, :data)
        ret = f(k, v)
        ret isa StructUtils.EarlyReturn && return ret
    end
    return
end

# Recursive helper to wrap nested Dicts and Arrays
wrap(x::Vector{<:Pair}) = StringPairs([k => wrap(v) for (k, v) in x])
wrap(x::Vector) = [wrap(v) for v in x]
wrap(x) = x

# Sample test data
raw = [
    "a" => 1,
    "b" => ["c" => 2, "d" => 3],
    "e" => [4, 5, 6],
    "f" => [["g" => 7], ["g" => 8]],
]
x = wrap(raw)

@testset "Selectors on StringPairs" begin
    # property names at root
    @test propertynames(x) == [:a, :b, :e, :f]

    # direct selection
    @test x.a == 1
    @test x["b"]["c"] == 2
    @test x.b.d == 3

    # select all direct values
    vals = collect(x[:])
    @test length(vals) == 4
    @test vals[1] == 1
    @test isa(vals[2], StringPairs)
    @test vals[3] == [4, 5, 6]
    @test isa(vals[4], Vector)

    # recursive flatten of all values
    flat = collect(x[~, :])
    @test flat == [1, 2, 3, 4, 5, 6, 7, 8]

    # recursive selection by key
    gvals = collect(x[~, "g"])
    @test gvals == [7, 8]

    # selection by array of keys
    be = collect(x[["b", "e"]])
    @test length(be) == 2
    @test isa(be[1], StringPairs)
    @test be[2] == [4, 5, 6]

    # filter selection at root by key predicate
    eonly = collect(x[:, (k, v) -> k == "e"])
    @test eonly == [[4, 5, 6]]

    # missing key throws KeyError
    @test_throws KeyError x["missing"]
end