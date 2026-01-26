using Test
import ClimaUtilities.DataStructures

@testset "Construct LRUCache" begin
    cache = DataStructures.LRUCache{Int, Int}(max_size = 10)
    @test isempty(cache.cache)
    @test isempty(cache.priority)
    @test cache.max_size == 10
end

@testset "setindex! for LRUCache" begin
    cache = DataStructures.LRUCache{String, Int}(max_size = 3)

    # Test adding new values (cache misses)
    @test setindex!(cache, 1, "a") == cache
    @test cache.priority == ["a"]
    cache["b"] = 2
    @test cache.priority == ["a", "b"]

    # Test updating existing key value
    cache["a"] = 3
    @test cache.cache == Dict{String, Int64}("a" => 3, "b" => 2)
    @test cache.priority == ["b", "a"]

    # Test enforcing size limit
    cache["c"] = 4
    cache["d"] = 5
    @test cache.cache == Dict{String, Int64}("a" => 3, "c" => 4, "d" => 5)
    @test cache.priority == ["a", "c", "d"]
    @test length(cache.priority) == length(keys(cache.cache)) == cache.max_size
end

@testset "getindex for LRUCache" begin
    cache = DataStructures.LRUCache{String, Int}(max_size = 3)

    # Test getting value for key that is in cache
    cache["a"] = 1
    @test getindex(cache, "a") == 1

    # Add "b", then access "a" to ensure getindex updates priority
    cache["b"] = 2
    cache["a"]
    @test cache.priority == ["b", "a"]

    # Test behavior if key not in cache
    @test_throws KeyError cache["c"]
end

@testset "length for LRUCache" begin
    cache = DataStructures.LRUCache{String, Int}(max_size = 3)

    @test length(cache) == 0

    cache["a"] = 1
    @test length(cache) == 1
end

@testset "== for LRUCache" begin
    cache1 = DataStructures.LRUCache{String, Int}(max_size = 3)
    cache2 = DataStructures.LRUCache{String, Int}(max_size = 4)

    @test cache1 != cache2

    cache3 = DataStructures.LRUCache{String, Int}(max_size = 3)

    @test cache1 == cache3

    cache1["a"] = 1

    @test cache1 != cache3
end

@testset "get! with default value" begin
    cache = DataStructures.LRUCache{String, Int}(max_size = 3)

    # Test adding new values (cache misses)
    @test get!(cache, "a", 1) == 1
    @test cache.priority == ["a"]
    @test get!(cache, "b", 2) == 2
    @test cache.priority == ["a", "b"]
    get!(cache, "c", 3)
    @test cache.priority == ["a", "b", "c"]

    # Test updating existing values (cache hits)
    @test get!(cache, "a", 4) == 1
    @test cache.priority == ["b", "c", "a"]

    # Test enforcing size limit
    get!(cache, "d", 5)
    @test cache.priority == ["c", "a", "d"]
    @test length(cache.priority) == length(keys(cache.cache)) == cache.max_size
end

@testset "get! with default callable" begin
    cache = DataStructures.LRUCache{String, Int}(max_size = 3)

    # Test adding new values (cache misses)
    @test get!(() -> 1, cache, "a") == 1
    @test cache.priority == ["a"]
    @test get!(() -> 2, cache, "b") == 2
    @test cache.priority == ["a", "b"]
    get!(() -> 3, cache, "c")
    @test cache.priority == ["a", "b", "c"]

    # Test updating existing values (cache hits)
    @test get!(() -> 4, cache, "a") == 1
    @test cache.priority == ["b", "c", "a"]

    # Test enforcing size limit
    get!(() -> 5, cache, "d")
    @test cache.priority == ["c", "a", "d"]
    @test length(cache.priority) == length(keys(cache.cache)) == cache.max_size
end

@testset "get with default value" begin
    cache = DataStructures.LRUCache{String, Int}(max_size = 3)

    # Test cache hits (key in dict)
    cache["a"] = 1
    cache["b"] = 2
    @test get(cache, "a", 3) == 1
    @test cache.priority == ["b", "a"]

    # Test cache miss (key not in dict)
    @test get(cache, "c", 3) == 3
    @test get(cache, "c", 2) == 2
    @test cache.priority == ["b", "a"]
end

@testset "get with default callable" begin
    cache = DataStructures.LRUCache{String, Int}(max_size = 3)

    # Test cache hits (key in dict)
    cache["a"] = 1
    cache["b"] = 2
    @test get(() -> 3, cache, "a") == 1
    @test cache.priority == ["b", "a"]

    # Test cache miss (key not in dict)
    @test get(() -> 3, cache, "c") == 3
    @test get(() -> 2, cache, "c") == 2
    @test cache.priority == ["b", "a"]
end

@testset "haskey LRU" begin
    cache = DataStructures.LRUCache{String, Int}(max_size = 3)

    # Test haskey if key exists
    cache["a"] = 1
    cache["b"] = 2
    @test haskey(cache, "a")

    # Test haskey doesnt mutate priority
    @test cache.priority == ["a", "b"]

    # Test if key does not exist
    @test !haskey(cache, "c")
    @test cache.priority == ["a", "b"]

    # Test if key added
    cache["c"] = 5
    @test haskey(cache, "c")
end

@testset "copy LRU" begin
    cache = DataStructures.LRUCache{String, Vector{Int64}}(max_size = 3)
    cache["a"] = [1]
    cache["b"] = [2]

    # Test copy for equivalence
    @test copy(cache).priority == cache.priority
    @test copy(cache).cache == cache.cache
    @test copy(cache).max_size == cache.max_size

    # Test that copying object not just reference by mutating a value in the copy
    copy(cache)["c"] = [3]
    @test cache.cache == Dict{String, Vector{Int64}}("a" => [1], "b" => [2])

    # Test that it is shallow copying
    value_for_a = copy(cache)["a"]
    push!(value_for_a, 1)
    @test cache.cache == Dict{String, Vector{Int64}}("a" => [1, 1], "b" => [2])
    @test cache.priority == ["a", "b"]
end

@testset "deepcopy" begin
    cache = DataStructures.LRUCache{String, Vector{Int64}}(max_size = 3)
    cache["a"] = [1]
    cache["b"] = [2]

    # Test copy for equivalence
    @test deepcopy(cache).priority == cache.priority
    @test deepcopy(cache).cache == cache.cache
    @test deepcopy(cache).max_size == cache.max_size

    # Test that it is deep copying
    value_for_a = deepcopy(cache)["a"]
    push!(value_for_a, 1)
    @test cache.cache == Dict{String, Vector{Int64}}("a" => [1], "b" => [2])
    @test cache.priority == ["a", "b"]
end

@testset "empty LRU" begin
    cache = DataStructures.LRUCache{String, Int64}(max_size = 3)
    cache["a"] = 1
    cache["b"] = 2

    # Test that it is returning a new empty LRUcache
    empty_cache = empty(cache)
    @test empty_cache != cache

    # Test that if not specified, key and value types are copied
    @test empty_cache.cache == Dict{String, Int64}()
    @test empty_cache.cache isa Dict{String, Int64}
    @test empty_cache.priority == []
    @test empty_cache.max_size == cache.max_size

    # Test empty if a second arg is given
    empty_cache = empty(cache, Char)
    @test empty_cache.cache == Dict{String, Char}()
    @test empty_cache.cache isa Dict{String, Char}
    @test empty_cache.priority isa Vector{String}

    # Test empty if two args are given
    empty_cache = empty(cache, Int64, String)
    @test empty_cache.cache == Dict{Int64, String}()
    @test empty_cache.cache isa Dict{Int64, String}
    @test empty_cache.priority isa Vector{Int64}
end

@testset "empty! LRU" begin
    cache = DataStructures.LRUCache{String, Int64}(max_size = 3)
    cache["a"] = 1
    cache["b"] = 2

    # Test that it is emptying LRUcache
    @test empty!(cache).cache == Dict{String, Int64}()
    @test cache.cache == Dict{String, Int64}()
    @test cache.cache isa Dict{String, Int64}
    @test cache.priority == []

    # Test that it is not returning a copy
    empty!(cache)["c"] = 3
    @test cache.cache == Dict{String, Int64}("c" => 3)
end

@testset "pop! LRU" begin
    cache = DataStructures.LRUCache{String, Int64}(max_size = 3)
    cache["a"] = 1
    cache["b"] = 2
    cache["c"] = 3

    # Test that it deletes and returns mapping for key if in cache
    @test pop!(cache, "b") == 2
    @test cache.cache == Dict{String, Int64}("a" => 1, "c" => 3)
    @test cache.priority == ["a", "c"]

    # Test popping key that doesn't exist without specifying default throws error
    @test_throws KeyError pop!(cache, "b")

    # Test that default returned if key not in cache
    @test pop!(cache, "b", 5) == 5
end

@testset "delete! LRU" begin
    cache = DataStructures.LRUCache{String, Int64}(max_size = 3)
    cache["a"] = 1
    cache["b"] = 2
    cache["c"] = 3

    # Test when key in cache, that it deletes and returns the cache
    @test delete!(cache, "b") isa DataStructures.LRUCache
    @test cache.cache == Dict{String, Int64}("a" => 1, "c" => 3)
    @test cache.priority == ["a", "c"]

    # Test that deleting key that isn't in cache just returns the cache
    @test delete!(cache, "b").cache == Dict{String, Int64}("a" => 1, "c" => 3)
    @test cache.priority == ["a", "c"]
end

@testset "merge LRU" begin
    cache1 = DataStructures.LRUCache{String, Int64}(max_size = 3)
    cache1["a"] = 1
    cache1["b"] = 2
    cache2 = DataStructures.LRUCache{String, Int64}(max_size = 3)
    cache2["c"] = 3
    cache2["a"] = 4

    # Test merging two LRUcaches throws error
    @test_throws ErrorException merge(cache1, cache2)
end
