using Atomix
using Atomix.Internal: referenceable
using Atomix: @atomic, @atomicreplace, @atomicswap
using Test


@testset "Aqua.jl" begin
    using Aqua
    Aqua.test_all(Atomix)
end


@testset "test_indexableref" begin
    A = ones(Int, 3)
    ref = Atomix.IndexableRef(A, (1,))
    @test eltype(ref) === Int
    @test Atomix.get(ref) === 1
    Atomix.set!(ref, 123)
    @test Atomix.get(ref) === 123
    @test Atomix.modify!(ref, +, 1) === (123 => 124)
    @test Atomix.get(ref) === 124
    @test Atomix.swap!(ref, 345) == 124
    @test Atomix.get(ref) === 345
    @test Atomix.replace!(ref, 345, 567) === (old = 345, success = true)
    @test Atomix.replace!(ref, 345, 567) === (old = 567, success = false)
end


@testset "test_referenceablearray" begin
    for a in Any[
        ones(Int, 3),
        view(ones(Int, 3), 1:2),
        view(ones(Int, 2, 3), 1:1, 1:2),
        ones(Int, 3)',
    ]
        ra = referenceable(a)
        @test size(ra) == size(a)
        @test IndexStyle(ra) == IndexStyle(a)
    end
end


@testset "test_fallback" begin

    mutable struct Atomic{T}
        @atomic x::T
    end
    
    a = Atomic(123)
    @test (@atomic a.x) == 123
    @test (@atomic :monotonic a.x) == 123
    @atomic a.x = 456
    @test (@atomic a.x) == 456
    @atomic :monotonic a.x = 123
    @test (@atomic a.x) == 123
    @test (@atomic a.x += 111) == 234
    @test (@atomic :monotonic a.x += 111) == 345
    @test (@atomic a.x + 111) == (345 => 456)
    @test (@atomic :monotonic a.x + 111) == (456 => 567)
    @test (@atomicswap a.x = 123) == 567
    @test (@atomicswap :monotonic a.x = 234) == 123
    @test (@atomicreplace a.x 234 => 123) == (old = 234, success = true)
    @test (@atomicreplace a.x 234 => 123) == (old = 123, success = false)
    @test (@atomicreplace :monotonic a.x 123 => 234) == (old = 123, success = true)
    @test (@atomicreplace :monotonic a.x 123 => 234) == (old = 234, success = false)
    @test (@atomicreplace :monotonic :monotonic a.x 234 => 123) ==
          (old = 234, success = true)
    @test (@atomicreplace :monotonic :monotonic a.x 234 => 123) ==
          (old = 123, success = false)
end


@testset "test_get" begin
    A = [42]
    @test (@atomic A[1]) === 42
    @test (@atomic A[end]) === 42
    @test (@atomic :monotonic A[begin]) === 42
    order = :monotonic
    @test (@atomic order A[begin]) === 42
end


@testset "test_get_2d" begin
    A = view([11 12; 21 22], 1:2, 1:2)
    @test IndexStyle(A) isa IndexCartesian
    @test (@atomic A[1]) === 11
    @test (@atomic A[2]) === 21
    @test (@atomic A[end]) === 22
    @test (@atomic A[2, 1]) === 21
    @test (@atomic A[end, 1]) === 21
    @test (@atomic A[1, end]) === 12
end


@testset "test_set" begin
    A = [42, 43]
    @atomic A[1] = 123
    @atomic A[end] = 124
    @test A[1] === 123
    @test A[end] === 124
end


@testset "test_inc" begin
    A = [1, 1]
    @test (@atomic A[1] += 123) === 124
    @test A[1] === 124
    @test (@atomic A[end] += 123) === 124
    @test A[end] === 124
end


@testset "test_swap" begin
    A = [1, 1]
    @test (@atomicswap A[1] = 123) === 1
    @test A[1] === 123
    @test (@atomicswap A[end] = 456) === 1
    @test A[end] === 456
end


@testset "test_cas" begin
    A = [1, 1]
    @test (@atomicreplace A[1] 1 => 123) == (old = 1, success = true)
    @test A[1] === 123
    @test (@atomicreplace A[1] 1 => 456) == (old = 123, success = false)
    @test A[1] === 123
    @test (@atomicreplace A[end] 1 => 789) == (old = 1, success = true)
    @test A[end] === 789
    update = 789 => 123
    @test (@atomicreplace A[end] update) == (old = 789, success = true)
    @test A[end] === 123
end


# KernelAbstractions backend tests
# Pass command-line argument to test suite to install the right backend, e.g.
#   julia> import Pkg
#   julia> Pkg.test("Atomix", test_args=["--Metal"])
if "--Metal" in ARGS
    import Pkg
    Pkg.add("Metal")
    include("test_atomix_metal.jl")
elseif "--CUDA" in ARGS
    import Pkg
    Pkg.add("CUDA")
    include("test_atomix_cuda.jl")
elseif "--oneAPI" in ARGS
    import Pkg
    Pkg.add("oneAPI")
    include("test_atomix_oneapi.jl")
elseif "--OpenCL" in ARGS
    import Pkg
    Pkg.add(["OpenCL", "pocl_jll"])
    include("test_atomix_opencl.jl")
end
