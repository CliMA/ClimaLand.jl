using CUDA
using CUDA: @allowscalar


@testset "AtomixCUDAExt:extension_found" begin
    @test !isnothing(Base.get_extension(Atomix, :AtomixCUDAExt))
end


function cuda(f)
    function g()
        f()
        nothing
    end
    CUDA.@cuda g()
end


# Not implemented:
#=
function test_get_set()
    A = CUDA.ones(Int, 3)
    cuda() do
        GC.@preserve A begin
            ref = Atomix.IndexableRef(A, (1,))
            x = Atomix.get(ref)
            Atomix.set!(ref, -x)
        end
    end
    @test collect(A) == [-1, 1, 1]
end
=#


@testset "AtomixCUDAExt:test_cas" begin
    idx = (
        data = 1,
        cas1_ok = 2,
        cas2_ok = 3,
        # ...
    )
    @assert minimum(idx) >= 1
    @assert maximum(idx) == length(idx)

    A = CUDA.zeros(Int, length(idx))
    cuda() do
        GC.@preserve A begin
            ref = Atomix.IndexableRef(A, (1,))
            (old, success) = Atomix.replace!(ref, 0, 42)
            A[idx.cas1_ok] = old == 0 && success
            (old, success) = Atomix.replace!(ref, 0, 43)
            A[idx.cas2_ok] = old == 42 && !success
        end
    end
    @test collect(A) == [42, 1, 1]
end


@testset "AtomixCUDAExt:test_inc" begin
    A = CUDA.CuVector(1:3)
    cuda() do
        GC.@preserve A begin
            ref = Atomix.IndexableRef(A, (1,))
            pre, post = Atomix.modify!(ref, +, 1)
            A[2] = pre
            A[3] = post
        end
    end
    @test collect(A) == [2, 1, 2]
end


@testset "AtomixCUDAExt:test_inc_sugar" begin
    A = CUDA.ones(Int, 3)
    cuda() do
        GC.@preserve A begin
            @atomic A[begin] += 1
        end
    end
    @test collect(A) == [2, 1, 1]
end
