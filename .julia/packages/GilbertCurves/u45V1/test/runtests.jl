using GilbertCurves
using Test

@testset "size $m,$n" for m = 1:20, n = 1:20
    list = gilbertindices((m,n))
    @test length(list) == m*n
    @test list[1] == CartesianIndex(1,1)
    if m >= n
        if !(n == 2 && isodd(m))
            @test list[end] == CartesianIndex(m,1)
        end
    else
        if !(m == 2 && isodd(n))
            @test list[end] == CartesianIndex(1,n)
        end
    end

    ndiag = 0
    for i = 1:m*n-1
        Δ = map(abs, (list[i+1] - list[i]).I)
        @test Δ[1] <= 1
        @test Δ[2] <= 1
        if Δ[1] + Δ[2] > 1
            ndiag += 1
        end        
    end
    if m > n && isodd(m) && iseven(n) || m < n && isodd(n) && iseven(m)
        @test ndiag <= 1
    else
        @test ndiag == 0
    end
    
    L = GilbertCurves.linearindices(list)
    @test !any(iszero, L)
    @test sum(L) == sum(1:m*n)
end
