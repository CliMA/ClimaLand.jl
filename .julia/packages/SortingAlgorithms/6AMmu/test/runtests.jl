using SortingAlgorithms
using Test
using Random

stable_algorithms = [TimSort, RadixSort, PagedMergeSort]
unstable_algorithms = [HeapSort, CombSort]

a = rand(1:10000, 1000)
am = [rand() < .9 ? i : missing for i in a]

for alg in [stable_algorithms; unstable_algorithms; SortingAlgorithms.TimSortAlg()]
    b = sort(a, alg=alg)
    @test issorted(b)
    ix = sortperm(a, alg=alg)
    b = a[ix]
    @test issorted(b)
    @test a[ix] == b

    # legacy 3-argument calling convention
    @test b == sort!(copy(a), alg, Base.Order.Forward)

    b = sort(a, alg=alg, rev=true)
    @test issorted(b, rev=true)
    ix = sortperm(a, alg=alg, rev=true)
    b = a[ix]
    @test issorted(b, rev=true)
    @test a[ix] == b

    b = sort(a, alg=alg, by=x->1/x)
    @test issorted(b, by=x->1/x)
    ix = sortperm(a, alg=alg, by=x->1/x)
    b = a[ix]
    @test issorted(b, by=x->1/x)
    @test a[ix] == b

    c = copy(a)
    permute!(c, ix)
    @test c == b

    invpermute!(c, ix)
    @test c == a

    if alg != RadixSort  # RadixSort does not work with Lt orderings or missing
        c = sort(a, alg=alg, lt=(>))
        @test b == c

        # Issue https://github.com/JuliaData/DataFrames.jl/issues/3340
        bm1 = sort(am, alg=alg)
        @test issorted(bm1)
        @test count(ismissing, bm1) == count(ismissing, am)

        bm2 = am[sortperm(am, alg=alg)]
        @test issorted(bm2)
        @test count(ismissing, bm2) == count(ismissing, am)

        bm3 = am[sortperm!(collect(eachindex(am)), am, alg=alg)]
        @test issorted(bm3)
        @test count(ismissing, bm3) == count(ismissing, am)

        if alg == TimSort # Stable
            @test all(bm1 .=== bm2 .=== bm3)
        end
    end

    c = sort(a, alg=alg, by=x->1/x)
    @test b == c
end

randnans(n) = reinterpret(Float64,[rand(UInt64)|0x7ff8000000000000 for i=1:n])

function randn_with_nans(n,p)
    v = randn(n)
    x = findall(rand(n).<p)
    v[x] = randnans(length(x))
    return v
end

Random.seed!(0xdeadbeef)

for n in [0:10..., 100, 101, 1000, 1001]
    r = 1:10
    v = rand(1:10,n)

    for ord in [Base.Order.Forward, Base.Order.Reverse]
        # insertion sort (stable) as reference
        pi = sortperm(v, alg=InsertionSort, order=ord)
        @test isperm(pi)
        si = v[pi]
        @test issorted(si, order=ord)
        @test all(issorted,[pi[si.==x] for x in r])
        c = copy(v)
        permute!(c, pi)
        @test c == si
        invpermute!(c, pi)
        @test c == v

        for alg in stable_algorithms
            p = sortperm(v, alg=alg, order=ord)
            @test p == pi
            s = copy(v)
            permute!(s, p)
            @test s == si
            invpermute!(s, p)
            @test s == v
        end

        for alg in unstable_algorithms
            p = sortperm(v, alg=alg, order=ord)
            @test isperm(p)
            @test v[p] == si
            s = copy(v)
            permute!(s, p)
            @test s == si
            invpermute!(s, p)
            @test s == v
        end
    end

    v = randn_with_nans(n,0.1)
    for ord in [Base.Order.Forward, Base.Order.Reverse],
        alg in [stable_algorithms; unstable_algorithms]
        # test float sorting with NaNs
        s = sort(v, alg=alg, order=ord)
        @test issorted(s, order=ord)

        # This tests that NaNs (which compare equivalent) are treated stably
        # even when the underlying algorithm is unstable. That it happens to
        # pass is not a part of the public API:
        @test reinterpret(UInt64, v[map(isnan, v)]) == reinterpret(UInt64, s[map(isnan, s)])

        # test float permutation with NaNs
        p = sortperm(v, alg=alg, order=ord)
        @test isperm(p)
        vp = v[p]
        @test isequal(vp,s)
        @test reinterpret(UInt64,vp) == reinterpret(UInt64,s)
    end
end

for T in (Float64, Int, UInt8)
    for alg in stable_algorithms
        for ord in [Base.Order.By(identity), Base.Order.By(_ -> 0), Base.Order.By(Base.Fix2(÷, 100))]
            for n in vcat(0:31, 40:11:100, 110:51:1000)
                v = rand(T, n)
                # use MergeSort to guarantee stable sorting in Julia 1.0
                @test sort(v, alg=alg, order=ord) == sort(v, alg=MergeSort, order=ord)
            end
        end
    end
end

# PagedMergeSort with small input without InitialOptimizations
# (https://github.com/JuliaCollections/SortingAlgorithms.jl/pull/71#discussion_r1292774352)
if isdefined(Base.Sort, :InitialOptimizations)
    v = [0,1]
    @test sort(v, alg=SortingAlgorithms.PagedMergeSortAlg()) == sort(v, alg=MergeSort)
end

include("aqua.jl")
