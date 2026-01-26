using ImageCore, Colors, ColorVectorSpace
using Test, Statistics

# Different access patterns (getindex)
function mysum_elt_boundscheck(A)
    s = zero(eltype(A))
    for a in A
        s += a
    end
    s
end
function mysum_index_boundscheck(A)
    s = zero(eltype(A))
    for I in eachindex(A)
        s += A[I]
    end
    s
end
function mysum_elt_inbounds(A)
    s = zero(eltype(A))
    @inbounds for a in A
        s += a
    end
    s
end
function mysum_index_inbounds_simd(A)
    s = zero(eltype(A))
    @inbounds @simd for I in eachindex(A)
        s += A[I]
    end
    s
end
# setindex!
function myfill1!(A, val)
    f = convert(eltype(A), val)
    for I in eachindex(A)
        A[I] = f
    end
    A
end
function myfill2!(A, val)
    f = convert(eltype(A), val)
    @inbounds @simd for I in eachindex(A)
        A[I] = f
    end
    A
end

# Rather than using BenchmarkTools (and thus run one test repeatedly,
# accumulating timings), we run the same test interleaving the two
# array types. This is designed to reduce the risk of spurious
# failure, particularly on shared machines like Travis where they may
# get "distracted" by other tasks
function test_getindex(f, ar, cv, n)
    t_ar = Array{Float64}(undef, n)
    t_cv = Array{Float64}(undef, n)
    # Store the results to prevent the compiler from eliding the call
    f_ar = Ref(f(ar))
    f_cv = Ref(f(cv))
    @test f_ar[] ≈ f_cv[]  # but this also gives us a chance to test correctness
    for i = 1:n
        t_ar[i] = (tstart = time(); f_ar[] = f(ar); time()-tstart)
        t_cv[i] = (tstart = time(); f_cv[] = f(cv); time()-tstart)
    end
    median(t_ar), median(t_cv)
end
function test_setindex(f, ar, cv, n)
    t_ar = Array{Float64}(undef, n)
    t_cv = Array{Float64}(undef, n)
    for i = 1:n
        t_ar[i] = @elapsed f(ar, zero(eltype(ar)))
        t_cv[i] = @elapsed f(cv, zero(eltype(cv)))
    end
    median(t_ar), median(t_cv)
end

cc_getindex_funcs = (mysum_elt_boundscheck,
                     mysum_index_boundscheck,
                     mysum_elt_inbounds,
                     mysum_index_inbounds_simd)
cc_setindex_funcs = (myfill1!,
                     myfill2!)

# Performance tolerances
isfast = VERSION >= v"1.6.0-DEV.1083"
chanvtol = Dict{Any,Int}(mysum_index_inbounds_simd => isfast ? 3 : 20,
                         mysum_elt_boundscheck => isfast ? 3 : 20,
                         myfill1! => 20,
                         myfill2! => isfast ? 3 : 20)
chanvdefault = isfast ? 3 : 10
colvtol = Dict{Any,Int}(mysum_elt_boundscheck=>isfast ? 3 : 5,
                        mysum_index_boundscheck=>isfast ? 3 : 5)
colvdefault = 3

ssz = (1000,300)

@info "Benchmark tests are warnings for now"
# @testset "benchmarks" begin
for T in (Float32, Float64)
    c = rand(RGB{T}, ssz...)
    a = copy(reinterpretc(T, c))
    vchan = channelview(c)
    vcol = colorview(RGB, a)

    # view versions
    rview = 2:ssz[1]-1
    csub = view(c, rview, :)
    asub = view(a, :, rview, :)
    vchansub = channelview(csub)
    vcolsub = colorview(RGB, asub)

    for (suite, testf) in ((cc_getindex_funcs, test_getindex),
                           (cc_setindex_funcs, test_setindex))
        for f in suite
            # channelview
            t_ar, t_cv = testf(f, a, vchan, 30)
            tol = haskey(chanvtol, f) ? chanvtol[f] : chanvdefault
            if t_cv >= tol*t_ar
                @warn "channelview1: failed on $f with eltype $T, time ratio $(t_cv/t_ar), tol $tol"
            end

            t_ar, t_cv = testf(f, asub, vchansub, 30)
            tol = haskey(chanvtol, f) ? chanvtol[f] : chanvdefault
            if t_cv >= tol*t_ar
                @warn "channelview2: failed on $f with eltype $T, time ratio $(t_cv/t_ar), tol $tol"
            end

            # colorview
            t_ar, t_cv = testf(f, c, vcol, 30)
            tol = haskey(colvtol, f) ? colvtol[f] : colvdefault
            if t_cv >= tol*t_ar
                @warn "colorview1: failed on $f with eltype $T, time ratio $(t_cv/t_ar), tol $tol"
            end

            t_ar, t_cv = testf(f, csub, vcolsub, 30)
            tol = haskey(colvtol, f) ? colvtol[f] : colvdefault
            if t_cv >= tol*t_ar
                @warn "colorview2: failed on $f with eltype $T, time ratio $(t_cv/t_ar), tol $tol"
            end
        end
    end
end
# end
