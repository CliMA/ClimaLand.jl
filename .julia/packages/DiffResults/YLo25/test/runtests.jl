using DiffResults, StaticArrays, Test

using DiffResults: DiffResult, GradientResult, JacobianResult, HessianResult,
                   value, value!, derivative, derivative!, gradient, gradient!,
                   jacobian, jacobian!, hessian, hessian!

@testset "DiffResult" begin

    k = 4
    n0, n1, n2 = rand(), rand(), rand()
    x0, x1, x2 = rand(k), rand(k, k), rand(k, k, k)
    s0, s1, s2 = SVector{k}(rand(k)), SMatrix{k,k}(rand(k, k)), SArray{Tuple{k,k,k}}(rand(k, k, k))
    rn = DiffResult(n0, n1, n2)
    rx = DiffResult(x0, x1, x2)
    rs = DiffResult(s0, s1, s2)
    rsmix = DiffResult(n0, s0, s1)

    issimilar(x, y) = typeof(x) == typeof(y) && size(x) == size(y)
    issimilar(x::DiffResult, y::DiffResult) = issimilar(value(x), value(y)) && all(issimilar, zip(x.derivs, y.derivs))
    issimilar(t::Tuple) = issimilar(t...)

    @test rn === DiffResult(n0, n1, n2)
    @test rx == DiffResult(x0, x1, x2)
    @test rs === DiffResult(s0, s1, s2)
    @test rsmix === DiffResult(n0, s0, s1)

    @test issimilar(GradientResult(x0), DiffResult(first(x0), x0))
    @test issimilar(JacobianResult(x0), DiffResult(x0, similar(x0, k, k)))
    @test issimilar(JacobianResult(similar(x0, k + 1), x0), DiffResult(similar(x0, k + 1), similar(x0, k + 1, k)))
    @test issimilar(HessianResult(x0), DiffResult(first(x0), x0, similar(x0, k, k)))

    @test GradientResult(s0) === DiffResult(first(s0), s0)
    @test JacobianResult(s0) === DiffResult(s0, zeros(SMatrix{k,k,Float64}))
    @test JacobianResult(SVector{k+1}(vcat(s0, 0.0)), s0) === DiffResult(SVector{k+1}(vcat(s0, 0.0)), zeros(SMatrix{k+1,k,Float64}))
    @test HessianResult(s0) === DiffResult(first(s0), s0, zeros(SMatrix{k,k,Float64}))

    @test eltype(rn) === typeof(n0)
    @test eltype(rx) === eltype(x0)
    @test eltype(rs) === eltype(s0)

    rn_copy = copy(rn)
    @test rn == rn_copy
    @test rn === rn_copy

    rx_copy = copy(rx)
    @test rx == rx_copy
    @test rx !== rx_copy

    rs_copy = copy(rs)
    @test rs == rs_copy
    @test rs === rs_copy

    rsmix_copy = copy(rsmix)
    @test rsmix == rsmix_copy
    @test rsmix === rsmix_copy

    @testset "value/value!" begin

        @test value(rn) === n0
        @test value(rx) === x0
        @test value(rs) === s0
        @test value(rsmix) === n0

        rn = value!(rn, n1)
        @test value(rn) === n1
        rn = value!(rn, n0)

        x0_new, x0_copy = rand(k), copy(x0)
        rx = value!(rx, x0_new)
        @test value(rx) === x0 == x0_new
        rx = value!(rx, x0_copy)

        s0_new = rand(k)
        rs = value!(rs, s0_new)
        @test value(rs) == s0_new
        @test typeof(value(rs)) === typeof(s0)
        rs = value!(rs, s0)

        rsmix = value!(rsmix, n1)
        @test value(rsmix) === n1
        rsmix = value!(rsmix, n0)

        rn = value!(exp, rn, n1)
        @test value(rn) === exp(n1)
        rn = value!(rn, n0)

        x0_new, x0_copy = rand(k), copy(x0)
        rx = value!(exp, rx, x0_new)
        @test value(rx) === x0 == exp.(x0_new)
        rx = value!(rx, x0_copy)

        s0_new = rand(k)
        rs = value!(exp, rs, s0_new)
        @test value(rs) == exp.(s0_new)
        @test typeof(value(rs)) === typeof(s0)
        rs = value!(rs, s0)

        rsmix = value!(exp, rsmix, n1)
        @test value(rsmix) === exp(n1)
        rsmix = value!(rsmix, n0)

        ksqrt = Int(sqrt(k))
        T = typeof(SMatrix{ksqrt,ksqrt}(rand(ksqrt,ksqrt)))
        rs_new = value!(rs, convert(T, value(rs)))
        @test rs_new === rs
    end

    @testset "derivative/derivative!" begin

        @test derivative(rn) === n1
        @test derivative(rn, Val{2}) === n2

        @test derivative(rx) === x1
        @test derivative(rx, Val{2}) === x2

        @test derivative(rs) === s1
        @test derivative(rs, Val{2}) === s2

        @test derivative(rsmix) === s0
        @test derivative(rsmix, Val{2}) === s1

        rn = derivative!(rn, n0)
        @test derivative(rn) === n0
        rn = derivative!(rn, n1)

        x1_new, x1_copy = rand(k, k), copy(x1)
        rx = derivative!(rx, x1_new)
        @test derivative(rx) === x1 == x1_new
        rx = derivative!(rx, x1_copy)

        s1_new = rand(k, k)
        rs = derivative!(rs, s1_new)
        @test derivative(rs) == s1_new
        @test typeof(derivative(rs)) === typeof(s1)
        rs = derivative!(rs, s1)

        s0_new = rand(k)
        rsmix = derivative!(rsmix, s0_new)
        @test derivative(rsmix) == s0_new
        @test typeof(derivative(rsmix)) === typeof(s0)
        rsmix = derivative!(rsmix, s0)

        rn = derivative!(rn, n1, Val{2})
        @test derivative(rn, Val{2}) === n1
        rn = derivative!(rn, n2, Val{2})

        x2_new, x2_copy = rand(k, k, k), copy(x2)
        rx = derivative!(rx, x2_new, Val{2})
        @test derivative(rx, Val{2}) === x2 == x2_new
        rx = derivative!(rx, x2_copy, Val{2})

        s2_new = rand(k, k, k)
        rs = derivative!(rs, s2_new, Val{2})
        @test derivative(rs, Val{2}) == s2_new
        @test typeof(derivative(rs, Val{2})) === typeof(s2)
        rs = derivative!(rs, s2, Val{2})

        s1_new = rand(k, k)
        rsmix = derivative!(rsmix, s1_new, Val{2})
        @test derivative(rsmix, Val{2}) == s1_new
        @test typeof(derivative(rsmix, Val{2})) === typeof(s1)
        rsmix = derivative!(rsmix, s1, Val{2})

        rn = derivative!(exp, rn, n0)
        @test derivative(rn) === exp(n0)
        rn = derivative!(rn, n1)

        x1_new, x1_copy = rand(k, k), copy(x1)
        rx = derivative!(exp, rx, x1_new)
        @test derivative(rx) === x1 == exp.(x1_new)
        rx = derivative!(exp, rx, x1_copy)

        s1_new = rand(k, k)
        rs = derivative!(exp, rs, s1_new)
        @test derivative(rs) == exp.(s1_new)
        @test typeof(derivative(rs)) === typeof(s1)
        rs = derivative!(exp, rs, s1)

        s0_new = rand(k)
        rsmix = derivative!(exp, rsmix, s0_new)
        @test derivative(rsmix) == exp.(s0_new)
        @test typeof(derivative(rsmix)) === typeof(s0)
        rsmix = derivative!(exp, rsmix, s0)

        rn = derivative!(exp, rn, n1, Val{2})
        @test derivative(rn, Val{2}) === exp(n1)
        rn = derivative!(rn, n2, Val{2})

        x2_new, x2_copy = rand(k, k, k), copy(x2)
        rx = derivative!(exp, rx, x2_new, Val{2})
        @test derivative(rx, Val{2}) === x2 == exp.(x2_new)
        rx = derivative!(exp, rx, x2_copy, Val{2})

        s2_new = rand(k, k, k)
        rs = derivative!(exp, rs, s2_new, Val{2})
        @test derivative(rs, Val{2}) == exp.(s2_new)
        @test typeof(derivative(rs, Val{2})) === typeof(s2)
        rs = derivative!(exp, rs, s2, Val{2})

        s1_new = rand(k, k)
        rsmix = derivative!(exp, rsmix, s1_new, Val{2})
        @test derivative(rsmix, Val{2}) == exp.(s1_new)
        @test typeof(derivative(rsmix, Val{2})) === typeof(s1)
        rsmix = derivative!(exp, rsmix, s1, Val{2})
    end

    @testset "gradient/gradient!" begin

        x1_new, x1_copy = rand(k, k), copy(x1)
        rx = gradient!(rx, x1_new)
        @test gradient(rx) === x1 == x1_new
        rx = gradient!(rx, x1_copy)

        s1_new = rand(k, k)
        rs = gradient!(rs, s1_new)
        @test gradient(rs) == s1_new
        @test typeof(gradient(rs)) === typeof(s1)
        rs = gradient!(rs, s1)

        x1_new, x1_copy = rand(k, k), copy(x1)
        rx = gradient!(exp, rx, x1_new)
        @test gradient(rx) === x1 == exp.(x1_new)
        rx = gradient!(exp, rx, x1_copy)

        s0_new = rand(k)
        rsmix = gradient!(exp, rsmix, s0_new)
        @test gradient(rsmix) == exp.(s0_new)
        @test typeof(gradient(rsmix)) === typeof(s0)
        rsmix = gradient!(exp, rsmix, s0)

        T = typeof(SVector{k*k}(rand(k*k)))
        rs_new = gradient!(rs, convert(T, gradient(rs)))
        @test rs_new === rs
    end

    @testset "jacobian/jacobian!"  begin

        x1_new, x1_copy = rand(k, k), copy(x1)
        rx = jacobian!(rx, x1_new)
        @test jacobian(rx) === x1 == x1_new
        rx = jacobian!(rx, x1_copy)

        s1_new = rand(k, k)
        rs = jacobian!(rs, s1_new)
        @test jacobian(rs) == s1_new
        @test typeof(jacobian(rs)) === typeof(s1)
        rs = jacobian!(rs, s1)

        x1_new, x1_copy = rand(k, k), copy(x1)
        rx = jacobian!(exp, rx, x1_new)
        @test jacobian(rx) === x1 == exp.(x1_new)
        rx = jacobian!(exp, rx, x1_copy)

        s0_new = rand(k)
        rsmix = jacobian!(exp, rsmix, s0_new)
        @test jacobian(rsmix) == exp.(s0_new)
        @test typeof(jacobian(rsmix)) === typeof(s0)
        rsmix = jacobian!(exp, rsmix, s0)

        T = typeof(SVector{k*k}(rand(k*k)))
        rs_new = jacobian!(rs, convert(T, jacobian(rs)))
        @test rs_new === rs
    end

    @testset "hessian/hessian!" begin

        x2_new, x2_copy = rand(k, k, k), copy(x2)
        rx = hessian!(rx, x2_new)
        @test hessian(rx) === x2 == x2_new
        rx = hessian!(rx, x2_copy)

        s2_new = rand(k, k, k)
        rs = hessian!(rs, s2_new)
        @test hessian(rs) == s2_new
        @test typeof(hessian(rs)) === typeof(s2)
        rs = hessian!(rs, s2)

        x2_new, x2_copy = rand(k, k, k), copy(x2)
        rx = hessian!(exp, rx, x2_new)
        @test hessian(rx) === x2 == exp.(x2_new)
        rx = hessian!(exp, rx, x2_copy)

        s1_new = rand(k, k)
        rsmix = hessian!(exp, rsmix, s1_new)
        @test hessian(rsmix) == exp.(s1_new)
        @test typeof(hessian(rsmix)) === typeof(s1)
        rsmix = hessian!(exp, rsmix, s1)

        T = typeof(SVector{k*k*k}(rand(k*k*k)))
        rs_new = hessian!(rs, convert(T, hessian(rs)))
        @test rs_new === rs

        @test size(gradient(HessianResult(x0))) == size(x0)
        @test size(gradient(HessianResult(x1))) == size(x1)
        @test size(gradient(HessianResult(x2))) == size(x2)

        @test HessianResult(Float32.(x0)).derivs[1] isa Vector{Float32}
    end
end
