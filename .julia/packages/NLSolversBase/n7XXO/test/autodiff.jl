@testset "autodiff" begin
    Random.seed!(0)

    # Should throw, as :wah is not a proper autodiff choice
    @test_throws ErrorException OnceDifferentiable(x->x, rand(10); autodiff=:wah)
    @test_throws ErrorException OnceDifferentiable(x->x, rand(10), 0.0; autodiff=:wah)
    @test_throws ErrorException TwiceDifferentiable(x->x, rand(10); autodiff=:wah)
    @test_throws ErrorException TwiceDifferentiable(x->x, rand(10), 0.0; autodiff=:wah)
    #@test_throws ErrorException TwiceDifferentiable(x->x, rand(10), 0.0, rand(10); autodiff=:wah)
    @test_throws ErrorException TwiceDifferentiable(x->x, x->x, rand(10); autodiff=:wah)
    @test_throws ErrorException TwiceDifferentiable(x->x, x->x, rand(10), 0.0; autodiff=:wah)
    #@test_throws ErrorException TwiceDifferentiable(x->x, x->x, rand(10), 0.0, rand(10); autodiff=:wah)

    for T in (OnceDifferentiable, TwiceDifferentiable)
        odad1 = T(x->5.0, rand(1); autodiff = :finite)
        odad2 = T(x->5.0, rand(1); autodiff = :forward)
        gradient!(odad1, rand(1))
        gradient!(odad2, rand(1))
        #    odad3 = T(x->5., rand(1); autodiff = :reverse)
        @test gradient(odad1) == [0.0]
        @test gradient(odad2) == [0.0]
        #    @test odad3.g == [0.0]
    end

    for a in (1.0, 5.0)
        xa = rand(1)
        odad1 = OnceDifferentiable(x->a*x[1], xa; autodiff = :finite)
        odad2 = OnceDifferentiable(x->a*x[1], xa; autodiff = :forward)
    #    odad3 = OnceDifferentiable(x->a*x[1], xa; autodiff = :reverse)
        gradient!(odad1, xa)
        gradient!(odad2, xa)
        @test gradient(odad1) ≈ [a]
        @test gradient(odad2) == [a]
    #    @test odad3.g == [a]
    end
    for a in (1.0, 5.0)
        xa = rand(1)
        odad1 = OnceDifferentiable(x->a*x[1]^2, xa; autodiff = :finite)
        odad2 = OnceDifferentiable(x->a*x[1]^2, xa; autodiff = :forward)
    #    odad3 = OnceDifferentiable(x->a*x[1]^2, xa; autodiff = :reverse)
        gradient!(odad1, xa)
        gradient!(odad2, xa)
        @test gradient(odad1) ≈ 2.0*a*xa
        @test gradient(odad2) == 2.0*a*xa
    #    @test odad3.g == 2.0*a*xa
    end
    nx = 2
    x = rand(nx)
    f(x) = sum(x.^3)
    fx = f(x)
    g(G, x) = copyto!(G, 3 .* x.^2)
    gx = g(NLSolversBase.alloc_DF(x, 0.0), x)
    h(H, x) = copyto!(H, Diagonal(6 .* x))
    hx = h(fill(0.0, nx, nx), x)
    @testset for dtype in (OnceDifferentiable, TwiceDifferentiable)
        @testset for autodiff in (:finite, :forward, AutoForwardDiff())
            # :forward should be exact, but :finite will not be
            differentiable = dtype(f, copy(x); autodiff = autodiff)
            value!(differentiable, copy(x))
            @test isapprox(value(differentiable), fx)
            clear!(differentiable)
            value_gradient!(differentiable, x)
            @test isapprox(value(differentiable), fx)
            @test isapprox(gradient(differentiable), gx)
            clear!(differentiable)
            gradient!(differentiable, x)
            @test isapprox(gradient(differentiable), gx)
            clear!(differentiable)
            if dtype == TwiceDifferentiable
                hessian!(differentiable, x)
                if autodiff == :finite
                    # we have to increase the tolerance here, as the hessian is
                    # not very accurate
                    @test isapprox(hessian(differentiable), hx; atol = 1e-6)
                else
                    @test hessian(differentiable) == hx
                end
            end
        end
    end
    @testset for autodiff in (:finite, :forward, AutoForwardDiff())
        td = TwiceDifferentiable(x->sum(x), (G, x)->copyto!(G, fill!(copy(x),1)), copy(x); autodiff = autodiff)
        value(td)
        value!(td, x)
        value_gradient!(td, x)
        gradient!(td, x)
        hessian!(td, x)
    end
    @testset for autodiff in (:finite, :forward, AutoForwardDiff())
        for nd = (NonDifferentiable(x->sum(x), copy(x)), NonDifferentiable(x->sum(x), copy(x), 0.0))
            td = TwiceDifferentiable(nd; autodiff = autodiff)
            value(td)
            value!(td, x)
            value_gradient!(td, x)
            gradient!(td, x)
            hessian!(td, x)
        end
        for od = (OnceDifferentiable(x->sum(x), (G, x)->copyto!(G, fill!(copy(x),1)), copy(x)), OnceDifferentiable(x->sum(x), (G, x)->copyto!(G, fill!(copy(x),1)), copy(x), 0.0))
            td = TwiceDifferentiable(od; autodiff = autodiff)
            value(td)
            value!(td, x)
            value_gradient!(td, x)
            gradient!(td, x)
            hessian!(td, x)
        end
    end
    @testset "autodiff ℝᴺ → ℝᴺ" begin
        function f!(F, x)
            F[1] = 1.0 - x[1]
            F[2] = 10.0*(x[2] - x[1]^2)
            F
        end
        function j!(J, x)
            J[1,1] = -1.0
            J[1,2] = 0.0
            J[2,1] = -20.0*x[1]
            J[2,2] = 10.0
            J
        end
        # Some random x
        x = rand(2)
        # A type and shape-correct F vector
        F = similar(x)
        # A type and shape-correct J matrix
        J = NLSolversBase.alloc_DF(x, F)
        # Default, should be :central
        od_fd = OnceDifferentiable(f!, x, F)
        value_jacobian!(od_fd, x)
        @test value(od_fd) == f!(F, x)
        # Can't test equality here
        @test jacobian(od_fd) ≈ j!(J, x)
        # Specifically :central
        od_fd_2 = OnceDifferentiable(f!, x, F, :central)
        @test value!(od_fd_2, x) == f!(F, x)
        # Can't test equality here
        @test jacobian!(od_fd_2, x) ≈ j!(J, x)
        # Test that they're identical -> they used the same scheme
        @test jacobian(od_fd) == jacobian(od_fd_2)

        od_ad = OnceDifferentiable(f!, x, F, :forward)
        @test value!(od_ad, x) == f!(F, x)
        # Can't test equality here
        @test jacobian!(od_ad, x) ≈ j!(J, x)

        od_ad_2 = OnceDifferentiable(f!, x, F, :forward)
        @test value!(od_ad_2, x) == f!(F, x)
        # Can't test equality here
        @test jacobian!(od_ad_2, x) ≈ j!(J, x)
        # Test that they're identical -> they used the same scheme
        @test jacobian(od_ad) == jacobian(od_ad_2)
        @testset "error handling" begin
            x = rand(2)
            F = similar(x)
            # Wrong symbol
            @test_throws ErrorException OnceDifferentiable(f!, x, F, :foo)
            # Wrong bool
            @test_throws ErrorException OnceDifferentiable(f!, x, F, false)
        end
    end
end
@testset "value/grad" begin
    a = 3.0
    x_seed = rand(1)
    odad1 = OnceDifferentiable(x->a*x[1]^2, x_seed)
    value_gradient!(odad1, x_seed)
    @test gradient(odad1) ≈ 2 .* a .* (x_seed)
    @testset "call counters" begin
        @test f_calls(odad1) == 1
        @test g_calls(odad1) == 1
        @test h_calls(odad1) == 0
        value_gradient!(odad1, x_seed .+ 1.0)
        @test f_calls(odad1) == 2
        @test g_calls(odad1) == 2
        @test h_calls(odad1) == 0
    end
    @test gradient(odad1) ≈ 2 .* a .* (x_seed .+ 1.0)
end

@testset "residual function" begin
    function f!(F, x)
        F[1] = 2x[1]+x[2]
        F[2] = x[1]+x[2]^2
        F[3] = x[1]^2+x[2]^2
        F
    end

    function j!(J, x)
        J[1, 1] = 2
        J[1, 2] = 1
        J[2, 1] = 1
        J[2, 2] = 2 * x[2]
        J[3, 1] = 2 * x[1]
        J[3, 2] = 2 * x[2]
        J
    end

    F = zeros(3)
    J = zeros(3, 2)
    x_init = [1., 1.]
    od = OnceDifferentiable(f!, x_init,  F)
    value_jacobian!(od, x_init)
    @test length(value(od)) == length(F)
    @test value(od) ≈ f!(F, x_init)
    @test jacobian(od) ≈ j!(J, x_init)
end
