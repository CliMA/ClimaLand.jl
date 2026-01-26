@testset "objective types" begin
    x_seed = [0.0, 0.0]
    g_seed = [0.0, 0.0]
    h_seed = [0.0 0.0; 0.0 0.0]
    f_x_seed = 8157.682077608529

    nd = NonDifferentiable(exponential, x_seed)
    @test nd.f == exponential
    @test value(nd) == 0.0
    @test nd.f_calls == [0]

    od = OnceDifferentiable(exponential, exponential_gradient!, nothing, x_seed, 0.0, g_seed)
    @test od.f == exponential
    @test od.df == exponential_gradient!
    @test value(od) == 0.0
    @test od.f_calls == [0]
    @test od.df_calls == [0]
    od.x_df .= x_seed
    gold = copy(od.DF)
    xnew = rand(eltype(x_seed), size(x_seed))
    gnew = gradient(od, xnew)
    @test od.x_df == x_seed
    @test od.DF == gold
    @test gnew == gradient(od, xnew)

    td = TwiceDifferentiable(exponential, exponential_gradient!, nothing, exponential_hessian!, x_seed, 0.0, g_seed, h_seed)
    @test td.f == exponential
    @test td.df == exponential_gradient!
    @test value(td) == 0.0
    @test td.f_calls == [0]
    @test td.df_calls == [0]
    @test td.h_calls == [0]

    @testset "no fg!" begin
        Random.seed!(324)
        od = OnceDifferentiable(exponential, exponential_gradient!, x_seed, 0.0, g_seed)
        xrand = rand(2)
        value_gradient!(od, xrand)
        ndod = NonDifferentiable(od, xrand)
        @test value(ndod, xrand) === value(od, xrand)
        fcache = value(od)
        gcache = copy(gradient(od))
        value_gradient!(od, zeros(2))
        gradient!(od, xrand)
        @test value(od, zeros(2)) == od.F
        @test value(od, zeros(2)) == value(od)
        @test gradient(od) == gcache

        od = OnceDifferentiable(exponential, exponential_gradient!, x_seed)
        xrand = rand(2)
        value_gradient!(od, xrand)
        fcache = value(od)
        gcache = copy(gradient(od))
        value_gradient!(od, zeros(2))
        gradient!(od, xrand)
        @test value(od, zeros(2)) == od.F
        @test value(od, zeros(2)) == value(od)
        @test gradient(od) == gcache

        td = TwiceDifferentiable(exponential, exponential_gradient!, exponential_hessian!, x_seed, 0.0, g_seed)
        xrand = rand(2)
        value_gradient!(td, xrand)
        ndtd = NonDifferentiable(td, xrand)
        @test value(ndtd, xrand) === value(td, xrand)

        fcache = value(td)
        gcache = copy(gradient(td))
        value_gradient!(td, zeros(2))
        gradient!(td, xrand)
        @test value(td, zeros(2)) == td.F
        @test value(td, zeros(2)) == value(td)
        @test gradient(td) == gcache
    end
end
