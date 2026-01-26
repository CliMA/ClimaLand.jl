@testset "incomplete objectives" begin
    import NLSolversBase: df, fdf, make_f, make_df, make_fdf
    function f(x)
        sum(x->x^2,x)
    end

    g!(G, x) = copyto!(G, 2 .* x)
    g(x) = 2 .* x

    h!(H, x) = copy!(H, Diagonal(fill(2, length(n))))
    h(x) = Diagonal(fill(2, length(n)))

    function fg!(G, x)
        copyto!(G, 2 .* x)
        sum(x->x^2,x)
    end
    function just_fg!(F, G, x)
        !(G == nothing) && copyto!(G, 2 .* x)
        !(F == nothing) && sum(x->x^2,x)
    end
    fg(x) = f(x), g(x)
    function just_fgh!(F, G, H, x)
        !(H == nothing) && copy!(H, Diagonal(fill(2, length(n))))
        !(G == nothing) && copyto!(G, 2 .* x)
        !(F == nothing) && return sum(x->x^2,x)
    end
    fgh(x) = f(x), g(x), h(x)
    fdf!_real = only_fg!(just_fg!)
    fdf_real = only_fg(fg)

    function just_hv!(Hv, x, v)
        copyto!(Hv, 2 .* v)
    end
    function just_fghv!(F, G, Hv, x, v)
        G  === nothing || copyto!(G, 2 .* x)
        Hv === nothing || copyto!(Hv, 2 .* v)
        F  === nothing || return sum(x->x^2, x)
    end

    df_fdf_real = only_g_and_fg(g, fg)
    Random.seed!(3259)
    x = rand(10)
    G_cache = similar(x)
    G_cache2 = similar(G_cache)

    @test df(fdf!_real) === nothing
    @test df(fdf_real) === nothing
    @test df(df_fdf_real) === g
    @test df(df_fdf_real)(x) == g(x)

    @test fdf(fdf!_real) === just_fg!
    @test fdf(fdf_real) === fg
    @test df(df_fdf_real) == g
    @test fdf(df_fdf_real) == fg
    @test df(df_fdf_real)(x) == g(x)
    @test fdf(df_fdf_real)(x) == fg(x)

    for FDF in (fdf_real, fdf!_real)
        @test make_f(FDF, x, x[1])(x) == f(x)
        make_df(FDF, x, x[1])(G_cache, x)
        g!(G_cache2, x)
        @test G_cache == G_cache2
        f1 = make_fdf(FDF, x, x[1])(G_cache, x.*2)
        f2 = fg!(G_cache2, x.*2)
        @test G_cache == G_cache2
        @test f1 == f2
    end

    nd_fg = NonDifferentiable(only_fg(fg), x)
    nd_fg! = NonDifferentiable(only_fg!(just_fg!), x)
    for ND in (nd_fg, nd_fg!)
        value!(ND, x)
        value(ND) == f(x)
    end
    od_fg = OnceDifferentiable(only_fg(fg), x)
    od_fg! = OnceDifferentiable(only_fg!(just_fg!), x)
    od_fgh! = TwiceDifferentiable(only_fgh!(just_fgh!), x)
    _F = zero(eltype(x))
    od_fgh! = TwiceDifferentiable(only_fgh!(just_fgh!), x, _F)
    od_fgh! = TwiceDifferentiable(only_fgh!(just_fgh!), x, _F, similar(x))
    od_fgh! = TwiceDifferentiable(only_fgh!(just_fgh!), x, _F, similar(x), NLSolversBase.alloc_DF(x, _F))
#    od_fgh = TwiceDifferentiable(only_fgh(fgh), x)
    for OD in (od_fg, od_fg!, od_fgh!)#, od_fgh)
        value!(OD, x)
        @test value(OD) == f(x)
        gradient!(OD, x)
        @test gradient(OD) == g(x)
        @test gradient(OD, x) == g(x)
        value_gradient!(OD, 2 .* x)
        @test value(OD) == f(2 .* x)
        @test gradient(OD) == g(2 .* x)
    end

    # Incomplete TwiceDifferentiableHv
    v = randn(10)
    od_fg_and_hv = TwiceDifferentiableHV(only_fg_and_hv!(just_fg!, just_hv!), x)
    od_fghv      = TwiceDifferentiableHV(only_fghv!(just_fghv!), x)
    ndtdhv = NonDifferentiable(od_fghv, v)
    @test value(ndtdhv, v) === value(od_fghv, v)

    for OD in (od_fg_and_hv, od_fghv)
        gradient!(OD, x)
        @test gradient(OD) == g(x)
        hv_product!(OD, x, v)
        @test OD.Hv == 2v
        OD.f(x) == f(x)
        _g = similar(x)
        OD.fdf(_g, x)
        @test _g == g(x)
        @test OD.fdf(_g, x) == f(x)
    end
end
@testset "incomplete objectives vectors" begin
    function tf(x)
        x.^2
    end
    function tf(F, x)
        copyto!(F, tf(x))
    end

    tj!(J, x) = copyto!(J, Matrix(Diagonal(x)))
    tj(x) = Matrix(Diagonal(x))

    function tfj!(F, J, x)
        copyto!(J, Matrix(Diagonal(x)))
        copyto!(F, tf(x))
    end
    function just_tfj!(F, J, x)
        !(J == nothing) && copyto!(J, Matrix(Diagonal(x)))
        !(F == nothing) && copyto!(F, tf(x))
    end
    tfj(x) = tf(x), tj(x)

    fdf!_real = only_fj!(just_tfj!)
    fdf_real = only_fj(tfj)

    df_fdf_real = only_j_and_fj(tj, tfj)
    Random.seed!(3259)
    x = rand(10)
    J_cache = similar(Matrix(Diagonal(x)))
    J_cache2 = similar(Matrix(Diagonal(x)))
    F_cache = similar(x)
    F_cache2 = similar(x)

    @test df(fdf!_real) === nothing
    @test df(fdf_real) === nothing
    @test df(df_fdf_real) === tj
    @test df(df_fdf_real)(x) == tj(x)

    @test fdf(fdf!_real) === just_tfj!
    @test fdf(fdf_real) === tfj
    @test df(df_fdf_real) == tj
    @test fdf(df_fdf_real) == tfj
    @test df(df_fdf_real)(x) == tj(x)
    @test fdf(df_fdf_real)(x) == tfj(x)

    for FDF in (fdf_real, fdf!_real)
        @test make_f(FDF, x, x)(F_cache, x) == tf(x)
        make_df(FDF, x, x)(J_cache, x)
        tj!(J_cache2, x)
        @test J_cache == J_cache2
        make_fdf(FDF, x, x)(F_cache, J_cache, x.*2)
        tfj!(F_cache2, J_cache2, x.*2)
        @test F_cache == F_cache2
        @test J_cache == J_cache2
    end

    nd_fj = NonDifferentiable(only_fj(tfj), x, x)
    nd_fj! = NonDifferentiable(only_fj!(just_tfj!), x, x)
    for ND in (nd_fj, nd_fj!)
        value!(ND, x)
        value(ND) == tf(x)
    end
    od_fj = OnceDifferentiable(only_fj(tfj), x, x)
    od_fj! = OnceDifferentiable(only_fj!(just_tfj!), x, x)
    for OD in (od_fj, od_fj!)
        value!(OD, x)
        @test value(OD) == tf(x)
        jacobian!(OD, x)
        @test jacobian(OD) == tj(x)
        @test jacobian(OD, x) == tj(x)
        value_jacobian!(OD, 2 .* x)
        @test value(OD) == tf(2 .* x)
        @test jacobian(OD) == tj(2 .* x)
    end

end

@testset "https://github.com/JuliaNLSolvers/Optim.jl/issues/718" begin
    f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
    function g!(G, x)
      G[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
      G[2] = 200.0 * (x[2] - x[1]^2)
    end
    function h!(H, x)
      H[1, 1] = 2.0 - 400.0 * x[2] + 1200.0 * x[1]^2
      H[1, 2] = -400.0 * x[1]
      H[2, 1] = -400.0 * x[1]
      H[2, 2] = 200.0
    end

    function fg!(F,G,x)
      G == nothing || g!(G,x)
      F == nothing || return f(x)
      nothing
    end
    function fgh!(F,G,H,x)
      G == nothing || g!(G,x)
      H == nothing || h!(H,x)
      F == nothing || return f(x)
      nothing
    end
    
    gx = [0.0,0.0]
    x=[0.0,0.0]

    @test NLSolversBase.make_f(only_fgh!(fgh!),[0.0,0.0],0.0)(x) == 1.0
    @test NLSolversBase.make_df(only_fgh!(fgh!),[0.0,0.0],0.0)(gx, x) == nothing
    @test gx == [-2.0, 0.0]

    gx = [0.0,0.0]
    @test NLSolversBase.make_fdf(only_fgh!(fgh!),[0.0,0.0],0.0)(gx, x) == 1.0
    @test gx == [-2.0, 0.0]

end
