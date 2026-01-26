@testset "interface" begin
    @testset "single-valued" begin
        x_seed = [0.0, 0.0]
        g_seed = [0.0, 0.0]
        h_seed = [0.0 0.0; 0.0 0.0]
        f_x_seed = 8157.682077608529
        g_x_seed = [-218.39260013257694, -48618.50356545231]
        h_x_seed = [982.7667005965963 0.0; 0.0 307917.1892478646]

        x_alt = [1.0, 1.0]
        f_x_alt = 57.316431861603284
        g_x_alt = [-5.43656365691809, -218.39260013257694]
        h_x_alt = [16.30969097075427 0.; 0. 982.7667005965963]

        # Construct instances
        nd = NonDifferentiable(exponential, x_seed)
        od = OnceDifferentiable(exponential, exponential_gradient!, exponential_value_gradient!, x_seed, 0.0, g_seed)
        td = TwiceDifferentiable(exponential, exponential_gradient!, exponential_value_gradient!, exponential_hessian!, x_seed, 0.0, g_seed, h_seed)

        # Force evaluation
        value!!(nd, x_seed)
        value_gradient!!(od, x_seed)
        value_gradient!!(td, x_seed)
        hessian!!(td, x_seed)

        # Test that values are the same, and that values match those
        # calculated by the value(obj, x) methods
        @test value(nd) == value(od) == value(td) == f_x_seed
        @test value(nd, x_seed) == value(od, x_seed) == value(td, x_seed)

        # Test that the gradients match the intended values
        @test gradient(od) == gradient(td) == g_x_seed
        # Test that the Hessian matches the intended value
        @test hessian(td) == h_x_seed
        # Test hv_product! for TwiceDifferentiable
        v = [0.111, -1234]
        @test hv_product!(td, x_seed, v) == h_x_seed * v

        # Test that the call counters got incremented
        @test nd.f_calls == od.f_calls == td.f_calls == [2]
        @test od.df_calls == td.df_calls == [1]
        @test td.h_calls == [1]

        # Test that the call counters do not get incremented
        # with single-"bang" methods...
        value!(nd, x_seed)
        value_gradient!(od, x_seed)
        value_gradient!(td, x_seed)
        hessian!(td, x_seed)

        @test nd.f_calls == od.f_calls == td.f_calls == [2]
        @test od.df_calls == td.df_calls == [1]
        @test td.h_calls == [1]

        # ... and that they do with double-"bang" methods
        value!!(nd, x_seed)
        value_gradient!!(od, x_seed)
        value_gradient!!(td, x_seed)
        hessian!!(td, x_seed)

        @test nd.f_calls == od.f_calls == td.f_calls == [3]
        @test od.df_calls == td.df_calls == [2]
        @test td.h_calls == [2]

        # Test that gradient doesn't work for NonDifferentiable, but does otherwise
        @test_throws ErrorException gradient!(nd, x_alt)
        gradient!(od, x_alt)
        gradient!(td, x_alt)

        @test value(nd) == value(od) == value(td) == f_x_seed
        @test gradient(td) == g_x_alt
        @test gradient(td) == [gradient(td, i) for i = 1:length(x_seed)]
        @test hessian(td) == h_x_seed
        @test nd.f_calls == od.f_calls == td.f_calls == [3]
        @test od.df_calls == td.df_calls == [3]
        @test td.h_calls == [2]

        @test_throws ErrorException hessian!(nd, x_alt)
        @test_throws ErrorException hessian!(od, x_alt)
        hessian!(td, x_alt)

        @test value(nd) == value(od) == value(td) == f_x_seed
        @test gradient(td) == g_x_alt
        @test hessian(td) == h_x_alt
        @test nd.f_calls == od.f_calls == td.f_calls == [3]
        @test od.df_calls == td.df_calls == [3]
        @test td.h_calls == [3]

        value!(nd, x_alt)
        value!(od, x_alt)
        value!(td, x_alt)
        @test value(nd) == value(od) == value(td) == f_x_alt
        @test gradient(td) == g_x_alt
        @test hessian(td) == h_x_alt
        @test nd.f_calls == od.f_calls == td.f_calls == [4]
        @test od.df_calls == td.df_calls == [3]
        @test td.h_calls == [3]

        @test_throws ErrorException value_gradient!(nd, x_seed)
        value_gradient!(od, x_seed)
        value_gradient!(td, x_seed)
        @test value(od) == value(td) == f_x_seed
        # change x_f manually to test branch
        od.x_f .*= 0
        td.x_f .*= 0
        value_gradient!(od, x_seed)
        value_gradient!(td, x_seed)
        @test value(od) == value(td) == f_x_seed
        # change x_df manually to test branch
        # only df is meant to be re-calculated as
        # d.x_f == x_seed
        # and that the hessian counter doesn't increment
        od.x_df .*= 0
        td.x_df .*= 0
        value_gradient!(od, x_seed)
        value_gradient!(td, x_seed)
        @test value(od) == value(td) == f_x_seed
        @test gradient(td) == g_x_seed
        @test hessian(td) == h_x_alt
        @test od.f_calls == td.f_calls == [5]
        @test od.df_calls == td.df_calls == [4]
        @test td.h_calls == [3]

        # test the non-mutating gradient() function
        gradient!(od, fill(1.0, 2))
        gradient!(td, fill(1.0, 2))
        od_df_old = copy(gradient(od))
        td_df_old = copy(gradient(td))
        gradient!(od, fill(-1.0, 2))
        gradient!(td, fill(-1.0, 2))
        @test od_df_old == gradient(od, fill(1.0, 2))
        @test td_df_old == gradient(td, fill(1.0, 2))
        # Mutate obj.DF directly to check that the gradient
        # is not recalculated if the same x is reused.
        # We obviously cannot just rerun for the same x twice
        # as we would then not be able to tell if it was just
        # calculated again, or the cache was simply returned
        # as intended.
        od.DF .= 0.0
        td.DF .= 0.0
        @test gradient(od) == zeros(2)
        @test gradient(td) == zeros(2)

        clear!(nd)
        clear!(od)
        clear!(td)
        @test nd.f_calls == [0,]
        @test od.f_calls == [0,]
        @test td.f_calls == [0,]
        @test od.df_calls == [0,]
        @test td.df_calls == [0,]
        @test td.h_calls == [0,]

        @testset "TwiceDifferentiableHV" begin
            for (name, prob) in MultivariateProblems.UnconstrainedProblems.examples
                if prob.istwicedifferentiable
                    hv!(storage::Vector, x::Vector, v::Vector) = begin
                        n = length(x)
                        H = Matrix{Float64}(undef, n, n)
                        MVP.hessian(prob)(H, x)
                        storage .= H * v
                    end
                    fg!(G::Vector, x::Vector) = begin
                        MVP.gradient(prob)(G, x)
                        MVP.objective(prob)(x), G
                    end
                    ddf = TwiceDifferentiableHV(MVP.objective(prob), fg!, hv!, prob.initial_x)
                    x = rand(prob.initial_x, length(prob.initial_x))
                    v = rand(prob.initial_x, length(prob.initial_x))
                    G = NLSolversBase.alloc_DF(x, 0.0)
                    H = NLSolversBase.alloc_H(x, 0.0)
                    MVP.hessian(prob)(H, x)
                    @test hv_product!(ddf, x, v) == H*v
                    @test hv_product(ddf) == H*v
                    @test hv_product(ddf) == ddf.Hv
                    F, G = fg!(G, x)
                    @test gradient!(ddf, x) == G
                    @test value!(ddf, x) == F
                    @test f_calls(ddf) == 1
                    @test g_calls(ddf) == 1
                    @test h_calls(ddf) == 0
                    @test hv_calls(ddf) == 1
                    clear!(ddf)
                    @test f_calls(ddf) == 0
                    @test g_calls(ddf) == 0
                    @test h_calls(ddf) == 0
                    @test hv_calls(ddf) == 0
                end
            end
        end
    end
    @testset "multivalued" begin
        # Test example: Rosenbrock MINPACK
        function f!(F::Vector, x::Vector)
            F[1] = 1 - x[1]
            F[2] = 10(x[2]-x[1]^2)
            F
        end
        function j!(J::Matrix, x::Vector)
            J[1,1] = -1
            J[1,2] = 0
            J[2,1] = -20x[1]
            J[2,2] = 10
            J
        end

        x_seed = [0.0, 0.0]
        F_seed = [0.0, 0.0]
        J_seed = [0.0 0.0; 0.0 0.0]
        F_x_seed = [1.0, 0.0]
        J_x_seed = [-1.0 0.0; -0.0 10.0]

        x_alt = [0.5, 0.5]
        F_x_alt = [0.5, 2.5]
        J_x_alt = [-1.0 0.0; -10.0 10.0]

        # Construct instances
        nd = NonDifferentiable(f!, x_seed, F_seed)
        od = OnceDifferentiable(f!, j!, x_seed, F_seed, J_seed)

        # Force evaluation
        value!!(nd, nd.F, x_seed)
        value_jacobian!!(od, od.F, od.DF, x_seed)

        # Test that values are the same, and that values match those
        # calculated by the value(obj, x) methods
        @test value(nd) == value(od) == F_x_seed
        @test value(nd, x_seed) == value(od, x_seed)
        @test value(nd, x_alt) == value(od, x_alt)

        # Test that the Jacobians match the intended values
        @test jacobian(od) == J_x_seed

        # Test that the call counters got incremented
        @test nd.f_calls == od.f_calls == [3]
        @test od.df_calls == [1]

        # Test that the call counters do not get incremented
        # with single-"bang" methods...
        value!(nd, x_seed)
        value_jacobian!(od, x_seed)

        @test nd.f_calls == od.f_calls == [3]
        @test od.df_calls == [1]

        # ... and that they do with double-"bang" methods
        value!!(nd, x_seed)
        value_jacobian!!(od, x_seed)

        @test nd.f_calls == od.f_calls == [4]
        @test od.df_calls == [2]

        # Test that jacobian doesn't work for NonDifferentiable, but does otherwise
        @test_throws ErrorException jacobian!(nd, x_alt)
        jacobian!(od, x_alt)

        @test value(nd) == value(od) == F_x_seed
        @test jacobian(od) == J_x_alt
        @test nd.f_calls == od.f_calls == [4]
        @test od.df_calls == [3]

        @test value(nd) == value(od) == F_x_seed
        @test jacobian(od) == J_x_alt
        @test nd.f_calls == od.f_calls == [4]
        @test od.df_calls == [3]

        value!(nd, x_alt)
        value!(od, x_alt)
        @test value(nd) == value(od) == F_x_alt
        @test jacobian(od) == J_x_alt
        @test nd.f_calls == od.f_calls == [5]
        @test od.df_calls == [3]

        @test_throws ErrorException value_jacobian!(nd, x_seed)
        value_jacobian!(od, x_seed)
        @test value(od) == F_x_seed
        # change x_f manually to test branch
        od.x_f .*= 0
        value_jacobian!(od, x_seed)
        @test value(od) == F_x_seed
        # change x_df manually to test branch
        # only df is meant to be re-calculated as
        # d.x_f == x_seed
        od.x_df .*= 0
        value_jacobian!(od, x_seed)
        @test value(od) == F_x_seed
        @test jacobian(od) == J_x_seed
        @test od.f_calls == [6]
        @test od.df_calls == [4]

        clear!(nd)
        clear!(od)
        @test nd.f_calls == [0,]
        @test od.f_calls == [0,]
        @test od.df_calls == [0,]

        # Test that the correct branch gets called if jacobian hasn't been
        # calculated yet
        xx = rand(2)
        @test jacobian!(od, xx) == j!(NLSolversBase.alloc_DF(xx, xx), xx)
        @test od.x_df == xx

        # Test the branching in value_jacobian! works as expected
        xxx = rand(2)
        value_jacobian!(od, xxx)
        @test xxx == od.x_f == od.x_df
        xxx2 = rand(2)
        value!(od, xxx2)
        @test xxx2 == od.x_f
        @test xxx == od.x_df
        xxx3 = rand(2)
        jacobian!(od, xxx3)
        @test xxx2 == od.x_f
        @test xxx3 == od.x_df

    end
end
