@testset "sparse" begin
    @testset "𝐑ⁿ → 𝐑" begin
        f(x) = sum(x->x.^2, x)
        g(G, x) = copy!(G, 2 .* x)
        h(H, x) = H .= sparse(2.0I, size(H)...)

        obj_dense = TwiceDifferentiable(f, g, h, rand(40))
        @test !issparse(obj_dense.H)

        obj_sparse = TwiceDifferentiable(f, g, h, rand(40), 0.0, rand(40), sparse(1.0I, 40, 40))
        @test typeof(obj_sparse.H) <: SparseMatrixCSC

        function fgh!(F, G, H, x)
            if !(F == nothing)
                fx = sum(x->x.^2, x)
            end

            if !(G == nothing)
                copy!(G, 2 .* x)
            end
            if !(H == nothing)
                H .= sparse(2.0I, size(H)...)
            end
            return fx
        end

        obj_fgh = TwiceDifferentiable(NLSolversBase.only_fgh!(fgh!), rand(40), 0.0, rand(40), sparse(1.0I, 40, 40))
        @test typeof(obj_fgh.H) <: SparseMatrixCSC
    end
    @testset "𝐑ⁿ → 𝐑ⁿ" begin
        f(F, x) = copy!(F, 2 .* x)
        j(J, x) = J .= sparse(2.0I, size(J)...)

        # Test that with no spec on the Jacobian cache it is dense
        obj_dense = OnceDifferentiable(f, j, rand(40), rand(40))
        @test !issparse(obj_dense.DF)

        obj_dense = OnceDifferentiable(f, j, rand(40), rand(40), rand(40))
        @test !issparse(obj_dense.DF)


        obj_sparse = OnceDifferentiable(f, j, rand(40), rand(40), sparse(1.0I, 40, 40))
        @test typeof(obj_sparse.DF) <: SparseMatrixCSC

        function fj!(F, J, x)
            if !(F == nothing)
                copy!(G, 2 .* x)
            end
            if !(J == nothing)
                J .= sparse(2.0I, size(J)...)
            end
            return fx
        end

        obj_fj = OnceDifferentiable(NLSolversBase.only_fj!(fj!), rand(40), rand(40), sparse(1.0I, 40, 40))
        @test typeof(obj_fj.DF) <: SparseMatrixCSC
    end
end
