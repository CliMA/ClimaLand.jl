@testset "Subtyping" begin
    for ad in every_ad()
        sparse_ad = AutoSparse(ad)
        @test sparse_ad isa AbstractADType
        if mode(ad) isa ForwardMode
            @test mode(sparse_ad) isa ForwardMode
        elseif mode(ad) isa ForwardOrReverseMode
            @test mode(sparse_ad) isa ForwardOrReverseMode
        elseif mode(ad) isa ReverseMode
            @test mode(sparse_ad) isa ReverseMode
        elseif mode(ad) isa SymbolicMode
            @test mode(sparse_ad) isa SymbolicMode
        end
        @test dense_ad(ad) == ad
        @test dense_ad(sparse_ad) == ad
        @test sparsity_detector(sparse_ad) == NoSparsityDetector()
        @test coloring_algorithm(sparse_ad) == NoColoringAlgorithm()
    end
end

@testset "Sparsity detector" begin
    f_jac1(x) = vcat(x, x)
    f_jac2(x) = hcat(x, x)
    function f_jac! end
    f_hess(x) = sum(x)

    @testset "NoSparsityDetector" begin
        sd = NoSparsityDetector()

        for x in (rand(2), rand(2, 3)), f in (f_jac1, f_jac2)

            y = f(x)
            Js = jacobian_sparsity(f, x, sd)
            @test Js isa AbstractMatrix
            @test size(Js) == (length(y), length(x))
            @test all(isone, Js)
        end

        for x in (rand(2), rand(2, 3)), y in (rand(5), rand(5, 6))

            Js = jacobian_sparsity(f_jac!, y, x, sd)
            @test Js isa AbstractMatrix
            @test size(Js) == (length(y), length(x))
            @test all(isone, Js)
        end

        for x in (rand(2), rand(2, 3))
            Hs = hessian_sparsity(f_hess, x, sd)
            @test Hs isa AbstractMatrix
            @test size(Hs) == (length(x), length(x))
            @test all(isone, Hs)
        end
    end

    @testset "KnownJacobianSparsityDetector" begin
        @testset "Jacobian sparsity detection" begin
            @testset "Out-of-place functions" begin
                for sx in ((2,), (2, 3)), f in (f_jac1, f_jac2)

                    x = rand(sx...)
                    nx = length(x)
                    Jref = rand(Bool, 2 * nx, nx)
                    sd = KnownJacobianSparsityDetector(Jref)
                    Js = jacobian_sparsity(f, x, sd)
                    @test Js isa AbstractMatrix
                    @test size(Js) == (2 * nx, nx)
                    @test Js === Jref
                end
            end
            @testset "In-place functions" begin
                for sx in ((2,), (2, 3)), sy in ((5,), (5, 6))

                    x, y = rand(sx...), rand(sy...)
                    nx, ny = length(x), length(y)
                    Jref = rand(Bool, ny, nx)
                    sd = KnownJacobianSparsityDetector(Jref)
                    Js = jacobian_sparsity(f_jac!, y, x, sd)
                    @test Js isa AbstractMatrix
                    @test size(Js) == (ny, nx)
                    @test Js === Jref
                end
            end
        end
        @testset "Exceptions: Hessian sparsity detection not supported" begin
            for sx in ((2,), (2, 3))
                x = rand(sx...)
                nx = length(x)
                Href = rand(Bool, nx, nx)
                sd = KnownJacobianSparsityDetector(Href)
                @test_throws ArgumentError hessian_sparsity(f_hess, x, sd)
            end
        end
        @testset "Exceptions: DimensionMismatch" begin
            sd = KnownJacobianSparsityDetector(rand(Bool, 6, 7)) # wrong Jacobian size
            for x in (rand(2), rand(2, 3)), f in (f_jac1, f_jac2)

                @test_throws DimensionMismatch jacobian_sparsity(f, x, sd)
            end
            for x in (rand(2), rand(2, 3)), y in (rand(5), rand(5, 6))

                @test_throws DimensionMismatch jacobian_sparsity(f_jac!, y, x, sd)
            end
        end
    end

    @testset "KnownHessianSparsityDetector" begin
        @testset "Hessian sparsity detection" begin
            for sx in ((2,), (2, 3))
                x = rand(sx...)
                nx = length(x)
                Href = rand(Bool, nx, nx)
                sd = KnownHessianSparsityDetector(Href)

                Hs = hessian_sparsity(f_hess, x, sd)
                @test Hs isa AbstractMatrix
                @test size(Hs) == (nx, nx)
                @test Hs === Href
            end
        end
        @testset "Exceptions: Jacobian sparsity detection not supported" begin
            @testset "Out-of-place functions" begin
                for sx in ((2,), (2, 3)), f in (f_jac1, f_jac2)

                    x = rand(sx...)
                    nx = length(x)
                    Jref = rand(Bool, 2 * nx, nx)
                    sd = KnownHessianSparsityDetector(Jref)
                    @test_throws ArgumentError jacobian_sparsity(f, x, sd)
                end
            end
            @testset "In-place functions" begin
                for sx in ((2,), (2, 3)), sy in ((5,), (5, 6))

                    x, y = rand(sx...), rand(sy...)
                    nx, ny = length(x), length(y)
                    Jref = rand(Bool, ny, nx)
                    sd = KnownHessianSparsityDetector(Jref)
                    @test_throws ArgumentError jacobian_sparsity(f_jac!, y, x, sd)
                end
            end
        end
        @testset "Exceptions: DimensionMismatch" begin
            sd = KnownHessianSparsityDetector(rand(Bool, 2, 3)) #wrong Hessian size
            for x in (rand(2), rand(2, 3))
                @test_throws DimensionMismatch hessian_sparsity(f_hess, x, sd)
            end
        end
    end
end

@testset "Coloring algorithm" begin
    ca = NoColoringAlgorithm()

    for M in (rand(2, 3), rand(3, 2))
        cv = column_coloring(M, ca)
        @test cv isa AbstractVector{<:Integer}
        @test length(cv) == size(M, 2)
        @test allunique(cv)

        rv = row_coloring(M, ca)
        @test rv isa AbstractVector{<:Integer}
        @test length(rv) == size(M, 1)
        @test allunique(rv)
    end

    M = rand(3, 3)
    sv = symmetric_coloring(M, ca)
    @test sv isa AbstractVector{<:Integer}
    @test length(sv) == size(M, 1) == size(M, 2)
    @test allunique(sv)
end
