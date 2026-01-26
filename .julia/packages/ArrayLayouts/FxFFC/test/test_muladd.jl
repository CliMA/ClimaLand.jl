module TestMulAdd

using ArrayLayouts, FillArrays, Random, StableRNGs, LinearAlgebra, Test, Quaternions, StaticArrays
using ArrayLayouts: DenseColumnMajor, AbstractStridedLayout, AbstractColumnMajor, DiagonalLayout, mul, Mul, zero!

Random.seed!(0)
@testset "Multiplication" begin
    @testset "zero!" begin
        for A in (randn(5,5), [randn(5,5),randn(4,4)], [SVector(2,3), SVector(3,4)])
            zero!(A)
            @test iszero(A)
        end
    end

    @testset "Matrix * Vector" begin
        @testset "eltype" begin
            v = MulAdd(1,zeros(Int,2,2), zeros(Float64,2),0,zeros(2))
            A = MulAdd(1,zeros(Int,2,2), zeros(Float64,2,2), 0,zeros(2,2))
            @test eltype(v) == eltype(A) == Float64
            @test @inferred(axes(v)) == (@inferred(axes(v,1)),) == (Base.OneTo(2),)
            @test @inferred(size(v)) == (@inferred(size(v,1)),) == (2,)
            @test @inferred(axes(A)) == (@inferred(axes(A,1)),@inferred(axes(A,2))) == (Base.OneTo(2),Base.OneTo(2))
            @test @inferred(size(A)) == (@inferred(size(A,1)),size(A,2)) == (2,2)

            Ã = materialize(A)
            @test Ã isa Matrix{Float64}
            @test Ã == zeros(2,2)

            c = similar(v)
            fill!(c,NaN)
            @test copyto!(c,v) == zeros(2)
            fill!(c,NaN)
            c .= v
            @test c == zeros(2)
        end

        @testset "gemv Float64" begin
            for A in (randn(5,5), view(randn(5,5),:,:), view(randn(5,5),1:5,:),
                    view(randn(5,5),1:5,1:5), view(randn(5,5),:,1:5)),
                b in (randn(5), view(randn(5),:), view(randn(5),1:5), view(randn(9),1:2:9))
                c = similar(b);

                muladd!(1.0,A,b,0.0,c)
                @test c == A*b == BLAS.gemv!('N', 1.0, A, b, 0.0, similar(c))
                copyto!(c, MulAdd(1.0,A,b,0.0,c))
                @test c == A*b == BLAS.gemv!('N', 1.0, A, b, 0.0, similar(c))
                c .= MulAdd(1.0,A,b,0.0,c)
                @test c == A*b == BLAS.gemv!('N', 1.0, A, b, 0.0, similar(c))

                M = MulAdd(2.0, A,b, 0.0, c)
                @test M isa MulAdd{<:AbstractColumnMajor,<:AbstractStridedLayout,<:AbstractColumnMajor}
                c .= M
                @test c == BLAS.gemv!('N', 2.0, A, b, 0.0, similar(c))
                copyto!(c, M)
                @test c == BLAS.gemv!('N', 2.0, A, b, 0.0, similar(c))

                c = copy(b)
                M = MulAdd(1.0, A,b, 1.0, c)
                @test M isa MulAdd{<:AbstractColumnMajor,<:AbstractStridedLayout,<:AbstractColumnMajor}
                copyto!(c, M)
                @test c == BLAS.gemv!('N', 1.0, A, b, 1.0, copy(b))
            end
        end

        @testset "gemv mixed array types" begin
            (A, b, c) = (randn(5,5), randn(5), 1.0:5.0)
            d = similar(b)
            d .= MulAdd(1.0,A,b,1.0,c)
            @test d == BLAS.gemv!('N', 1.0, A, b, 1.0, Vector{Float64}(c))
        end

        @testset "gemv Complex" begin
            for T in (ComplexF64,),
                    A in (randn(T,5,5), view(randn(T,5,5),:,:)),
                    b in (randn(T,5), view(randn(T,5),:))
                c = similar(b);

                c .= MulAdd(one(T), A, b, zero(T), c)
                @test c == BLAS.gemv!('N', one(T), A, b, zero(T), similar(c))

                c = copy(b)
                muladd!(2one(T),A,b,one(T),c)
                @test c == BLAS.gemv!('N', 2one(T), A, b, one(T), copy(b))
            end
        end

        @testset "Matrix{Int} * Vector{Vector{Int}}" begin
            A, x =  [1 2; 3 4] , [[1,2],[3,4]]
            X = reshape([[1 2],[3 4], [5 6], [7 8]],2,2)
            @test mul(A,x) == A*x
            @test mul(A,X) == A*X
            @test mul(X,A) == X*A
        end

        @testset "Diagonal Fill" begin
            for (A, B) in (([1:4;], [3:6;]), (reshape([1:16;],4,4), reshape(2 .* [1:16;],4,4)))
                D = Diagonal(Fill(3, 4))
                M = MulAdd(2, D, A, 3, B)
                @test copy(M) == mul!(B, D, A, 2, 3)
                M = MulAdd(1, D, A, 0, B)
                @test copy(M) == mul!(B, D, A)
            end

            A, B = [1:4;], reshape([3:6;], 4, 1)
            D = Diagonal(Fill(3, 1))
            M = MulAdd(2, A, D, 3, B)
            @test copy(M) == (VERSION >= v"1.9" ? mul!(B, A, D, 2, 3) : 2 * A * D + 3 * B)
            M = MulAdd(1, A, D, 0, B)
            @test copy(M) == (VERSION >= v"1.9" ? mul!(B, A, D) : A * D)
        end
    end

    @testset "Matrix * Matrix" begin
        @testset "gemm" begin
            for A in (randn(5,5), view(randn(5,5),:,:), view(randn(5,5),1:5,:),
                    view(randn(5,5),1:5,1:5), view(randn(5,5),:,1:5)),
                B in (randn(5,5), view(randn(5,5),:,:), view(randn(5,5),1:5,:),
                        view(randn(5,5),1:5,1:5), view(randn(5,5),:,1:5))
                C = similar(B);
                D = similar(C);

                C .= MulAdd(1.0,A,B,0.0,C)
                @test C == BLAS.gemm!('N', 'N', 1.0, A, B, 0.0, D)

                C .= MulAdd(2.0,A,B,0.0,C)
                @test C == BLAS.gemm!('N', 'N', 2.0, A, B, 0.0, D)

                C = copy(B)
                C .= MulAdd(2.0,A,B,1.0,C)
                @test C == BLAS.gemm!('N', 'N', 2.0, A, B, 1.0, copy(B))
            end

            A, B = ones(100, 100), ones(100, 100)
            C = ones(100, 100)
            C .= MulAdd(2,A,B,1,C)
            @test C ≈ BLAS.gemm!('N', 'N', 2.0, A, B, 1.0, copy(B))

            A, B = Float64[i+j for i in 1:100, j in 1:100], Float64[i+j for i in 1:100, j in 1:100]
            C = ones(100, 100)
            C .= MulAdd(2,A,B,1,C)
            @test_broken C ≈ BLAS.gemm!('N', 'N', 2.0, A, B, 1.0, copy(B))
        end

        @testset "gemm Complex" begin
            for T in (ComplexF64,),
                    A in (randn(T,5,5), view(randn(T,5,5),:,:)),
                    B in (randn(T,5,5), view(randn(T,5,5),:,:))
                C = similar(B);

                C .= MulAdd(one(T),A,B,zero(T),C)
                @test C == BLAS.gemm!('N', 'N', one(T), A, B, zero(T), similar(C))

                C .= MulAdd(2one(T),A,B,zero(T),C)
                @test C == BLAS.gemm!('N', 'N', 2one(T), A, B, zero(T), similar(C))

                C = copy(B)
                C .= MulAdd(one(T),A,B,one(T),C)
                @test C == BLAS.gemm!('N', 'N', one(T), A, B, one(T), copy(B))
            end
        end

        @testset "gemm mixed array types" begin
            (A, B, C) = (randn(5,5), randn(5,5), reshape(1.0:25.0,5,5))
            D = similar(B)
            D .= MulAdd(1.0,A,B,1.0, C)
            @test D == BLAS.gemm!('N', 'N', 1.0, A, B, 1.0, Matrix{Float64}(C))
        end
    end

    @testset "adjtrans" begin
        @testset "gemv adjtrans" begin
            for A in (randn(5,5), view(randn(5,5),:,:), view(randn(5,5),1:5,:),
                    view(randn(5,5),1:5,1:5), view(randn(5,5),:,1:5)),
                        b in (randn(5), view(randn(5),:), view(randn(5),1:5), view(randn(9),1:2:9)),
                        Ac in (transpose(A), A')
                c = similar(b);

                c .= MulAdd(1.0,Ac,b,0.0,c)
                @test c == BLAS.gemv!('T', 1.0, A, b, 0.0, similar(c))

                c .= MulAdd(2.0,Ac,b,0.0,c)
                @test c == BLAS.gemv!('T', 2.0, A, b, 0.0, similar(c))

                c = copy(b)
                c .= MulAdd(1.0,Ac,b,1.0,c)
                @test c == BLAS.gemv!('T', 1.0, A, b, 1.0, copy(b))
            end
        end

        @testset "gemv adjtrans mixed types" begin
            (A, b, c) = (randn(5,5), randn(5), 1.0:5.0)
            for Ac in (transpose(A), A')
                d = similar(b)
                d .= MulAdd(1.0,Ac,b,1.0, c)
                @test d == BLAS.gemv!('T', 1.0, A, b, 1.0, Vector{Float64}(c))

                d .= MulAdd(3.0,Ac,b,2.0, c)
                @test d == BLAS.gemv!('T', 3.0, A, b, 2.0, Vector{Float64}(c))
            end
        end

        @testset "gemv adjtrans Complex" begin
            for T in (ComplexF64,),
                    A in (randn(T,5,5), view(randn(T,5,5),:,:)),
                    b in (randn(T,5), view(randn(T,5),:)),
                    (Ac,trans) in ((transpose(A),'T'), (A','C'))
                c = similar(b);

                c .= MulAdd(one(T),Ac,b,zero(T),c)
                @test c == BLAS.gemv!(trans, one(T), A, b, zero(T), similar(c))

                c = copy(b)
                c .= MulAdd(3one(T),Ac,b,2one(T),c)
                @test c == BLAS.gemv!(trans, 3one(T), A, b, 2one(T), copy(b))
            end

            C = randn(6,6) + im*randn(6,6)
            V = view(C', 2:3, 3:4)
            c = randn(2) + im*randn(2)
            @test muladd!(1.0+0im,V,c,0.0+0im,similar(c,2)) == BLAS.gemv!('C', 1.0+0im, Matrix(V'), c, 0.0+0im, similar(c,2))
            @test muladd!(1.0+0im,V',c,0.0+0im,similar(c,2)) == BLAS.gemv!('N', 1.0+0im, Matrix(V'), c, 0.0+0im, similar(c,2))
        end

        @testset "gemm adjtrans" begin
            for A in (randn(5,5), view(randn(5,5),:,:)),
                B in (randn(5,5), view(randn(5,5),1:5,:))
                for Ac in (transpose(A), A')
                    C = copy(B)
                    C .= MulAdd(3.0, Ac, B, 2.0, C)
                    @test C == BLAS.gemm!('T', 'N', 3.0, A, B, 2.0, copy(B))
                end
                for Bc in (transpose(B), B')
                    C = copy(B)
                    C .= MulAdd(3.0, A, Bc, 2.0, C)
                    @test C == BLAS.gemm!('N', 'T', 3.0, A, B, 2.0, copy(B))
                end
                for Ac in (transpose(A), A'), Bc in (transpose(B), B')
                    C = copy(B)
                    d = similar(C)
                    d .= MulAdd(3.0, Ac, Bc, 2.0, C)
                    @test d == BLAS.gemm!('T', 'T', 3.0, A, B, 2.0, copy(B))
                end
            end
        end

        @testset "symv adjtrans" begin
            A = randn(100,100)
            x = randn(100)
            for As in (Symmetric(A), Hermitian(A), Symmetric(A)', Symmetric(view(A,:,:)',:L), view(Symmetric(A),:,:))
                @test muladd!(1.0,As,x,0.0,similar(x)) == BLAS.symv!('U', 1.0, A, x, 0.0, similar(x))
            end

            for As in (Symmetric(A,:L), Hermitian(A,:L), Symmetric(A,:L)', Symmetric(view(A,:,:)',:U), view(Hermitian(A,:L),:,:))
                @test muladd!(1.0,As,x,0.0,similar(x)) == BLAS.symv!('L', 1.0, A, x, 0.0, similar(x))
            end

            T = ComplexF64
            A = randn(T,100,100)
            x = randn(T,100)
            for As in (Symmetric(A), transpose(Symmetric(A)), Symmetric(transpose(view(A,:,:)),:L), view(Symmetric(A),:,:))
                @test muladd!(one(T),As,x,zero(T),similar(x)) == BLAS.symv!('U', one(T), A, x, zero(T), similar(x))
            end

            for As in (Symmetric(A,:L), transpose(Symmetric(A,:L)), Symmetric(transpose(view(A,:,:)),:U), view(Symmetric(A,:L),:,:))
                @test muladd!(one(T),As,x,zero(T),similar(x)) == BLAS.symv!('L', one(T), A, x, zero(T), similar(x))
            end

            for As in(Hermitian(A), Hermitian(A)', view(Hermitian(A),:,:))
                @test muladd!(one(T),As,x,zero(T),similar(x)) == BLAS.hemv!('U', one(T), A, x, zero(T), similar(x))
            end

            for As in(Hermitian(A,:L), Hermitian(A,:L)', view(Hermitian(A,:L),:,:))
                @test muladd!(one(T),As,x,zero(T),similar(x)) == BLAS.hemv!('L', one(T), A, x, zero(T), similar(x))
            end
        end
    end

    @testset "Mul" begin
        @testset "Mixed types" begin
            A = randn(5,6)
            b = rand(Int,6)
            c = Array{Float64}(undef, 5)
            d = similar(c)
            copyto!(d, MulAdd(1, A, b, 0.0, d))
            @test d == copyto!(similar(d), MulAdd(1, A, b, 0.0, d)) ≈ A*b
            @test copyto!(similar(d), MulAdd(1, A, b, 1.0, d)) ≈ A*b + d

            @test (similar(d) .= MulAdd(1, A, b, 1.0, d)) == copyto!(similar(d), MulAdd(1, A, b, 1.0, d))

            B = rand(Int,6,4)
            C = Array{Float64}(undef, 5, 4)

            @test_throws DimensionMismatch muladd!(1.0,A,B,0.0,similar(C,4,4))

            A = randn(Float64,20,22)
            B = randn(ComplexF64,22,24)
            C = similar(B,20,24)
            @test copyto!(similar(C), MulAdd(1.0, A, B, 0.0, C)) ≈ A*B
        end

        @testset "no allocation" begin
            A = randn(5,5); x = randn(5); y = randn(5); c = similar(y);
            muladd!(2.0, A, x, 3.0, y)
            @test @allocated(muladd!(2.0, A, x, 3.0, y)) == 0
            VA = view(A,:,1:2)
            vx = view(x,1:2)
            vy = view(y,:)
            muladd!(2.0, VA, vx, 3.0, vy)
            # spurious allocations in tests
            @test @allocated(muladd!(2.0, VA, vx, 3.0, vy)) < 100
        end

        @testset "BigFloat" begin
            A = BigFloat.(randn(5,5))
            x = BigFloat.(randn(5))
            @test A*x == muladd!(1.0,A,x,0.0,copy(x))
            @test_throws UndefRefError muladd!(1.0,A,x,0.0,similar(x))
        end

        @testset "Scalar * Vector" begin
            A, x =  [1 2; 3 4] , [[1,2],[3,4]]
            @test muladd!(1.0,A,x,0.0, 1.0*x) == A*x
        end

        @testset "adjoint" begin
            A = randn(5,5)
            V = view(A',1:3,1:3)
            b = randn(3)
            @test MemoryLayout(V) isa RowMajor
            @test MemoryLayout(V') isa ColumnMajor
            @test strides(V) == (5,1)
            @test strides(V') == (1,5)
            @test mul(V, b) ≈ V*b ≈ Matrix(V)*b
            @test mul(V', b) ≈ V'b ≈ Matrix(V)'*b
        end
    end

    @testset "Lmul/Rmul" begin
        @testset "tri Lmul" begin
            @testset "Float * Float vector" begin
                A = randn(Float64, 100, 100)
                x = randn(Float64, 100)

                L = Lmul(UpperTriangular(A),x)
                @test size(L) == (size(L,1),) == (100,)
                @test axes(L) == (axes(L,1),) == (Base.OneTo(100),)
                @test eltype(L) == Float64
                @test length(L) == 100

                @test similar(L) isa Vector{Float64}
                @test similar(L,Int) isa Vector{Int}

                @test mul(UpperTriangular(A),x) ==
                            copy(Lmul(UpperTriangular(A),x)) ==
                            ArrayLayouts.lmul(UpperTriangular(A),x) ==
                            ArrayLayouts.lmul!(UpperTriangular(A),copy(x)) ==
                            copyto!(similar(x),Lmul(UpperTriangular(A),x)) ==
                            UpperTriangular(A)*x ==
                            BLAS.trmv!('U', 'N', 'N', A, copy(x))
                @test copyto!(similar(x),Lmul(UnitUpperTriangular(A),x)) ==
                            UnitUpperTriangular(A)*x ==
                            BLAS.trmv!('U', 'N', 'U', A, copy(x))
                @test copyto!(similar(x),Lmul(LowerTriangular(A),x)) ==
                            LowerTriangular(A)*x ==
                            BLAS.trmv!('L', 'N', 'N', A, copy(x))
                @test ArrayLayouts.lmul!(UnitLowerTriangular(A),copy(x)) ==
                            UnitLowerTriangular(A)x ==
                            BLAS.trmv!('L', 'N', 'U', A, copy(x))

                @test ArrayLayouts.lmul!(UpperTriangular(A)',copy(x)) ==
                            ArrayLayouts.lmul!(LowerTriangular(A'),copy(x)) ==
                            UpperTriangular(A)'x ==
                            BLAS.trmv!('U', 'T', 'N', A, copy(x))
                @test ArrayLayouts.lmul!(UnitUpperTriangular(A)',copy(x)) ==
                            ArrayLayouts.lmul!(UnitLowerTriangular(A'),copy(x)) ==
                            UnitUpperTriangular(A)'*x ==
                            BLAS.trmv!('U', 'T', 'U', A, copy(x))
                @test ArrayLayouts.lmul!(LowerTriangular(A)',copy(x)) ==
                            ArrayLayouts.lmul!(UpperTriangular(A'),copy(x)) ==
                            LowerTriangular(A)'*x ==
                            BLAS.trmv!('L', 'T', 'N', A, copy(x))
                @test ArrayLayouts.lmul!(UnitLowerTriangular(A)',copy(x)) ==
                            ArrayLayouts.lmul!(UnitUpperTriangular(A'),copy(x)) ==
                            UnitLowerTriangular(A)'*x ==
                            BLAS.trmv!('L', 'T', 'U', A, copy(x))

                @test_throws DimensionMismatch mul(UpperTriangular(A), randn(5))
                @test_throws DimensionMismatch ArrayLayouts.lmul!(UpperTriangular(A), randn(5))
            end

            @testset "Float * Complex vector"  begin
                T = ComplexF64
                A = randn(T, 100, 100)
                x = randn(T, 100)

                @test (y = copy(x); y .= Lmul(UpperTriangular(A),y) ) ==
                            UpperTriangular(A)x ==
                            BLAS.trmv!('U', 'N', 'N', A, copy(x))
                @test (y = copy(x); y .= Lmul(UnitUpperTriangular(A),y) ) ==
                            UnitUpperTriangular(A)x ==
                            BLAS.trmv!('U', 'N', 'U', A, copy(x))
                @test (y = copy(x); y .= Lmul(LowerTriangular(A),y) ) ==
                            LowerTriangular(A)x ==
                            BLAS.trmv!('L', 'N', 'N', A, copy(x))
                @test (y = copy(x); y .= Lmul(UnitLowerTriangular(A),y) ) ==
                            UnitLowerTriangular(A)x ==
                            BLAS.trmv!('L', 'N', 'U', A, copy(x))
                @test LowerTriangular(A')  == UpperTriangular(A)'

                @test (y = copy(x); y .= Lmul(transpose(UpperTriangular(A)),y) ) ==
                            (similar(x) .= Lmul(LowerTriangular(transpose(A)),x)) ==
                            transpose(UpperTriangular(A))x ==
                            BLAS.trmv!('U', 'T', 'N', A, copy(x))
                @test (y = copy(x); y .= Lmul(transpose(UnitUpperTriangular(A)),y) ) ==
                            (similar(x) .= Lmul(UnitLowerTriangular(transpose(A)),x)) ==
                            transpose(UnitUpperTriangular(A))x ==
                            BLAS.trmv!('U', 'T', 'U', A, copy(x))
                @test (y = copy(x); y .= Lmul(transpose(LowerTriangular(A)),y) ) ==
                            (similar(x) .= Lmul(UpperTriangular(transpose(A)),x)) ==
                            transpose(LowerTriangular(A))x ==
                            BLAS.trmv!('L', 'T', 'N', A, copy(x))
                @test (y = copy(x); y .= Lmul(transpose(UnitLowerTriangular(A)),y) ) ==
                            (similar(x) .= Lmul(UnitUpperTriangular(transpose(A)),x)) ==
                            transpose(UnitLowerTriangular(A))x ==
                            BLAS.trmv!('L', 'T', 'U', A, copy(x))

                @test (y = copy(x); y .= Lmul(UpperTriangular(A)',y) ) ==
                            (similar(x) .= Lmul(LowerTriangular(A'),x)) ==
                            UpperTriangular(A)'x ==
                            BLAS.trmv!('U', 'C', 'N', A, copy(x))
                @test (y = copy(x); y .= Lmul(UnitUpperTriangular(A)',y) ) ==
                            (similar(x) .= Lmul(UnitLowerTriangular(A'),x)) ==
                            UnitUpperTriangular(A)'x ==
                            BLAS.trmv!('U', 'C', 'U', A, copy(x))
                @test (y = copy(x); y .= Lmul(LowerTriangular(A)',y) ) ==
                            (similar(x) .= Lmul(UpperTriangular(A'),x)) ==
                            LowerTriangular(A)'x ==
                            BLAS.trmv!('L', 'C', 'N', A, copy(x))
                @test (y = copy(x); y .= Lmul(UnitLowerTriangular(A)',y) ) ==
                            (similar(x) .= Lmul(UnitUpperTriangular(A'),x)) ==
                            UnitLowerTriangular(A)'x ==
                            BLAS.trmv!('L', 'C', 'U', A, copy(x))
            end

            @testset "Float * Float Matrix" begin
                A = randn(Float64, 100, 100)
                x = randn(Float64, 100, 100)

                @test UpperTriangular(A)*x ≈ (similar(x) .= Lmul(UpperTriangular(A), x))
            end

            @testset "adjtrans" begin
                for T in (Float64, ComplexF64)
                    A = randn(T,100,100)
                    b = randn(T,100)

                    @test copy(Lmul(UpperTriangular(A)', b)) == UpperTriangular(A)'b
                    @test copy(Lmul(UnitUpperTriangular(A)', b)) == UnitUpperTriangular(A)'b
                    @test copy(Lmul(LowerTriangular(A)', b)) == LowerTriangular(A)'b
                    @test copy(Lmul(UnitLowerTriangular(A)', b)) == UnitLowerTriangular(A)'b

                    @test copy(Lmul(transpose(UpperTriangular(A)), b)) == transpose(UpperTriangular(A))b
                    @test copy(Lmul(transpose(UnitUpperTriangular(A)), b)) == transpose(UnitUpperTriangular(A))b
                    @test copy(Lmul(transpose(LowerTriangular(A)), b)) == transpose(LowerTriangular(A))b
                    @test copy(Lmul(transpose(UnitLowerTriangular(A)), b)) == transpose(UnitLowerTriangular(A))b

                    B = randn(T,100,100)

                    @test copy(Lmul(UpperTriangular(A)', B)) == UpperTriangular(A)'B
                    @test copy(Lmul(UnitUpperTriangular(A)', B)) == UnitUpperTriangular(A)'B
                    @test copy(Lmul(LowerTriangular(A)', B)) == LowerTriangular(A)'B
                    @test copy(Lmul(UnitLowerTriangular(A)', B)) == UnitLowerTriangular(A)'B

                    @test copy(Lmul(transpose(UpperTriangular(A)), B)) == transpose(UpperTriangular(A))B
                    @test copy(Lmul(transpose(UnitUpperTriangular(A)), B)) == transpose(UnitUpperTriangular(A))B
                    @test copy(Lmul(transpose(LowerTriangular(A)), B)) == transpose(LowerTriangular(A))B
                    @test copy(Lmul(transpose(UnitLowerTriangular(A)), B)) == transpose(UnitLowerTriangular(A))B
                end

                for T in (Float64, ComplexF64)
                    A = big.(randn(T,100,100))
                    b = big.(randn(T,100))

                    @test copy(Lmul(UpperTriangular(A)', b)) ≈ UpperTriangular(A)'b
                    @test copy(Lmul(UnitUpperTriangular(A)', b)) ≈ UnitUpperTriangular(A)'b
                    @test copy(Lmul(LowerTriangular(A)', b)) ≈ LowerTriangular(A)'b
                    @test copy(Lmul(UnitLowerTriangular(A)', b)) ≈ UnitLowerTriangular(A)'b

                    @test copy(Lmul(transpose(UpperTriangular(A)), b)) ≈ transpose(UpperTriangular(A))b
                    @test copy(Lmul(transpose(UnitUpperTriangular(A)), b)) ≈ transpose(UnitUpperTriangular(A))b
                    @test copy(Lmul(transpose(LowerTriangular(A)), b)) ≈ transpose(LowerTriangular(A))b
                    @test copy(Lmul(transpose(UnitLowerTriangular(A)), b)) ≈ transpose(UnitLowerTriangular(A))b

                    B = big.(randn(T,100,100))

                    @test copy(Lmul(UpperTriangular(A)', B)) ≈ UpperTriangular(A)'B
                    @test copy(Lmul(UnitUpperTriangular(A)', B)) ≈ UnitUpperTriangular(A)'B
                    @test copy(Lmul(LowerTriangular(A)', B)) ≈ LowerTriangular(A)'B
                    @test copy(Lmul(UnitLowerTriangular(A)', B)) ≈ UnitLowerTriangular(A)'B

                    @test copy(Lmul(transpose(UpperTriangular(A)), B)) ≈ transpose(UpperTriangular(A))B
                    @test copy(Lmul(transpose(UnitUpperTriangular(A)), B)) ≈ transpose(UnitUpperTriangular(A))B
                    @test copy(Lmul(transpose(LowerTriangular(A)), B)) ≈ transpose(LowerTriangular(A))B
                    @test copy(Lmul(transpose(UnitLowerTriangular(A)), B)) ≈ transpose(UnitLowerTriangular(A))B
                end
            end

            @testset "BigFloat" begin
                A = BigFloat.(randn(10,10))
                b = BigFloat.(randn(10))
                @test mul(UpperTriangular(A), b) == UpperTriangular(A)*b
                @test mul(UnitUpperTriangular(A), b) == UnitUpperTriangular(A)*b
                @test mul(LowerTriangular(A), b) == LowerTriangular(A)*b
                @test mul(UnitLowerTriangular(A), b) == UnitLowerTriangular(A)*b
                @test_throws DimensionMismatch mul(UpperTriangular(A), randn(5))
                @test_throws DimensionMismatch ArrayLayouts.lmul!(UpperTriangular(A), randn(5))
                @test_throws DimensionMismatch ArrayLayouts.lmul!(UnitUpperTriangular(A), randn(5))
                @test_throws DimensionMismatch ArrayLayouts.lmul!(LowerTriangular(A), randn(5))
                @test_throws DimensionMismatch ArrayLayouts.lmul!(UnitLowerTriangular(A), randn(5))
            end
        end

        @testset "tri Rmul" begin
            for T in (Float64, ComplexF64)
                A = randn(T, 100,100)
                B = randn(T, 100,100)
                R = Rmul(copy(A), UpperTriangular(B))
                @test size(R) == (size(R,1),size(R,2)) == (100,100)
                @test axes(R) == (axes(R,1),axes(R,2)) == (Base.OneTo(100),Base.OneTo(100))
                @test eltype(R) == T
                @test length(R) == 100^2

                @test similar(R) isa Matrix{T}
                @test similar(R,Int) isa Matrix{Int}

                R2 = deepcopy(R)
                @test mul(A, UpperTriangular(B)) == BLAS.trmm('R', 'U', 'N', 'N', one(T), B, A) == copyto!(similar(R2), R2) == materialize!(R)
                @test R.A ≠ A
                @test BLAS.trmm('R', 'U', 'T', 'N', one(T), B, A) == copy(Rmul(A, transpose(UpperTriangular(B)))) == A*transpose(UpperTriangular(B))
                @test BLAS.trmm('R', 'U', 'N', 'U', one(T), B, A) == copy(Rmul(A, UnitUpperTriangular(B)))== A*UnitUpperTriangular(B)
                @test BLAS.trmm('R', 'U', 'T', 'U', one(T), B, A) == copy(Rmul(A, transpose(UnitUpperTriangular(B)))) == A*transpose(UnitUpperTriangular(B))
                @test BLAS.trmm('R', 'L', 'N', 'N', one(T), B, A) == copy(Rmul(A, LowerTriangular(B)))== A*LowerTriangular(B)
                @test BLAS.trmm('R', 'L', 'T', 'N', one(T), B, A) == copy(Rmul(A, transpose(LowerTriangular(B)))) == A*transpose(LowerTriangular(B))
                @test BLAS.trmm('R', 'L', 'N', 'U', one(T), B, A) == copy(Rmul(A, UnitLowerTriangular(B)))== A*UnitLowerTriangular(B)
                @test BLAS.trmm('R', 'L', 'T', 'U', one(T), B, A) == copy(Rmul(A, transpose(UnitLowerTriangular(B)))) == A*transpose(UnitLowerTriangular(B))

                if T <: Complex
                    @test BLAS.trmm('R', 'U', 'C', 'N', one(T), B, A) == copy(Rmul(A, UpperTriangular(B)')) == A*UpperTriangular(B)'
                    @test BLAS.trmm('R', 'U', 'C', 'U', one(T), B, A) == copy(Rmul(A, UnitUpperTriangular(B)')) == A*UnitUpperTriangular(B)'
                    @test BLAS.trmm('R', 'L', 'C', 'N', one(T), B, A) == copy(Rmul(A, LowerTriangular(B)')) == A*LowerTriangular(B)'
                    @test BLAS.trmm('R', 'L', 'C', 'U', one(T), B, A) == copy(Rmul(A, UnitLowerTriangular(B)')) == A*UnitLowerTriangular(B)'
                end
            end

            T = Float64
            A = big.(randn(100,100))
            B = big.(randn(100,100))
            @test ArrayLayouts.rmul!(copy(A),UpperTriangular(B)) == A*UpperTriangular(B)
        end

        @testset "Diagonal" begin
            A = randn(5,5)
            B = Diagonal(randn(5))
            @test MemoryLayout(B) == DiagonalLayout{DenseColumnMajor}()

            @test A*B == ArrayLayouts.rmul!(copy(A),B) == mul(A,B)
            @test B*A == ArrayLayouts.lmul!(B,copy(A)) == mul(B,A)
            @test B*B == ArrayLayouts.lmul!(B, copy(B)) == mul(B, B)
        end

        @testset "tri * tri" begin
            A = randn(5,5)
            B = randn(5,5)
            @test UpperTriangular(A) * UpperTriangular(B) ≈ mul(UpperTriangular(A) , UpperTriangular(B))
            @test LowerTriangular(A) * LowerTriangular(B) ≈ mul(LowerTriangular(A) , LowerTriangular(B))
            @test LowerTriangular(A) * UpperTriangular(B) ≈ mul(LowerTriangular(A) , UpperTriangular(B))
            @test UnitUpperTriangular(A) * UnitUpperTriangular(B) ≈ mul(UnitUpperTriangular(A) , UnitUpperTriangular(B))
            @test UnitLowerTriangular(A) * UnitLowerTriangular(B) ≈ mul(UnitLowerTriangular(A) , UnitLowerTriangular(B))
            @test UnitLowerTriangular(A) * UnitUpperTriangular(B) ≈ mul(UnitLowerTriangular(A) , UnitUpperTriangular(B))
            @test UnitUpperTriangular(A) * UpperTriangular(B) ≈ mul(UnitUpperTriangular(A) , UpperTriangular(B))
            @test @inferred(mul(UpperTriangular(A) , UpperTriangular(B))) isa UpperTriangular
            @test @inferred(mul(LowerTriangular(A) , LowerTriangular(B))) isa LowerTriangular
            @test @inferred(mul(UnitUpperTriangular(A) , UnitUpperTriangular(B))) isa UnitUpperTriangular
            @test @inferred(mul(UnitLowerTriangular(A) , UnitLowerTriangular(B))) isa UnitLowerTriangular
            @test @inferred(mul(UpperTriangular(A) , UnitUpperTriangular(B))) isa UpperTriangular
            @test @inferred(mul(LowerTriangular(A), UpperTriangular(B))) isa Matrix
            @test @inferred(mul(UnitUpperTriangular(A), LowerTriangular(B))) isa Matrix
        end

        @testset "diag * tri" begin
            D = Diagonal(randn(5))
            U = UpperTriangular(randn(5,5))
            @test mul(D,U) == D*U
            @test mul(U,D) == U*D
        end
    end

    @testset "MulAdd" begin
        A = randn(5,5)
        B = randn(5,4)
        C = randn(5,4)
        b = randn(5)
        c = randn(5)

        M = MulAdd(2.0,A,B,3.0,C)
        @test size(M) == size(C)
        @test size(M,1) == size(C,1)
        @test size(M,2) == size(C,2)
        @test_broken size(M,3) == size(C,3)
        @test length(M) == length(C)
        @test axes(M) == axes(C)
        @test eltype(M) == Float64
        @test materialize(M) ≈ 2.0A*B + 3.0C

        @test_throws DimensionMismatch materialize(MulAdd(2.0,A,randn(3),1.0,B))
        @test_throws DimensionMismatch materialize(MulAdd([1,2],A,B,[1,2],C))
        @test_throws DimensionMismatch materialize(MulAdd(2.0,A,B,3.0,randn(3,4)))
        @test_throws DimensionMismatch materialize(MulAdd(2.0,A,B,3.0,randn(5,5)))

        B = randn(5,5)
        C = randn(5,5)
        @test materialize(MulAdd(2.0,Diagonal(A),Diagonal(B),3.0,Diagonal(C))) ==
            muladd!(2.0,Diagonal(A),Diagonal(B),3.0,Diagonal(copy(C))) == 2.0Diagonal(A)*Diagonal(B) + 3.0*Diagonal(C)
        @test_broken materialize(MulAdd(2.0,Diagonal(A),Diagonal(B),3.0,Diagonal(C))) isa Diagonal

        @test materialize(MulAdd(1.0, Eye(5), A, 3.0, C)) == muladd!(1.0, Eye(5), A, 3.0, copy(C)) == A + 3.0C
        @test materialize(MulAdd(1.0, A, Eye(5), 3.0, C)) == muladd!(1.0, A, Eye(5), 3.0, copy(C)) == A + 3.0C

        @testset "Degenerate" begin
            C = BigFloat.(randn(3,3))
            @test muladd!(1.0, BigFloat.(randn(3,0)), randn(0,3), 2.0, copy(C) ) == 2C
            C = BigFloat.(randn(0,0))
            @test muladd!(1.0, BigFloat.(randn(0,3)), randn(3,0), 2.0, copy(C)) == C
        end

        @testset "Generic muladd" begin
            A = view(reshape(1:6,2,3),1:2,:)
            @test IndexStyle(typeof(A)) == IndexCartesian()
            b = randn(3)
            c = similar(b,2)
            @test muladd!(1.0, A, b, 0.0, c) == A*b
        end
    end

    @testset "Q" begin
        Q = qr(randn(5,5)).Q
        b = randn(5)
        B = randn(5,5)
        @test Q*1.0 ≈ ArrayLayouts.lmul!(Q, Matrix{Float64}(I, 5, 5))
        @test Q*b == ArrayLayouts.lmul!(Q, copy(b)) == mul(Q,b)
        @test Q*B == ArrayLayouts.lmul!(Q, copy(B)) == mul(Q,B)
        @test B*Q == ArrayLayouts.rmul!(copy(B), Q) == mul(B,Q)
        @test 1.0*Q ≈ ArrayLayouts.rmul!(Matrix{Float64}(I, 5, 5), Q)
        @test 1.0*Q' ≈ ArrayLayouts.rmul!(Matrix{Float64}(I, 5, 5), Q')
        @test Q*Q ≈ mul(Q,Q)
        @test Q'*b == ArrayLayouts.lmul!(Q', copy(b)) == mul(Q',b)
        @test Q'*B == ArrayLayouts.lmul!(Q', copy(B)) == mul(Q',B)
        @test B*Q' == ArrayLayouts.rmul!(copy(B), Q') == mul(B,Q')
        @test Q*Q' ≈ mul(Q,Q')
        @test Q'*Q' ≈ mul(Q',Q')
        @test Q'*Q ≈ mul(Q',Q)
        @test Q*UpperTriangular(B) ≈ mul(Q, UpperTriangular(B))
        @test UpperTriangular(B)*Q ≈ mul(UpperTriangular(B), Q)

        Q = qr(randn(7,5)).Q
        b = randn(5)
        B = randn(5,5)
        @test Q*1.0 ≈ ArrayLayouts.lmul!(Q, Matrix{Float64}(I, 7, 7))
        @test Q*b == mul(Q,b)
        @test Q*B == mul(Q,B)
        @test 1.0*Q ≈ ArrayLayouts.rmul!(Matrix{Float64}(I, 7, 7), Q)
        @test Q*Q ≈ mul(Q,Q)
        @test B*Q' == mul(B,Q')
        @test Q*Q' ≈ mul(Q,Q')
        @test Q'*Q' ≈ mul(Q',Q')
        @test Q'*Q ≈ mul(Q',Q)
        VERSION >= v"1.8-" && @test Q*UpperTriangular(B) ≈ mul(Q, UpperTriangular(B))
        @test UpperTriangular(B)*Q' ≈ mul(UpperTriangular(B), Q')
    end

    @testset "Mul" begin
        A = randn(5,5)
        b = randn(5)
        B = randn(5,5)

        M = Mul(A,b)
        @test size(M) == (size(M,1),) == (5,)
        @test length(M) == 5
        @test axes(M) == (axes(M,1),) == (Base.OneTo(5),)
        @test M[1] ≈ M[CartesianIndex(1)] ≈ (A*b)[1]
        @test ArrayLayouts.mul!(similar(b), A, b) ≈ A*b

        M = Mul(A,B)
        @test size(M) == (size(M,1),size(M,2)) == (5,5)
        @test length(M) == 25
        @test axes(M) == (axes(M,1),axes(M,2)) == (Base.OneTo(5),Base.OneTo(5))
        @test M[1,1] ≈ M[CartesianIndex(1,1)] ≈ M[1] ≈ (A*B)[1,1]
        @test similar(M) isa Matrix{Float64}
        @test similar(M, Int) isa Matrix{Int}
        @test similar(M, Int, (Base.OneTo(5),)) isa Vector{Int}

        M = Mul(b, b')
        @test size(M) == (size(M,1),size(M,2)) == (5,5)
        @test length(M) == 25
        @test axes(M) == (axes(M,1),axes(M,2)) == (Base.OneTo(5),Base.OneTo(5))
        @test M[1,1] ≈ M[CartesianIndex(1,1)] ≈ M[1] ≈ (b*b')[1,1]
    end

    @testset "Dot/Dotu" begin
        a = randn(5)
        b = randn(5)
        c = randn(5) + im*randn(5)
        d = randn(5) + im*randn(5)

        @test ArrayLayouts.dot(a,b) ≈ ArrayLayouts.dotu(a,b) ≈ mul(a',b)
        @test ArrayLayouts.dot(a,b) ≈ dot(a,b)
        @test eltype(Dot(a,1:5)) == Float64

        @test ArrayLayouts.dot(c,d) == mul(c',d)
        @test ArrayLayouts.dotu(c,d) == mul(transpose(c),d)
        @test ArrayLayouts.dot(c,d) ≈ dot(c,d)
        @test ArrayLayouts.dotu(c,d) ≈ BLAS.dotu(c,d)

        @test ArrayLayouts.dot(c,b) == mul(c',b)
        @test ArrayLayouts.dotu(c,b) == mul(transpose(c),b)
        @test ArrayLayouts.dot(c,b) ≈ dot(c,b)

        @test ArrayLayouts.dot(a,d) == mul(a',d)
        @test ArrayLayouts.dotu(a,d) == mul(transpose(a),d)
        @test ArrayLayouts.dot(a,d) ≈ dot(a,d)
    end

    @testset "adjtrans muladd" begin
        A,B = [1 2], [[1,2]',[3,4]']
        B̃ = [transpose([1,2]), transpose([3,4])]
        @test copy(MulAdd(A,B)) == A*B
        @test eltype(MulAdd(A,B)) == eltype(B)
        @test copy(MulAdd(A,B̃)) == A*B̃
        @test eltype(MulAdd(A,B̃)) == eltype(B̃)
    end

    @testset "Bidiagonal" begin
        BidiagU = Bidiagonal(randn(5), randn(4), :U)
        BidiagL = Bidiagonal(randn(5), randn(4), :L)
        Tridiag = Tridiagonal(randn(4), randn(5), randn(4))
        SymTri = SymTridiagonal(randn(5), randn(4))
        Diag = Diagonal(randn(5))
        @test typeof(mul(BidiagU,Diag)) <: Bidiagonal
        @test typeof(mul(BidiagL,Diag)) <: Bidiagonal
        @test typeof(mul(Tridiag,Diag)) <: Tridiagonal
        @test typeof(mul(SymTri,Diag))  <: Tridiagonal

        @test typeof(mul(BidiagU,Diag)) <: Bidiagonal
        @test typeof(mul(Diag,BidiagL)) <: Bidiagonal
        @test typeof(mul(Diag,Tridiag)) <: Tridiagonal
        @test typeof(mul(Diag,SymTri))  <: Tridiagonal
    end

    @testset "tiled_blasmul!" begin
        rng = StableRNG(1)
        X = randn(rng, ComplexF64, 8, 4)
        Y = randn(rng, 8, 2)
        @test mul(Y',X) ≈ Y'X

        for A in (randn(5,5), view(randn(5,5),:,:), view(randn(5,5),1:5,:),
                    view(randn(5,5),1:5,1:5), view(randn(5,5),:,1:5)),
            B in (randn(5,5), view(randn(5,5),:,:), view(randn(5,5),1:5,:),
                    view(randn(5,5),1:5,1:5), view(randn(5,5),:,1:5))
            C = similar(B);
            D = similar(C);

            C .= MulAdd(1,A,B,0,C)
            @test C ≈ BLAS.gemm!('N', 'N', 1.0, A, B, 0.0, D)

            C = copy(B)
            C .= MulAdd(2,A,B,1,C)
            @test C ≈ BLAS.gemm!('N', 'N', 2.0, A, B, 1.0, copy(B))
        end
    end

    @testset "Vec * Adj" begin
        @test ArrayLayouts.mul(1:5, (1:4)') == (1:5) * (1:4)'
    end

    @testset "Fill" begin
        mutable struct MFillMat{T} <: FillArrays.AbstractFill{T,2,NTuple{2,Base.OneTo{Int}}}
            x :: T
            sz :: NTuple{2,Int}
        end
        MFillMat(x::T, sz::Vararg{Int,2}) where {T} = MFillMat{T}(x, sz)
        Base.size(M::MFillMat) = M.sz
        FillArrays.getindex_value(M::MFillMat) = M.x
        Base.copyto!(M::MFillMat, A::Broadcast.Broadcasted) = (M.x = only(unique(A)); M)
        Base.copyto!(M::MFillMat, A::Broadcast.Broadcasted{<:Base.Broadcast.AbstractArrayStyle{0}}) = (M.x = only(unique(A)); M)

        M = MulAdd(1, Fill(2,4,4), Fill(3,4,4), 2, MFillMat(2,4,4))
        X = copy(M)
        @test X == Fill(28,4,4)

        M = MulAdd(1, Fill(2,4,4), Fill(3,4,4), 0, MFillMat(2,4,4))
        X = copy(M)
        @test X == Fill(24,4,4)
    end

    @testset "non-commutative" begin
        A = [quat(rand(4)...) for i in 1:4, j in 1:4]
        B = [quat(rand(4)...) for i in 1:4, j in 1:4]
        C = [quat(rand(4)...) for i in 1:4, j in 1:4]
        α, β = quat(0,0,0,1), quat(0,1,0,0)
        M = MulAdd(α, A, B, β, C)
        @test copy(M) ≈ mul!(copy(C), A, B, α, β) ≈ A * B * α + C * β

        SA = Symmetric(A)
        M = MulAdd(α, SA, B, β, C)
        @test copy(M) ≈ mul!(copy(C), SA, B, α, β) ≈ SA * B * α + C * β

        B = [quat(rand(4)...) for i in 1:4]
        C = [quat(rand(4)...) for i in 1:4]
        M = MulAdd(α, A, B, β, C)
        @test copy(M) ≈ mul!(copy(C), A, B, α, β) ≈ A * B * α + C * β

        M = MulAdd(α, SA, B, β, C)
        @test copy(M) ≈ mul!(copy(C), SA, B, α, β) ≈ SA * B * α + C * β

        A = [quat(rand(4)...) for i in 1:4]
        B = [quat(rand(4)...) for i in 1:1, j in 1:1]
        C = [quat(rand(4)...) for i in 1:4, j in 1:1]
        M = MulAdd(α, A, B, β, C)
        @test copy(M) ≈ mul!(copy(C), A, B, α, β) ≈ A * B * α + C * β

        D = Diagonal(Fill(quat(rand(4)...), 4))
        b = [quat(rand(4)...) for i in 1:4]
        c = [quat(rand(4)...) for i in 1:4]
        M = MulAdd(α, D, b, β, c)
        @test copy(M) ≈ mul!(copy(c), D, b, α, β) ≈ D * b * α + c * β

        D = Diagonal(Fill(quat(rand(4)...), 1))
        b = [quat(rand(4)...) for i in 1:4]
        c = [quat(rand(4)...) for i in 1:4, j in 1:1]
        M = MulAdd(α, b, D, β, c)
        if VERSION >= v"1.9"
            @test copy(M) ≈ mul!(copy(c), b, D, α, β) ≈ b * D * α + c * β
        else
            @test copy(M) ≈ b * D * α + c * β
        end
    end

    @testset "Error paths" begin
        if VERSION >= v"1.10.0"
            Q = qr(rand(2,2), ColumnNorm()).Q
        else
            Q = qr(rand(2,2), Val(true)).Q
        end
        v = rand(Float32, 3)
        @test_throws DimensionMismatch ArrayLayouts.materialize!(ArrayLayouts.Rmul(v, Q))
        @test_throws DimensionMismatch ArrayLayouts.materialize!(ArrayLayouts.Rmul(v, Q'))
        @test_throws DimensionMismatch ArrayLayouts.materialize!(ArrayLayouts.Lmul(Q, v))
        @test_throws DimensionMismatch ArrayLayouts.materialize!(ArrayLayouts.Lmul(Q', v))
    end

    @testset "dual" begin
        a = randn(5)
        X = randn(5,6)
        @test copyto!(similar(a,6)', MulAdd(2.0, a', X, 3.0, Zeros(6)')) ≈ 2a'*X
        @test copyto!(transpose(similar(a,6)), MulAdd(2.0, a', X, 3.0, Zeros(6)')) ≈ 2a'*X
        @test copyto!(transpose(similar(a,6)), MulAdd(2.0, transpose(a), X, 3.0, transpose(Zeros(6)))) ≈ 2a'*X
    end
end

end
