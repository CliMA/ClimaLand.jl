module TestLdiv

using ArrayLayouts, LinearAlgebra, FillArrays, Test
import ArrayLayouts: ApplyBroadcastStyle, QRCompactWYQLayout, QRCompactWYLayout, QRPackedQLayout, QRPackedLayout
using StaticArrays

@testset "Ldiv" begin
    @testset "Float64 \\ *" begin
        A = randn(5,5)
        b = randn(5)
        B = randn(5,5)
        M = Ldiv(A,b)

        @test M[1] ≈ (A\b)[1]
        @test Ldiv(A,B)[1,1] ≈ (A\B)[1,1]

        @test ndims(M) == 1
        @test size(M) == (5,)
        @test similar(M) isa Vector{Float64}
        @test materialize(M) isa Vector{Float64}
        @test materialize(M) == (A\b) == ldiv(A,b)

        @test Base.BroadcastStyle(typeof(Ldiv(A,b))) isa ApplyBroadcastStyle

        @test copyto!(similar(b), Ldiv(A,b)) ==
                    (similar(b) .= Ldiv(A,b)) == materialize(Ldiv(A,b)) ==
                  (A\b) == (b̃ =  copy(b); LAPACK.gesv!(copy(A), b̃); b̃)

        @test copyto!(similar(b), Ldiv(UpperTriangular(A) , b)) ≈ UpperTriangular(A) \ b
        @test copyto!(similar(b), Ldiv(UpperTriangular(A) , b)) ==
                    (similar(b) .= Ldiv(UpperTriangular(A),b)) ==
                    BLAS.trsv('U', 'N', 'N', A, b)

        @test copyto!(similar(b), Ldiv(UpperTriangular(A)' , b)) ≈ UpperTriangular(A)' \ b
        @test copyto!(similar(b), Ldiv(UpperTriangular(A)' , b)) ==
                    (similar(b) .= Ldiv(UpperTriangular(A)',b)) ==
                    copyto!(similar(b), Ldiv(transpose(UpperTriangular(A)) , b)) ==
                            (similar(b) .= Ldiv(transpose(UpperTriangular(A)),b)) ==
                    BLAS.trsv('U', 'T', 'N', A, b)
    end

    @testset "ComplexF64 \\ *" begin
        T = ComplexF64
        A = randn(T,5,5)
        b = randn(T,5)
        @test copyto!(similar(b), Ldiv(A,b)) ==
                    (similar(b) .= Ldiv(A,b)) ==
                  (A\b) ≈ (b̃ =  copy(b); LAPACK.gesv!(copy(A), b̃); b̃)

        @test copyto!(similar(b), Ldiv(UpperTriangular(A) , b)) ≈ UpperTriangular(A) \ b
        @test copyto!(similar(b), Ldiv(UpperTriangular(A) , b)) ==
                    (similar(b) .= Ldiv(UpperTriangular(A),b)) ==
                    BLAS.trsv('U', 'N', 'N', A, b)

        @test copyto!(similar(b), Ldiv(UpperTriangular(A)' , b)) ≈ UpperTriangular(A)' \ b
        @test copyto!(similar(b), Ldiv(UpperTriangular(A)' , b)) ==
                    (similar(b) .= Ldiv(UpperTriangular(A)',b)) ==
                    BLAS.trsv('U', 'C', 'N', A, b)

        @test copyto!(similar(b), Ldiv(transpose(UpperTriangular(A)) , b)) ≈ transpose(UpperTriangular(A)) \ b
        @test copyto!(similar(b), Ldiv(transpose(UpperTriangular(A)) , b)) ==
                            (similar(b) .= Ldiv(transpose(UpperTriangular(A)),b)) ==
                    BLAS.trsv('U', 'T', 'N', A, b)
    end

    @testset "BigFloat Triangular \\" begin
        for T in (Float64, ComplexF64)
            A = big.(randn(T,5,5))
            b = big.(randn(T,5))
            @test ArrayLayouts.ldiv!(UpperTriangular(A),copy(b)) ≈ UpperTriangular(A)\b
            @test ArrayLayouts.ldiv!(UnitUpperTriangular(A),copy(b)) ≈ UnitUpperTriangular(A)\b
            @test ArrayLayouts.ldiv!(LowerTriangular(A),copy(b)) ≈ LowerTriangular(A)\b
            @test ArrayLayouts.ldiv!(UnitLowerTriangular(A),copy(b)) ≈ UnitLowerTriangular(A)\b

            @test ArrayLayouts.ldiv!(UpperTriangular(A)',copy(b)) ≈ UpperTriangular(A)'\b
            @test ArrayLayouts.ldiv!(UnitUpperTriangular(A)',copy(b)) ≈ UnitUpperTriangular(A)'\b
            @test ArrayLayouts.ldiv!(LowerTriangular(A)',copy(b)) ≈ LowerTriangular(A)'\b
            @test ArrayLayouts.ldiv!(UnitLowerTriangular(A)',copy(b)) ≈ UnitLowerTriangular(A)'\b

            @test ArrayLayouts.ldiv!(transpose(UpperTriangular(A)),copy(b)) ≈ transpose(UpperTriangular(A))\b
            @test ArrayLayouts.ldiv!(transpose(UnitUpperTriangular(A)),copy(b)) ≈ transpose(UnitUpperTriangular(A))\b
            @test ArrayLayouts.ldiv!(transpose(LowerTriangular(A)),copy(b)) ≈ transpose(LowerTriangular(A))\b
            @test ArrayLayouts.ldiv!(transpose(UnitLowerTriangular(A)),copy(b)) ≈ transpose(UnitLowerTriangular(A))\b

            B = big.(randn(T,5,5))
            @test ArrayLayouts.ldiv!(UpperTriangular(A),copy(B)) ≈ UpperTriangular(A)\B
        end
    end

    @testset "Triangular \\ matrix" begin
        A = randn(5,5)
        B = randn(5,5)
        M =  Ldiv(UpperTriangular(A), B)
        @test Base.Broadcast.broadcastable(M) === M
        @test UpperTriangular(A) \ B ≈ copyto!(similar(B) , Ldiv(UpperTriangular(A), B)) ≈ (B .= Ldiv(UpperTriangular(A), B))
        @test UpperTriangular(A') \ B ≈ copyto!(similar(B) , Ldiv(UpperTriangular(A'), B)) ≈ (B .= Ldiv(UpperTriangular(A'), B))
        @test LowerTriangular(A') \ B ≈ copyto!(similar(B) , Ldiv(LowerTriangular(A'), B)) ≈ (B .= Ldiv(LowerTriangular(A'), B))
        @test UnitUpperTriangular(A) \ B ≈ copyto!(similar(B) , Ldiv(UnitUpperTriangular(A), B)) ≈ (B .= Ldiv(UnitUpperTriangular(A), B))
        @test UnitUpperTriangular(A') \ B ≈ copyto!(similar(B) , Ldiv(UnitUpperTriangular(A'), B)) ≈ (B .= Ldiv(UnitUpperTriangular(A'), B))
        @test UnitLowerTriangular(A') \ B ≈ copyto!(similar(B) , Ldiv(UnitLowerTriangular(A'), B)) ≈ (B .= Ldiv(UnitLowerTriangular(A'), B))

        C = randn(5,5) + im * randn(5,5)
        D = randn(5,5) + im * randn(5,5)
        @test UpperTriangular(C) \ D ≈ copyto!(similar(D) , Ldiv(UpperTriangular(C), D)) ≈ (D .= Ldiv(UpperTriangular(C), D))
        @test UpperTriangular(C') \ D ≈ copyto!(similar(D) , Ldiv(UpperTriangular(C'), D)) ≈ (D .= Ldiv(UpperTriangular(C'), D))
        @test LowerTriangular(C') \ D ≈ copyto!(similar(D) , Ldiv(LowerTriangular(C'), D)) ≈ (D .= Ldiv(LowerTriangular(C'), D))
        @test UnitUpperTriangular(C) \ D ≈ copyto!(similar(D) , Ldiv(UnitUpperTriangular(C), D)) ≈ (D .= Ldiv(UnitUpperTriangular(C), D))
        @test UnitUpperTriangular(C') \ D ≈ copyto!(similar(D) , Ldiv(UnitUpperTriangular(C'), D)) ≈ (D .= Ldiv(UnitUpperTriangular(C'), D))
        @test UnitLowerTriangular(C') \ D ≈ copyto!(similar(D) , Ldiv(UnitLowerTriangular(C'), D)) ≈ (D .= Ldiv(UnitLowerTriangular(C'), D))
    end

    @testset "Int" begin
        A = [1 2 ; 3 4]; b = [5,6];
        @test eltype(Ldiv(A, b)) == Float64
    end

    @testset "Rdiv" begin
        @testset "Float64 \\ *" begin
            A = randn(3,5)
            B = randn(5,5)
            M = Rdiv(A,B)

            @test eltype(M) == Float64
            @test ndims(M) == 2
            @test size(M) == (3,5)
            @test axes(M) == (Base.OneTo(3),Base.OneTo(5))
            @test similar(M) isa Matrix{Float64}
        end
    end

    @testset "Diagonal" begin
        D = Diagonal(randn(5))
        F = Eye(5)
        A = randn(5,5)
        @test copy(Ldiv(D,F)) == ldiv(D,F) == D \ F
        @test ldiv(F,D) == F \ D
        @test ldiv(D,A) ≈ D \ A
        @test ldiv(A,D) == A \ D
        @test ldiv(F,A) ≈ F \ A
        @test ldiv(A,F) == A \ F

        @test copy(Rdiv(D,F)) == rdiv(D,F) == D / F
        @test rdiv(F,D) == F / D
        @test rdiv(A,D) ≈ A / D
        @test rdiv(A,F) == A / F
    end

    @testset "QR" begin
        @testset "QRCompactWYQ" begin
            for T in (Float64, ComplexF64)
                A = randn(T,10,10)
                b = randn(T,10)
                B = randn(T,10,10)
                F = qr(A)
                Q = F.Q
                @test MemoryLayout(F) isa QRCompactWYLayout
                @test MemoryLayout(Q) isa QRCompactWYQLayout
                @test ArrayLayouts.lmul!(Q,copy(b)) == lmul!(Q,copy(b))
                @test ArrayLayouts.lmul!(Q',copy(b)) == ArrayLayouts.ldiv!(Q,copy(b)) == lmul!(Q',copy(b))
                @test ArrayLayouts.lmul!(Q,copy(B)) == lmul!(Q,copy(B))
                @test ArrayLayouts.lmul!(Q',copy(B)) == ArrayLayouts.ldiv!(Q,copy(B)) == lmul!(Q',copy(B))
                @test ArrayLayouts.rmul!(copy(B),Q) == rmul!(copy(B),Q)
                @test ArrayLayouts.rmul!(copy(B),Q') == ArrayLayouts.rdiv!(copy(B),Q) == rmul!(copy(B),Q')
                @test ArrayLayouts.ldiv!(F,copy(b)) ≈ ldiv!(F,copy(b)) # only approx since we use BLAS.trsv!
                @test ArrayLayouts.ldiv!(F,copy(B)) ≈ ldiv!(F,copy(B)) # only approx since we use BLAS.trsv!

                @test copyto!(similar(b), Lmul(Q,b)) == Q*b
                @test copyto!(similar(B), Lmul(Q,B)) == Q*B
            end
        end

        @testset "QRPacked" begin
            @testset "Blas square" begin
                for T in (Float64, ComplexF64)
                    A = randn(T,10,10)
                    b = randn(T,10)
                    B = randn(T,10,10)
                    F = LinearAlgebra.qrfactUnblocked!(copy(A))
                    Q = F.Q
                    @test MemoryLayout(F) isa QRPackedLayout
                    @test MemoryLayout(Q) isa QRPackedQLayout
                    @test ArrayLayouts.lmul!(Q,copy(b)) == lmul!(Q,copy(b))
                    @test ArrayLayouts.lmul!(Q',copy(b)) == ArrayLayouts.ldiv!(Q,copy(b)) == lmul!(Q',copy(b))
                    @test ArrayLayouts.lmul!(Q,copy(B)) == lmul!(Q,copy(B))
                    @test ArrayLayouts.lmul!(Q',copy(B)) == ArrayLayouts.ldiv!(Q,copy(B)) == lmul!(Q',copy(B))
                    @test ArrayLayouts.rmul!(copy(B),Q) == rmul!(copy(B),Q)
                    @test ArrayLayouts.rmul!(copy(B),Q') ≈ ArrayLayouts.rdiv!(copy(B),Q) ≈ rmul!(copy(B),Q')
                    @test ArrayLayouts.ldiv!(F,copy(b)) ≈ ldiv!(F,copy(b)) # only approx since we use BLAS.trsv!
                    @test ArrayLayouts.ldiv!(F,copy(B)) ≈ ldiv!(F,copy(B)) # only approx since we use BLAS.trsv!

                    @test copyto!(similar(b), Lmul(Q,b)) == Q*b
                    @test copyto!(similar(B), Lmul(Q,B)) == Q*B
                end
            end
            @testset "BigFloat" begin
                A = BigFloat.(randn(10,10))
                b = BigFloat.(randn(10))
                B = BigFloat.(randn(10,10))
                F = LinearAlgebra.qrfactUnblocked!(copy(A))
                Q = F.Q
                @test MemoryLayout(F) isa QRPackedLayout
                @test MemoryLayout(Q) isa QRPackedQLayout
                @test ArrayLayouts.lmul!(Q,copy(b)) == lmul!(Q,copy(b))
                @test ArrayLayouts.lmul!(Q',copy(b)) == ArrayLayouts.ldiv!(Q,copy(b)) == lmul!(Q',copy(b))
                @test ArrayLayouts.lmul!(Q,copy(B)) == lmul!(Q,copy(B))
                @test ArrayLayouts.lmul!(Q',copy(B)) == ArrayLayouts.ldiv!(Q,copy(B)) == lmul!(Q',copy(B))
                @test ArrayLayouts.rmul!(copy(B),Q) == rmul!(copy(B),Q)
                @test ArrayLayouts.rmul!(copy(B),Q') ≈ ArrayLayouts.rdiv!(copy(B),Q) ≈ rmul!(copy(B),Q')
                @test ArrayLayouts.ldiv!(F,copy(b)) ≈ ldiv!(F,copy(b)) # only approx since we use BLAS.trsv!
                @test ArrayLayouts.ldiv!(F,copy(B)) ≈ ldiv!(F,copy(B)) # only approx since we use BLAS.trsv!

                @test copyto!(similar(b), Lmul(Q,b)) == Q*b
                @test copyto!(similar(B), Lmul(Q,B)) == Q*B
            end
            @testset "rectangular" begin
                for T in (Float64,BigFloat)
                    A = T.(randn(12,10))
                    b = T.(randn(12))
                    B = T.(randn(12,12))
                    F = LinearAlgebra.qrfactUnblocked!(copy(A))
                    Q = F.Q
                    @test ArrayLayouts.lmul!(Q,copy(b)) == lmul!(Q,copy(b))
                    @test ArrayLayouts.lmul!(Q',copy(b)) == ArrayLayouts.ldiv!(Q,copy(b)) == lmul!(Q',copy(b))
                    @test ArrayLayouts.lmul!(Q,copy(B)) == lmul!(Q,copy(B))
                    @test ArrayLayouts.lmul!(Q',copy(B)) == ArrayLayouts.ldiv!(Q,copy(B)) == lmul!(Q',copy(B))
                    @test ArrayLayouts.rmul!(copy(B),Q) == rmul!(copy(B),Q)
                    @test ArrayLayouts.rmul!(copy(B),Q') ≈ ArrayLayouts.rdiv!(copy(B),Q) ≈ rmul!(copy(B),Q')
                    @test ArrayLayouts.ldiv!(F,copy(b)) ≈ ldiv!(F,copy(b)) # only approx since we use BLAS.trsv!
                    @test ArrayLayouts.ldiv!(F,copy(B)) ≈ ldiv!(F,copy(B)) # only approx since we use BLAS.trsv!

                    @test_throws BoundsError ArrayLayouts.ldiv!(F,randn(10))

                    A = T.(randn(10,12))
                    b = T.(randn(10))
                    B = T.(randn(10,10))
                    F = LinearAlgebra.qrfactUnblocked!(copy(A))
                    Q = F.Q
                    @test ArrayLayouts.lmul!(Q,copy(b)) == lmul!(Q,copy(b))
                    @test ArrayLayouts.lmul!(Q',copy(b)) == ArrayLayouts.ldiv!(Q,copy(b)) == lmul!(Q',copy(b))
                    @test ArrayLayouts.lmul!(Q,copy(B)) == lmul!(Q,copy(B))
                    @test ArrayLayouts.lmul!(Q',copy(B)) == ArrayLayouts.ldiv!(Q,copy(B)) == lmul!(Q',copy(B))
                    @test ArrayLayouts.rmul!(copy(B),Q) == rmul!(copy(B),Q)
                    @test ArrayLayouts.rmul!(copy(B),Q') ≈ ArrayLayouts.rdiv!(copy(B),Q) ≈ rmul!(copy(B),Q')

                    @test_throws DimensionMismatch ArrayLayouts.ldiv!(F,copy(b))
                    @test_throws DimensionMismatch ArrayLayouts.ldiv!(F,copy(B))

                    b = T.(randn(12))
                    B = T.(randn(12,12))
                    @test ArrayLayouts.ldiv!(F,copy(b)) ≈ A\b[1:10]
                    @test ArrayLayouts.ldiv!(F,copy(B)) ≈ A\B[1:10,:]
                end
            end
            @testset "ldiv!" begin
                A = randn(5,5)
                F = LinearAlgebra.qrfactUnblocked!(copy(A))
                b = randn(5)
                v = view(copy(b),:)
                @test ArrayLayouts.ldiv!(F, v) === v
                @test ArrayLayouts.ldiv!(F, view(copy(b),:)) ≈ A \ b
            end

            @testset "row/colsupport" begin
                n = 10
                B = Bidiagonal(randn(n), randn(n-1), :L)
                Q = LinearAlgebra.qrfactUnblocked!(Matrix(B))
                Q = LinearAlgebra.QRPackedQ(Tridiagonal(Q.factors), Q.τ)
                @test rowsupport(Q, 4) ≡ colsupport(Q', 4) ≡ 3:10
                @test colsupport(Q, 4) ≡ rowsupport(Q', 4) ≡ Base.OneTo(5)
            end
        end
    end

    @testset "Bidiagonal" begin
        n = 10
        L = Bidiagonal(randn(n), randn(n-1), :L)
        U = Bidiagonal(randn(n), randn(n-1), :U)
        b = randn(n)
        B = randn(n,2)

        @test ArrayLayouts.ldiv(L,b) == L \ b
        @test ArrayLayouts.ldiv(U,b) == U \ b
        @test ArrayLayouts.ldiv(L,B) == L \ B
        @test ArrayLayouts.ldiv(U,B) == U \ B

        @test_throws DimensionMismatch ArrayLayouts.ldiv(L,randn(3))
        @test_throws DimensionMismatch materialize!(Ldiv(L,randn(3)))
    end

    @testset "Diagonal Ldiv bug (BandedMatrices #188)" begin
        n = 5
        B = randn(n,n)
        D = Diagonal(-collect(1.0:n))
        @test ArrayLayouts.ldiv!(similar(B), D, B) == D \ B
        @test ArrayLayouts.ldiv!(similar(B), B, D) == B \ D
        @test ArrayLayouts.rdiv!(similar(B), B, D) ==  B / D
        @test_broken ArrayLayouts.rdiv!(similar(B), D, B) ==  D / B
    end

    @testset "error paths" begin
        v = rand(Int,3)
        A = rand(2,2)
        S = SMatrix{2,2}(A)
        errA(U) = DimensionMismatch("second dimension of A, $(size(U,2)), does not match length of x, $(length(v))")
        errS(U) = DimensionMismatch("second axis of A, $(axes(U,2)), and first axis of B, $(axes(v,1)), must match")
        for (M, errf) in ((A, errA), (S, errS))
            U = UpperTriangular(M)
            err = errf(U)
            @test_throws err ArrayLayouts.materialize!(ArrayLayouts.Ldiv(U, v))
            UU = UnitUpperTriangular(M)
            @test_throws err ArrayLayouts.materialize!(ArrayLayouts.Ldiv(UU, v))
            L = LowerTriangular(M)
            @test_throws err ArrayLayouts.materialize!(ArrayLayouts.Ldiv(L, v))
            UL = UnitLowerTriangular(M)
            @test_throws err ArrayLayouts.materialize!(ArrayLayouts.Ldiv(UL, v))
        end
    end
end

end
