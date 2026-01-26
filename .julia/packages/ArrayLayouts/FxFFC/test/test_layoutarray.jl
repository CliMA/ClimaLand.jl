module TestLayoutArray

using ArrayLayouts, LinearAlgebra, FillArrays, Test, SparseArrays, Random
using ArrayLayouts: sub_materialize, ColumnNorm, RowMaximum, CRowMaximum, @_layoutlmul, Mul
import ArrayLayouts: triangulardata, MemoryLayout
import LinearAlgebra: Diagonal, Bidiagonal, Tridiagonal, SymTridiagonal

struct MyMatrix{T,M<:AbstractMatrix{T}} <: LayoutMatrix{T}
    A::M
end

Base.getindex(A::MyMatrix, k::Int, j::Int) = A.A[k,j]
Base.setindex!(A::MyMatrix, v, k::Int, j::Int) = setindex!(A.A, v, k, j)
Base.size(A::MyMatrix) = size(A.A)
Base.strides(A::MyMatrix) = strides(A.A)
Base.elsize(::Type{<:MyMatrix{T}}) where {T} = sizeof(T)
Base.cconvert(::Type{Ptr{T}}, A::MyMatrix{T}) where {T} = Base.cconvert(Ptr{T}, A.A)
Base.unsafe_convert(::Type{Ptr{T}}, A::MyMatrix{T}) where {T} = Base.unsafe_convert(Ptr{T}, A.A)
MemoryLayout(::Type{MyMatrix{T,M}}) where {T,M} = MemoryLayout(M)
Base.copy(A::MyMatrix) = MyMatrix(copy(A.A))
ArrayLayouts.bidiagonaluplo(M::MyMatrix) = ArrayLayouts.bidiagonaluplo(M.A)
for MT in (:Diagonal, :Bidiagonal, :Tridiagonal, :SymTridiagonal)
    @eval $MT(M::MyMatrix) = $MT(M.A)
end

struct MyVector{T,V<:AbstractVector{T}} <: LayoutVector{T}
    A::V
end

MyVector(M::MyVector) = MyVector(M.A)
Base.getindex(A::MyVector, k::Int) = A.A[k]
Base.setindex!(A::MyVector, v, k::Int) = setindex!(A.A, v, k)
Base.size(A::MyVector) = size(A.A)
Base.strides(A::MyVector) = strides(A.A)
Base.elsize(::Type{<:MyVector{T}}) where {T} = sizeof(T)
Base.cconvert(::Type{Ptr{T}}, A::MyVector{T}) where {T} = Base.cconvert(Ptr{T}, A.A)
Base.unsafe_convert(::Type{Ptr{T}}, A::MyVector{T}) where T = Base.unsafe_convert(Ptr{T}, A.A)
MemoryLayout(::Type{MyVector{T,V}}) where {T,V} = MemoryLayout(V)
Base.copy(A::MyVector) = MyVector(copy(A.A))

# These need to test dispatch reduces to ArrayLayouts.mul, etc.
@testset "LayoutArray" begin
    @testset "LayoutVector" begin
        a = MyVector([1.,2,3])
        B = randn(3,3)
        b = randn(3)

        @test a == a.A == Vector(a)
        @test a[1:3] == a.A[1:3]
        @test a[:] == a
        @test (a')[1,:] == (a')[1,1:3] == a
        @test sprint(show, "text/plain", a) == "$(summary(a)):\n 1.0\n 2.0\n 3.0"
        @test B*a ≈ B*a.A
        @test B'*a ≈ B'*a.A
        @test transpose(B)*a ≈ transpose(B)*a.A
        @test b'a ≈ transpose(b)a ≈ a'b ≈ transpose(a)b ≈ b'a.A
        @test qr(B).Q*a ≈ qr(B).Q*a.A

        @test a'a == transpose(a)a == dot(a,a) == dot(a,a.A) == dot(a.A,a) == 14
        v = view(a,1:3)
        @test copy(v) == sub_materialize(v) == a[1:3]
        @test dot(v,a) == dot(v,a.A) == dot(a,v) == dot(a.A,v) == dot(v,v) == 14
        @test norm(v) == norm(a) == norm([1,2,3])

        V = view(a',:,1:3)
        @test copy(V) == sub_materialize(V) == (a')[:,1:3]

        s = SparseVector(3, [1], [2])
        @test a's == s'a == dot(a,s) == dot(s,a) == dot(s,a.A)
    end

    @testset "LayoutMatrix" begin
        A = MyMatrix(randn(5,5))
        for (kr,jr) in ((1:2,2:3), (:,:), (:,1:2), (2:3,:), ([1,2],3:4), (:,[1,2]), ([2,3],:),
                        (2,:), (:,2), (2,1:3), (1:3,2))
            @test A[kr,jr] == A.A[kr,jr]
        end
        b = randn(5)
        B = randn(5,5)
        for Tri in (UpperTriangular, UnitUpperTriangular, LowerTriangular, UnitLowerTriangular)
            @test ldiv!(Tri(A), copy(b)) ≈ ldiv!(Tri(A.A), copy(b)) ≈ Tri(A.A) \ MyVector(b)
            @test ldiv!(Tri(A), copy(B)) ≈ ldiv!(Tri(A.A), copy(B)) ≈ Tri(A.A) \ MyMatrix(B)
            if VERSION ≥ v"1.9"
                @test rdiv!(copy(b)', Tri(A)) ≈ rdiv!(copy(b)', Tri(A.A)) ≈ MyVector(b)' / Tri(A.A)
                @test rdiv!(copy(B), Tri(A)) ≈ rdiv!(copy(B), Tri(A.A)) ≈ B / Tri(A.A)
            end
            @test lmul!(Tri(A), copy(b)) ≈ lmul!(Tri(A.A), copy(b)) ≈ Tri(A.A) * MyVector(b)
        end

        @test copyto!(MyMatrix(Array{Float64}(undef,5,5)), A) == A
        @test copyto!(Array{Float64}(undef,5,5), A) == A
        @test copyto!(MyMatrix(Array{Float64}(undef,5,5)), A.A) == A
        @test copyto!(view(MyMatrix(Array{Float64}(undef,5,5)),1:3,1:3), view(A,1:3,1:3)) == A[1:3,1:3]
        @test copyto!(view(MyMatrix(Array{Float64}(undef,5,5)),:,:), A) == A
        @test copyto!(MyMatrix(Array{Float64}(undef,3,3)), view(A,1:3,1:3)) == A[1:3,1:3]
        @test copyto!(view(MyMatrix(Array{Float64}(undef,5,5)),:,:), A.A) == A
        @test copyto!(Array{Float64}(undef,3,3), view(A,1:3,1:3)) == A[1:3,1:3]
        @test copyto!(MyMatrix(Array{Float64}(undef,5,5)), A') == A'
        @test copyto!(MyMatrix(Array{Float64}(undef,5,5)), view(A',:,:)) == A'
        @test copyto!(Array{Float64}(undef,5,5), A') == A'
        @test copyto!(Array{Float64}(undef,5,5), view(A',:,:)) == A'
        @test copyto!(view(MyMatrix(Array{Float64}(undef,5,5)),:,:), A') == A'

        @test copyto!(view(MyMatrix(Array{Float64}(undef,5,5)),:,:), view(A',:,:)) == A'

        @testset "factorizations" begin
            @test qr(A).factors ≈ qr(A.A).factors
            @test qr(A,ColumnNorm()).factors ≈ qr(A.A,ColumnNorm()).factors
            @test lu(A).factors ≈ lu(A.A).factors
            @test lu(A,RowMaximum()).factors ≈ lu(A.A,RowMaximum()).factors
            @test_throws ErrorException qr!(A)
            @test lu!(copy(A)).factors ≈ lu(A.A).factors
            b = randn(5)
            @test A \ b == A.A \ b == A.A \ MyVector(b) == ldiv!(lu(A.A), copy(b))
            @test A \ b == ldiv!(lu(A), copy(b))
            @test lu(A).L == lu(A.A).L
            @test lu(A).U == lu(A.A).U
            @test lu(A).p == lu(A.A).p
            @test lu(A).P == lu(A.A).P

            @test qr(A) isa LinearAlgebra.QRCompactWY
            @test inv(A) ≈ inv(A.A)

            S = Symmetric(MyMatrix(reshape(inv.(1:25),5,5) + 10I))
            @test cholesky(S).U ≈ @inferred(cholesky!(deepcopy(S))).U
            @test cholesky(S, CRowMaximum()).U ≈ cholesky(Matrix(S), CRowMaximum()).U
            @test cholesky!(deepcopy(S), CRowMaximum()).U ≈ cholesky(Matrix(S), CRowMaximum()).U
            @test cholesky(S) \ b ≈ cholesky(Matrix(S)) \ b ≈ cholesky(Matrix(S)) \ MyVector(b)
            @test cholesky(S, CRowMaximum()) \ b ≈ cholesky(Matrix(S), CRowMaximum()) \ b
            @test cholesky(S, CRowMaximum()) \ b ≈ ldiv!(cholesky(Matrix(S), CRowMaximum()), copy(b))
            @test cholesky(S) \ b ≈ Matrix(S) \ b ≈ Symmetric(Matrix(S)) \ b
            @test cholesky(S) \ b ≈ Symmetric(Matrix(S)) \ MyVector(b)
            if VERSION >= v"1.9"
                @test S \ b ≈ Matrix(S) \ b ≈ Symmetric(Matrix(S)) \ b
                @test S \ b ≈ Symmetric(Matrix(S)) \ MyVector(b)
            end

            S = Symmetric(MyMatrix(reshape(inv.(1:25),5,5) + 10I), :L)
            @test cholesky(S).U ≈ @inferred(cholesky!(deepcopy(S))).U
            @test cholesky(S,CRowMaximum()).U ≈ cholesky(Matrix(S),CRowMaximum()).U
            @test cholesky(S) \ b ≈ Matrix(S) \ b ≈ Symmetric(Matrix(S), :L) \ b
            @test cholesky(S) \ b ≈ Symmetric(Matrix(S), :L) \ MyVector(b)
            if VERSION >= v"1.9"
                @test S \ b ≈ Matrix(S) \ b ≈ Symmetric(Matrix(S), :L) \ b
                @test S \ b ≈ Symmetric(Matrix(S), :L) \ MyVector(b)
            end

            @testset "ldiv!" begin
                c = MyVector(randn(5))
                @test ldiv!(lu(A), MyVector(copy(c))) ≈ A \ c
                @test_throws ErrorException ldiv!(eigen(randn(5,5)), c)
                @test ArrayLayouts.ldiv!(svd(A.A), Vector(c)) ≈ ArrayLayouts.ldiv!(similar(c), svd(A.A), c) ≈ A \ c
                if VERSION ≥ v"1.8"
                    @test ArrayLayouts.ldiv!(similar(c), transpose(lu(A.A)), copy(c)) ≈ A'\c
                end

                B = Bidiagonal(randn(5), randn(4), :U)
                @test ldiv!(B, MyVector(copy(c))) ≈ B \ c
                @test ldiv!(Transpose(B), MyVector(copy(c))) ≈ transpose(B) \ c
                @test ldiv!(B', MyVector(copy(c))) ≈ B' \ c

                @test ldiv!(cholesky(S), MyVector(copy(c))) ≈ S \ c

                @test B \ MyVector(fill([1.,2],5)) ≈ B \ fill([1.,2],5)
                @test Transpose(B) \ MyVector(fill([1.,2],5)) ≈ transpose(B) \ fill([1.,2],5)
                @test Adjoint(B) \ MyVector(fill([1.,2],5)) ≈ transpose(B) \ fill([1.,2],5)
            end
        end
        Bin = randn(5,5)
        B = MyMatrix(copy(Bin))
        muladd!(1.0, A, A, 2.0, B)
        @test B == A.A^2 + 2Bin

        @testset "tiled_blasmul!" begin
            B = MyMatrix(copy(Bin))
            muladd!(1.0, Ones(5,5), A, 2.0, B)
        end

        @testset "mul!" begin
            @test mul!(B, A, A) ≈ A*A
            @test mul!(B, A', A) ≈ A'*A
            @test mul!(B, A, A') ≈ A*A'
            @test mul!(B, A', A') ≈ A'*A'
            @test mul!(B, A, Bin) ≈ A*Bin

            @test mul!(copy(B), A, A, 2, 3) ≈ 2A*A + 3B
            @test mul!(copy(B), A', A, 2, 3) ≈ 2A'*A + 3B
            @test mul!(copy(B), A, A', 2, 3) ≈ 2A*A' + 3B
            @test mul!(copy(B), A', A', 2, 3) ≈ 2A'*A' + 3B
            @test mul!(copy(B), Bin, A, 2, 3) ≈ 2Bin*A + 3B
            @test mul!(copy(B), Bin, A', 2, 3) ≈ 2Bin*A' + 3B
            @test mul!(copy(B), A, Bin, 2, 3) ≈ 2A*Bin + 3B
            @test mul!(copy(B), A', Bin, 2, 3) ≈ 2A'*Bin + 3B
            @test mul!(copy(B), A, Bin, 2, 3) ≈ 2A*Bin + 3B
            @test mul!(copy(B), A, Bin', 2, 3) ≈ 2A*Bin' + 3B
            @test mul!(copy(B), Bin', A, 2, 3) ≈ 2Bin'*A + 3B

            @test mul!(copy(B), A, Diagonal(Bin), 2, 3) ≈ 2A*Diagonal(Bin) + 3B
            @test mul!(copy(B), Diagonal(Bin), A, 2, 3) ≈ 2Diagonal(Bin)*A + 3B

            @test mul!(view(copy(B), 1:3, 1:3), view(A, 1:3, 1:3), view(B, 1:3, 1:3)) ≈ A[1:3,1:3]*B[1:3,1:3]
            @test mul!(Matrix{Float64}(undef, 3, 3), view(A, 1:3, 1:3), view(B, 1:3, 1:3)) ≈ A[1:3,1:3]*B[1:3,1:3]
            @test mul!(MyMatrix(Matrix{Float64}(undef, 3, 3)), view(A, 1:3, 1:3), view(B, 1:3, 1:3)) ≈ A[1:3,1:3]*B[1:3,1:3]
        end

        @testset "generic_blasmul!" begin
            A = BigFloat.(randn(5,5))
            Bin = BigFloat.(randn(5,5))
            B = copy(Bin)
            muladd!(1.0, Ones(5,5), A, 2.0, B)
            @test B == Ones(5,5)*A + 2.0Bin
        end

        C = MyMatrix(Float64[1 2; 3 4])
        @test sprint(show, "text/plain", C) == "$(summary(C)):\n 1.0  2.0\n 3.0  4.0"

        @testset "layoutldiv" begin
            A = MyMatrix(randn(5,5))
            x = randn(5)
            X = randn(5,5)
            t = view(randn(10),[1,3,4,6,7])
            T = view(randn(10,5),[1,3,4,6,7],:)
            t̃ = copy(t)
            T̃ = copy(T)
            B = Bidiagonal(randn(5),randn(4),:U)
            D = Diagonal(randn(5))
            @test ldiv!(A, copy(x)) ≈ A\x
            @test A\t ≈ A\t̃
            # QR is not general enough
            @test_broken ldiv!(A, t) ≈ A\t
            @test ldiv!(A, copy(X)) ≈ A\X
            @test A\T ≈ A\T̃
            VERSION >= v"1.9" && @test A/T ≈ A/T̃
            @test_broken ldiv!(A, T) ≈ A\T
            @test B\A ≈ B\Matrix(A)
            @test D \ A ≈ D \ Matrix(A)
            @test transpose(B)\A ≈ transpose(B)\Matrix(A) ≈ Transpose(B)\A ≈ Adjoint(B)\A
            @test B'\A ≈ B'\Matrix(A)
            @test A\A ≈ I
            VERSION >= v"1.9" && @test A/A ≈ I
            @test A\MyVector(x) ≈ A\x
            @test A\MyMatrix(X) ≈ A\X

            if VERSION >= v"1.9"
                @test A/A ≈ A.A / A.A
                @test x' / A ≈ x' / A.A
                @test transpose(x) / A ≈ transpose(x) / A.A 
                @test transpose(x) / A isa Transpose
                @test x' / A isa Adjoint
            end

            @test D \ UpperTriangular(A) ≈ D \ UpperTriangular(A.A)
            @test UpperTriangular(A) \ D ≈ UpperTriangular(A.A) \ D
        end

        @testset "dot" begin
            A = MyMatrix(randn(5,5))
            b = randn(5)
            @test dot(b, A, b) ≈ b'*(A*b) ≈ b'A*b
            @test dot(b, A, b) ≈ transpose(b)*(A*b) ≈ transpose(b)A*b
        end

        @testset "dual vector * symmetric (#40)" begin
            A = randn(3,3)
            x = rand(3)
            @test x' * Symmetric(MyMatrix(A)) ≈ x'Symmetric(A)
            @test transpose(x) * Symmetric(MyMatrix(A)) ≈ transpose(x)Symmetric(A)
        end

        @testset "map(copy, ::Diagonal)" begin
            # this is needed in BlockArrays
            D = Diagonal([MyMatrix(randn(2,2)), MyMatrix(randn(2,2))])
            @test map(copy, D) == D
        end

        @testset "permutedims(::Diagonal)" begin
            D = Diagonal(MyVector(randn(5)))
            @test permutedims(D) ≡ D
        end
    end

    @testset "l/rmul!" begin
        b = MyVector(randn(5))
        A = MyMatrix(randn(5,5))
        @test lmul!(2, deepcopy(b)) == rmul!(deepcopy(b), 2) == 2b
        @test lmul!(2, deepcopy(A)) == rmul!(deepcopy(A), 2) == 2A
        @test lmul!(2, deepcopy(A)') == rmul!(deepcopy(A)', 2) == 2A'
        @test lmul!(2, transpose(deepcopy(A))) == rmul!(transpose(deepcopy(A)), 2) == 2transpose(A)
        @test lmul!(2, Symmetric(deepcopy(A))) == rmul!(Symmetric(deepcopy(A)), 2) == 2Symmetric(A)
        @test lmul!(2, Hermitian(deepcopy(A))) == rmul!(Hermitian(deepcopy(A)), 2) == 2Hermitian(A)

        C = randn(ComplexF64,5,5)
        @test ArrayLayouts.lmul!(2, Hermitian(copy(C))) == ArrayLayouts.rmul!(Hermitian(copy(C)), 2) == 2Hermitian(C)


        @test ldiv!(2, deepcopy(b)) == rdiv!(deepcopy(b), 2) == 2\b
        @test ldiv!(2, deepcopy(A)) == rdiv!(deepcopy(A), 2) == 2\A
        @test ldiv!(2, deepcopy(A)') == rdiv!(deepcopy(A)', 2) == 2\A'
        @test ldiv!(2, transpose(deepcopy(A))) == rdiv!(transpose(deepcopy(A)), 2) == 2\transpose(A)
        @test ldiv!(2, Symmetric(deepcopy(A))) == rdiv!(Symmetric(deepcopy(A)), 2) == 2\Symmetric(A)
        @test ldiv!(2, Hermitian(deepcopy(A))) == rdiv!(Hermitian(deepcopy(A)), 2) == 2\Hermitian(A)
        @test ArrayLayouts.ldiv!(2, Hermitian(copy(C))) == ArrayLayouts.rdiv!(Hermitian(copy(C)), 2) == 2\Hermitian(C)
    end

    @testset "pow/I" begin
        A = randn(2,2)
        B = MyMatrix(A)
        @test B^2 ≈ A^2
        @test B^2.3 ≈ A^2.3
        @test B^(-1) ≈ inv(A)
        @test B + I ≈ I + B ≈ A + I
        @test B - I ≈ A - I
        @test I - B ≈ I - B
    end

    @testset "Diagonal" begin
        D = Diagonal(MyVector(randn(5)))
        D̃ = Diagonal(Vector(D.diag))
        B = randn(5,5)
        B̃ = MyMatrix(B)
        @test D*D ≈ Matrix(D)^2
        if VERSION ≥ v"1.7-"
            @test D^2 ≈ D*D
        end
        @test D*B ≈ Matrix(D)*B
        @test B*D ≈ B*Matrix(D)
        @test D*B̃ ≈ Matrix(D)*B̃
        @test B̃*D ≈ B̃*Matrix(D)
        @test D*D̃ ≈ D̃*D

        @test D\D ≈ I
        @test D\B ≈ Matrix(D)\B
        @test B\D ≈ B\Matrix(D)
        @test D\B̃ ≈ Matrix(D)\B̃
        @test B̃\D ≈ B̃\Matrix(D)
        @test D\D̃ ≈ D̃\D
        @test B̃/D ≈ B̃/Matrix(D)

        @testset "Diagonal * Bidiagonal/Tridiagonal with structured diags" begin
            n = size(D,1)
            BU = Bidiagonal(map(MyVector, (rand(n), rand(n-1)))..., :U)
            MBU = MyMatrix(BU)
            BL = Bidiagonal(map(MyVector, (rand(n), rand(n-1)))..., :L)
            MBL = MyMatrix(BL)
            S = SymTridiagonal(map(MyVector, (rand(n), rand(n-1)))...)
            MS = MyMatrix(S)
            T = Tridiagonal(map(MyVector, (rand(n-1), rand(n), rand(n-1)))...)
            MT = MyMatrix(T)
            DA, BUA, BLA, SA, TA = map(Array, (D, BU, BL, S, T))
            if VERSION >= v"1.11"
                @test D * BU ≈ DA * BUA
                @test BU * D ≈ BUA * DA
                @test D * MBU ≈ DA * BUA
                @test MBU * D ≈ BUA * DA
                @test D * BL ≈ DA * BLA
                @test BL * D ≈ BLA * DA
                @test D * MBL ≈ DA * BLA
                @test MBL * D ≈ BLA * DA
            end
            if VERSION >= v"1.12.0-DEV.824"
                @test D * S ≈ DA * SA
                @test D * MS ≈ DA * SA
                @test D * T ≈ DA * TA
                @test D * MT ≈ DA * TA
                @test S * D ≈ SA * DA
                @test MS * D ≈ SA * DA
                @test T * D ≈ TA * DA
                @test MT * D ≈ TA * DA
            end
        end
    end

    @testset "Adj/Trans" begin
        A = MyMatrix(randn(5,5))
        T = UpperTriangular(randn(5,5))
        D = Diagonal(MyVector(randn(5)))

        @test D * A' ≈ D * A.A'
        @test A' * D ≈ A.A' * D

        @test T * D ≈ T .* D.diag'
        @test D * T ≈ D.diag .* T
        @test A * Adjoint(T) ≈ A.A * Adjoint(T)
        @test A * Transpose(T) ≈ A.A * Transpose(T)
        @test Adjoint(T) * A ≈ Adjoint(T) * A.A
        @test Transpose(T) * A ≈ Transpose(T) * A.A
        @test Transpose(T)A' ≈ Adjoint(T)A' ≈ Adjoint(T)Transpose(A) ≈ Transpose(T)Transpose(A)
        @test Transpose(A)Adjoint(T) ≈ A'Adjoint(T) ≈ A'Transpose(T) ≈ Transpose(A)Transpose(T)

        @test Zeros(5)' * A ≡ Zeros(5)'
        @test transpose(Zeros(5)) * A ≡ transpose(Zeros(5))

        @test A' * Zeros(5) ≡ Zeros(5)
        @test Zeros(3,5) * A ≡ Zeros(3,5)
        @test A * Zeros(5,3) ≡ Zeros(5,3)
        @test A' * Zeros(5,3) ≡ Zeros(5,3)
        @test transpose(A) * Zeros(5,3) ≡ Zeros(5,3)
        @test A' * Zeros(5) ≡ Zeros(5)
        @test transpose(A) * Zeros(5) ≡ Zeros(5)

        b = MyVector(randn(5))
        @test A' * b ≈ A' * b.A

        @test b'*Zeros(5) == 0
        @test transpose(b)*Zeros(5) == 0
        @test_throws DimensionMismatch b'*Zeros(6)
        @test_throws DimensionMismatch transpose(b)*Zeros(6)
    end

    @testset "AbstractQ" begin
        A = MyMatrix(randn(5,5))
        Q = qr(randn(5,5)).Q
        @test Q'*A ≈ Q'*A.A
        @test Q*A ≈ Q*A.A
        @test A*Q ≈ A.A*Q
        @test A*Q' ≈ A.A*Q'
    end


    @testset "concat" begin
        a = MyVector(randn(5))
        A = MyMatrix(randn(5,6))

        @test [a; a] == [Array(a); a] == [a; Array(a)]
        @test [A; A] == [Array(A); A] == [A; Array(A)]

        @test [Array(a) A A] == [Array(a) Array(A) Array(A)]
        @test [a A A] == [Array(a) Array(A) Array(A)]
        @test [a Array(A) A] == [Array(a) Array(A) Array(A)]
    end

    @testset "dot" begin
        a = MyVector(randn(5))
        @test dot(a, Zeros(5)) ≡ dot(Zeros(5), a) == 0.0
    end

    @testset "layout_getindex scalar" begin
        A = MyMatrix(rand(5,4))
        @test A[6] == ArrayLayouts.layout_getindex(A,6) == A[1,2]
        @test A[1,2] == A[CartesianIndex(1,2)] == ArrayLayouts.layout_getindex(A,CartesianIndex(1,2))
    end

    @testset "structured axes" begin
        A = MyMatrix(rand(5,5))
        x = MyVector(rand(5))
        @test axes(Symmetric(A)) == axes(Symmetric(view(A,1:5,1:5))) == axes(A)
        @test axes(UpperTriangular(A)) == axes(UpperTriangular(view(A,1:5,1:5))) == axes(A)
        @test axes(Diagonal(x)) == axes(Diagonal(Vector(x)))
    end

    @testset "adjtrans *" begin
        A = MyMatrix(rand(5,5))
        x = MyVector(rand(5))

        @test x'A ≈ transpose(x)A ≈ x.A'A.A
        @test x'A' ≈ x'transpose(A) ≈ transpose(x)A' ≈ transpose(x)transpose(A) ≈ x.A'A.A'
    end

    @testset "Zeros mul" begin
        A = MyMatrix(randn(5,5))
        @test A*Zeros(5) ≡ Zeros(5)
        @test A*Zeros(5,2) ≡ Zeros(5,2)
        @test Zeros(1,5) * A ≡ Zeros(1,5)
        @test Zeros(5)' * A ≡ Zeros(5)'
        @test transpose(Zeros(5)) * A ≡ transpose(Zeros(5))

        @test Zeros(5)' * A * (1:5) == 
            transpose(Zeros(5)) * A * (1:5) ==
            (1:5)' * A * Zeros(5) ==
            transpose(1:5) * A * Zeros(5) ==
            Zeros(5)' * A * Zeros(5) ==
            transpose(Zeros(5)) * A * Zeros(5) == 0.0
    end

    @testset "triple *" begin
        D = Diagonal(1:5)
        y = MyVector(randn(5))
        @test (1:5)' *  D * y ≈ transpose(1:5) *  D * y ≈ (1:5)' * D * y.A
        @test y' * D * y ≈ transpose(y) * D * y ≈ y.A' * D * y.A
        @test y' * D * (1:5) ≈ y.A' * D * (1:5)
        @test y' * Diagonal(y) isa Adjoint
        @test transpose(y) * Diagonal(y) isa Transpose

        @test y' * D isa Adjoint
        @test transpose(y) * D isa Transpose

        @test Zeros(5)' * D * y == transpose(Zeros(5)) * D * y == 0.0
        @test y' * D * Zeros(5) == transpose(y) * D * Zeros(5) == 0.0
        @test Zeros(5)' * Diagonal(y) ≡ Zeros(5)'
        @test transpose(Zeros(5)) * Diagonal(y) ≡ transpose(Zeros(5))
        @test Zeros(5)' * Diagonal(y) * y == 0.0
        @test transpose(Zeros(5)) * Diagonal(y) * y == 0.0
        @test y' * Diagonal(y) * Zeros(5) == 0.0
        @test transpose(y) * Diagonal(y) * Zeros(5) == 0.0
        @test Zeros(5)' * Diagonal(y) * Zeros(5) == 0.0
        @test transpose(Zeros(5)) * Diagonal(y) * Zeros(5) == 0.0
    end

    @testset "rmul with lazy and Diagonal" begin
        D = Diagonal(1:5)
        y = MyVector(randn(5))
        @test mul(view(y', :, 1:5), D) isa Adjoint
        @test mul(view(transpose(y), :, 1:5), D) isa Transpose
    end

    @testset "Tri * Tri" begin
        A = MyMatrix(randn(3,3))
        @test UpperTriangular(A) * LowerTriangular(A) ≈ UpperTriangular(A.A) * LowerTriangular(A.A)
        @test UpperTriangular(A) * UnitUpperTriangular(A) ≈ UpperTriangular(A.A) * UnitUpperTriangular(A.A)
        @test UpperTriangular(A') * UnitUpperTriangular(A) ≈ UpperTriangular(A.A') * UnitUpperTriangular(A.A)
        @test UpperTriangular(A) * UnitUpperTriangular(A') ≈ UpperTriangular(A.A) * UnitUpperTriangular(A.A')
        @test UpperTriangular(A') * UnitUpperTriangular(A') ≈ UpperTriangular(A.A') * UnitUpperTriangular(A.A')
    end

    @testset "mul! involving a triangular" begin
        A = MyMatrix(rand(4,4))
        UA = UpperTriangular(A)
        MA = Matrix(A)
        MUA = Matrix(UA)
        B = rand(4,4)
        UB = UpperTriangular(B)
        @test mul!(zeros(4,4), A, UB) ≈ MA * UB
        @test mul!(ones(4,4), A, UB, 2, 2) ≈ 2 * MA * UB .+ 2
        @test mul!(zeros(4,4), UA, B) ≈ MUA * B
        @test mul!(ones(4,4), UA, B, 2, 2) ≈ 2 * MUA * B .+ 2
        @test mul!(zeros(4,4), UA, UB) ≈ MUA * UB
        @test mul!(ones(4,4), UA, UB, 2, 2) ≈ 2 * MUA * UB .+ 2
        @test mul!(zeros(4,4), UB, A) ≈ UB * MA
        @test mul!(ones(4,4), UB, A, 2, 2) ≈ 2 * UB * MA .+ 2
        @test mul!(zeros(4,4), UB, UA) ≈ UB * MUA
        @test mul!(ones(4,4), UB, UA, 2, 2) ≈ 2 * UB * MUA .+ 2
        @test mul!(zeros(4,4), B, UA) ≈ B * MUA
        @test mul!(ones(4,4), B, UA, 2, 2) ≈ 2 * B * MUA .+ 2
        @test mul!(zeros(4,4), A, UA) ≈ MA * MUA
        @test mul!(ones(4,4), A, UA, 2, 2) ≈ 2 * MA * MUA .+ 2
        @test mul!(zeros(4,4), UA, A) ≈ MUA * MA
        @test mul!(ones(4,4), UA, A, 2, 2) ≈ 2 * MUA * MA .+ 2
        @test mul!(zeros(4,4), UA, UA) ≈ MUA * MUA
        @test mul!(ones(4,4), UA, UA, 2, 2) ≈ 2 * MUA * MUA .+ 2


        v = rand(4)
        @test mul!(zeros(4), UA, v) ≈ MUA * v
        @test mul!(ones(4), UA, v, 2, 2) ≈ 2 * MUA * v .+ 2
    end

    if isdefined(LinearAlgebra, :copymutable_oftype)
        @testset "copymutable_oftype" begin
            A = MyMatrix(randn(3,3))
            @test LinearAlgebra.copymutable_oftype(A, BigFloat) == A
        end
    end

    @testset "sparse" begin
        @testset "MyVector" begin
            V = MyVector([1:4;])
            V2 = MyVector(2*[1:4;])
            S = 2*sparse(V)
            copyto!(V, S)
            @test S == V2
            V = MyVector([1:4;])
            copyto!(view(V, :), S)
            @test S == V2
        end
        @testset "MyMatrix" begin
            M = MyMatrix(reshape([1:4;], 2, 2))
            M2 = MyMatrix(reshape(2*[1:4;], 2, 2))
            S = 2*sparse(M)
            copyto!(M, S)
            @test S == M2
            M = MyMatrix(reshape([1:4;], 2, 2))
            copyto!(view(M, :, :), S)
            @test S == M2
        end
    end

    @testset "mul! with subarrays" begin
            A = MyMatrix(randn(3,3))
            V = view(A, 1:3, 1:3)
            B = randn(3,3)
            x = randn(3)
            D = Diagonal(1:3)

            @test mul!(similar(B), V, B) ≈ A * B
            @test mul!(similar(B), B, V) ≈ B * A
            @test mul!(similar(B), V, V) ≈ A^2
            @test mul!(similar(B), V, A) ≈ A * A
            @test mul!(similar(B), A, V) ≈ A * A
            @test mul!(MyMatrix(randn(3,3)), A, V) ≈ A * A
            @test mul!(similar(x), V, x) ≈ V * x

            @test mul!(copy(B), V, B, 2.0, 3.0) ≈ 2A * B + 3B
            @test mul!(copy(B), B, V, 2.0, 3.0) ≈ 2B * A + 3B
            @test mul!(copy(B), V, V, 2.0, 3.0) ≈ 2A^2 + 3B
            @test mul!(copy(B), V, A, 2.0, 3.0) ≈ 2A * A + 3B
            @test mul!(copy(B), A, V, 2.0, 3.0) ≈ 2A * A + 3B
            @test mul!(MyMatrix(copy(B)), A, V, 2.0, 3.0) ≈ 2A * A + 3B
            @test mul!(copy(x), V, x, 2.0, 3.0) ≈ 2A * x + 3x

            @test D * V == D * A == D * A.A
            @test V * D == A * D == A.A  * D
            @test mul!(copy(B), D, A, 2.0, 3.0) ≈ mul!(copy(B), D, V, 2.0, 3.0) ≈ 2D * A + 3B
            @test mul!(copy(B), A, D, 2.0, 3.0) ≈ mul!(copy(B), V, D, 2.0, 3.0) ≈ 2A * D + 3B
    end

    @testset "Diagonal * Q" begin
        A = randn(10,5)
        Q̃ = qr(A).Q
        Q = LinearAlgebra.QRCompactWYQ(MyMatrix(Q̃.factors), Q̃.T)
        D = Diagonal(1:10)
        @test ArrayLayouts.mul(D, Q) ≈ ArrayLayouts.mul(D, Q̃) ≈ D*Q
        @test ArrayLayouts.mul(Q, D) ≈ ArrayLayouts.mul(Q̃, D) ≈ Q*D
    end

    @testset "QR" begin
        A = randn(10,5)
        Q̃ = qr(A).Q
        Q = LinearAlgebra.QRCompactWYQ(MyMatrix(Q̃.factors), Q̃.T)
        b = randn(10)
        B = randn(10,3)
        @test Q*b ≈ Q̃*b ≈ lmul!(Q, copy(b)) ≈ lmul!(Q, MyVector(copy(b)))
        @test Q'b ≈ Q̃'b ≈ lmul!(Q', copy(b)) ≈ lmul!(Q', MyVector(copy(b)))
        @test Q*B ≈ Q̃*B ≈ lmul!(Q, copy(B)) ≈ lmul!(Q, MyMatrix(copy(B)))
        @test Q'B ≈ Q̃'B ≈ lmul!(Q', copy(B)) ≈ lmul!(Q', MyMatrix(copy(B)))

        @test B'Q ≈ B'Q̃ ≈ rmul!(copy(B'),Q) ≈ rmul!(MyMatrix(copy(B')), Q)
        @test B'Q' ≈ B'Q̃' ≈ rmul!(copy(B'),Q') ≈ rmul!(MyMatrix(copy(B')), Q')

        @test_throws ErrorException lmul!(Q, big.(B))
        @test_throws ErrorException lmul!(Q', big.(B))
        @test_throws ErrorException rmul!(big.(B'),Q')
        @test_throws ErrorException rmul!(big.(B'),Q)

        @test_broken Q*b ≈ lmul!(Q̃, MyVector(copy(b))) # broken due to commented out code for ambiguity with MatrixFactorizations.jl
    end
end

struct MyUpperTriangular{T} <: AbstractMatrix{T}
    A::UpperTriangular{T,Matrix{T}}
end

MyUpperTriangular{T}(::UndefInitializer, n::Int, m::Int) where T = MyUpperTriangular{T}(UpperTriangular(Array{T}(undef, n, m)))
MyUpperTriangular(A::AbstractMatrix{T}) where T = MyUpperTriangular{T}(UpperTriangular(Matrix{T}(A)))
Base.convert(::Type{MyUpperTriangular{T}}, A::MyUpperTriangular{T}) where T = A
Base.convert(::Type{MyUpperTriangular{T}}, A::MyUpperTriangular) where T = MyUpperTriangular(convert(AbstractArray{T}, A.A))
Base.convert(::Type{MyUpperTriangular}, A::MyUpperTriangular)= A
Base.convert(::Type{AbstractArray{T}}, A::MyUpperTriangular) where T = MyUpperTriangular(convert(AbstractArray{T}, A.A))
Base.convert(::Type{AbstractMatrix{T}}, A::MyUpperTriangular) where T = MyUpperTriangular(convert(AbstractArray{T}, A.A))
Base.convert(::Type{MyUpperTriangular{T}}, A::AbstractArray{T}) where T = MyUpperTriangular{T}(A)
Base.convert(::Type{MyUpperTriangular{T}}, A::AbstractArray) where T = MyUpperTriangular{T}(convert(AbstractArray{T}, A))
Base.convert(::Type{MyUpperTriangular}, A::AbstractArray{T}) where T = MyUpperTriangular{T}(A)
Base.getindex(A::MyUpperTriangular, kj...) = A.A[kj...]
Base.getindex(A::MyUpperTriangular, ::Colon, j::AbstractVector) = MyUpperTriangular(A.A[:,j])
Base.setindex!(A::MyUpperTriangular, v, kj...) = setindex!(A.A, v, kj...)
Base.size(A::MyUpperTriangular) = size(A.A)
Base.similar(::Type{MyUpperTriangular{T}}, m::Int, n::Int) where T = MyUpperTriangular{T}(undef, m, n)
Base.similar(::MyUpperTriangular{T}, m::Int, n::Int) where T = MyUpperTriangular{T}(undef, m, n)
Base.similar(::MyUpperTriangular, ::Type{T}, m::Int, n::Int) where T = MyUpperTriangular{T}(undef, m, n)
LinearAlgebra.factorize(A::MyUpperTriangular) = factorize(A.A)

MemoryLayout(::Type{MyUpperTriangular{T}}) where T = MemoryLayout(UpperTriangular{T,Matrix{T}})
triangulardata(A::MyUpperTriangular) = triangulardata(A.A)

@_layoutlmul MyUpperTriangular


@testset "MyUpperTriangular" begin
    A = randn(5,5)
    B = randn(5,5)
    x = randn(5)
    U = MyUpperTriangular(A)

    @test lmul!(U, copy(x)) ≈ U * x
    @test lmul!(U, copy(B)) ≈ U * B

    @test_skip lmul!(U,view(copy(B),collect(1:5),1:5)) ≈ U * B

    @test MyMatrix(A) / U ≈ A / U
    VERSION >= v"1.9" && @test U / MyMatrix(A) ≈ U / A
end

# Tests needed for InfiniteRandomArrays.jl (see https://github.com/DanielVandH/InfiniteRandomArrays.jl/issues/5) 
using ..InfiniteArrays

@testset "* for infinite layouts" begin
    tup = InfSymTridiagonal(), InfTridiagonal(), InfBidiagonal('U'),
        InfBidiagonal('L'),
        InfUnitUpperTriangular(), InfUnitLowerTriangular(),
        InfUpperTriangular(), InfLowerTriangular(), InfDiagonal();
    for (i, A) in enumerate(tup)
        A_up, A_lo = A isa Union{UpperTriangular,UnitUpperTriangular}, A isa Union{LowerTriangular,UnitLowerTriangular}
        for (j, B) in enumerate(tup)
            B_up, B_lo = B isa Union{UpperTriangular,UnitUpperTriangular}, B isa Union{LowerTriangular,UnitLowerTriangular}
            ((A_up && B_lo) || (A_lo && B_up)) && continue
            C = A * B 
            _C = [C[i, j] for i in 1:100, j in 1:100] # else we need to fix the C[1:100, 1:100] MethodError from _getindex(::Mul, ...). This is easier
            @test _C ≈ Matrix(A[1:100, 1:102]) * Matrix(B[1:102, 1:100])
        end
    end
end

@testset "disambiguation with FillArrays" begin
    v = [1,2,3]
    lv = MyVector(v)
    F = Fill(2, 3, 3)
    @test F * lv == F * v
    @test lv' * F == v' * F
    @test transpose(lv) * F == transpose(v) * F
end

end
