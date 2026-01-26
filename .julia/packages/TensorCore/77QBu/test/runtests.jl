using TensorCore
using LinearAlgebra
using Test

@testset "Ambiguities" begin
    if VERSION >= v"1.6-"
        @test isempty(detect_ambiguities(TensorCore))
    else
        @test isempty(detect_ambiguities(TensorCore, Base, Core, LinearAlgebra))
    end
end

@testset "tensor and hadamard" begin
    for T in (Int, Float32, Float64, BigFloat)
        a = [T[1, 2], T[-3, 7]]
        b = [T[5, 11], T[-13, 17]]
        @test map(⋅, a, b) == map(dot, a, b) == [27, 158]
        @test map(⊙, a, b) == map(hadamard, a, b) == [a[1].*b[1], a[2].*b[2]]
        @test map(⊗, a, b) == map(tensor, a, b) == [a[1]*transpose(b[1]), a[2]*transpose(b[2])]
        @test hadamard!(fill(typemax(Int), 2), T[1, 2], T[-3, 7]) == [-3, 14]
        @test tensor!(fill(typemax(Int), 2, 2), T[1, 2], T[-3, 7]) == [-3 7; -6 14]
    end

    @test_throws DimensionMismatch [1,2] ⊙ [3]
    @test_throws DimensionMismatch hadamard!([0, 0, 0], [1,2], [-3,7])
    @test_throws DimensionMismatch hadamard!([0, 0], [1,2], [-3])
    @test_throws DimensionMismatch hadamard!([0, 0], [1], [-3,7])
    @test_throws DimensionMismatch tensor!(Matrix{Int}(undef, 2, 2), [1], [-3,7])
    @test_throws DimensionMismatch tensor!(Matrix{Int}(undef, 2, 2), [1,2], [-3])

    u, v = [2+2im, 3+5im], [1-3im, 7+3im]
    @test u ⋅ v == conj(u[1])*v[1] + conj(u[2])*v[2]
    @test u ⊙ v == [u[1]*v[1], u[2]*v[2]]
    @test u ⊗ v == [u[1]*v[1] u[1]*v[2]; u[2]*v[1] u[2]*v[2]]
    @test hadamard(u, v) == u ⊙ v
    @test tensor(u, v)   == u ⊗ v
    dest = similar(u)
    @test hadamard!(dest, u, v) == u ⊙ v
    dest = Matrix{Complex{Int}}(undef, 2, 2)
    @test tensor!(dest, u, v) == u ⊗ v

    for (A, B, b) in (([1 2; 3 4], [5 6; 7 8], [5,6]),
                      ([1+0.8im 2+0.7im; 3+0.6im 4+0.5im],
                       [5+0.4im 6+0.3im; 7+0.2im 8+0.1im],
                       [5+0.6im,6+0.3im]))
        @test A ⊗ b == cat(A*b[1], A*b[2]; dims=3)
        @test A ⊗ B == cat(cat(A*B[1,1], A*B[2,1]; dims=3),
                           cat(A*B[1,2], A*B[2,2]; dims=3); dims=4)
    end

    A, B = reshape(1:27, 3, 3, 3), reshape(1:4, 2, 2)
    @test A ⊗ B == [a*b for a in A, b in B]

    # Adjoint/transpose is a dual vector, not an AbstractMatrix
    v = [1,2]
    @test_throws ErrorException v ⊗ v'
    @test_throws ErrorException v ⊗ transpose(v)
    @test_throws ErrorException v' ⊗ v
    @test_throws ErrorException transpose(v) ⊗ v
    @test_throws ErrorException v' ⊗ v'
    @test_throws ErrorException transpose(v) ⊗ transpose(v)
    @test_throws ErrorException v' ⊗ transpose(v)
    @test_throws ErrorException transpose(v) ⊗ v'
    @test_throws ErrorException A ⊗ v'
    @test_throws ErrorException A ⊗ transpose(v)

    # Docs comparison to `kron`
    v, w = [1,2,3], [5,7]
    @test kron(v,w) == vec(w ⊗ v)
    @test w ⊗ v == reshape(kron(v,w), (length(w), length(v)))
end

@testset "arrays of arrays" begin
    @testset "pseudo-complex" begin

        # This is a matrix representation of the complex numbers:
        jm = [0 1; -1 0]
        om = [1 0; 0 1]
        @test jm * jm == -om
        @test jm * om == om * jm == jm
        @test om * om == om

        @test jm' == -jm

        mreal(x::Matrix) = tr(x)/2
        jmag(x::Matrix) = tr(-x * jm)/2
        @test mreal((2om+3jm)^2) == real((2+3im)^2)
        @test jmag((2om+3jm)^2) == imag((2+3im)^2)

        v = [im, 2, 3+im]
        vm = [jm, 2om, 3om + jm]

        w = [4 + im, 5, 6+7im]
        wm = [4om + jm, 5om, 6om+7jm]
        @test real.(w) == mreal.(wm)
        @test imag.(w) == jmag.(wm)

        # hadamard & tensor
        @test real.(v ⊙ w) == mreal.(vm ⊙ wm)
        @test imag.(v ⊙ w) == jmag.(vm ⊙ wm)

        @test real.(v ⊗ w) == mreal.(vm ⊗ wm)
        @test imag.(v ⊗ w) == jmag.(vm ⊗ wm)

        # kron
        @test mreal.(kron(vm, wm)) == real.(kron(v, w))
        @test jmag.(kron(vm, wm)) == imag.(kron(v, w))

        @test vec(kron(v', w')') == kron(v, w)
        @test vec(kron(vm', wm')') == kron(vm, wm)

        # times, cross & dot
        outervw = v * w'
        @test outervw == mreal.(vm * wm') + im * jmag.(vm * wm')

        vcrossw = v × w
        @test vcrossw == mreal.(vm × wm) + im * jmag.(vm × wm)

        dvw = dot(v, w) # 36 + 11im
        @test dvw == v' * w
        @test dvw == sum(v[i]' * w[i] for i in 1:3)

        m_star = vm' * wm  # recursive adjoint is important here!
        if !(VERSION < v"1.5-")
            @test dvw == mreal(m_star) + im * jmag(m_star) # fails on 1.4
        end

        m_wrong = first(permutedims(vm) * wm)
        @test mreal(m_wrong) + im * jmag(m_wrong) == 20 + 31im

        m_sum = sum(vm[i]' * wm[i] for i in 1:3)
        @test dvw == mreal(m_sum) + im * jmag(m_sum)

        @test_broken dvw == dot(vm, wm) # 72, wtf?
        sum(dot(vm[i], wm[i]) for i in 1:3)
        # I think this behaviour was introduced here:
        # https://github.com/JuliaLang/julia/pull/27401
        # I don't actually see an argument for this, and there also seem to be no tests
        # of this in LinearAlgebra (although I have not changed & run them to check).

        # dot with kron
        kdot = kron(v, v)' * kron(w, w)
        kmdot = kron(vm, vm)' * kron(wm, wm)
        if !(VERSION < v"1.5-")
            @test kdot == mreal(kmdot) + im * jmag(kmdot)
        end

    end
    @testset "block matrix & vector" begin

        # A matrix of matrices, and a vector of vectors:
        M = [rand(Int8,2,2).+3im for i in 1:3, j in 1:4]
        N = [rand(Int8,2,2).+5im for i in 1:3, j in 1:4]

        V = [rand(Int8,2).+7im for i in 1:3]
        W = [rand(Int8,2).+11im for j in 1:4]

        # ordinary LinearAlgebra
        @test M * N' == (N * M')'

        @test M * W == (W' * M')' # recursive adjoint is important here!
        @test V' * M == (M' * V)'

        @test V' * V == dot(V,V) # here the behaviour of dot is not wrong
        if VERSION >= v"1.4"
            @test V' * M * W == dot(V,M,W)
        end

        # kron & TensorCore
        @test M ⊙ M == M .* M != map(⊙, M, M)

        @test first(M ⊗ N) == first(M) * first(N)
        @test size(M ⊗ N) == (3,4,3,4)

        @test_throws MethodError kron(W, W) # what should this do, if anything?
        @test_throws MethodError W ⊗ V

    end

end

@testset "boxdot" begin

    # Matrices and vectors
    A = [1 2+im; 3 4im]
    B = [5im 6; 7+im 8]
    c = [1, 2+3im]
    d = [3im, 4-5im]

    @test A ⊡ B == A * B
    @test A ⊡ c == A ⊡ c
    @test c ⊡ A == vec(transpose(c) * A)
    @test c ⊡ d == sum(c .* d)

    @test A' ⊡ B == A' * B
    @test A ⊡ B' == A * B'
    @test A' ⊡ B' == A' * B'

    # Dual vectors
    @test c' ⊡ d == dot(c, d)
    @test c ⊡ d' == conj(dot(c, d))
    @test c' ⊡ d' == conj(dot(c', d))

    @test transpose(c) ⊡ d == sum(c .* d)
    @test c ⊡ transpose(d) == sum(c .* d)
    @test transpose(c) ⊡ transpose(d) == sum(c .* d)
    @test transpose(c) ⊡ adjoint(d) == sum(c .* conj(d))
    @test adjoint(c) ⊡ transpose(d) == dot(c,d)

    @test A ⊡ c' isa Adjoint
    @test A ⊡ c' == (c ⊡ A')'
    @test c' ⊡ A isa Adjoint
    @test c' ⊡ A == (A' * c)'

    @test B' ⊡ c' isa Adjoint
    @test B' ⊡ c' == (c ⊡ B)'
    @test c' ⊡ B' isa Adjoint
    @test c' ⊡ B' == (B * c)'

    @test B' ⊡ transpose(c) isa Transpose
    @test B ⊡ transpose(c) == transpose(c ⊡ transpose(B))
    @test transpose(c) ⊡ B' isa Transpose
    @test transpose(c) ⊡ B == transpose(transpose(B) * c)

    # Higher-dimensional arrays
    E3 = cat(A, -B, dims=3)
    F4 = cat(E3, conj(E3 .+ 1), dims=4)
    E3adjoint = conj(permutedims(E3, (3,2,1)))
    E3lazy = PermutedDimsArray(permutedims(E3, (3,2,1)), (3,2,1))
    F4lazy = PermutedDimsArray(permutedims(F4, (4,3,2,1)), (4,3,2,1))
    @test E3lazy == E3
    @test F4lazy == F4

    @test E3 ⊡ A == reshape(reshape(E3, 4,2) * A, 2,2,2)
    @test size(A ⊡ E3) == (2, 2, 2)
    @test size(B ⊡ F4) == (2, 2, 2, 2)
    @test size(E3 ⊡ F4) == (2, 2, 2, 2, 2)
    @test E3lazy ⊡ A == E3 ⊡ A
    @test E3lazy ⊡ F4lazy == E3 ⊡ F4
    @test F4lazy ⊡ E3lazy == F4 ⊡ E3

    @test c ⊡ E3 == reshape(transpose(c) * reshape(E3, 2,4), 2,2)
    @test size(F4 ⊡ c) == (2, 2, 2)

    @test c' ⊡ E3 isa Matrix
    @test c' ⊡ E3 == conj(c) ⊡ E3
    @test c' ⊡ E3 == (E3adjoint ⊡ c)'
    @test transpose(c) ⊡ E3 == c ⊡ E3
    @test c ⊡ E3lazy == c ⊡ E3

    @test E3 ⊡ c' isa Matrix
    @test E3 ⊡ c' == E3 ⊡ conj(c)
    @test E3 ⊡ c' == (c ⊡ E3adjoint)'
    @test E3 ⊡ transpose(c) == E3 ⊡ c
    @test E3 ⊡ transpose(c) == transpose(c ⊡ PermutedDimsArray(E3, (3,2,1)))

    # Errors
    @test_throws DimensionMismatch ones(3) ⊡ ones(4)
    @test_throws DimensionMismatch ones(3) ⊡ ones(4)'
    @test_throws DimensionMismatch ones(3)' ⊡ ones(4)
    @test_throws DimensionMismatch ones(3)' ⊡ ones(4)'
    @test_throws DimensionMismatch ones(2,3) ⊡ ones(4)
    @test_throws DimensionMismatch ones(2,3) ⊡ ones(4)'
    @test_throws DimensionMismatch ones(2,3) ⊡ ones(4,5)

    # In-place
    @test boxdot!(similar(c), A, c) == A * c
    if VERSION >= v"1.3"
        @test boxdot!(similar(c), A, c, 100) == A * c * 100
        @test boxdot!(copy(c), B, d, 100, -5) == B * d * 100 .- 5 .* c
    end

    @test boxdot!(similar(c), A, c') == A * conj(c)
    @test boxdot!(similar(c,1), c, d') == [sum(c .* conj(d))]

    @test boxdot!(similar(c)', c', A) == c' * A
    @test boxdot!(similar(c,1,2), c', A) == c' * A

    @test boxdot!(similar(c,1), c', d) == [dot(c, d)]
end

@testset "_adjoint" begin
    A = [1 2+im; 3 4im]
    E3 = cat(A, -A, dims=3)
    F4 = cat(E3, conj(E3 .+ 1), dims=4)
    M3 = [rand(ComplexF32, 2,2) for _ in 1:3, _ in 1:4, _ in 1:5];

    @test TensorCore._adjoint(E3) == conj(permutedims(E3, (3,2,1)))
    @test TensorCore._transpose(F4) == permutedims(F4, (4,3,2,1))
    @test TensorCore._adjoint(A) == A'
    @test TensorCore._transpose(A) == transpose(A)

    @test TensorCore._adjoint(M3)[5,4,3] == adjoint(M3[3,4,5])
    @test TensorCore._transpose(M3)[1,2,3] == transpose(M3[3,2,1])

    @test TensorCore._adjoint(E3 ⊗ F4) == TensorCore._adjoint(F4) ⊗ TensorCore._adjoint(E3)
    @test TensorCore._adjoint(E3 ⊡ F4) == TensorCore._adjoint(F4) ⊡ TensorCore._adjoint(E3)
    @test TensorCore._adjoint(M3 ⊡ TensorCore._adjoint(M3)) ≈ M3 ⊡ TensorCore._adjoint(M3)
end
