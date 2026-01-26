#
using SciMLOperators, Zygote, LinearAlgebra
using Random

using SciMLOperators
using SciMLOperators: AbstractSciMLOperator,
                      IdentityOperator, NullOperator,
                      AdjointOperator, TransposedOperator,
                      InvertedOperator, InvertibleOperator,
                      BatchedDiagonalOperator, AddedOperator, ComposedOperator,
                      AddedScalarOperator, ComposedScalarOperator, ScaledOperator,
                      has_mul, has_ldiv

Random.seed!(0)
n = 3
N = n * n
K = 12

t = rand()
u0 = rand(N, K)
ps = rand(N)

s = rand()
v = rand(N, K)
M = rand(N, N)
Mi = inv(M)

sca_update_func = (a, u, p, t) -> sum(p) * s
vec_update_func = (b, u, p, t) -> Diagonal(p) * v
mat_update_func = (A, u, p, t) -> Diagonal(p) * M
inv_update_func = (A, u, p, t) -> Mi * inv(Diagonal(p))
tsr_update_func = (A, u, p, t) -> reshape(p, n, n) |> copy

α = ScalarOperator(zero(Float32), update_func = sca_update_func)
L_dia = DiagonalOperator(zeros(N, K); update_func = vec_update_func)
L_mat = MatrixOperator(zeros(N, N); update_func = mat_update_func)
L_mi = MatrixOperator(zeros(N, N); update_func = inv_update_func)
L_aff = AffineOperator(L_mat, L_mat, zeros(N, K); update_func = vec_update_func)
L_sca = α * L_mat
L_inv = InvertibleOperator(L_mat, L_mi)
L_fun = FunctionOperator((v, u, p, t) -> Diagonal(p) * v, u0, u0; batch = true,
    op_inverse = (v, u, p, t) -> inv(Diagonal(p)) * v)

Ti = MatrixOperator(zeros(n, n); update_func = tsr_update_func)
To = deepcopy(Ti)
L_tsr = TensorProductOperator(To, Ti)

for (LType, L) in ((IdentityOperator, IdentityOperator(N)),
    (NullOperator, NullOperator(N)),
    (MatrixOperator, L_mat),
    (AffineOperator, L_aff),
    (ScaledOperator, L_sca),
    (InvertedOperator, InvertedOperator(L_mat)),
    (InvertibleOperator, L_inv),
    (BatchedDiagonalOperator, L_dia),
    (AddedOperator, L_mat + L_dia),
    (ComposedOperator, L_mat * L_dia),
    (TensorProductOperator, L_tsr),
    (FunctionOperator, L_fun),

    ## ignore wrappers
    # (AdjointOperator, AdjointOperator(rand(N,N) |> MatrixOperator) |> adjoint),
    # (TransposedOperator, TransposedOperator(rand(N,N) |> MatrixOperator) |> transpose),

    (ScalarOperator, α),
    (AddedScalarOperator, α + α),
    (ComposedScalarOperator, α * α))
    @assert L isa LType

    # Cache the operator for efficient application
    L_cached = cache_operator(L, u0)

    # Updated loss function using the new interface:
    # v is the action vector, u0 is the update vector
    loss_mul = function (p)
        v = Diagonal(p) * u0
        # Use new interface: L(v, u, p, t)
        w = L_cached(v, u0, p, t)
        l = sum(w)
    end

    loss_div = function (p)
        v = Diagonal(p) * u0

        # Update coefficients first, then apply inverse
        L_updated = update_coefficients(L_cached, u0, p, t)
        w = L_updated \ v

        l = sum(w)
    end

    @testset "$LType" begin
        l_mul = loss_mul(ps)
        g_mul = Zygote.gradient(loss_mul, ps)[1]

        if L isa NullOperator
            @test isa(g_mul, Nothing)
        else
            @test !isa(g_mul, Nothing)
        end

        if has_ldiv(L)
            l_div = loss_div(ps)
            g_div = Zygote.gradient(loss_div, ps)[1]

            @test !isa(g_div, Nothing)
        end
    end
end
