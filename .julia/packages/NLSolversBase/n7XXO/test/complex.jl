@testset "complex" begin
	# Just catch that it doesn't error for now
	x0 = randn(ComplexF64, 3)
	@test_nowarn OnceDifferentiable(x -> sum(abs2, x), x0)
	@test_nowarn  NLSolversBase.value!(OnceDifferentiable(x -> sum(abs2, x), x0), x0)
	@test_nowarn  NLSolversBase.gradient!(OnceDifferentiable(x -> sum(abs2, x), x0), x0)
	@test_nowarn  NLSolversBase.value_gradient!(OnceDifferentiable(x -> sum(abs2, x), x0), x0)
end