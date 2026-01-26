correct = sqrt(π) * exp(-1/4) / 2
f(x) = sin(sqrt(x)) / sqrt(x)
α = 1/2

let
    x, w = laguerre(10, α)
    Q = sum(w .* f.(x))
    @test Q ≈ correct
end

let
    x, w = laguerre(10, α, left)
    X = x[2:end]
    W = w[2:end]
    Q = w[1] + sum(W .* f.(X))
    @test Q ≈ correct
end
