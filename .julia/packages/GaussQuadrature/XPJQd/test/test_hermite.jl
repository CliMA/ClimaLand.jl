correct = sqrt(π/sqrt(ℯ)) / 2
let
    x, w = hermite(10)
    Q = sum(w .* ( x .- 1.0) .* sin.(x))
    @test Q ≈ correct
end
