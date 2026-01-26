correct = ℯ + 1/ℯ
for endpt in [neither, left, right, both]
    x, w = jacobi(10, 0.0, 1.0, endpt)
    Q = sum(w .* exp.(x))
    @test Q ≈ correct
end
