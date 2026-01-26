
I = π * besseli(0.0, 1.0)
for endpt in [neither, left, right, both]
    x, w = chebyshev(10, 1, endpt)
    Q = sum(w .* exp.(x))
    @test I ≈ Q
end

I = π * besseli(1.0, 1.0)
for endpt in [neither, left, right, both]
    x, w = chebyshev(10, 2, endpt)
    Q = sum(w .* exp.(x))
    @test I ≈ Q
end
