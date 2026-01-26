I1 = 1 - exp(-1)
f1(x) = (1-x)*exp(-x) # Gradshteyn 4.351(1)

for endpt in [neither, both, left, right]
    x, w = logweight(10, 0, endpt)
    @test I1 ≈ sum(w .* f1.(x))
    x, w = logweight(10, 0.0, endpt)
    @test I1 ≈ sum(w .* f1.(x))
end

I2 = (π/3)^2 * 2 * sqrt(3)
f2(x) = (1-x^2) / (1+x^3) # Gradshteyn 4.255(1), p=3/2, q=3
for endpt in [neither, both, left, right]
    x, w = logweight(10, -1/2, endpt)
    @test I2 ≈ sum(w .* f2.(x))
end
