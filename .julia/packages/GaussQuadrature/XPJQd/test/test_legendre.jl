correct = ( pi/4 + atan(3.0) ) / 2
for endpt in [neither, both, left, right]
    x, w = legendre(40, endpt)
    if endpt in [left, both]
        @test x[1] == -1
    end
    if endpt in [right, both]
        @test x[end] == 1
    end
    Q = sum(w .* 1 ./ ( 1 .+ (2x .- 1).^2))
    @test Q ≈ correct
end
