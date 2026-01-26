module InterfaceTests

using Contour, Test

function setup()
    nx, ny = 10, 10

    xs = sort!(rand(nx))
    ys = sort!(rand(ny))
    zs = rand(nx, ny)

    xs, ys, zs
end

xs, ys, zs = setup()

cs = @inferred contours(xs, ys, zs)
for c in levels(cs)
    for l in lines(c)
        x, y = coordinates(l)
        @assert typeof(x) == typeof(y) == Vector{Float64}
        xy = vertices(l)
        @test xy isa Vector{Tuple{Float64,Float64}}
    end
end

end # InterfaceTests
