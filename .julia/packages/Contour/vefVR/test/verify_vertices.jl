module VerticesTests

using Contour, Test
using Base.MathConstants: π, φ
using LinearAlgebra: Diagonal
using OffsetArrays
using StaticArrays

# Setup test axes that will be shared among the tests

# Shift the axes so that they do not line up with
# integer values
Δ = 0.01
X = collect(0:Δ:4) .+ π
Y = collect(0:Δ:3) .+ φ

# TEST CASE 1
#
# f(x,y) = x^2 + y^2
#
Z = [(x^2 + y^2)::Float64 for x in X, y in Y]
h = rand() * (maximum(Z) - minimum(Z)) + minimum(Z)

contourlevels = Contour.contour(X, Y, Z, h; VT=SVector{2,Float64})
for line in contourlevels.lines
    # Contour vertices lie on a circle around the origin
    for v in line.vertices
        @test isapprox(v[1]^2 + v[2]^2, h, atol=0.1Δ)
    end

    # coordinates() returns the correct values
    xs, ys = coordinates(line)
    xs .== [v[1] for v in line.vertices]
    ys .== [v[2] for v in line.vertices]
end

# TEST CASE 1.5
#
# Same as test case 1, but with a shifted center
#
h = 1
x0, y0 = (.5, -1.5)
Z = Float64[(x - x0)^2 + (y - y0)^2 for x in X, y in Y]

contourlevels = Contour.contour(X, Y, Z, h, VT=SVector{2,Float64})
for line in contourlevels.lines
    for v in line.vertices
        @test isapprox((v[1] - x0)^2 + (v[2] - y0)^2, h, atol=0.01Δ)
    end
end

# TEST CASE 2
#
# Check that ambigious cells (5, 10) are handled correctly
# Case 5: z_center > h
X = [1.0, 2.0]
Y = [2.0, 1.0]

Z = float([1 0;
           0 1])
h = 0.1

lines = Contour.contour(X, Y, Z, h).lines
@test length(lines) == 2

for line in lines
    @test length(line.vertices) == 2
    d = line.vertices[2] .- line.vertices[1]
    @test d[2] / d[1] ≈ -1.0
end

# Case 5: z_center < h
Z = float([1 0;
           0 1])
h = 0.9

lines = Contour.contour(X, Y, Z, h).lines

@test length(lines) == 2

for line in lines
    @test length(line.vertices) == 2
    d = line.vertices[2] .- line.vertices[1]
    @test d[2] / d[1] ≈ 1.0
end

# Case 10: z_center > h
Z = float([0 1;
           1 0])
h = 0.1

lines = Contour.contour(X, Y, Z, h).lines
@test length(lines) == 2

for line in lines
    @test length(line.vertices) == 2
    d = line.vertices[2] .- line.vertices[1]
    @test d[2] / d[1] ≈ 1.0
end

# Case 10: z_center < h

Z = float([0 1;
           1 0])
h = 0.9

lines = Contour.contour(X, Y, Z, h, VT=SVector{2, Float64}).lines
@test length(lines) == 2

for line in lines
    @test length(line.vertices) == 2
    d = line.vertices[2] .- line.vertices[1]
    @test d[2] / d[1] ≈ -1.0
end

# Test curvilinear coordinates
θ = range(0.0, stop=2π,length=100)
R = range(1.0, stop=2.0, length=100)
ζ = ComplexF64[r*exp(im*ϕ) for ϕ in θ, r in R]
x, y, z = real.(ζ), imag.(ζ), abs.(ζ)

h = 1 + rand()
xs, ys = coordinates(contour(x, y, z, h, VT=SVector{2, Float64}).lines[1])
@test all(xs.^2 + ys.^2 .≈ h^2)


# Test offset arrays
offset_x, offset_y = -10, 27
z = cumsum(cumsum(randn(20,20); dims=1); dims=2)
zoff = OffsetArray(z, offset_x, offset_y)

x, y = axes(z)
xoff, yoff = axes(zoff)
curves = Contour.contour(x,y,z,0.5, VT=SVector{2, Float64})
curves_off = Contour.contour(xoff, yoff, zoff, 0.5,VT=SVector{2, Float64})

# sort offset and non-offset curves to the same order
offset = SVector(offset_x, offset_y)
lns = sort(Contour.lines(curves); by=c->sum(sum.(c.vertices.+[offset])))
lns_off = sort(Contour.lines(curves_off); by=c->sum(sum.(c.vertices)))

# verify that each line matches a possibly circularly shifted or reversed
# line from the offset array
opencurve(a) = first(a) == last(a) ? a[1:end-1] : a
for (c1, c2) in zip(lns, lns_off)
    o1 = opencurve(c1.vertices) .+ [offset]
    o2 = opencurve(c2.vertices)
    m = length(o1)
    @test m == length(o2)
    cshifts = [circshift(o2,i) for i=1:m]
    @test any(o1 ≈ c || o1 ≈ reverse(c) for c in cshifts)
end



# Test Known Bugs

# Issue #12
x = float(collect(1:3));
y = copy(x);
z = Diagonal(ones(3));
contours(x, y, z)

# Test handling of saddle points
#
# f(x,y) = x^2 - y^2
#
Δ = 0.01
X = collect(-3:Δ:3)
Y = collect(-3:Δ:3)

Z = [(x^2 - y^2)::Float64 for x in X, y in Y]
h = rand() * (maximum(Z) - minimum(Z)) + minimum(Z)

contourlevels = Contour.contour(X, Y, Z, h) #, VT=SVector{2, Float64})
for line in contourlevels.lines
    # Contour vertices lie on a circle around the origin
    for v in line.vertices
        @test isapprox(v[1]^2 - v[2]^2, h, atol=0.1Δ)
    end

    # coordinates() returns the correct values
    xs, ys = coordinates(line)
    xs .== [v[1] for v in line.vertices]
    ys .== [v[2] for v in line.vertices]
end


# Test range API
#
# f(x,y) = x^2 - y^2
#
Δ = 0.01
X = -3:Δ:3
Y = -3:Δ:3

Z = [(x^2 - y^2)::Float64 for x in X, y in Y]
h = rand() * (maximum(Z) - minimum(Z)) + minimum(Z)

contourlevels = Contour.contour(X, Y, Z, h)
for line in contourlevels.lines
    # Contour vertices lie on a circle around the origin
    for v in line.vertices
        @test isapprox(v[1]^2 - v[2]^2, h, atol=0.1Δ)
    end

    # coordinates() returns the correct values
    xs, ys = coordinates(line)
    xs .== [v[1] for v in line.vertices]
    ys .== [v[2] for v in line.vertices]
end

# Test that closed contours are identified correctly
# when ambiguous cells are involved

Z = float([1 1 1 1 1 1
           1 0 0 1 1 1
           1 0 0 1 1 1
           1 1 1 0 0 1
           1 1 1 0 0 1
           1 1 1 1 1 1])

X = Y = collect(0:0.2:1)
h = 0.75
contourlevels = Contour.contour(X, Y, Z, h)

# There should be only one closed contour
@test length(contourlevels.lines) == 1

# Test contour location on a realistic dataset
include("testdata.jl")
cts = Contour.contours(x, y, z)
@test length(cts.contours) == 10
cts_ct = (8, 8, 8, 8, 126, 7, 5, 5, 5, 4)

# Length need to be sorted as hashing might change see #62
lines_ct = sort!.([[138, 220, 469, 138, 469, 208, 143, 143],
            [220, 475, 140, 210, 146, 475, 140, 146],
            [222, 481, 140, 214, 481, 145, 140, 145],
            [228, 485, 214, 142, 485, 142, 146, 146],
            [9, 7, 515, 501, 9, 9, 9, 7, 7, 9, 7, 15, 7, 18, 39, 9, 7, 7, 9, 9, 7, 7, 7, 7, 9, 7, 18, 9, 7, 7, 9, 9, 35, 9, 9, 9, 7, 7, 9, 35, 9, 7, 9, 7, 9, 9, 7, 7, 7, 9, 7, 7, 7, 7, 9, 7, 9, 7, 9, 7, 7, 9, 9, 9, 7, 9, 7, 7, 9, 7, 7, 9, 7, 9, 7, 7, 7, 7, 7, 9, 7, 7, 7, 9, 7, 7, 7, 7, 9, 7, 7, 7, 7, 7, 7, 7, 9, 7, 9, 9, 7, 7, 7, 7, 7, 7, 9, 7, 7, 7, 7, 7, 7, 7, 9, 7, 9, 7, 7, 7, 7, 9, 7, 7, 7, 7],
            [5, 17, 14, 34, 34, 14, 5],
            [29, 12, 29, 15, 12],
            [26, 10, 13, 26, 10],
            [11, 11, 23, 11, 21],
            [18, 19, 7, 7]])
for i in eachindex(cts_ct)
    @test length(cts.contours[i].lines) == cts_ct[i]
    @test all(lines_ct[i] .== sort!([length(c.vertices) for c in cts.contours[i].lines]))
end


# support non-float z
using StatsBase

N = 10000
x = randn(N)
y = randn(N)
H = fit(Histogram, (x, y), closed = :left)
contours(midpoints.(H.edges)..., H.weights)

# Integer support/Conversion
contours(-5:5, -5:5, (-5:5)*(-5:5)')
contours(-5:5, -5:5, (-5:5)*(-5:5)', VT=NTuple{2,Float16})

end
