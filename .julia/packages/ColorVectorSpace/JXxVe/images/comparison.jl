using Colors

Base.mean(c1::C,c2::C) where {C<:Color} = convert(C, (convert(XYZ, c1) + convert(XYZ, c2))/2)  # colorimetrically "proper"

# RGB linear vector space averaging
function meanrgb(c1::C, c2::C) where C<:AbstractRGB
    C((red(c1)+red(c2))/2,
      (green(c1)+green(c2))/2,
      (blue(c1)+blue(c2))/2)
end

n = 10
colors = distinguishable_colors(n, Array(RGB{Float32},0))
m = 50 # number of repeated pixels in each block
grid = Array(eltype(colors), n*m, n*m)
z = zero(RGB{Float32})
for i = 1:n
    for j = i:n   # i == j is a sanity check
        prop = mean(colors[i], colors[j])
        vs = meanrgb(colors[i], colors[j])
        for l = 1:m, k = 1:m
            grid[(i-1)*m+k,(j-1)*m+l] = grid[(j-1)*m+k,(i-1)*m+l] = k > l ? prop : vs
        end
    end
end
# Put a boundary around the diagonal elements
for i = 1:n
    for k = 1:m
        grid[(i-1)*m+k,(i-1)*m+1] = grid[(i-1)*m+k,(i-1)*m+m] = z
        grid[(i-1)*m+1,(i-1)*m+k] = grid[(i-1)*m+m,(i-1)*m+k] = z
    end
end

using Images

gridu8 = map(Clamp(RGB{N0f8}), grid)
img = Image(gridu8, spatialorder=["x","y"])
imwrite(img, "comparison.png")
