"""
    area(poly)

Calculate the signed area via the shoelace formula. If the points are ordered
counterclockwise, the result will be positive. If the points are ordered clockwise,
the result will be negative.
"""
function area(poly::T) where T
    a = zero(eltype(eltype(T)))
    @inbounds for i in 1:length(poly)-1
        p1 = poly[i]
        p2 = poly[i+1]
        a += p1[1]*p2[2]-p2[1]*p1[2]
    end
    return a/2
end

function centroid(poly::T) where T
    ZT = zero(eltype(eltype(T)))
    a = ZT
    c = zero(eltype(poly))

    @inbounds for i in 1:length(poly)-1
        p = poly[i]
        n = poly[i+1]
        k = p[1]*n[2]-n[1]*p[2]
        a += k
        c = c .+ (p.+n).*k
    end
    a /= 2
    return c./(6*a)
end
