# SPDX-License-Identifier: MIT
# Author: Pierre Lairez



@genpredicate function orient(u :: 2, v :: 2, w :: 2)
    u = u - w
    v = v - w

    group!(u...)
    group!(v...)

    ext(u, v)
end
@doc """

    orient(p :: 2, q :: 2, r :: 2) -> Int

* Return 1 if `r` is on the left of the oriented line defined by `p` and `q`.
* Return –1 if `r` is on the right.
* Return 0 if `r` is on the line or if `p == q`.

""" orient(::NTuple{2, Float64}, ::NTuple{2, Float64}, ::NTuple{2, Float64})


@genpredicate function incircle(p :: 2, q :: 2, r :: 2, a :: 2)

    qp = q - p
    rp = r - p
    ap = a - p
    aq = a - q
    rq = r - q

    group!(qp...)
    group!(ap...)
    group!(rp...)
    group!(rq..., aq...)

    ext(qp, ap)*inp(rp, rq) - ext(qp, rp)*inp(ap, aq)
end
@doc """
    incircle(a :: 2, b :: 2, c :: 2, p :: 2) -> Int


Assume that `a`, `b` and `c` define a counterclockwise triangle.

* Return 1 if `p` is strictly inside the circumcircle of this triangle.
* Return –1 if `p` is outside.
* Return 0 if `p` is on the circle.

If the triangle is oriented clockwise, the signs are reversed.
If `a`, `b` and `c` are collinear, this degenerate to an orientation test.

If two of the four arguments are equal, return 0.

""" incircle(::NTuple{2, Float64}, ::NTuple{2, Float64}, ::NTuple{2, Float64}, ::NTuple{2, Float64})



@genpredicate function closest(p :: 2, q :: 2, a :: 2)
    qp = q - p
    pa = p - a
    qa = q - a

    group!(qp...)
    group!(pa..., qa...)

    inp(qp, pa+qa)
end
@doc """
    closest(p :: 2, q :: 2, a :: 2) -> Int

* Return 1 if `a` is closer to `p` than to `q`.
* Return –1 if `a` is closer to `q` than to `p`.
* Return 0 is `a` is equally close to both.
""" closest(::NTuple{2, Float64}, ::NTuple{2, Float64}, ::NTuple{2, Float64})



"""
    sameside(p :: 2, a :: 2, b :: 2)

Assume that the three arguments are collinear, on some line L

* Return 1 if `a` and `b` are on the same side of `p` on L
* Return -1 if `a` and `b` are on different sides
* Return 0 if `a == p` or `b == p`.
"""
function sameside(p :: Tuple{Float64, Float64}, a :: Tuple{Float64, Float64}, b :: Tuple{Float64, Float64})
    if a < p && b < p || a > p && b > p
        return 1
    elseif a < p && b > p || a > p && b < p
        return -1
    else
        return 0
    end
end

function sameside(p, a, b)
    sameside(coord(p), coord(a), coord(b))
end

"""
    meet(p :: 2, q :: 2, a :: 2, b :: 2)

* Return 1 if the open line segments `(p, q)` and `(a, b)` meet in a single point.
* Return 0 if the the closed line segments `[p, q]` and `[a, b]` meet in one or several points.
* Return –1 otherwise.
"""
function meet(p, q, a, b)
    pqa = orient(p, q, a)
    pqb = orient(p, q, b)
    abp = orient(a, b, p)
    abq = orient(a, b, q)

    if opposite_signs(pqa, pqb) && opposite_signs(abp, abq)
        return 1
    elseif (pqa != 0 && pqa == pqb)  || (abq != 0  && abp == abq)
        return -1
    elseif pqa == 0 && pqb == 0
        # all four points are collinear
        if sameside(p, a, b) == 1 && sameside(q, a, b) == 1 && sameside(a, p, q) == 1 && sameside(b, p, q) == 1
            return -1
        else
            return 0
        end
    else
        return 0
    end
end


@genpredicate function parallelorder(a :: 2, b :: 2, p :: 2, q :: 2)
    δ = b - a
    qp = q - p

    group!(δ...)
    group!(qp...)
    ext(δ, qp)
end
@doc """
    parallelorder(a :: 2, b :: 2, p :: 2, q :: 2) -> Int

Consider the oriented line defined by `a` and `b`
and the parallel lines passing through `p` and `q` respectively, with the same orientation.

* return 1 if the line passing through `p` is left of the line passing through `q`.
* return -1 in the reverse situation.
* return 0 if `a` and `b` are equal or if the parallel lines passing through `p` and `q` are equal.

This is a robust version of to `orient(b-a, q-p, 0)`.
Note also that `orient(a, b, c) == parallelorder(a, b, a, c)`.
""" parallelorder(:: NTuple{2, Float64}, :: NTuple{2, Float64}, :: NTuple{2, Float64}, :: NTuple{2, Float64})




"""
    rotation(pts :: AbstractVector)

Gives the [rotation number](https://en.wikipedia.org/wiki/Winding_number#Turning_number) of the polygonal path defined by the elements of `pts`.
"""
function rotation(pts :: AbstractVector{T}) where T
    u = rand(SVector{2, Float64})
    origin = (0.0, 0.0)

    n = length(pts)
    @assert n >= 3

    pts = copy(pts)
    push!(pts, pts[1], pts[2])

    r = 0
    o1 = parallelorder(origin, u, pts[1], pts[2])
    for i in 1:n
        o2 = parallelorder(origin, u, pts[i+1], pts[i+2])
        if opposite_signs(o1, o2)
            ro = parallelorder(pts[i], pts[i+1], pts[i+1], pts[i+2])
            @assert ro != 0
            if ro > 0 && o1 < 0
                r += 1
            elseif ro < 0 && o1 > 0
                r -= 1
            end
        end
        o1 = o2
    end

    return r
end



@genpredicate function intersectorder(a :: 2, b :: 2, pa :: 2, pb :: 2, qa :: 2, qb :: 2)
    pp = pb - pa
    qq = qb - qa
    δ = b - a
    p0 = pa - a
    q0 = qa - a

    group!(pp...)
    group!(qq...)
    group!(δ...)
    group!(p0..., q0...)

    det(ext(p0, pp), ext(q0, qq),
        ext(δ, pp), ext(δ, qq))
end
@doc """
    intersectorder(a :: 2, b :: 2, pa :: 2, pb :: 2, qa :: 2, qb :: 2) -> Int

Consider the oriented line *L* defined by `a` and `b`, the line *P* defined by
`pa` and `pb` and the line *Q* defined by `qa` and `qb`.

Assumes that `parallelorder(a, b, pa, pb)` and `parallelorder(a, b, qa, qb)` have the same sign.
    Otherwise, the result has the opposite sign.

* return -1 if the intersection of *P* with *L* comes before the intersection of *Q* with *L*, following the orientation of *L*.
* return 1 in the reverse situation
* return 0 in case of equality or degeneracy.
""" intersectorder(:: NTuple{2, Float64}, :: NTuple{2, Float64}, :: NTuple{2, Float64}, :: NTuple{2, Float64}, :: NTuple{2, Float64}, :: NTuple{2, Float64})



@genpredicate function lengthcompare(a :: 2, b :: 2, c :: 2, d :: 2)
    group!(a..., b..., c..., d...)
    return inp(a - b, a - b) - inp(c - d, c - d)
end

@doc """
    lengthcompare(a :: 2, b :: 2, c :: 2, d :: 2) -> Int

* return -1 if the distance between `a` and `b` is smaller than the distance between `c` and `d`
* return 1 in the reverse situation
* return 0 in case the distance are equal
"""
lengthcompare(:: NTuple{2, Float64}, :: NTuple{2, Float64}, :: NTuple{2, Float64}, :: NTuple{2, Float64})
