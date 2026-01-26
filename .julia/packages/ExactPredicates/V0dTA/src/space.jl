# SPDX-License-Identifier: MIT
# Author: Pierre Lairez


@genpredicate nogeneric function closest(p :: 3, q :: 3, a :: 3)
    qp = q - p
    pa = p - a
    qa = q - a

    Codegen.group!(qp...)
    Codegen.group!(pa..., qa...)

    inp(qp, pa+qa)
end
@doc """
    closest(p :: 3, q :: 3, a :: 3) -> Int

Return 1 if `a` is closer to `p` than to `q`.
Return –1 if `a` is closer to `q` than to `p`.
Return 0 is `a` is equally close to both.
""" closest(::NTuple{3, Float64},::NTuple{3, Float64},::NTuple{3, Float64})



@genpredicate function orient(p :: 3, q :: 3, r :: 3, a :: 3)
    pa = p - a
    qa = q - a
    ra = r - a

    Codegen.group!(pa...)
    Codegen.group!(qa...)
    Codegen.group!(ra...)

    det(pa..., qa..., ra...)
end
@doc """
    orient(p :: 3, q :: 3, r :: 3, a :: 3) -> Int

Consider the oriented plane on which the triangle `pqr` is positively oriented.

* Return 1 if `a` is below this plane.
* Return –1 if `a` is above this plane.
* Return 0 if `a` lies on this plane.

""" orient(::NTuple{3, Float64},::NTuple{3, Float64},::NTuple{3, Float64},::NTuple{3, Float64})



@genpredicate function insphere(p :: 3, q :: 3, r :: 3, s :: 3, a :: 3)
    p = p - a
    q = q - a
    r = r - a
    s = s - a

    Codegen.group!(p..., q..., r..., s...)

    det(p..., abs2(p), q..., abs2(q), r..., abs2(r), s..., abs2(s))
end
@doc """
    insphere(p :: 3, q :: 3, r :: 3, s :: 3, a :: 3)

* Return 1 if `a` is inside the circumscribed sphere defined by the four points `p`, `q`, `r` and `s`.
* Return –1 if `a` is outside.
* Return 0 is `a` lies on the sphere or if the four points are coplanar.

""" insphere(::NTuple{3, Float64},::NTuple{3, Float64},::NTuple{3, Float64},::NTuple{3, Float64},::NTuple{3, Float64})

