# SPDX-License-Identifier: MIT
# Author: Pierre Lairez

module ExactPredicates

export incircle, orient, closest, insphere, meet, sameside, rotation, parallelorder, intersectorder, lengthcompare

include("Codegen.jl")

using StaticArrays
using IntervalArithmetic

using .Codegen: group!, @genpredicate
import .Codegen: coord

export coord


# A few helper functions

const R2 = SVector{2}
const R3 = SVector{3}

ext(u :: R2, v :: R2) = u[1] * v[2] - u[2] * v[1]
inp(u :: R2, v :: R2) = u[1] * v[1] + u[2] * v[2]
inp(u :: R3, v :: R3) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3]

det(a,b,c,d) = a*d-b*c
det(a,b,c,d,e,f,g,h,i) = a*det(e,f,h,i) - d*det(b,c,h,i) + g*det(b,c,e,f)

@inline function det(m1, m5, m9, m13, m2, m6, m10, m14, m3, m7, m11, m15, m4, m8, m12, m16)
    p34 = m11*m16 - m12*m15
    p23 = m10*m15 - m11*m14
    p12 = m9 *m14 - m10*m13
    p13 = m9 *m15 - m11*m13
    p14 = m9 *m16 - m12*m13
    p24 = m10*m16 - m12*m14
    return (
          m1*(m6*p34 - m7*p24 + m8*p23)
        - m2*(m5*p34 - m7*p14 + m8*p13)
        + m3*(m5*p24 - m6*p14 + m8*p12)
        - m4*(m5*p23 - m6*p13 + m7*p12)
    )
end

abs2(u :: SVector) = inp(u, u)


# if x and y are integer in {-1, 0, 1}, then xor(x,y)==-2 is equivalent to x*y == -1.
opposite_signs(x :: Int, y :: Int) = (xor(x, y) == -2)



include("plane.jl")
include("space.jl")


end
