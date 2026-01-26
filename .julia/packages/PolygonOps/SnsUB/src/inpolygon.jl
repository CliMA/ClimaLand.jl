

abstract type MembershipCheckAlgorithm end

"""
Algorithm by Hao and Sun (2018):
https://doi.org/10.3390/sym10100477
"""
struct HaoSun <: MembershipCheckAlgorithm end

"""
Hormann-Agathos (2001) Point in Polygon algorithm:
https://doi.org/10.1016/S0925-7721(01)00012-8
"""
struct HormannAgathos <: MembershipCheckAlgorithm end

"""
    inpolygon(p, poly)
    inpolygon(p, poly, [::MembershipCheckAlgorithm])
    inpolygon(p, poly, [::MembershipCheckAlgorithm]; in = 1, on = -1, out = 0)

check the membership of `p` in `poly` where `poly` is an `AbstractVector` of `AbstractVector`s.
`poly` should have the first and last elements equal.

*Returns:*
- in = 1
- on = -1
- out = 0

*MembershipCheckAlgorithm:*

- `HaoSun()`
- `HormannAgathos()`

Default is `HaoSun()` as it has the best performance and works invariant of winding order and self-intersections.
However the HaoSun algorithm is new and bugs may be possible. The classic HormannAgathos algorithm
is provided, however it is sensitive to self-intersections and winding order so may produce different results.

*Custom return values:*

The return logic can be customized to return alternate values and types by specify the `in`, `on`, and `out` keywords.
For example, to treat the 'on' and 'in' cases the same and return a `Bool`:
`inpolygon(p, poly, in=true, on=true, out=false)`

Algorithm by Hao and Sun (2018):
https://doi.org/10.3390/sym10100477

Hormann-Agathos (2001) Point in Polygon algorithm:
https://doi.org/10.1016/S0925-7721(01)00012-8
"""
function inpolygon(p, poly; kwargs...)
    inpolygon(p, poly, HaoSun(); kwargs...)
end

function inpolygon(p, poly, ::Union{HaoSun,Type{HaoSun}}; in::T=1, on::T=-1, out::T=0) where {T}
    k = 0

    xp = p[1]
    yp = p[2]

    validate_poly(poly)

    PT = eltype(p)

    @inbounds for i in UnitRange(firstindex(poly),lastindex(poly)-1)
        v1 = poly[i][2] - yp
        v2 = poly[i+1][2] - yp

        if v1 < zero(PT) && v2 < zero(PT) || v1 > zero(PT) && v2 > zero(PT)
            continue
        end

        u1 = poly[i][1] - xp
        u2 = poly[i+1][1] - xp

        f = (u1 * v2) - (u2 * v1)

        if v2 > zero(PT) && v1 <= zero(PT)
            if f > zero(PT)
                k += 1
            elseif iszero(f)
                return on
            end
        elseif v1 > zero(PT) && v2 <= zero(PT)
            if f < zero(PT)
                k += 1
            elseif iszero(f)
                return on
            end
        elseif iszero(v2) && v1 < zero(PT)
            iszero(f) && return on
        elseif iszero(v1) && v2 < zero(PT)
            iszero(f) && return on
        elseif iszero(v1) && iszero(v2)
            if u2 <= zero(PT) && u1 >= zero(PT)
                return on
            elseif u1 <= zero(PT) && u2 >= zero(PT)
                return on
            end
        end
    end

    iszero(k % 2) && return out
    return in
end

function detq(q1,q2,r)
    @inbounds (q1[1]-r[1])*(q2[2]-r[2])-(q2[1]-r[1])*(q1[2]-r[2])
end


function inpolygon(r, poly, ::Union{HormannAgathos,Type{HormannAgathos}};
                   in::T=1, on::T=-1, out::T=0) where{T}
    c = false

    validate_poly(poly)

    @inbounds for i in UnitRange(firstindex(poly),lastindex(poly)-1)
        q1 = poly[i]
        q2 = poly[i+1]
        if q1 == r
            # throw(VertexException())
            return on
        end
        if q2[2] == r[2]
            if q2[1] == r[1]
                #throw(VertexException())
                return on
            elseif (q1[2] == r[2]) && ((q2[1] > r[1]) == (q1[1] < r[1]))
                #throw(EdgeException())
                return on
            end
        end
        if (q1[2] < r[2]) != (q2[2] < r[2]) # crossing
            if q1[1] >= r[1]
                if q2[1] > r[1]
                    c = !c
                elseif (detq(q1,q2,r) > 0) == (q2[2] > q1[2]) # right crossing
                    c = !c
                end
            elseif q2[1] > r[1]
                if (detq(q1,q2,r) > 0) == (q2[2] > q1[2]) # right crossing
                    c = !c
                end
            end
        end
    end
    return c ? in : out
end
