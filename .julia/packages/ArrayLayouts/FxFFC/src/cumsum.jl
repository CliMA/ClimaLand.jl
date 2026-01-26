"""
    RangeCumsum(range)

represents the cumsum of a `AbstractRange`.
"""
struct RangeCumsum{T, RR<:AbstractRange{T}} <: LayoutVector{T}
    range::RR
end

size(c::RangeCumsum) = size(c.range)
axes(c::RangeCumsum) = axes(c.range)

Base.parent(r::RangeCumsum) = r.range

function Base.promote_rule(::Type{RangeCumsum{T1, R1}}, ::Type{RangeCumsum{T2, R2}}) where {T1,T2,R1,R2}
    R = promote_type(R1, R2)
    RangeCumsum{promote_type(T1, T2), R}
end

==(a::RangeCumsum, b::RangeCumsum) = a.range == b.range
BroadcastStyle(::Type{<:RangeCumsum{<:Any,RR}}) where RR = BroadcastStyle(RR)

function _half_prod(a::Integer, b::Integer)
    iseven(a) ? (a÷2) * b : a * (b÷2)
end
function _onethird_prod(a::Integer, b::Integer)
    mod(a, 3) == 0 ? (a÷3) * b : a * (b÷3)
end

function _getindex(r::AbstractRange{<:Real}, k)
    v = first(r)
    s = step(r)
    # avoid overflow, if possible
    k * v + s * _half_prod(k, k-1)
end
Base.@propagate_inbounds _getindex(r::AbstractRange, k) = sum(r[range(firstindex(r), length=k)])

Base.@propagate_inbounds function getindex(c::RangeCumsum{<:Any,<:AbstractRange}, k::Integer)
    @boundscheck checkbounds(c, k)
    r = c.range
    _getindex(r, k-firstindex(r)+1)
end

Base.@propagate_inbounds getindex(c::RangeCumsum, kr::OneTo) = RangeCumsum(c.range[kr])

Base.@propagate_inbounds view(c::RangeCumsum, kr::OneTo) = c[kr]

first(r::RangeCumsum) = first(r.range)
last(r::RangeCumsum) = sum(r.range)
diff(r::RangeCumsum) = r.range[firstindex(r)+1:end]
isempty(r::RangeCumsum) = isempty(r.range)

function Base.sum(r::RangeCumsum{<:Real})
    N = length(r)
    v = first(r)
    s = step(r.range)
    # avoid overflow, if possible
    halfnnp1 = _half_prod(N, N+1)
    v * halfnnp1 + s * _onethird_prod(halfnnp1, N-1)
end

union(a::RangeCumsum{<:Any,<:OneTo}, b::RangeCumsum{<:Any,<:OneTo}) =
    RangeCumsum(OneTo(max(last(a.range), last(b.range))))

sort!(a::RangeCumsum{<:Any,<:Base.OneTo}) = a
sort(a::RangeCumsum{<:Any,<:Base.OneTo}) = a
function sort(a::RangeCumsum{<:Any,<:AbstractUnitRange{<:Integer}})
    r = parent(a)
    #= A RangeCumsum with a range (-a:b) for positive (a,b) may be viewed as a concatenation
    of two components: (-a:a) and (a+1:b). The second is strictly increasing.
    We only sort the first component, and concatenate the second to the result.
    =#
    a1 = RangeCumsum(r[firstindex(r):searchsortedlast(r, -first(r))])
    a2 = RangeCumsum(r[searchsortedfirst(r, -first(r)+1):end])
    vcat(sort!(Vector(a1)), a2)
end
Base.issorted(a::RangeCumsum{<:Any,<:Base.OneTo}) = true
function Base.issorted(a::RangeCumsum{<:Any,<:AbstractUnitRange{<:Integer}})
    r = parent(a)
    r2 = r[firstindex(r):searchsortedlast(r, zero(eltype(r)))]
    # at max one negative value is allowed
    length(r2) <= 1 + (last(r) >= 0)
end

struct _InitialValue end

_reduce_empty(init) = init
_reduce_empty(::_InitialValue) = throw(ArgumentError("RangeCumsum must be non-empty"))

function Base.minimum(a::RangeCumsum{<:Any, <:OneTo}; init = _InitialValue())
    isempty(a) && return _reduce_empty(init)
    first(a)
end
function Base.maximum(a::RangeCumsum{<:Any, <:OneTo}; init = _InitialValue())
    isempty(a) && return _reduce_empty(init)
    last(a)
end
function Base.maximum(a::RangeCumsum{<:Any, <:AbstractUnitRange{<:Integer}}; init = _InitialValue())
    isempty(a) && return _reduce_empty(init)
    r = parent(a)
    if -first(r) in r
        r2 = r[searchsortedfirst(r, -first(r)+1):end]
        max(zero(eltype(r)), sum(r2))
    else
        max(first(r), sum(r))
    end
end
function Base.minimum(a::RangeCumsum{<:Any, <:AbstractUnitRange{<:Integer}}; init = _InitialValue())
    isempty(a) && return _reduce_empty(init)
    r = parent(a)
    if zero(eltype(r)) in r
        r2 = r[firstindex(r):searchsortedlast(r, zero(eltype(r)))]
        min(sum(r2), zero(eltype(r)))
    else
        min(first(r), sum(r))
    end
end

convert(::Type{RangeCumsum{T,R}}, r::RangeCumsum) where {T,R} = RangeCumsum{T,R}(convert(R, r.range))

function Broadcast.broadcasted(::typeof(-), r::RangeCumsum)
    RangeCumsum(.-r.range)
end
function Broadcast.broadcasted(::typeof(*), x::Number, r::RangeCumsum)
    RangeCumsum(x * r.range)
end
function Broadcast.broadcasted(::typeof(*), r::RangeCumsum, x::Number)
    RangeCumsum(r.range * x)
end

Base.show(io::IO, r::RangeCumsum) = print(io, RangeCumsum, "(", r.range, ")")
