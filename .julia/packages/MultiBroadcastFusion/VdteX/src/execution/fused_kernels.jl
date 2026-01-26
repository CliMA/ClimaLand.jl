@make_type FusedMultiBroadcast
@make_fused fused_direct FusedMultiBroadcast fused_direct
@make_fused fused_assemble FusedMultiBroadcast fused_assemble
@make_get_fused fused_direct FusedMultiBroadcast get_fused_direct

struct MBF_CPU end
struct MBF_CUDA end
device(x::AbstractArray) = MBF_CPU()

function Base.copyto!(fmb::FusedMultiBroadcast)
    # Since we intercept Base.copyto!, we have not yet
    # called Base.Broadcast.instantiate (as this is done
    # in materialize, which has been stripped away), so,
    # let's call it here.
    fmb′ = FusedMultiBroadcast(
        map(fmb.pairs) do p
            Pair(p.first, Base.Broadcast.instantiate(p.second))
        end,
    )
    (; pairs) = fmb′ # (Pair(dest1, bc1),Pair(dest2, bc2),...)
    dest = first(pairs).first
    fused_copyto!(fmb′, device(dest))
end

Base.@propagate_inbounds function rcopyto_at!(
    pair::Pair,
    i::Vararg{T},
) where {T}
    dest, src = pair.first, pair.second
    @inbounds dest[i...] = src[i...]
    return nothing
end
Base.@propagate_inbounds function rcopyto_at!(
    pairs::Tuple,
    i::Vararg{T},
) where {T}
    rcopyto_at!(first(pairs), i...)
    rcopyto_at!(Base.tail(pairs), i...)
end
Base.@propagate_inbounds rcopyto_at!(
    pairs::Tuple{<:Any},
    i::Vararg{T},
) where {T} = rcopyto_at!(first(pairs), i...)
@inline rcopyto_at!(pairs::Tuple{}, i::Vararg{T}) where {T} = nothing

# This is better than the baseline.
function fused_copyto!(fmb::FusedMultiBroadcast, ::MBF_CPU)
    (; pairs) = fmb
    destinations = map(x -> x.first, pairs)
    ei = if eltype(destinations) <: Vector
        eachindex(destinations...)
    else
        eachindex(IndexCartesian(), destinations...)
    end
    for (dest, bc) in pairs
        @inbounds @simd ivdep for i in ei
            dest[i] = bc[i]
        end
    end
    return destinations
end


# This should, in theory be better, but it seems like inlining is
# failing somewhere.
# function fused_copyto!(fmb::FusedMultiBroadcast, ::MBF_CPU)
#     (; pairs) = fmb
#     destinations = map(x -> x.first, pairs)
#     ei = if eltype(destinations) <: Vector
#         eachindex(destinations...)
#     else
#         eachindex(IndexCartesian(), destinations...)
#     end
#     @inbounds @simd ivdep for i in ei
#         MBF.rcopyto_at!(pairs, i)
#     end
# end
