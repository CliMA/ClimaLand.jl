# TODO: respect ordering
module AtomixMetalExt

using Atomix: Atomix, IndexableRef
using Metal: Metal, MtlDeviceArray

const MtlIndexableRef{Indexable<:MtlDeviceArray} = IndexableRef{Indexable}

function Atomix.get(ref::MtlIndexableRef, order)
    error("not implemented")
end

function Atomix.set!(ref::MtlIndexableRef, v, order)
    error("not implemented")
end

@inline function Atomix.replace!(
    ref::MtlIndexableRef,
    expected,
    desired,
    success_ordering,
    failure_ordering,
)
    ptr = Atomix.pointer(ref)
    expected = convert(eltype(ref), expected)
    desired = convert(eltype(ref), desired)
    begin
        old = Metal.atomic_compare_exchange_weak_explicit(ptr, expected, desired)
    end
    return (; old = old, success = old === expected)
end


# CAS is needed for FP ops on ThreadGroup memory
@inline function Atomix.modify!(ref::IndexableRef{<:MtlDeviceArray{<:AbstractFloat, <:Any, Metal.AS.ThreadGroup}} , op::OP, x, order) where {OP}
    x = convert(eltype(ref), x)
    ptr = Atomix.pointer(ref)
    old = Metal.atomic_fetch_op_explicit(ptr, op, x)
    return old => op(old, x)
end

@inline function Atomix.modify!(ref::MtlIndexableRef, op::OP, x, order) where {OP}
    x = convert(eltype(ref), x)
    ptr = Atomix.pointer(ref)
    begin
        old = if op === (+)
            Metal.atomic_fetch_add_explicit(ptr, x)
        elseif op === (-)
            Metal.atomic_fetch_sub_explicit(ptr, x)
        elseif op === (&)
            Metal.atomic_fetch_and_explicit(ptr, x)
        elseif op === (|)
            Metal.atomic_fetch_or_explicit(ptr, x)
        elseif op === xor
            Metal.atomic_fetch_xor_explicit(ptr, x)
        elseif op === min
            Metal.atomic_fetch_min_explicit(ptr, x)
        elseif op === max
            Metal.atomic_fetch_max_explicit(ptr, x)
        else
            error("not implemented")
        end
    end
    return old => op(old, x)
end

end  # module AtomixMetalExt
