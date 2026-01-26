# TODO: respect ordering
module AtomixoneAPIExt

using Atomix: Atomix, IndexableRef
using oneAPI: oneAPI, oneDeviceArray

const oneIndexableRef{Indexable<:oneDeviceArray} = IndexableRef{Indexable}

function Atomix.get(ref::oneIndexableRef, order)
    error("not implemented")
end

function Atomix.set!(ref::oneIndexableRef, v, order)
    error("not implemented")
end

@inline function Atomix.replace!(
    ref::oneIndexableRef,
    expected,
    desired,
    success_ordering,
    failure_ordering,
)
    ptr = Atomix.pointer(ref)
    expected = convert(eltype(ref), expected)
    desired = convert(eltype(ref), desired)
    begin
        old = oneAPI.atomic_cmpxchg!(ptr, expected, desired)
    end
    return (; old = old, success = old === expected)
end

@inline function Atomix.modify!(ref::oneIndexableRef, op::OP, x, order) where {OP}
    x = convert(eltype(ref), x)
    ptr = Atomix.pointer(ref)
    begin
        old = if op === (+)
            oneAPI.atomic_add!(ptr, x)
        elseif op === (-)
            oneAPI.atomic_sub!(ptr, x)
        elseif op === (&)
            oneAPI.atomic_and!(ptr, x)
        elseif op === (|)
            oneAPI.atomic_or!(ptr, x)
        elseif op === xor
            oneAPI.atomic_xor!(ptr, x)
        elseif op === min
            oneAPI.atomic_min!(ptr, x)
        elseif op === max
            oneAPI.atomic_max!(ptr, x)
        else
            error("not implemented")
        end
    end
    return old => op(old, x)
end

end  # module AtomixoneAPIExt
