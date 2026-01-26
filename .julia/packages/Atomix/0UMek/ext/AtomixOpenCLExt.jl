# TODO: respect ordering
module AtomixOpenCLExt

using Atomix: Atomix, IndexableRef
using OpenCL: SPIRVIntrinsics, CLDeviceArray

const CLIndexableRef{Indexable<:CLDeviceArray} = IndexableRef{Indexable}

function Atomix.get(ref::CLIndexableRef, order)
    error("not implemented")
end

function Atomix.set!(ref::CLIndexableRef, v, order)
    error("not implemented")
end

@inline function Atomix.replace!(
    ref::CLIndexableRef,
    expected,
    desired,
    success_ordering,
    failure_ordering,
)
    ptr = Atomix.pointer(ref)
    expected = convert(eltype(ref), expected)
    desired = convert(eltype(ref), desired)
    begin
        old = SPIRVIntrinsics.atomic_cmpxchg!(ptr, expected, desired)
    end
    return (; old = old, success = old === expected)
end

@inline function Atomix.modify!(ref::CLIndexableRef, op::OP, x, order) where {OP}
    x = convert(eltype(ref), x)
    ptr = Atomix.pointer(ref)
    begin
        old = if op === (+)
            SPIRVIntrinsics.atomic_add!(ptr, x)
        elseif op === (-)
            SPIRVIntrinsics.atomic_sub!(ptr, x)
        elseif op === (&)
            SPIRVIntrinsics.atomic_and!(ptr, x)
        elseif op === (|)
            SPIRVIntrinsics.atomic_or!(ptr, x)
        elseif op === xor
            SPIRVIntrinsics.atomic_xor!(ptr, x)
        elseif op === min
            SPIRVIntrinsics.atomic_min!(ptr, x)
        elseif op === max
            SPIRVIntrinsics.atomic_max!(ptr, x)
        else
            error("not implemented")
        end
    end
    return old => op(old, x)
end

end  # module AtomixOpenCLExt

