# TODO: respect ordering
module AtomixCUDAExt

using Atomix: Atomix, IndexableRef
using CUDA: CUDA, CuDeviceArray

const CuIndexableRef{Indexable<:CuDeviceArray} = IndexableRef{Indexable}

function Atomix.get(ref::CuIndexableRef, order)
    error("not implemented")
end

function Atomix.set!(ref::CuIndexableRef, v, order)
    error("not implemented")
end

@inline function Atomix.replace!(
    ref::CuIndexableRef,
    expected,
    desired,
    success_ordering,
    failure_ordering,
)
    ptr = Atomix.pointer(ref)
    expected = convert(eltype(ref), expected)
    desired = convert(eltype(ref), desired)
    begin
        old = CUDA.atomic_cas!(ptr, expected, desired)
    end
    return (; old = old, success = old === expected)
end

@inline function Atomix.modify!(ref::CuIndexableRef, op::OP, x, order) where {OP}
    x = convert(eltype(ref), x)
    ptr = Atomix.pointer(ref)
    begin
        old = if op === (+)
            CUDA.atomic_add!(ptr, x)
        elseif op === (-)
            CUDA.atomic_sub!(ptr, x)
        elseif op === (&)
            CUDA.atomic_and!(ptr, x)
        elseif op === (|)
            CUDA.atomic_or!(ptr, x)
        elseif op === xor
            CUDA.atomic_xor!(ptr, x)
        elseif op === min
            CUDA.atomic_min!(ptr, x)
        elseif op === max
            CUDA.atomic_max!(ptr, x)
        else
            error("not implemented")
        end
    end
    return old => op(old, x)
end

end  # module AtomixCUDAExt
