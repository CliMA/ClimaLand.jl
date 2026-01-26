# TODO: move this to UnsafeAtomics
julia_ordering_name(::UnsafeAtomics.Internal.LLVMOrdering{name}) where {name} = name
julia_ordering_name(::typeof(UnsafeAtomics.acquire_release)) = :acquire_release
julia_ordering_name(::typeof(UnsafeAtomics.sequentially_consistent)) =
    :sequentially_consistent

julia_syncscope_name(::UnsafeAtomics.Internal.LLVMSyncScope{name}) where {name} = name
julia_syncscope_name(::typeof(UnsafeAtomics.none)) = :system

include("atomics.jl")

@inline UnsafeAtomics.load(ptr::LLVMPtr, order::Ordering, sync::UnsafeAtomics.Internal.SyncScope) =
    atomic_pointerref(ptr, Val{julia_ordering_name(order)}(), Val{julia_syncscope_name(sync)}())

@inline function UnsafeAtomics.store!(ptr::LLVMPtr, x, order::Ordering, sync::UnsafeAtomics.Internal.SyncScope)
    atomic_pointerset(ptr, x, Val{julia_ordering_name(order)}(), Val{julia_syncscope_name(sync)}())
    return
end

mapop(op::OP) where {OP} = op
mapop(::typeof(UnsafeAtomics.right)) = right

@inline UnsafeAtomics.modify!(ptr::LLVMPtr, op::OP, x, order::Ordering, sync::UnsafeAtomics.Internal.SyncScope) where {OP} =
    atomic_pointermodify(ptr, mapop(op), x, Val{julia_ordering_name(order)}(), Val{julia_syncscope_name(sync)}())

@inline UnsafeAtomics.cas!(
    ptr::LLVMPtr,
    expected,
    desired,
    success_order::Ordering,
    failure_order::Ordering,
    sync::UnsafeAtomics.Internal.SyncScope,
) = atomic_pointerreplace(
    ptr,
    expected,
    desired,
    Val{julia_ordering_name(success_order)}(),
    Val{julia_ordering_name(failure_order)}(),
    Val{julia_syncscope_name(sync)}(),
)
