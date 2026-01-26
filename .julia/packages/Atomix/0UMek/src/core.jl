@inline function Atomix.get(ref, order)
    ptr = Atomix.pointer(ref)
    root = Atomix.gcroot(ref)
    GC.@preserve root begin
        UnsafeAtomics.load(ptr, order)
    end
end

@inline function Atomix.set!(ref, v, order)
    v = Atomix.asstorable(ref, v)
    ptr = Atomix.pointer(ref)
    root = Atomix.gcroot(ref)
    GC.@preserve root begin
        UnsafeAtomics.store!(ptr, v, order)
    end
end

@inline function Atomix.replace!(ref, expected, desired, success_ordering, failure_ordering)
    expected = Atomix.asstorable(ref, expected)
    desired = Atomix.asstorable(ref, desired)
    ptr = Atomix.pointer(ref)
    root = Atomix.gcroot(ref)
    GC.@preserve root begin
        UnsafeAtomics.cas!(ptr, expected, desired, success_ordering, failure_ordering)
    end
end

@inline function Atomix.modify!(ref, op::OP, x, ord) where {OP}
    x = Atomix.asstorable(ref, x)
    ptr = Atomix.pointer(ref)
    root = Atomix.gcroot(ref)
    GC.@preserve root begin
        UnsafeAtomics.modify!(ptr, op, x, ord)
    end
end

Atomix.asstorable(ref, v) = convert(eltype(ref), v)
