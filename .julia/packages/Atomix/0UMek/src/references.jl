"""
    Atomix.IndexableRef{Indexable}

An object `IndexableRef(data, indices)` represents a *reference* to the location
`data[indices...]`.

A reference object supports `Atomix.pointer`, `Atomix.asstorable`, and `eltype`:

```jldoctest
julia> using Atomix

julia> a = [111, 222, 333];

julia> ref = Atomix.IndexableRef(a, (1,));

julia> Atomix.pointer(ref) === pointer(a, 1)
true

julia> Atomix.asstorable(ref, 1.0)
1

julia> eltype(ref)
$Int
```

To customize the behavior of atomic updates of an `Indexable`, define the following methods:

* `Atomix.get(ref::Atomix.IndexableRef{Indexable}, order) -> v::eltype(ref)`
* `Atomix.set!(ref::Atomix.IndexableRef{Indexable}, v, order)`
* `Atomix.replace!(ref::Atomix.IndexableRef{Indexable}, expected, desired, success_order,
  failure_order) -> (; old, success)`.
* `Atomix.modify!(ref::Atomix.IndexableRef{Indexable}, op, x, order) -> (old => new)`

The ordering arguments (`order`, `success_order`, and `failure_order`) are one of:

* `Atomix.monotonic`
* `Atomix.acquire`
* `Atomix.release`
* `Atomix.acquire_release` (`Atomix.acq_rel`)
* `Atomix.sequentially_consistent` (`Internal.seq_cst`)
"""
Atomix.IndexableRef

Atomix.IndexableRef(data, indices) = IndexableRef(data, indices)

struct IndexableRef{Indexable,Indices} <: Atomix.IndexableRef{Indexable}
    data::Indexable
    indices::Indices
end

"""
    referenceable(mutable) -> r

Return a *referenceable* object `r` whose index and field access yield a *reference* to the
corresponding memory location instead of the value stored.

TODO: record field access

TODO: record nested field and index access to update nested element (e.g.,
array-of-large-structs) atomically?

Note: The idea is identical to Referenceables.jl but Atomix.jl has a minimal implementation
in its internal.
"""
function referenceable end

referenceable(xs::AbstractArray) = ReferenceableArray(xs)

Base.eltype(::Type{<:IndexableRef{Indexable,<:Any}}) where {Indexable} = eltype(Indexable)

const IntIndexableRef{N,Indexable} = IndexableRef{Indexable,NTuple{N,Int}}

# TODO: Use Referenceables.jl?
# TODO: Don't subtype `AbstractArray`?  It is crazy that `xs[1]` and `xs[1, 1]` where
# `xs::ReferenceableArray{2}` return different values.  Or maybe just define `==` on the
# refs?
struct ReferenceableArray{N,Data<:AbstractArray{<:Any,N}} <:
       AbstractArray{Union{IntIndexableRef{1,Data},IntIndexableRef{N,Data}},N}
    data::Data
end

Base.size(a::ReferenceableArray) = size(a.data)
Base.IndexStyle(::Type{<:ReferenceableArray{<:Any,Data}}) where {Data} =
    Base.IndexStyle(Data)

@propagate_inbounds function Base.getindex(a::ReferenceableArray, i::Int)
    @boundscheck checkbounds(a.data, i)
    return IndexableRef(a.data, (i,))::IntIndexableRef{1}
end

@propagate_inbounds function Base.getindex(
    a::ReferenceableArray{N},
    I::Vararg{Int,N},
) where {N}
    @boundscheck checkbounds(a.data, I...)
    return IndexableRef(a.data, I)::IntIndexableRef{N}
end

@inline Atomix.pointer(ref::IntIndexableRef{1}) = pointer(ref.data, ref.indices[1])

@inline function Atomix.pointer(ref::IntIndexableRef)
    i = LinearIndices(ref.data)[ref.indices...]
    return pointer(ref.data, i)
end

Atomix.gcroot(ref) = ref
Atomix.gcroot(ref::IndexableRef) = ref.data
