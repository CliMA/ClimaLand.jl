# TODO: Support Symbol ordering

"""
    Atomix.get(ref, order) -> x
    Atomix.get(ref) -> x

Atomically load the value `x` stored in `ref` with ordering `order`.  The default ordering
`Atomix.sequentially_consistent` is used when not specified.

# Examples
```jldoctest
julia> using Atomix

julia> a = [111, 222, 333];

julia> ref = Atomix.IndexableRef(a, (1,));

julia> Atomix.get(ref)
111
```
"""
Atomix.get

"""
    Atomix.set!(ref, x, order)
    Atomix.set!(ref, x)

Atomically store the value `x` in `ref` with ordering `order`.  The default ordering
`Atomix.sequentially_consistent` is used when not specified.

# Examples
```jldoctest
julia> using Atomix

julia> a = [111, 222, 333];

julia> ref = Atomix.IndexableRef(a, (1,));

julia> Atomix.set!(ref, 123);

julia> a[1]
123
```
"""
Atomix.set!

"""
    Atomix.modify!(ref, op, x, order) -> (old => new)
    Atomix.modify!(ref, op, x) -> (old => new)

Atomically update `ref` from stored value `old` to `new = op(old, x)` with ordering `order`
(default: `Atomix.sequentially_consistent`).  Return a pair `old => new`.

# Examples
```jldoctest
julia> using Atomix

julia> a = [111, 222, 333];

julia> ref = Atomix.IndexableRef(a, (1,));

julia> Atomix.modify!(ref, +, 123)
111 => 234
```
"""
Atomix.modify!

"""
    Atomix.swap!(ref, new, order) -> old
    Atomix.swap!(ref, new) -> old

Swap the `old` stored in `ref` with the `new` value and establish the memory ordering
`order` (default: `Atomix.sequentially_consistent`).

Notes for implementers: `Atomix.swap!(ref, new, order)` is defined as `Atomix.modify!(ref,
Atomix.right, x, order)`.  Thus, only `Atomix.modify!` has to be defined.

# Examples
```jldoctest
julia> using Atomix

julia> a = [111, 222, 333];

julia> ref = Atomix.IndexableRef(a, (1,));

julia> Atomix.swap!(ref, 123)
111

julia> a[1]
123
```
"""
Atomix.swap!

"""
    Atomix.replace!(ref, expected, desired, success_order, fail_order) -> (; old, success)
    Atomix.replace!(ref, expected, desired, order) -> (; old, success)
    Atomix.replace!(ref, expected, desired) -> (; old, success)

Atomically replace the value stored in `ref` to `desired` if `expected` is stored.  A named
tuple `(; old::eltype(ref), success::Bool)` is returned.

# Examples
```jldoctest
julia> using Atomix

julia> a = [111, 222, 333];

julia> ref = Atomix.IndexableRef(a, (1,));

julia> Atomix.replace!(ref, 111, 123)
(old = 111, success = true)
```
"""
Atomix.replace!

@inline Atomix.get(ref) = Atomix.get(ref, seq_cst)
@inline Atomix.set!(ref, x) = Atomix.set!(ref, x, seq_cst)

@inline Atomix.modify!(ref, op::OP, x) where {OP} = Atomix.modify!(ref, op, x, seq_cst)

@inline Atomix.replace!(ref, expected, desired, order::Ordering = seq_cst) =
    Atomix.replace!(ref, expected, desired, order, order)

@inline Atomix.swap!(ref, x, order::Ordering = seq_cst) =
    first(Atomix.modify!(ref, right, x, order))
