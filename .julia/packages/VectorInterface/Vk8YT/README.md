# VectorInterface


| **Documentation** | **Build Status** | **License** |
|:-----------------:|:----------------:|:-----------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][aqua-img]][aqua-url] [![CI][github-img]][github-url] [![][codecov-img]][codecov-url] | [![license][license-img]][license-url] |

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://jutho.github.io/VectorInterface.jl/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://jutho.github.io/VectorInterface.jl/stable

[github-img]: https://github.com/Jutho/VectorInterface.jl/workflows/CI/badge.svg
[github-url]: https://github.com/Jutho/VectorInterface.jl/actions?query=workflow%3ACI

[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

[codecov-img]: https://codecov.io/gh/Jutho/VectorInterface.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/Jutho/VectorInterface.jl

[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat
[license-url]: LICENSE.md

---

An Julia interface proposal for vector-like objects.

### Motivation

This package proposes a Julia interface for vector-like objects, i.e. a small set of methods that should be supported by types representing objects that can be interpreted as living in a vector space.

Having such a unified interface is useful, because this small collection of methods can be used to write a large collection of algorithms in a universal and agnostic manner. Unfortunately, the current interfaces existing in the Julia ecosystem are not adequate for this purpose.

### Context

Recall the basic properties of vectors in a given vector space. A vector space is a set of objects, called vectors, with the basic properties that:

* Vectors can be added (which is commutative, there is a neutral element and inverses)
* The neutral element is called the zero vector
* Vectors can be rescaled with a scalar coefficient, taking values in an underlying scalar field. (and there are all kind of relations with vector addition, such as distributivity, ...)
Together, these two operations give rise to the construction of linear combinations. On top of those properties, many vector spaces used in practical applications have more structure.
* Often there is an inner product, or at least a norm.
* Finally, many useful vector spaces admit a finite basis, such that any vector can be written as a linear combination of a finite number of basis vectors.

Many quantities in science and engineering (and thus in scientific computing) behave as vectors, typically with real or complex numbers as underlying field. Even for quantities that do not (but rather belong to a manifold), their derivatives (tangents and cotangents) do, which is important for automatic differentiation.

More importantly, many algorithms can be formulated using the basic operations (linear combinations and inner products), in combination with some recipe for spitting out new vectors (e.g. applying a linear map): gradient optimization, ODE integrators, Krylov methods, ...

### Current situation and problems in Julia

The most elementary Julia type that acts as a vector is `Vector{T<:Number}`, and by extension (almost *) any subtype of `AbstractArray{T<:Number}`. However, many other types which are not subtypes of `AbstractArray` also are vectors conceptually, e.g. types for representing functions (ApproxFun ecoystem), `Tangent` types in the AD ecosystem, and many more.

However, I have several gripes with the Julia interface (or lack thereof) to access the basic vector operations, which makes writing generic algorithms hard. In particular:

1.  Vector addition and vector rescaling almost always works simply as `v + w` and `v * α`. However, for efficiency, we often want to be able to do these operations in place. For instances of `AbstractArray`, we then need to be `using LinearAlgebra` (which also includes a lot of stuff that we might not need). For vector rescaling within the `AbstractArray` family, we can use `rmul!(v, α)` or `lmul!(α, v)` or `mul!(w, v, α)`, but that interface is conflated with the concept of matrix multiplication etc, whereas these are two very different concepts. There used to be a `scale!` method (in analogy with a corresponding method in BLAS level 1) in older Julia 0.x versions (though it was also not perfect, and was also used for matrix times diagonal matrix). For vector addition on the other hand, the options for  in-place methods are the cryptic `axpy!` and `axpby!` methods (referring to their BLAS level 1 analogues), with thus a very un-Julian interface (the vector that is modified is the final rather than the first argument).

2.  In programming, rather than the scalar field (reals or complex), we of course want to know the specific scalar type, with which the vectors can be natively rescaled. For an instance `v` of `AbstractArray{T<:Number}`, this is the type `T` and it can be obtained as `eltype(v)`. However, because `eltype` is also used by the iteration interface, other types which might have an iteration behaviour that is distinct from their vector-like behaviour, cannot overload `eltype` for both purposes. An example from `Base` would be a nested array, e.g. the type `Vector{Vector{T<:Number}}` still constitutes a set of vectors with scalar type `T`, but `eltype` equal to `Vector{T}`

3.  To get the zero vector associated with some vector `v`, we can use `zero(v)`, which is fine as an interface, as it is defined to be the neutral element with respect to `+`. However, `zero` on e.g. a nested array will fail because of how it is implemented. Furthermore, to make a zero vector in place, you could use `fill!(v, 0)` for `v::AbstractArray{T<:Number}`, but that is a very `AbstractArray` specific interface. The only more general solution is to restort to scaling by zero, e.g. using `rmul!(v, 0)`, if available.

4.  Closely related to the previous two points, we often want to be able to create equivalent vectors but with a modified scalar type. For vectors in `AbstractArray{T<:Number}`, there is `similar(v, T′)` for this, but that is again a very array-specific method, and fails once more for e.g. nested arrays.

5.  Most (but not all) vector-like objects in programming belong to a finite-dimensional vector space. For `v::AbstractArray{T<:Number}`, this dimension is given by `length(v)`, but again this interface is also used for the iteration length, and so new types might face an incompatibility as with `eltype`. And for structured arrays, `length(v)` might also not express the vector space dimension, for e.g. `UpperTriangular{T<:Number}`, the natural vector space dimension is `n*(n+1)/2`, not `n*n`.

6.  The inner product and norm corresponds to the Julia methods `LinearAlgebra.dot` and `LinearAlgebra.norm`. Unlike in some of the previous points, `dot` and `norm` natively support nested arrays. However, `dot` is so loose in its implementation, that it also happily computes an inner product between things which are probably not vectors from the same vector space, such as `dot( (1, 2, [3, 4]), [[1], (2,), (3,4)])`. In particular, `dot` and `norm` also accept tuples, whereas tuples do not behave as vectors with respect to the previous methods (`+`, `*`, `zero`).

In summary, the main problem is that there actually is no formal standardized vector interface in Julia, despite its broad potential applicability for writing very generic algorithms. There are standardized interfaces for containers (`AbstractArray`) and for iterators, which have become conflated with a hypothetical vector interface.

### Existing solutions

Different ecosystems have responded to this hiatus in different ways. Several Krylov and optimization packages merely restrict their applicability to instances of `(Abstract)Array{T<:Number}` or even simply `Vector{T<:Number}` (like their Fortran and C analogues would probably do). The DifferentialEquations.jl ecosystem does more or less the same, restricting to `AbstractArray` (if I remember correctly), but provides a bunch of packages such as `RecursiveArrayTools.jl` and `ArrayInterface.jl` to accommodate for more complex use cases. Finally, the AD ecosystem (Zygote.jl and ChainRules.jl) use custom `Tangent` types for which they define the necessary operations, using a lot of internal machinery to destructure custom types.

Forcing everything to be a subtype of `AbstractArray` is an unsatisfactory solution. Some vector like objects might not be naturally represented with respect to a basis, and thus have no notion of indexing, and might not even be finite-dimensional. The `AbstractArray` or container interface is and should be distinct from a general vector (space) interface.

### New solution
With VectorInterface.jl, I have tried to create a simple package to resolve my gripes. As I hope that I am not alone with those, I would like this to be useful for the community and could eventually evolve into a standardized interface. Therefore, I would very much value comments. Everything is up for bikeshedding. I tried to come up with a design which is compatible with `LinearAlgebra` (e.g. not stealing names) and does not commit type piracy. Currently, VectorInterface.jl provides the following minimalistic interface:

*   `scalartype(v)`: accesses the native scalar type `T<:Number` of a vector-like object; also works in the type domain (i.e. `scalartype(typeof(v))`).
*   `zerovector(v)`, `zerovector!(v)` and `zerovector!!(v)`: produce a zero vector of the same type as `v`; the second method tries to do this in-place for mutable types, whereas the third method is inspired by BangBang.jl and tries to do it in-place when possible, and out of place otherwise.

    *Comment: Ideally, the `zerovector` functionality would be provided by `Base.zero`.*

*   `zerovector(v, S<:Number)` creates a zero vector of similar type, but with a modified scalar type that is now given by `S`. `zerovector!(v)` sets `v` equal to the zero vector; in that case the scalar type cannot be changed. Finally `zerovector!!(v, S)` works, but reduces to `zerovector(v, S)` if `scalartype(v)` is not `S`.

    *Comment: Given that there is a tendency to zero out uninitialized memory, I think it is fine to merge the concept of constructing a new vector with different scalar type with that of constructing the zero vector.*

*   `scale(v, α)`, `scale!(v, α)` and `scale!!(v, α)` rescale the vector `v` with the scalar coefficient `α`. The second method tries to do this in place, but will fail if `α` cannot be converted to `scalartype(v)` (or if `v` contains immutable contents), whereas the third method is the BangBang-style unified solution. There is also `scale!(w, v, α)` and `scale!!(w, v, α)` to rescale `v` out of place, but storing the result in `w`.

*   `add(w, v, [, α = 1, β = 1])`, `add!(w, v, [, α = 1, β = 1])` and `add!!(w, v, [, α = 1, β = 1]))` compute `w * β + v * α`, where (by now self-explanatory) the second method stores the result in `w` (and errors if not possible), and the third method tries to store in `w` but doesn't error.

*   `inner(v, w)` works almost equivalent to `dot(v, w)`, is sometimes a bit more efficient, but also more strict in what arguments it allows.

*   `norm(v)` simply reexports `LinearAlgebra.norm`

These methods are implemented for instances `v` of type `<:Number` (scalars are also vectors over themselves) and `<:AbstractArray` (both `<:AbstractArray{<:Number}` and nested array).

In addition, the interface is currently also defined for tuples and named tuples, again with arbitrary nesting. So instances of e.g. `Vector{NTuple{3,Matrix{Float64}}}` are currently also supported.

Furthermore, in-place methods will work recursively so as to try to maximally work in place. What I mean by that is, if you have nested vectors, say `v = [[1., 2.], [3., 4.]]`, then `rmul!(v, 0.5)` will keep the outer array, but will work on the inner arrays using regular `*` multiplication, and will thus allocate two new inner arrays in this example. In contrast, `scale!(v, 0.5)` does work in-place on the inner arrays and is free of allocations.

Similarly, for `v = ((1., 2.), [3., 4.])`, `scale!!(v, 0.5)` could of course not work in-place on the outer tuple or inner tuple, but would still work in-place on the inner array. Hence, the return value of `scale!!(v, 0.5)` is of course `((0.5, 1.), [1.5, 2.])`, but after his operation, `v` would be `((1., 2.), [1.5, 2.])`.

### Open questions

Upon development, there were several decisions to be made, about which I was not sure.
* Should tuples actually be supported? In Base, they are not treated as vectors, though they are supported by `dot` and `norm`.
  -> Currently, tuples and even named tuples are supported (which was a request)
* Should there be some `vectordim` function (name up to debate) to probe the vector space dimension?
  -> No such function is currently defined
* Should `LinearAlgebra.dot` be exported? Or should `inner` just become `dot` and should and is the looseness of `dot` of no concern?
  -> `dot` is currently not exported; `inner` is the preferred method for inner products
* Should there be fallbacks in place for user types that did already implement `rmul!`, `mul!`, `axpy!`, `axpby!` and `dot` (the latter relating to the previous question)?
  -> Currently, there is a fallback to these methods, but with the intention to remove these fallbacks in future versions.

All thoughts and comments are very welcome.

(*) There is one (actually two) subtype of `AbstractArray` in `LinearAlgebra` that does not behave as a vector, in the sense that its instances cannot represent the zero vector or cannot be rescaled without changing its type, namely `UnitUpperTriangular` and `UnitLowerTriangular`. The fixed unit diagonal prevents these types from constituting a vector space. It seems like the unit diagonal also poses issues for broadcasting, as operations that are preserving ones are much more rare than operations preserving zeros (which is necessary for any structured or unstructured sparseness).
