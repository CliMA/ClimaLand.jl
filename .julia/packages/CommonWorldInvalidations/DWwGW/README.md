# CommonWorldInvalidations.jl

CommonWorldInvalidations.jl is a package that fixes the common unfixable invalidators by simply forcing
the invalidations. It's made to be an unchanged dependency that is then reused by downstream libraries
with `PrecompileTools.@recompile_invalidations` so the common invalidators are all handled during the
precompilation stage a single time.

## What is this package doing and why is it necessary?

The reason for this package essentially comes down to over-aggressive world-splitting optimizations. There's
multiple sources on this optimization:

* https://discourse.julialang.org/t/avoiding-vectors-of-abstract-types/61883/20
* https://discourse.julialang.org/t/how-is-it-that-new-julia-programmers-tend-to-abuse-type-annotations/108465/19
* https://discourse.julialang.org/t/does-julia-create-a-1-5-language-problem/107984/110
* https://discourse.julialang.org/t/static-jl-vs-staticnumbers-jl/87228/21

Basically what happens is that Julia's base image specializes on all of the dispatches it sees in the
world that it builds. So for example, in Julia's Base image, you see that `<(x,y)` always returns a Bool.
What the world-splitting optimization does then is that if `x::Any` and `y::Any` are uninferred variables
in some function, the optimization produces a code which does `<(x,y)::Bool`, i.e. it assumes that the
`<` operation will always return a `Bool` and specializes all generated code on this "fact".

However, enter a library like Symbolics.jl which adds a method `<(x::Num, y)::Num`, i.e. `x < y` is
represented symbolically rather than eagerly evaluated into a Bool. This breaks the world-splitting
assumptions and thus every single code that assumed `<` would output a Bool has to be recompiled.
However, you can see that it's not only Symbolics.jl that does this, but Static.jl, TaylorModels.jl,
tracers defined in JuMP, ..., there's a huge list of libraries that can trigger this invalidation. 

**The issue isn't that these packages are abnormal, it's that Julia's Base image is really abnormal!**
**It's a world with effectively no standard Julia code and no standard Julia package**
**and the Base image is built to hyper-specialize on that exact set of code**

So what can we do about this? Enter CommonWorldInvalidations.jl. CommonWorldInvalidations.jl collects the common
invalidation paths seen throughout the Julia ecosystem, which includes `<` assuming Bool, but also
`||`, `>`, etc. every other Boolean operator, and forces the de-specialization of this optimization
by simply defining a few more methods. Packages can then reply on CommonWorldInvalidations.jl and use
`@recompile_invalidations` which forces the Julia package image builder to recompile all code
invalidated through these changes, which then causes the world post-import to be consistent.

## What happens if you try to do this without CommonWorldInvalidations.jl?

Say for example you don't rely on CommonWorldInvalidations.jl and you do this invalidation in a library
like Symbolics. First of all, if you don't do anything about invalidations you will simply get
bad "time to first x" (TTFX) issues, i.e. after `using Symbolics` you will notice your REPL is
slower (since it relies on uninferred Bool evaluations in some places), Pkg is slower, ... the
entire Julia process is a mess. So what you can do is use PrecompileTools.@recompile_invalidations
on the methods which define these operations. This will force the Julia package image builder
to recompile all code that gets invalidated by this operation, thus making it all happen at
precompilation time, and it will get stored into your Symbolics package image. This means after
a bit more precompilation time, your Pkg, REPL, etc. is back to being swift again. Neat, problem
solved.

However, if you do a change to the code in Symbolics, say improve the `simplify` function
and do a release, the Julia package image system will see that Symbolics updated and thus
redo precompilation. In this precompilation redo you can force this de-specialization, which will
thus cause Pkg, REPL, etc. code to all recompile again. But why is it recompiling if you didn't
actually change any of the code that is related to this invalidation fix? This is because the package
image builder works at the package level, not the method level, and thus it will always recompile
the whole package, not just the parts that changed.

Another issue is that if you do this to Symbolics.jl but someone else separately handles this in
Static.jl, then `using Symbolics, Static` vs `using Static, Symbolics` can cause the way invalidations
are handled to trigger slightly differently, which can force more precompilation if both packages
`@recompile_invalidations`. I am unsure if this is a bug or a necessary requirement in the package
compilation system in its current design.

But whatever the matter is, it's very clear that `@recompile_invalidations` is simply being used at
too high of a level in those scenarios. What we really want is to simply trigger the invalidations
we know we will have in Pkg, the REPL, etc., force the construction of a new image, and keep that
around for the future. And since we know "most" (according to JuliaHub statistics, around 1300 out of
6000 for Static.jl alone) packages will hit these same invalidators, we might as well put it as the 
very bottom of the dependency chain so that this only happens once, and if this package never changes 
(or almost never changes), this precompilation only happens the first time you install Julia packages. 
All subsequent updates can then reuse the world that comes after this invalidation fix.

That of course leads directly to CommonWorldInvalidations.jl.

## Why not just de-specialize these functions in Base?

That is another possible solution. That solution requires convincing everyone that whatever these
1000 packages want is something that should always be de-specialized, which may make the REPL 10% slower
in cases where these packages aren't used. Of course I throw out 10% as a number I haven't actually measured.
The point is, you'd have to go measure a ton of things and convince a ton of people on each method
despecialization, which may take a lot of time and effort. We would like to see this change happen in Base,
but until that happens, CommonWorldInvalidations.jl effectively renders this a non-issue and we can 
move on with our lives with this simple little hack in place.