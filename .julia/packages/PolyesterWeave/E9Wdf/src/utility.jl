"""
    disable_polyester_threads(f::F)

A context manager function that disables Polyester threads without affecting the scheduling
of `Base.Treads.@threads`. Particularly useful for cases when Polyester has been used to
multithread an inner small problem that is now to be used in an outer embarassingly parallel
problem (in such cases it is best to multithread only at the outermost level).

This call will disable all PolyesterWeave threads, including those in
`LoopVectorization.@tturbo` and `Octavian.matmul`.

Avoid calling it as `@threads for i in 1:n disable_polyester_threads(f) end` as that creates
unnecessary per-thread overhead. Rather call it in the outermost scope, e.g. as given here:

```
disable_polyester_threads() do
    @threads for i in 1:n f()
end
```
"""
function disable_polyester_threads(f::F) where {F}
  _, r = request_threads(Threads.nthreads())
  try
    f()
  finally
    foreach(free_threads!, r)
  end
end
