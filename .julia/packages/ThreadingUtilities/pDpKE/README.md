# ThreadingUtilities

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSIMD.github.io/ThreadingUtilities.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSIMD.github.io/ThreadingUtilities.jl/dev)
[![Continuous Integration](https://github.com/JuliaSIMD/ThreadingUtilities.jl/workflows/CI/badge.svg)](https://github.com/JuliaSIMD/ThreadingUtilities.jl/actions?query=workflow%3ACI)
[![Continuous Integration (Julia nightly)](https://github.com/JuliaSIMD/ThreadingUtilities.jl/workflows/CI%20(Julia%20nightly)/badge.svg)](https://github.com/JuliaSIMD/ThreadingUtilities.jl/actions?query=workflow%3A%22CI+%28Julia+nightly%29%22)
[![Coverage](https://codecov.io/gh/JuliaSIMD/ThreadingUtilities.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSIMD/ThreadingUtilities.jl)

Utilities for low overhead threading in Julia.

Please see the [documentation](https://JuliaSIMD.github.io/ThreadingUtilities.jl/stable).

If you're using Windows, please note that Windows often allocates memory when neither Mac or Linux do. I do not know why. If you can help diagnose/fix the problem, please take a look at `count_allocated()` in `/test/staticarrays.jl`.

If you are trying to mix code relying on `ThreadingUtitlies.jl` (e.g. `LoopVectorization.jl`'s threaded code, `Polyester.jl`, or `Octavian.jl`) with other threaded code, you can use `ThreadingUtilities.sleep_all_tasks()` to put `ThreadingUtilities`'s tasks to sleep before switching from `ThreadingUtilities` to the other threaded code. Manually putting `ThreadingUtilities`'s tasks to sleep in this way should stop them from interfearing.

