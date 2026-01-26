include("testsetup.jl")

include("test-suite-preamble.jl")

include("internals.jl")
include("threadingutilities.jl")
if Threads.nthreads() > 3
  include("staticarrays.jl")
end
include("threadpool.jl")
include("warnings.jl")

include("aqua.jl") # run the Aqua.jl tests last

