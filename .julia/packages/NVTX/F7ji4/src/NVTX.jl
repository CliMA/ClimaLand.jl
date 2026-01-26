module NVTX

import Colors, Libdl
using NVTX_jll, JuliaNVTXCallbacks_jll

const NVTX_ACTIVE = Ref{Bool}(false)

"""
    NVTX.isactive()

Determine if Nsight Systems profiling is currently active.
"""
isactive() = NVTX_ACTIVE[]

"""
    NVTX.activate()

Activate the NVTX APIs. This is called automatically when importing NVTX.jl while
running under Nsight profilers. It can also be called manually to enable NVTX when
running under a different profiler that is compatible with NVTX. Note that in such
cases, you may have to manually set the `NVTX_INJECTION64_PATH` environment variable
and have it point to a library that can handle the NVTX APIs.
"""
function activate()
    NVTX_ACTIVE[] = true
    initialize()
    name_threads_julia()
    callbacks = split(get(ENV, "JULIA_NVTX_CALLBACKS", ""), [',','|'])
    enable_gc_hooks(;
        gc="gc" in callbacks,
        alloc="alloc" in callbacks,
        free="free" in callbacks
    )
    enable_inference_hook("inference" in callbacks)
end

function __init__()
    if haskey(ENV, "NVTX_INJECTION32_PATH") || haskey(ENV, "NVTX_INJECTION64_PATH")
        activate()
    end
end

include("api.jl")
include("julia.jl")
include("macro.jl")

end # module
