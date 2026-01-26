# convenience functions for using Julia

"""
    gettid()

Get the system thread ID of the current Julia thread. This is compatible with
[`name_os_thread`](@ref).
"""
function gettid end

@static if Sys.islinux()
    if Sys.ARCH == :x86_64
        gettid() = ccall(:syscall, UInt32, (Clong, Clong...), 186)
    elseif Sys.ARCH == :aarch64
        gettid() = ccall(:syscall, UInt32, (Clong, Clong...), 178)
    elseif Sys.ARCH == :powerpc64le || Sys.ARCH == :ppc64le
        gettid() = ccall(:syscall, UInt32, (Clong, Clong...), 207)
    end
elseif Sys.iswindows()
    gettid() = ccall(:GetCurrentThreadId, UInt32,())
end

"""
    name_threads_julia([namefn])

Name the threads owned by the Julia process using `namefn()`. The default is
`namefn() = "julia thread \$(Threads.threadid())"`.

This function is called at module initialization if the profiler is active, so
should not generally need to be called manually unless a custom name is required.
"""
function name_threads_julia(fn = () -> "julia thread $(Threads.threadid())")
    Threads.@threads :static for _ = 1:Threads.nthreads()
        name_os_thread(gettid(), fn())
    end
end

# domain used for instrumenting Julia runtime
const JULIA_DOMAIN = Domain("Julia")
const GC_MESSAGE = StringHandle(JULIA_DOMAIN, "GC")
const GC_ALLOC_MESSAGE = StringHandle(JULIA_DOMAIN, "alloc")
const GC_FREE_MESSAGE = StringHandle(JULIA_DOMAIN, "free")
const GC_COLOR = Ref{UInt32}(Colors.ARGB32(Colors.colorant"brown").color)
const GC_ALLOC_COLOR = Ref{UInt32}(Colors.ARGB32(Colors.colorant"goldenrod1").color)
const GC_FREE_COLOR = Ref{UInt32}(Colors.ARGB32(Colors.colorant"dodgerblue").color)
const INFERENCE_COLOR = Ref{UInt32}(Colors.ARGB32(Colors.colorant"turquoise3").color)

"""
    NVTX.enable_gc_hooks(;gc=true, alloc=false, free=false)

Add NVTX hooks for the Julia garbage collector:
 - `gc`: instrument GC invocations as ranges,
 - `alloc`: instrument calls to alloc as marks (payload will contain allocation size), and
 - `free`: instrument calls to free as marks.

If the `JULIA_NVTX_CALLBACKS` environment variable is set, this will be
automatically called at module initialization. `JULIA_NVTX_CALLBACKS` should be
either a comma (`,`) or bar (`|`) separated list of callbacks to enable. For
example, setting it to `gc|alloc|free` will enable all hooks. The 
`--env-var` argument can be helpful for setting this variable, e.g.
```sh
nsys profile --env-var=JULIA_NVTX_CALLBACKS=gc|alloc|free julia --project script.jl
```
"""
function enable_gc_hooks(;gc::Bool=true, alloc::Bool=false, free::Bool=false)
    if gc || alloc || free
        init!(JULIA_DOMAIN)
        unsafe_store!(cglobal((:julia_domain,libjulia_nvtx_callbacks),Ptr{Cvoid}), JULIA_DOMAIN.ptr)
    end
    if gc
        init!(GC_MESSAGE)
        unsafe_store!(cglobal((:gc_message,libjulia_nvtx_callbacks),Ptr{Cvoid}), GC_MESSAGE.ptr)
        unsafe_store!(cglobal((:gc_color,libjulia_nvtx_callbacks),UInt32), GC_COLOR[])
        # https://github.com/JuliaLang/julia/blob/v1.8.3/src/julia.h#L879-L883
        name_category(JULIA_DOMAIN, 1+0, "GC auto")
        name_category(JULIA_DOMAIN, 1+1, "GC full")
        name_category(JULIA_DOMAIN, 1+2, "GC incremental")
    end
    if alloc
        init!(GC_ALLOC_MESSAGE)
        unsafe_store!(cglobal((:gc_alloc_message,libjulia_nvtx_callbacks),Ptr{Cvoid}), GC_ALLOC_MESSAGE.ptr)
        unsafe_store!(cglobal((:gc_alloc_color,libjulia_nvtx_callbacks),UInt32), GC_ALLOC_COLOR[])
    end
    if free
        init!(GC_FREE_MESSAGE)
        unsafe_store!(cglobal((:gc_free_message,libjulia_nvtx_callbacks),Ptr{Cvoid}), GC_FREE_MESSAGE.ptr)
        unsafe_store!(cglobal((:gc_free_color,libjulia_nvtx_callbacks),UInt32), GC_FREE_COLOR[])
    end

    ccall(:jl_gc_set_cb_pre_gc, Cvoid, (Ptr{Cvoid}, Cint),
        cglobal((:nvtx_julia_gc_cb_pre,libjulia_nvtx_callbacks)), gc)
    ccall(:jl_gc_set_cb_post_gc, Cvoid, (Ptr{Cvoid}, Cint),
        cglobal((:nvtx_julia_gc_cb_post,libjulia_nvtx_callbacks)), gc)
    ccall(:jl_gc_set_cb_notify_external_alloc, Cvoid, (Ptr{Cvoid}, Cint),
        cglobal((:nvtx_julia_gc_cb_alloc,libjulia_nvtx_callbacks)), alloc)
    ccall(:jl_gc_set_cb_notify_external_free, Cvoid, (Ptr{Cvoid}, Cint),
        cglobal((:nvtx_julia_gc_cb_free,libjulia_nvtx_callbacks)), free)
    return nothing
end

typeinf_ext_nvtx(mi::Base.Core.MethodInstance, world::UInt) = typeinf_ext_nvtx(Base.Core.Compiler.NativeInterpreter(world), mi)
function typeinf_ext_nvtx(interp::Base.Core.Compiler.AbstractInterpreter, linfo::Base.Core.MethodInstance)
    method = linfo.def
    types = linfo.specTypes.parameters[2:end]
    message = "$(method.name)($(join([string("::", t) for t in types], ", "))) @ $(method.module) $(method.file):$(method.line)"
    id = range_start(JULIA_DOMAIN; message, color = INFERENCE_COLOR[], category = 11)
    ret = Core.Compiler.typeinf_ext_toplevel(interp, linfo)
    range_end(id)
    return ret
end
precompile(typeinf_ext_nvtx, (Base.Core.Compiler.NativeInterpreter, Base.Core.MethodInstance))
precompile(typeinf_ext_nvtx, (Base.Core.MethodInstance, UInt))

"""
    NVTX.enable_inference_hook(active::Bool=true)

Add hooks for method inference. Can also be activated by adding `inference` to
the `JULIA_NVTX_CALLBACKS` environment variable.
"""
function enable_inference_hook(enable::Bool=true)
    if enable
        name_category(JULIA_DOMAIN, 11, "compiler inference")
        ccall(:jl_set_typeinf_func, Cvoid, (Any,), typeinf_ext_nvtx)
    else
        ccall(:jl_set_typeinf_func, Cvoid, (Any,), Core.Compiler.typeinf_ext_toplevel)
    end
end
