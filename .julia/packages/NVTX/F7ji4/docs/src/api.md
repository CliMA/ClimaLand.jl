# API

```@meta
CurrentModule = NVTX
```

## Macros

The following macros provide the most convenient usage of this package

```@docs
@mark
@range
@annotate
@category
```

## Julia interactions

```@docs
enable_gc_hooks()
enable_inference_hook
name_threads_julia()
```

## Low-level API

These closely map to the C API

```@docs
isactive
initialize
activate
Domain
StringHandle
mark
range_start
range_end
range_push
range_pop
name_category
gettid
name_os_thread
```