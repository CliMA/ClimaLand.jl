function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    eltypes = (N0f8, N0f16, Float32, Float64)        # eltypes of parametric colors
    pctypes = (Gray, RGB, AGray, GrayA, ARGB, RGBA)  # parametric colors
    cctypes = (Gray24, AGray32, RGB24, ARGB32)       # non-parametric colors
    realtypes = (Float16, Float32, Float64, Int)     # types for mixed Normed/Real operations
    for R in realtypes, T in eltypes, C in pctypes
        precompile(Tuple{typeof(+),C{T},C{T}})
        precompile(Tuple{typeof(-),C{T},C{T}})
        precompile(Tuple{typeof(*),R,C{T}})
        precompile(Tuple{typeof(*),C{T},R})
        precompile(Tuple{typeof(/),C{T},R})
    end
end
