using PrecompileTools

let
@setup_workload begin
    function pcarray(f::F, ::Type{A}, sz) where {F,A}
        a = f(A(undef, sz))
        fill!(a, zero(eltype(a)))
        return first(a)
    end
    function pcm(a1, a2; fillvalue = zero(eltype(a1)))
        v = mosaic(a1, a2; fillvalue=fillvalue)
        return first(v)
    end
    function pcmv(a; fillvalue = zero(eltype(a)))
        v = mosaicview(a; fillvalue=fillvalue)
        return first(v)
    end

    eltypes = (N0f8, N0f16, Float32, Float64)    # eltypes of parametric colors
    pctypes = (Gray, RGB)                        # parametric colors
    cctypes = ()                                 # non-parametric colors (e.g., Gray24)
    dims  = (1, 2, 3, 4)
    szs   = ((2,), (2, 2), (2, 2, 2), (2, 2, 2, 2))
    @compile_workload begin
        for T in eltypes
            clamp01(zero(T))
            clamp01nan(zero(T))
            scaleminmax(zero(T), oneunit(T))
            scalesigned(oneunit(T))
            scalesigned(zero(T), oneunit(T) / 2, oneunit(T))
            for C in pctypes
                clamp01(zero(C{T}))
                clamp01nan(zero(C{T}))
                colorsigned(zero(C{T}), oneunit(C{T}))
            end
        end
        for C in cctypes
            clamp01(zero(C))
            clamp01nan(zero(C))
            C === AGray32 && continue
            colorsigned(zero(C), oneunit(C))
        end
        # For the arrays, it's better to make them and exercise them so we get the getindex/setindex!
        # methods precompiled too.
        for sz in szs
            for T in eltypes
                T <: FixedPoint || continue
                pcarray(rawview, Array{T,length(sz)}, sz)
            end
            pcarray(normedview, Array{UInt8,length(sz)}, sz)
            pcarray(a->normedview(N0f8, a),  Array{UInt8, length(sz)}, sz)
            pcarray(a->normedview(N0f16, a), Array{UInt16,length(sz)}, sz)
        end
        for sz in szs
            for T in eltypes
                for C in pctypes
                    nc = sizeof(C{T}) ÷ sizeof(T)
                    if nc == 1
                        pcarray(a->colorview(C,    a), Array{T,length(sz)}, sz)
                        pcarray(a->colorview(C{T}, a), Array{T,length(sz)}, sz)
                    else
                        pcarray(a->colorview(C,    a), Array{T,length(sz)+1}, (nc, sz...))
                        pcarray(a->colorview(C{T}, a), Array{T,length(sz)+1}, (nc, sz...))
                    end
                    pcarray(channelview, Array{C{T},length(sz)}, sz)
                    if T<:FixedPoint
                        R = FixedPointNumbers.rawtype(T)
                        if nc == 1
                            pcarray(a->colorview(C,    normedview(T, a)), Array{R,length(sz)}, sz)
                            pcarray(a->colorview(C{T}, normedview(T, a)), Array{R,length(sz)}, sz)
                        else
                            pcarray(a->colorview(C,    normedview(T, a)), Array{R,length(sz)+1}, (nc, sz...))
                            pcarray(a->colorview(C{T}, normedview(T, a)), Array{R,length(sz)+1}, (nc, sz...))
                        end
                    end
                end
            end
            T, C = Bool, Gray
            pcarray(a->colorview(C,    a), Array{T,length(sz)}, sz)
            pcarray(a->colorview(C{T}, a), Array{T,length(sz)}, sz)
            pcarray(channelview, Array{C{T},length(sz)}, sz)
        end
        for T in eltypes
            a = zeros(T, (2, 2))
            a3 = zeros(T, (2, 2, 2))
            pcm(a, a)
            pcmv(a3)
            for C in (Gray, RGB, GrayA, RGBA)
                a = zeros(C{T}, (2, 2))
                a3 = zeros(C{T}, (2, 2, 3))
                pcm(a, a)
                pcmv(a3)
                if C === RGB
                    pcm(a, a; fillvalue=zero(Gray{T}))
                elseif C === RGBA
                    pcm(a, a; fillvalue=zero(GrayA{T}))
                end
            end
        end
        pcm(zeros(Float32, 2, 2), zeros(Float64, 2, 2))  # heterogeneous
    end
end
end
