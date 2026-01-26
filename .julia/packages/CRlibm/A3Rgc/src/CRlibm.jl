module CRlibm

    import CRlibm_jll

# All functions in CRlibm except the power function, according to Section 0.4 of
# the PDF manual (page 8)
# Source: ./deps/src/crlibm1.0beta4/docs/latex/0_getting-started.tex
# Section "Currently available functions"

const functions = [:exp, :expm1,
                   :log, :log1p, :log2, :log10,
                   :sin, :cos, :tan, :asin, :acos, :atan,
                   :sinpi, :cospi, :tanpi, :atanpi,
                   :sinh, :cosh]

if VERSION ≥ v"1.10"
    # `sinpi`, `cospi`, `tanpi`, `atanpi` are only available since Julia v1.10
    const MPFR_functions = [:exp, :expm1,
                            :log, :log1p, :log2, :log10,
                            :sin, :cos, :tan, :asin, :acos, :atan,
                            :sinpi, :cospi, :tanpi, :atanpi,
                            :sinh, :cosh]
else
    const MPFR_functions = [:exp, :expm1,
                            :log, :log1p, :log2, :log10,
                            :sin, :cos, :tan, :asin, :acos, :atan,
                            :sinh, :cosh]
end

# The underlying CRlibm tests fail on 32 bit build, use MPFR
is_32_bit = Int == Int32

for f ∈ functions
    crlibm_f_rd = string(f, "_rd")
    crlibm_f_ru = string(f, "_ru")
    crlibm_f_rn = string(f, "_rn")
    crlibm_f_rz = string(f, "_rz")
    crlibm_f = Symbol(:crlibm_, f)
    mpfr_f = Symbol(:mpfr_, f)
    str_mpfr_f = string(:mpfr_, f)

    if !is_32_bit
        @eval $crlibm_f(x::Float16, r::RoundingMode) = Float16($crlibm_f(Float64(x), r), r)
        @eval $crlibm_f(x::Float32, r::RoundingMode) = Float32($crlibm_f(Float64(x), r), r)
        @eval $crlibm_f(x::Float64, ::RoundingMode{:Down})    = ccall(($crlibm_f_rd, CRlibm_jll.libcrlibm), Float64, (Float64,), x)
        @eval $crlibm_f(x::Float64, ::RoundingMode{:Up})      = ccall(($crlibm_f_ru, CRlibm_jll.libcrlibm), Float64, (Float64,), x)
        @eval $crlibm_f(x::Float64, ::RoundingMode{:Nearest}) = ccall(($crlibm_f_rn, CRlibm_jll.libcrlibm), Float64, (Float64,), x)
        @eval $crlibm_f(x::Float64, ::RoundingMode{:ToZero})  = ccall(($crlibm_f_rz, CRlibm_jll.libcrlibm), Float64, (Float64,), x)
    end

    if f ∈ MPFR_functions
        @eval $mpfr_f(x::Float16, r::RoundingMode) = Float16($mpfr_f(Float64(x), r), r)
        @eval $mpfr_f(x::Float32, r::RoundingMode) = Float32($mpfr_f(Float64(x), r), r)
        @eval $mpfr_f(x::Float64, r::RoundingMode) = Float64($mpfr_f(BigFloat(x; precision = precision(x)), r), r)
        @eval function $mpfr_f(x::BigFloat, r::RoundingMode)
            z = BigFloat(; precision = precision(x))
            ccall(($str_mpfr_f, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), z, x, r)
            return z
        end
    end
end

"""
    setup(use_MPFR=false)

Define correctly-rounded standard mathematical functions. See `CRlibm.functions`
for the complete list.

The functions are not exported; for instance, use `CRlibm.sin(0.5, RoundDown)`.

If `use_MPFR` is `true`, then the functions wrap the corresponding MPFR
functionality (i.e., `BigFloat`).
"""
function setup(use_MPFR::Bool=false)
    for f ∈ functions
        mpfr_f = Symbol(:mpfr_, f)
        @eval $f(x::Real) = $f(x, RoundNearest)
        @eval $f(x::BigFloat, r::RoundingMode) = $mpfr_f(x, r)
    end

    if use_MPFR
        @info "CRlibm is shadowing MPFR"
        for f ∈ functions
            mpfr_f = Symbol(:mpfr_, f)
            @eval $f(x::Real, r::RoundingMode) = $mpfr_f(float(x), r)
        end
    else
        for f ∈ functions
            crlibm_f = Symbol(:crlibm_, f)
            @eval $f(x::Real, r::RoundingMode) = $crlibm_f(float(x), r)
        end
    end

    return use_MPFR
end

setup(is_32_bit)

end
