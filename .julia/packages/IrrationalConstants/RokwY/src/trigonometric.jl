# Functions return `Float64`, consistent with Base
# https://github.com/JuliaLang/julia/pull/42595
# Values at poles are defined to be consistent with `cot(0)` and `cot(π)`
# https://github.com/JuliaLang/julia/issues/7123
# https://github.com/JuliaLang/julia/blob/e3d366f1966595ba737220df49e220610823b331/base/mathconstants.jl#L130

# `sin`
Base.sin(::Twoπ) = 0.0
Base.sin(::Fourπ) = 0.0
Base.sin(::Halfπ) = 1.0
Base.sin(::Quartπ) = Float64(invsqrt2)

# `cos`
Base.cos(::Twoπ) = 1.0
Base.cos(::Fourπ) = 1.0
Base.cos(::Halfπ) = 0.0
Base.cos(::Quartπ) = Float64(invsqrt2)

# `sincos`
Base.sincos(::Twoπ) = (0.0, 1.0)
Base.sincos(::Fourπ) = (0.0, 1.0)
Base.sincos(::Halfπ) = (1.0, 0.0)
Base.sincos(::Quartπ) = (Float64(invsqrt2), Float64(invsqrt2))

# `tan`
Base.tan(::Twoπ) = 0.0
Base.tan(::Fourπ) = 0.0
Base.tan(::Halfπ) = 1/0
Base.tan(::Quartπ) = 1.0

# `csc`, `sec`, and `cot` are defined automatically, so we do not define them
# there is one exception where we can improve accuracy:
Base.csc(::Quartπ) = Float64(sqrt2)
Base.sec(::Quartπ) = Float64(sqrt2)
