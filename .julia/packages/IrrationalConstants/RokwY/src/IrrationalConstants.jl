module IrrationalConstants

export
    twoπ,       # 2π
    fourπ,      # 4π
    halfπ,      # π / 2
    quartπ,     # π / 4
    invπ,       # 1 / π
    twoinvπ,    # 2 / π
    fourinvπ,   # 4 / π
    inv2π,      # 1 / (2π)
    inv4π,      # 1 / (4π)
    sqrt2,      # √2
    sqrt3,      # √3
    sqrtπ,      # √π
    sqrt2π,     # √2π
    sqrt4π,     # √4π
    sqrthalfπ,  # √(π / 2)
    invsqrt2,   # 1 / √2
    invsqrtπ,   # 1 / √π
    invsqrt2π,  # 1 / √2π
    loghalf,    # log(1 / 2)
    logtwo,     # log(2)
    logten,     # log(10)
    logπ,       # log(π)
    log2π,      # log(2π)
    log4π       # log(4π)

include("macro.jl")
include("stats.jl")
include("trigonometric.jl")

end # module
