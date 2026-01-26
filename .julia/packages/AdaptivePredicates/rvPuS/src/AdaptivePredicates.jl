module AdaptivePredicates

@static if isdefined(Base, :Memory)
    const Vec{T} = Memory{T}
else
    const Vec{T} = Vector{T}
end

include("init.jl")
include("macros.jl")
include("arithmetic.jl")
include("predicates.jl")

export orient2, orient3, incircle, insphere
export orient2p, orient3p, incirclep, inspherep
export orient2fast, orient3fast, incirclefast, inspherefast

if VERSION ≥ v"1.11.0-DEV.469"
    eval(Meta.parse("public orient3adapt_cache, incircleadapt_cache, insphereexact_cache"))
end

end # module
